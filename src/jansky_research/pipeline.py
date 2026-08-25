"""Pipeline orchestration: fetch -> analyze -> metrics -> paper inputs.

The single entry point shared by ``make pipeline``, the notebooks, and the Airflow DAG. Keeping
one code path is what makes the result reproducible three ways. Offline-first: with ``offline=True``
(or when the real CHIME catalogue can't be downloaded) it runs on the synthetic fixture, so tests
and CI never touch the network.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from . import data, frbstats, report

__all__ = ["analyze", "build_catalog", "load_catalog_csv", "metrics_dict", "run"]

# Candidate CHIME/FRB Catalog 1 column names -> our canonical keys. The loader is tolerant
# because the public CSV schema has varied across releases; unmatched columns are ignored and
# the real-data run (GATE 2) confirms the mapping.
_COLUMN_ALIASES = {
    "mjd": ("mjd_400", "mjd_inf", "mjd", "bary_mjd_400"),
    "fluence": ("fluence", "fluence_fitb", "fluence_jy_ms"),
    "dm": ("dm_fitb", "dm_exc_ne2001", "dm", "dm_obs"),
    "width": ("width_fitb", "bc_width", "width", "width_ms"),
    "repeater_name": ("repeater_name", "repeater_of", "previous_name"),
    "sub_num": ("sub_num", "sub_burst"),
    "excluded_flag": ("excluded_flag",),
}


def _first_present(header: list[str], names: tuple[str, ...]) -> str | None:
    lower = {h.lower(): h for h in header}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def load_catalog_csv(path: str | Path) -> dict[str, np.ndarray]:
    """Load a CHIME/FRB catalogue CSV into the canonical catalogue dict.

    Returns arrays ``mjd``, ``fluence``, ``dm``, ``width`` and a boolean ``repeater`` mask (a row
    is a repeater burst when its repeater-name field is set and not a null marker like ``-``).
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"empty catalogue: {path}")
    header = list(rows[0].keys())
    cols = {key: _first_present(header, names) for key, names in _COLUMN_ALIASES.items()}

    def _floats(colname: str | None) -> np.ndarray:
        if colname is None:
            return np.full(len(rows), np.nan)
        out = np.empty(len(rows))
        for i, r in enumerate(rows):
            try:
                out[i] = float(r[colname])
            except (TypeError, ValueError):
                out[i] = np.nan
        return out

    # Non-repeaters carry a null marker in the repeater-name field; in CHIME/FRB Catalog 1 that
    # marker is "-9999" (other catalogues use "-" or blank).
    sentinels = {"-9999", "-9999.0", "-", "--", "", "nan"}
    rep_col = cols["repeater_name"]
    names = (
        np.array([(r.get(rep_col, "") or "").strip() for r in rows])
        if rep_col
        else np.array([""] * len(rows))
    )
    repeater = np.array([nm not in sentinels for nm in names])
    # Width upper limits ("<0.0001" strings) parse to NaN and are silently dropped by the
    # KS comparison; count them per class so the cut is disclosed (all 26 in Cat 1 are
    # non-repeaters, so dropping them is conservative for the width claim).
    wcol = cols["width"]
    is_limit = (
        np.array([str(r.get(wcol, "")).strip().startswith("<") for r in rows])
        if wcol
        else np.zeros(len(rows), bool)
    )
    cat = {
        "mjd": _floats(cols["mjd"]),
        "fluence": _floats(cols["fluence"]),
        "dm": _floats(cols["dm"]),
        "width": _floats(cols["width"]),
        "repeater": repeater,
        "repeater_name": names,
        "width_is_limit": is_limit,
        "excluded_flag": _floats(cols["excluded_flag"]),
    }
    # One row per *event*: the CHIME catalogue stores each multi-component burst as several
    # sub_num rows (600 rows = 536 events). Treating sub-bursts as independent would
    # pseudo-replicate near-identical DMs and inflate KS significance, so keep sub_num == 0.
    if cols["sub_num"] is not None:
        sub = _floats(cols["sub_num"])
        keep = ~(sub > 0)  # sub_num == 0 or missing
        cat = {k: v[keep] for k, v in cat.items()}
    return cat


def build_catalog(*, offline: bool = False) -> tuple[dict[str, np.ndarray], str]:
    """Return ``(catalog, source)`` — the real CHIME catalogue, or the synthetic fixture offline.

    A fetch failure RAISES rather than silently substituting the synthetic fixture: the old
    fallback could rewrite the paper's macros with synthetic values on any network hiccup
    while every guard stayed green (the round-9 referee's split-brain finding).
    """
    if not offline:
        try:
            path = data.fetch("chime-frb-catalog")
        except Exception as exc:
            raise RuntimeError(
                "could not fetch the CHIME/FRB catalogue; pass --offline explicitly to run "
                "on the synthetic fixture"
            ) from exc
        return load_catalog_csv(path), "chime-frb-catalog"
    return frbstats.synthetic_catalog(), "synthetic"


def metrics_dict(stats: frbstats.BurstStats, source: str) -> dict:
    """A flat, JSON-serialisable summary of the analysis (read by the paper's macros)."""
    return {
        "source": source,
        "n_bursts": stats.n_bursts,
        "n_repeater_bursts": stats.n_repeater_bursts,
        "n_repeater_sources": stats.n_repeater_sources,
        "weibull": asdict(stats.weibull),
        "weibull_clustered": stats.weibull.clustered,
        "energy": asdict(stats.energy),
        "ks": stats.ks,
    }


def analyze(catalog: dict[str, np.ndarray], source: str = "unknown") -> dict:
    """Run the burst-statistics analysis and return the metrics dict.

    Beyond the burst-level summary this commits the round-9 statistics: the source-level KS
    restatement (the source, not the burst, is the unit of analysis), the joint
    gamma/f_min bootstrap, the Cat 1 comparand, the source-cluster Weibull CI, the wait-time
    cadence structure, and the disclosed sample cuts.
    """
    stats_ = frbstats.summarise(catalog)
    metrics = metrics_dict(stats_, source)

    rep_mask = np.asarray(catalog["repeater"], dtype=bool)
    rep: dict[str, np.ndarray] = {
        k: np.asarray(catalog[k])[rep_mask] for k in ("dm", "fluence", "width") if k in catalog
    }
    one: dict[str, np.ndarray] = {
        k: np.asarray(catalog[k])[~rep_mask] for k in ("dm", "fluence", "width") if k in catalog
    }
    if "repeater_name" in catalog:
        labels = np.asarray(catalog["repeater_name"])[rep_mask]
        metrics["ks_source_level"] = frbstats.compare_populations_source_level(rep, one, labels)
        rep_mjd = np.asarray(catalog["mjd"])[rep_mask]
        lo, hi = frbstats.weibull_cluster_ci(rep_mjd, labels)
        metrics["weibull_k_cluster_ci"] = [round(lo, 2), round(hi, 2)]
        # cadence structure of the within-source waits, and the sources that contribute any
        waits = frbstats.grouped_wait_times(rep_mjd, labels)
        metrics["n_wait_sources"] = int(
            sum(1 for g in np.unique(labels) if frbstats.wait_times(rep_mjd[labels == g]).size)
        )
        metrics["n_waits_same_transit"] = int(np.sum(waits <= 0.02))
        metrics["n_waits_sidereal_day"] = int(np.sum((waits >= 0.95) & (waits <= 1.05)))

    # honest gamma error: joint bootstrap over (gamma, f_min)
    metrics["energy_bootstrap"] = frbstats.bootstrap_power_law(np.asarray(catalog["fluence"]))
    # the comparand: Cat 1's own cumulative slope, transformed to a differential index
    gamma_ref = 1.0 + abs(frbstats.CHIME_CAT1_ALPHA_CUM)
    gamma_ref_err = frbstats.CHIME_CAT1_ALPHA_CUM_ERR
    g = metrics["energy"]["gamma"]
    g_sd = max(metrics["energy_bootstrap"]["gamma_boot_sd"], metrics["energy"]["gamma_err"])
    metrics["gamma_ref"] = gamma_ref
    metrics["gamma_ref_err"] = gamma_ref_err
    metrics["gamma_agreement_sigma"] = round(
        abs(g - gamma_ref) / float(np.hypot(g_sd, gamma_ref_err)), 2
    )
    # disclosed sample cuts and variants
    if "width_is_limit" in catalog:
        wl = np.asarray(catalog["width_is_limit"], bool)
        metrics["n_width_limits_rep"] = int(np.sum(wl & rep_mask))
        metrics["n_width_limits_one"] = int(np.sum(wl & ~rep_mask))
    if "excluded_flag" in catalog:
        ex = np.asarray(catalog["excluded_flag"], float)
        ok = ~(ex == 1.0)
        if (~ok).any():
            fit_noexcl = frbstats.fit_power_law(np.asarray(catalog["fluence"])[ok], auto_xmin=True)
            metrics["n_excluded_flag"] = int((~ok).sum())
            metrics["gamma_without_excluded"] = round(fit_noexcl.gamma, 2)
            metrics["f_min_without_excluded"] = round(fit_noexcl.f_min, 1)
    return metrics


def run(out_dir: str | Path = ".", *, offline: bool = False) -> dict:
    """Full pipeline: build the catalogue, analyse it, and write the paper inputs.

    Writes ``results/metrics.json``, the figures under ``papers/frbstats/figures/``, and
    ``papers/frbstats/generated/macros.tex``. Returns the metrics dict.
    """
    out = Path(out_dir)
    catalog, source = build_catalog(offline=offline)
    stats = frbstats.summarise(catalog)
    metrics = analyze(catalog, source)

    paper = out / "papers" / "frbstats"
    (out / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, out / "results" / "metrics.json")
    report.make_figures(catalog, stats, paper / "figures")
    report.write_macros(metrics, paper / "generated" / "macros.tex")
    return metrics


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the FRB burst-statistics pipeline.")
    parser.add_argument("--out", default=".", help="output root (default: repo root)")
    parser.add_argument("--offline", action="store_true", help="use the synthetic fixture")
    args = parser.parse_args(argv)
    metrics = run(args.out, offline=args.offline)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
