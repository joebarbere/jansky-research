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
    cat = {
        "mjd": _floats(cols["mjd"]),
        "fluence": _floats(cols["fluence"]),
        "dm": _floats(cols["dm"]),
        "width": _floats(cols["width"]),
        "repeater": repeater,
        "repeater_name": names,
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
    """Return ``(catalog, source)`` — the real CHIME catalogue, or the synthetic fixture offline."""
    if not offline:
        try:
            path = data.fetch("chime-frb-catalog")
            return load_catalog_csv(path), "chime-frb-catalog"
        except Exception:  # noqa: BLE001 - any failure falls back to the offline fixture
            pass
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
    """Run the burst-statistics analysis and return the metrics dict."""
    return metrics_dict(frbstats.summarise(catalog), source)


def run(out_dir: str | Path = ".", *, offline: bool = False) -> dict:
    """Full pipeline: build the catalogue, analyse it, and write the paper inputs.

    Writes ``results/metrics.json``, the figures under ``paper/figures/``, and
    ``paper/generated/macros.tex``. Returns the metrics dict.
    """
    out = Path(out_dir)
    catalog, source = build_catalog(offline=offline)
    stats = frbstats.summarise(catalog)
    metrics = metrics_dict(stats, source)

    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "results" / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    report.make_figures(catalog, stats, out / "paper" / "figures")
    report.write_macros(metrics, out / "paper" / "generated" / "macros.tex")
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
