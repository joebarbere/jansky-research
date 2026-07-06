"""CHIME/FRB Catalog 2: the first uniform repeater wait-time & duty-cycle census (plan 39).

Catalog 2 (arXiv:2601.09399, ApJS doi:10.3847/1538-4365/ae3828; 4,539 bursts, 83 repeaters)
defers all repeater timing statistics to its rate-census companion (Cook et al.,
arXiv:2605.08410) --- which computes rates, DM drifts, and the repeater fraction but **no
per-source Weibull shape, no periodograms, no duty cycles** (GATE-0 full-text pass, 2026-07-06).
All published periodicity results are single-source campaigns. This module computes ONE uniform
statistic across every Cat-2 repeater above a stated completeness cut: Weibull clustering shape
k (reusing `frbstats.fit_weibull_waits`), a Rayleigh-Z^2 periodogram (reusing
`frbperiod.rayleigh_z2`), and an activity-window/duty-cycle estimate at any significant period.

The Cat-1 lesson (`survey/period-findings.md`: exposure-blind FAPs are not rigorous) drives the
null: CHIME is a transit instrument, so burst TOAs live on a comb of period one SIDEREAL day
(0.99727 d). The released exposure product (CANFAR DOI 10.11570/25.0066,
`chimefrbcat2_exposure.h5`) is time-INTEGRATED (two nside-4096 HEALPix maps, upper/lower
transit), so a per-epoch exposure correction is not possible from public data --- stated
honestly. Instead the periodogram FAP comes from a **transit-comb-preserving scramble**: each
burst keeps its sidereal phase and gets a uniformly redrawn sidereal-day number across the
source's observed span. That preserves the daily-visibility comb exactly and destroys
multi-day structure; per-source total exposure (the catalogue's own `exp_up` hours) normalises
rates. Anchor: FRB 20180916B's 16.35-d period (CHIME/FRB 2020, Nature 582, 351) must be
re-found. NOTE: FRB 20240114A and 20240209A post-date the Cat-2 cutoff (2023-09-15) and are NOT
in the table --- the plan's second anchor is out of scope, documented in the findings.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .frbstats import fit_weibull_waits

__all__ = [
    "load_catalog2",
    "repeater_trains",
    "sidereal_scramble",
    "scramble_fap_periodogram",
    "activity_window",
    "census",
    "synthetic_repeater_set",
    "run",
]

CAT2_LOCAL = Path("data/chimefrbcat2.csv")
CAT2_DOI = "10.11570/25.0066"  # CANFAR; table also ships with ApJS 10.3847/1538-4365/ae3828
SIDEREAL_DAY = 0.9972695663  # days
SENTINEL = -9999.0


def load_catalog2(path: str | Path = CAT2_LOCAL) -> dict:
    """Load the Cat-2 burst table (CSV mirror of the CANFAR/ApJS release).

    One row per (event, sub-burst); this keeps the FIRST sub-burst row per event_id (the
    catalogue's TOA/DM/fluence for multi-component bursts differ per component; the first
    component's `mjd_400` is the event TOA convention used throughout). The -9999 sentinel
    becomes NaN.
    """
    import csv

    with open(path) as f:
        rows = list(csv.DictReader(f))
    seen: set[str] = set()
    keep = []
    for r in rows:
        if r["event_id"] not in seen:
            seen.add(r["event_id"])
            keep.append(r)

    def col(name: str, numeric: bool = True) -> np.ndarray:
        if not numeric:
            return np.array([r[name] for r in keep])
        v = np.array([float(r[name]) if r[name] not in ("", "None") else np.nan for r in keep])
        v[v == SENTINEL] = np.nan
        return v

    return {
        "tns_name": col("tns_name", numeric=False),
        "repeater_name": col("repeater_name", numeric=False),
        "mjd": col("mjd_400"),
        "mjd_err": col("mjd_400_err"),
        "dm": col("dm_fitb"),
        "dm_err": col("dm_fitb_err"),
        "fluence": col("fluence"),
        "width": col("width_fitb"),
        "scat_time": col("scat_time"),
        "exp_up_hr": col("exp_up"),
        "ra": col("ra"),
        "dec": col("dec"),
    }


def repeater_trains(cat: dict, *, min_bursts: int = 3) -> dict:
    """Per-repeater burst trains: sorted TOAs + DM/fluence/width arrays + exposure hours."""
    is_rep = np.array([n not in ("-9999", "", "None") for n in cat["repeater_name"]])
    out: dict[str, dict] = {}
    for name in np.unique(cat["repeater_name"][is_rep]):
        m = cat["repeater_name"] == name
        order = np.argsort(cat["mjd"][m])
        if m.sum() < min_bursts:
            continue
        out[str(name)] = {
            "mjd": cat["mjd"][m][order],
            "dm": cat["dm"][m][order],
            "dm_err": cat["dm_err"][m][order],
            "fluence": cat["fluence"][m][order],
            "width": cat["width"][m][order],
            "exp_up_hr": float(np.nanmedian(cat["exp_up_hr"][m])),
            "ra": float(np.nanmedian(cat["ra"][m])) if "ra" in cat else float("nan"),
            "dec": float(np.nanmedian(cat["dec"][m])) if "dec" in cat else float("nan"),
            "n_bursts": int(m.sum()),
        }
    return out


def sidereal_scramble(mjds: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Transit-comb-preserving scramble: keep sidereal phase, redraw the sidereal-day number.

    Each TOA is decomposed as mjd = n*SIDEREAL_DAY + phase; the integer n is redrawn uniformly
    over the source's observed span while the phase (the position within the daily transit
    window) is kept. The null therefore contains exactly the real data's daily-visibility comb
    --- the structure that makes exposure-blind FAPs non-rigorous on a transit survey --- but no
    multi-day clustering or periodicity.
    """
    t = np.asarray(mjds, float)
    n = np.floor(t / SIDEREAL_DAY)
    phase = t - n * SIDEREAL_DAY
    n_new = rng.integers(int(n.min()), int(n.max()) + 1, t.size)
    return np.sort(n_new * SIDEREAL_DAY + phase)


def _z2_grid(times: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Vectorised Rayleigh Z^2_1 over a period grid (same statistic as frbperiod.rayleigh_z2)."""
    t = np.asarray(times, float)
    phi = 2.0 * np.pi * ((t[None, :] / periods[:, None]) % 1.0)
    c = np.cos(phi).sum(axis=1)
    s = np.sin(phi).sum(axis=1)
    return (2.0 / t.size) * (c * c + s * s)


def default_period_grid(span_days: float) -> np.ndarray:
    """Trial periods from 1.5 d to span/3, uniform in frequency (2x oversampled)."""
    p_min, p_max = 1.5, max(span_days / 3.0, 4.5)
    n_f = max(int((1.0 / p_min - 1.0 / p_max) * span_days * 2), 16)
    freqs = np.linspace(1.0 / p_max, 1.0 / p_min, n_f)
    return 1.0 / freqs[::-1]


def scramble_fap_periodogram(
    mjds: np.ndarray,
    *,
    periods: np.ndarray | None = None,
    n_scramble: int = 200,
    seed: int = 0,
) -> dict:
    """Rayleigh periodogram whose peak FAP is calibrated by the sidereal scramble null.

    Returns the best period, its Z^2, and p = (k+1)/(n+1) where k counts scrambles whose
    maximum Z^2 over the SAME grid reaches the observed peak.
    """
    t = np.sort(np.asarray(mjds, float))
    span = float(t[-1] - t[0])
    if periods is None:
        periods = default_period_grid(span)
    z2 = _z2_grid(t, periods)
    k_best = int(np.argmax(z2))
    z2_obs = float(z2[k_best])
    rng = np.random.default_rng(seed)
    null_max = np.empty(n_scramble)
    for i in range(n_scramble):
        null_max[i] = _z2_grid(sidereal_scramble(t, rng), periods).max()
    p = float((null_max >= z2_obs).sum() + 1) / (n_scramble + 1)
    return {
        "periods": periods,
        "z2": z2,
        "best_period": float(periods[k_best]),
        "best_z2": z2_obs,
        "p_scramble": p,
        "null_z2max": null_max,
        "span_days": span,
    }


def activity_window(mjds: np.ndarray, period: float, *, containment: float = 0.9) -> dict:
    """Duty cycle at ``period``: the smallest phase arc containing ``containment`` of bursts.

    For FRB 20180916B the literature activity window is ~5 d of 16.35 d (duty ~0.31 for full
    containment). The arc is found by scanning circular windows over the sorted phases; the
    duty cycle is the arc length as a fraction of the period.
    """
    phases = np.sort((np.asarray(mjds, float) / period) % 1.0)
    n = phases.size
    k = max(int(np.ceil(containment * n)), 2)
    ext = np.concatenate([phases, phases + 1.0])
    arcs = ext[k - 1 : k - 1 + n] - ext[:n]
    i = int(np.argmin(arcs))
    return {
        "duty_cycle": float(arcs[i]),
        "phase_lo": float(ext[i] % 1.0),
        "containment": containment,
        "n_in_arc": k,
    }


def census(
    trains: dict,
    *,
    min_bursts_stats: int = 10,
    n_scramble: int = 200,
    fap_threshold: float = 0.01,
    seed: int = 0,
) -> list[dict]:
    """The uniform per-repeater census: Weibull k, scramble-calibrated periodogram, duty cycle.

    Sources below ``min_bursts_stats`` get rates only (their k posteriors are honestly
    unconstrained --- the stated completeness cut); duty cycles are quoted only where the
    periodogram peak beats ``fap_threshold`` against the sidereal-scramble null.
    """
    rows = []
    for j, (name, tr) in enumerate(sorted(trains.items())):
        t = tr["mjd"][np.isfinite(tr["mjd"])]
        row: dict = {
            "name": name,
            "n_bursts": int(t.size),
            "span_days": float(t.max() - t.min()) if t.size > 1 else 0.0,
            "exp_up_hr": tr["exp_up_hr"],
            "rate_per_hr": float(t.size / tr["exp_up_hr"]) if tr["exp_up_hr"] > 0 else np.nan,
        }
        if t.size >= min_bursts_stats and row["span_days"] > 30.0:
            fit = fit_weibull_waits(t, n_boot=200, seed=seed + j)
            pg = scramble_fap_periodogram(t, n_scramble=n_scramble, seed=seed + j)
            row.update(
                {
                    "weibull_k": fit.k,
                    "k_ci_low": fit.k_ci_low,
                    "k_ci_high": fit.k_ci_high,
                    "clustered": bool(fit.clustered),
                    "best_period": pg["best_period"],
                    # span/period: peaks with few cycles are activity-EPOCH degeneracies,
                    # not established periods -- the paper separates them by this number
                    "n_cycles": float(row["span_days"] / pg["best_period"]),
                    "best_z2": pg["best_z2"],
                    "p_scramble": pg["p_scramble"],
                }
            )
            if pg["p_scramble"] <= fap_threshold:
                row["duty_cycle"] = activity_window(t, pg["best_period"])["duty_cycle"]
        rows.append(row)
    return rows


def synthetic_repeater_set(
    *,
    k: float = 0.4,
    mean_wait: float = 8.0,
    period: float = 16.35,
    duty: float = 0.3,
    span: float = 1200.0,
    transit_window_min: float = 15.0,
    seed: int = 0,
) -> dict:
    """A synthetic repeater with KNOWN Weibull shape, period, and duty cycle, transit-sampled.

    Arrival times come from a Weibull renewal process (shape ``k``, mean wait ``mean_wait``
    days), are thinned to a periodic activity window (``period``, ``duty``), and are then
    snapped to the daily transit comb: each surviving burst lands at its day's transit, jittered
    within a ``transit_window_min``-minute window --- the same selection CHIME imposes. The
    recover-a-known must return k (clustered), the injected period, and a compatible duty cycle.
    """
    import math

    rng = np.random.default_rng(seed)
    scale = mean_wait / math.gamma(1.0 + 1.0 / k)  # Weibull mean = scale * Gamma(1 + 1/k)
    t, times = 0.0, []
    while t < span * 3:  # oversample; thinning + transit selection prune it
        t += scale * rng.weibull(k)
        times.append(t)
    times_arr = np.array(times)
    times_arr = times_arr[times_arr < span * 3]
    phase = (times_arr / period) % 1.0
    times_arr = times_arr[(phase < duty)]  # periodic activity window
    n_day = np.floor(times_arr / SIDEREAL_DAY)
    jitter = rng.uniform(0.0, transit_window_min / (24.0 * 60.0), times_arr.size)
    toas = np.unique(n_day * SIDEREAL_DAY + 0.3 + jitter)  # one transit per sidereal day
    toas = toas[toas < span]
    return {
        "mjd": toas + 59000.0,
        "true_k": k,
        "true_period": period,
        "true_duty": duty,
    }


def run(out: str = ".", *, offline: bool = True, n_scramble: int = 200) -> dict:
    """Offline: k/period/duty recover-a-known on transit-sampled synthetics; real: the census."""
    import json

    anchor_name = "FRB20180916B"
    # the duty-cycle gate cannot be finer than the scramble p-value resolution (k+1)/(n+1)
    fap_threshold = max(0.01, 2.0 / (n_scramble + 1))
    if offline:
        syn = synthetic_repeater_set(seed=0)
        trains = {"SYN-INJECTED": _train_from_mjd(syn["mjd"])}
        # a Poisson (k=1, no period) control must NOT come out clustered/periodic
        ctrl = synthetic_repeater_set(k=1.0, duty=1.0, seed=1)
        trains["SYN-POISSON"] = _train_from_mjd(ctrl["mjd"])
        rows = census(trains, n_scramble=n_scramble, fap_threshold=fap_threshold)
        source = "synthetic transit-sampled repeaters"
        extra = {
            "true_k": syn["true_k"],
            "true_period": syn["true_period"],
            "true_duty": syn["true_duty"],
        }
    else:  # pragma: no cover - needs the local Cat-2 mirror
        cat = load_catalog2()
        trains = repeater_trains(cat, min_bursts=2)  # every Cat-2 repeater gets a rate row
        rows = census(trains, n_scramble=n_scramble, fap_threshold=fap_threshold)
        source = f"CHIME/FRB Catalog 2 (CANFAR DOI {CAT2_DOI})"
        extra = {"n_repeaters_total": len(trains)}

    with_stats = [r for r in rows if "weibull_k" in r]
    clustered = [r for r in with_stats if r["clustered"]]
    significant = [r for r in with_stats if r["p_scramble"] <= fap_threshold]
    anchor = next((r for r in rows if r["name"] == anchor_name), None)
    metrics = {
        "source": source,
        "is_real": not offline,
        "n_sources": len(rows),
        "n_with_stats": len(with_stats),
        "n_clustered": len(clustered),
        "n_periodic_p01": len(significant),
        "median_k": round(float(np.median([r["weibull_k"] for r in with_stats])), 3)
        if with_stats
        else None,
        "periodic_names": sorted(r["name"] for r in significant),
        "anchor_period": round(anchor["best_period"], 2) if anchor else None,
        "anchor_p": anchor["p_scramble"] if anchor else None,
        "anchor_duty": round(anchor.get("duty_cycle", float("nan")), 3) if anchor else None,
        "rows": rows,
        **extra,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "frbwait_metrics.json").write_text(
        json.dumps(_json_safe(metrics), indent=2) + "\n"
    )
    _figure(rows, op / "papers" / "frbwait" / "figures")
    _write_macros(metrics, op / "papers" / "frbwait" / "generated" / "macros.tex")
    _write_census_table(rows, op / "papers" / "frbwait" / "generated" / "census_table.tex")
    return metrics


def _train_from_mjd(mjd: np.ndarray) -> dict:
    return {
        "mjd": np.asarray(mjd, float),
        "dm": np.full(mjd.size, 500.0),
        "dm_err": np.full(mjd.size, 0.5),
        "fluence": np.full(mjd.size, 5.0),
        "width": np.full(mjd.size, 1.0),
        "exp_up_hr": 100.0,
        "n_bursts": int(mjd.size),
    }


def _json_safe(x):  # noqa: ANN001, ANN202 - small recursive coercion helper
    if isinstance(x, dict):
        return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, np.ndarray):
        return None  # arrays (grids, null distributions) stay out of the JSON
    if isinstance(x, (np.floating, float)):
        xf = float(x)
        return xf if np.isfinite(xf) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def _figure(rows: list[dict], out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ks = [r["weibull_k"] for r in rows if "weibull_k" in r]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    if ks:
        ax1.hist(ks, bins=15, color="C0")
    ax1.axvline(1.0, color="C3", ls="--", label="Poisson (k=1)")
    ax1.set(xlabel="Weibull shape k", ylabel="repeaters", title="Clustering census")
    ax1.legend(fontsize=8)
    ps = [r["p_scramble"] for r in rows if "p_scramble" in r]
    ns = [r["n_bursts"] for r in rows if "p_scramble" in r]
    if ps:
        ax2.scatter(ns, ps, s=14)
    ax2.axhline(0.01, color="C3", ls="--")
    ax2.set(
        xscale="log",
        yscale="log",
        xlabel="bursts",
        ylabel="periodogram p (sidereal scramble)",
        title="Periodicity census",
    )
    fig.tight_layout()
    fig.savefig(out / "frbwait.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "fwReal" if m.get("is_real") else "fwSyn"
    lines = [
        "% Auto-generated by jansky_research.frbwait._write_macros -- do not edit.",
        "% Synthetic (fwSyn*) and real (fwReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under fwReal* (an offline rebuild resets fwReal* to placeholders by design).",
        rf"\newcommand{{\fwSource}}{{{m['source']}}}",
    ]
    keys = (
        ("N", "n_sources"),
        ("NStats", "n_with_stats"),
        ("NClustered", "n_clustered"),
        ("NPeriodic", "n_periodic_p01"),
        ("MedianK", "median_k"),
        ("AnchorPeriod", "anchor_period"),
        ("AnchorP", "anchor_p"),
        ("AnchorDuty", "anchor_duty"),
    )
    for ns in ("fwSyn", "fwReal"):
        live = ns == pref
        for macro, key in keys:
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _write_census_table(rows: list[dict], path: str | Path, *, top_n: int = 20) -> None:
    """Per-repeater census rows as LaTeX (top ``top_n`` by burst count; full set in the JSON)."""
    out = [
        "% Auto-generated by jansky_research.frbwait._write_census_table -- do not edit.",
        "% Columns: name & N & k (CI) & best period (d) & cycles & p_scramble & duty cycle",
    ]
    ranked = sorted((r for r in rows if "weibull_k" in r), key=lambda r: -r["n_bursts"])
    for r in ranked[:top_n]:
        # duty cycles are only meaningful at established (many-cycle) periods; low-cycle
        # peaks are epoch degeneracies and quoting a "duty cycle" there invites misreading
        duty = (
            f"{r['duty_cycle']:.2f}" if "duty_cycle" in r and r.get("n_cycles", 0) >= 10 else "--"
        )
        cyc = f"{r['n_cycles']:.0f}" if "n_cycles" in r else "--"
        out.append(
            f"{r['name']} & {r['n_bursts']} & "
            f"${r['weibull_k']:.2f}^{{+{r['k_ci_high'] - r['weibull_k']:.2f}}}"
            f"_{{-{r['weibull_k'] - r['k_ci_low']:.2f}}}$ & "
            f"{r['best_period']:.2f} & {cyc} & ${r['p_scramble']:.4g}$ & {duty} \\\\"
        )
    out.append(r"\hline")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(out) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Uniform CHIME Cat-2 repeater wait-time census.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--n-scramble", type=int, default=200)
    args = p.parse_args(argv)
    m = run(args.out, offline=args.offline, n_scramble=args.n_scramble)
    m["rows"] = f"[{len(m['rows'])} rows in results/frbwait_metrics.json]"
    print(json.dumps(_json_safe(m), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
