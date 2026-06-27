"""FRB repeater activity-periodicity search.

Some repeating fast radio bursts show **periodic activity windows** — bursts cluster at a preferred
phase of a multi-day cycle. The archetype is FRB 20180916B, with a 16.35-day period (CHIME/FRB
Collaboration 2020, Nature 582, 351). This module searches a repeater's burst arrival times (MJDs)
for such a period with a phase-folding **Rayleigh ($Z^2_1$) periodogram** — pure NumPy, CPU-only.

Honest scope (and the write-up must say this): the CHIME catalogue is a **transit survey** that
sees each source roughly once per sidereal day with strongly non-uniform exposure. A catalogue-only
periodogram therefore (a) is heavily aliased near 1 day and its beats, and (b) cannot give a
rigorous false-alarm probability without the survey's exposure model. We can *recover* a known
period as a peak and set rough limits, but significance from the simple analytic FAP below is an
*upper bound on confidence* only — it ignores exposure and aliasing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "PeriodogramResult",
    "false_alarm_prob",
    "period_search",
    "rayleigh_z2",
    "run",
    "search_repeaters",
    "synthetic_periodic_arrivals",
]


@dataclass(frozen=True)
class PeriodogramResult:
    """Result of a Rayleigh periodogram search over a grid of trial periods (days)."""

    periods: np.ndarray
    z2: np.ndarray
    best_period: float
    best_z2: float
    n_bursts: int
    fap: float  # analytic false-alarm probability of the peak (exposure-blind; an upper bound)


def rayleigh_z2(times: np.ndarray, period: float) -> float:
    """Rayleigh $Z^2_1$ statistic for phase concentration at a trial ``period``.

    Folds the arrival times to phase $\\phi_i = 2\\pi\\,(t_i/P \\bmod 1)$ and returns
    $Z^2_1 = (2/n)\\,[(\\sum\\cos\\phi_i)^2 + (\\sum\\sin\\phi_i)^2]$. It is large when the bursts
    cluster at one phase (a periodic activity window) and $\\sim 2$ for random phases. Under the
    no-signal hypothesis $Z^2_1$ follows a $\\chi^2$ with 2 degrees of freedom.
    """
    t = np.asarray(times, dtype=float)
    phi = 2.0 * np.pi * ((t / period) % 1.0)
    c = np.cos(phi).sum()
    s = np.sin(phi).sum()
    return float((2.0 / t.size) * (c * c + s * s))


def false_alarm_prob(z2_max: float, n_indep: int) -> float:
    """Exposure-blind false-alarm probability of a peak $Z^2_{\\max}$ over ``n_indep`` trials.

    For a single trial, $P(Z^2_1 > z) = e^{-z/2}$ ($\\chi^2_2$ survival). With ``n_indep``
    independent trial periods, $\\mathrm{FAP} = 1 - (1 - e^{-z/2})^{n_\\mathrm{indep}}$. This
    **ignores the survey exposure and daily aliasing**, so it is only an upper bound on confidence.
    """
    p_single = np.exp(-z2_max / 2.0)
    return float(1.0 - (1.0 - p_single) ** max(n_indep, 1))


def period_search(
    times: np.ndarray, periods: np.ndarray, *, span: float | None = None
) -> PeriodogramResult:
    """Rayleigh periodogram over ``periods`` (days); returns the peak period and its FAP.

    ``n_indep`` (for the FAP) is estimated as the number of independent frequencies across the data
    span, $\\sim (1/P_\\min - 1/P_\\max)\\,T$, capped by the grid size.
    """
    t = np.asarray(times, dtype=float)
    z2 = np.array([rayleigh_z2(t, P) for P in periods])
    k = int(np.argmax(z2))
    if span is None:
        span = float(t.max() - t.min())
    n_indep = min(int((1.0 / periods.min() - 1.0 / periods.max()) * span) + 1, periods.size)
    return PeriodogramResult(
        periods=periods,
        z2=z2,
        best_period=float(periods[k]),
        best_z2=float(z2[k]),
        n_bursts=int(t.size),
        fap=false_alarm_prob(float(z2[k]), n_indep),
    )


def synthetic_periodic_arrivals(
    period: float = 16.0,
    n: int = 40,
    active_frac: float = 0.3,
    span: float = 350.0,
    seed: int | None = 0,
) -> np.ndarray:
    """Generate burst MJDs clustered in one phase window of ``period`` (offline test fixture).

    Each burst lands in a random cycle within ``span`` days, at a phase uniformly inside
    ``[0, active_frac]`` of the period — i.e. a periodic activity window. Used by the tests and as
    the offline fallback.
    """
    rng = np.random.default_rng(seed)
    n_cycles = max(int(span / period), 1)
    cycle = rng.integers(0, n_cycles, n)
    phase = rng.uniform(0.0, active_frac, n)
    return np.sort((cycle + phase) * period)


def _synthetic_repeaters() -> tuple[np.ndarray, np.ndarray]:
    """Two synthetic repeaters — one 16.35-day periodic, one random — for the offline run/tests."""
    a = synthetic_periodic_arrivals(16.35, n=40, active_frac=0.15, span=400, seed=0)
    rng = np.random.default_rng(1)
    b = np.sort(rng.uniform(0.0, 400.0, 15))
    mjd = 58000.0 + np.concatenate([a, b])
    names = np.array(["SYN-PER"] * a.size + ["SYN-RND"] * b.size)
    return mjd, names


def search_repeaters(
    mjd: np.ndarray, names: np.ndarray, *, min_bursts: int = 8, periods: np.ndarray | None = None
) -> tuple[list[dict], np.ndarray]:
    """Run the periodogram per repeater source; sources below ``min_bursts`` are skipped.

    Returns ``(rows, periods)`` where each row is ``{name, n, searched, best_period, z2, fap}``.
    """
    if periods is None:
        periods = np.linspace(2.0, 100.0, 12000)
    mjd = np.asarray(mjd, dtype=float)
    names = np.asarray(names)
    rows: list[dict] = []
    for nm in sorted(set(names.tolist())):
        t = np.sort(mjd[names == nm])
        if t.size < min_bursts:
            rows.append(
                {
                    "name": nm,
                    "n": int(t.size),
                    "searched": False,
                    "best_period": None,
                    "z2": None,
                    "fap": None,
                }
            )
            continue
        r = period_search(t, periods)
        rows.append(
            {
                "name": nm,
                "n": int(t.size),
                "searched": True,
                "best_period": r.best_period,
                "z2": r.best_z2,
                "fap": r.fap,
            }
        )
    return rows, periods


def run(out: str = ".", *, offline: bool = False, min_bursts: int = 8) -> dict:
    """Search every catalogue repeater for activity periodicity; write results + a periodogram.

    Writes ``results/period_metrics.json``, ``survey/period_results.csv`` (per source), and a
    periodogram figure for the most significant detection. Returns the metrics dict.
    """
    import csv
    import json
    from pathlib import Path

    if offline:
        mjd, names = _synthetic_repeaters()
        source = "synthetic"
    else:  # pragma: no cover - network
        from . import pipeline

        cat, source = pipeline.build_catalog(offline=False)
        rep = np.asarray(cat["repeater"], dtype=bool)
        mjd = np.asarray(cat["mjd"])[rep]
        names = np.asarray(cat["repeater_name"])[rep]

    rows, periods = search_repeaters(mjd, names, min_bursts=min_bursts)
    searched = [r for r in rows if r["searched"]]
    significant = sorted((r for r in searched if r["fap"] < 0.01), key=lambda r: r["fap"])
    metrics = {
        "source": source,
        "n_sources": len(rows),
        "n_searchable": len(searched),
        "n_significant": len(significant),
        "min_bursts": min_bursts,
        "detections": [
            {
                "name": r["name"],
                "period_days": round(r["best_period"], 3),
                "z2": round(r["z2"], 1),
                "fap": r["fap"],
                "n": r["n"],
            }
            for r in significant
        ],
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "period_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    (op / "survey").mkdir(parents=True, exist_ok=True)
    with open(op / "survey" / "period_results.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["repeater", "n_bursts", "searched", "best_period_days", "z2", "fap"])
        for r in sorted(rows, key=lambda r: -r["n"]):
            w.writerow(
                [
                    r["name"],
                    r["n"],
                    r["searched"],
                    "" if r["best_period"] is None else f"{r['best_period']:.3f}",
                    "" if r["z2"] is None else f"{r['z2']:.1f}",
                    "" if r["fap"] is None else f"{r['fap']:.2e}",
                ]
            )
    if significant:
        _periodogram_figure(mjd, names, significant[0], periods, op / "paper" / "figures")
    return metrics


def _periodogram_figure(mjd, names, det, periods, out_dir):
    """Plot the Rayleigh periodogram of the most significant detection, marking its peak."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t = np.sort(np.asarray(mjd)[np.asarray(names) == det["name"]])
    z2 = np.array([rayleigh_z2(t, P) for P in periods])
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(periods, z2, lw=0.8)
    ax.axvline(det["best_period"], color="r", ls="--", label=f"peak {det['best_period']:.2f} d")
    ax.set(
        xlabel="trial period (days)",
        ylabel=r"Rayleigh $Z^2_1$",
        title=f"{det['name']} activity periodogram ({det['n']} bursts)",
    )
    ax.legend()
    p = out / "periodogram.pdf"
    fig.tight_layout()
    fig.savefig(p)
    plt.close(fig)
    return p
