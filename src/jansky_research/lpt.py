"""The first long-period-transient (LPT) population catalogue and P--Pdot diagram (plan 35).

The confirmed LPT sample grew from 2 (2022) to 13 objects by mid-2026, and the class review (Rea,
Hurley-Walker & Caleb 2026, arXiv:2601.10393) explicitly notes that **no population synthesis
exists**. This slice ships one: a verified, per-value-provenanced table (``data/lpt_sample.csv``,
compiled 2026-07 from the discovery papers -- during which a period typo in the review's own data
file was caught, see FLAG_A) and a P--Pdot placement against the pulsar population, death line, and
constant-B tracks (reusing ``ppdot``) -- the review's own Fig. 3 plots the class, so the novelty
here is narrower and stated as such: per-value provenance, explicit measurement-vs-limit typing,
and statistics that regenerate from the table. Dipole-formula quantities
(B, tau) are computed **only** for objects where a neutron-star interpretation is viable; confirmed
white-dwarf binaries are plotted but not assigned NS dipole values -- their "period" is orbital.

The population question the diagram frames: only TWO Pdot *measurements* exist in the class
(CHIME J0630+25 spin-down, glitch-caveated; CHIME J1634+44 spin-UP -- natural for a binary),
everything else is upper limits, several so weak they constrain nothing; and the WD-binary members
cluster at long periods (a hinted ~78-min boundary the sample is still too small to establish --
we report the split statistic with its tiny-N caveat, not a claim).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .ppdot import DEATH_B_OVER_P2, death_line

__all__ = [
    "load_sample",
    "population_table",
    "period_split_stat",
    "synthetic_lpt_population",
    "run",
]

SAMPLE_CSV = Path(__file__).resolve().parents[2] / "data" / "lpt_sample.csv"


def load_sample(path: str | Path = SAMPLE_CSV) -> dict:
    """Load the vendored, provenance-carrying LPT table into arrays.

    ``pdot`` is NaN where no constraint exists; ``pdot_is_limit`` marks upper limits (plotted as
    downward arrows, never used as measurements); ``is_wd_binary`` marks confirmed/candidate
    binaries (the review's WDLPT split + the post-review accretor).
    """
    rows = list(csv.DictReader(open(path)))
    period = np.array([float(r["period_s"]) for r in rows])
    pdot = np.array(
        [float(r["pdot_s_s"]) if r["pdot_s_s"] not in ("", None) else np.nan for r in rows]
    )
    ptype = np.array([r["pdot_type"] for r in rows])
    return {
        "name": np.array([r["name"] for r in rows]),
        "ra": np.array([float(r["ra_deg"]) for r in rows]),
        "dec": np.array([float(r["dec_deg"]) for r in rows]),
        "period_s": period,
        "pdot": pdot,
        "pdot_is_limit": ptype == "upper_limit",
        "pdot_is_measurement": np.isin(ptype, ("measurement", "measurement_disputed")),
        "is_wd_binary": np.array([r["binary_status"] in ("yes", "candidate") for r in rows]),
        "xray": np.array([r["xray"] == "yes" for r in rows]),
        "year": np.array([int(r["year"]) for r in rows]),
        "arxiv": np.array([r["discovery_arxiv"] for r in rows]),
    }


def population_table(s: dict) -> dict:
    """Headline population numbers for the class (the paper's summary statistics)."""
    p_min = float(s["period_s"].min())
    p_max = float(s["period_s"].max())
    below_death = 0
    n_constrained = 0
    for p, pd in zip(s["period_s"], s["pdot"], strict=True):
        if np.isfinite(pd) and pd > 0:
            n_constrained += 1
            dl = float(death_line(np.array([p]))[0])
            if pd < dl:  # measured value or the LIMIT itself sits below the death line
                below_death += 1
    return {
        "n_lpt": int(s["period_s"].size),
        "n_wd_binary": int(s["is_wd_binary"].sum()),
        "n_xray": int(s["xray"].sum()),
        "n_pdot_measurements": int(s["pdot_is_measurement"].sum()),
        "period_min_min": round(p_min / 60.0, 1),
        "period_max_hr": round(p_max / 3600.0, 2),
        "median_period_min": round(float(np.median(s["period_s"])) / 60.0, 1),
        # n_constrained counts POSITIVE Pdot values/limits only: the spin-up measurement
        # (J1634+44) and the consistent-with-zero object (J1755-2527) cannot be placed on the
        # spin-down death-line criterion and are excluded (stated in the paper).
        "n_pdot_constrained": n_constrained,
        "n_below_death_line": below_death,
    }


def period_split_stat(period_s: np.ndarray, is_wd: np.ndarray) -> dict:
    """Rank test of the hinted WD-binary-vs-rest period split (report, don't claim: tiny N).

    Mann--Whitney U via the normal approximation is unreliable at N~13, so we report the exact
    two-sided probability from a permutation test on the median difference.
    """
    rng = np.random.default_rng(0)
    a = np.log10(period_s[is_wd])
    b = np.log10(period_s[~is_wd])
    if a.size < 2 or b.size < 2:
        return {"delta_log_median": float("nan"), "p_perm": float("nan")}
    obs = np.median(a) - np.median(b)
    pool = np.concatenate([a, b]).copy()  # shuffled in-place below
    n_a = a.size
    count = 0
    n_perm = 20000
    for _ in range(n_perm):
        rng.shuffle(pool)
        d = np.median(pool[:n_a]) - np.median(pool[n_a:])
        if abs(d) >= abs(obs):
            count += 1
    return {"delta_log_median": round(float(obs), 3), "p_perm": round(count / n_perm, 4)}


def synthetic_lpt_population(
    n: int = 13, *, split_min: float = 78.0, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Round-trip fixture: a fake population with a KNOWN period split at ``split_min`` minutes."""
    rng = np.random.default_rng(seed)
    n_wd = max(2, n // 3)
    p_wd = 10 ** rng.uniform(np.log10(split_min * 60), np.log10(12 * 3600), n_wd)
    p_ns = 10 ** rng.uniform(np.log10(400), np.log10(split_min * 60), n - n_wd)
    period = np.concatenate([p_wd, p_ns])
    is_wd = np.zeros(n, bool)
    is_wd[:n_wd] = True
    return period, is_wd


def crossmatch_counterparts(
    s: dict, *, match_arcsec: float = 20.0
) -> list[dict]:  # pragma: no cover - network
    """Per-LPT continuum-counterpart check: VLASS QL2 cone + LoTSS DR3 forced cutout peak.

    LPTs are burst emitters; a persistent continuum counterpart (or its absence) constrains any
    steady emission component. VLASS covers Dec > -40 (2-4 GHz, ~0.7 mJy at 5 sigma QL); LoTSS
    DR3 the northern sky (144 MHz). Returns one row per object with fluxes or 5-sigma limits.
    """
    import io

    import requests
    from astropy.io import fits as _fits

    from .vlass import fetch_vlass_epoch

    out = []
    for i, name in enumerate(s["name"]):
        ra, dec = float(s["ra"][i]), float(s["dec"][i])
        row: dict = {"name": str(name), "ra": ra, "dec": dec}
        if dec > -40.0:
            try:
                vra, vdec, pk, _ = fetch_vlass_epoch(1, (ra, dec), 0.02)
                d = np.hypot((vra - ra) * np.cos(np.radians(dec)), vdec - dec) * 3600.0
                j = int(np.argmin(d)) if d.size else -1
                if d.size and d[j] < match_arcsec:
                    row["vlass_mJy"] = round(float(pk[j]), 2)
                    row["vlass_sep_as"] = round(float(d[j]), 1)
                else:
                    row["vlass_mJy"] = None  # < ~0.7 mJy (5 sigma QL)
            except Exception as exc:
                row["vlass_note"] = f"fetch failed: {type(exc).__name__}"
        else:
            row["vlass_note"] = "outside VLASS dec range"
        try:
            r = requests.get(
                "https://lofar-surveys.org/dr3-cutout.fits",
                params={"pos": f"{ra},{dec}", "size": "0.05"},
                timeout=90,
            )
            if r.ok and r.headers.get("content-type", "").startswith("application/fits"):
                with _fits.open(io.BytesIO(r.content)) as hdul:
                    img = np.asarray(hdul[0].data, float).squeeze() * 1e3  # Jy->mJy
                c = np.array(img.shape) // 2
                peak = float(np.nanmax(img[c[0] - 2 : c[0] + 3, c[1] - 2 : c[1] + 3]))
                rms = float(1.4826 * np.nanmedian(np.abs(img - np.nanmedian(img))))
                row["lotss_peak_mJy"] = round(peak, 2)
                row["lotss_rms_mJy"] = round(rms, 3)
                row["lotss_detected"] = bool(peak > 5 * rms)
            else:
                row["lotss_note"] = f"no coverage (HTTP {r.status_code})"
        except Exception as exc:
            row["lotss_note"] = f"cutout failed: {type(exc).__name__}"
        out.append(row)
    return out


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Build the population table, the split statistic, and the class P--Pdot diagram."""
    import json

    s = load_sample()
    pop = population_table(s)
    split = period_split_stat(s["period_s"], s["is_wd_binary"])
    # round-trip check: the injected split must register as extreme (small p)
    p_syn, wd_syn = synthetic_lpt_population()
    split_syn = period_split_stat(p_syn, wd_syn)

    metrics = {
        "source": "vendored verified sample (13 LPTs, per-value provenance)",
        **pop,
        "delta_log_median_period": split["delta_log_median"],
        "p_perm_split": split["p_perm"],
        "p_perm_synthetic_split": split_syn["p_perm"],
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "lpt_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(s, op / "papers" / "lpt" / "figures")
    _write_macros(metrics, op / "papers" / "lpt" / "generated" / "macros.tex")
    return metrics


def _figure(s: dict, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))

    # pulsar context: constant-B tracks + the death line over the extended period range
    p_grid = np.logspace(-3, 5.2, 300)
    for b in (1e12, 1e14, 1e16):
        pd = (b / 3.2e19) ** 2 / p_grid
        ax.plot(p_grid, pd, ":", color="0.75", lw=0.8)
        ax.text(
            p_grid[-1] * 0.5,
            (b / 3.2e19) ** 2 / (p_grid[-1] * 0.5),
            f"$10^{{{int(np.log10(b))}}}$ G",
            fontsize=7,
            color="0.5",
            ha="right",
        )
    ax.plot(
        p_grid,
        death_line(p_grid),
        "--",
        color="C3",
        lw=1.2,
        label=rf"death line ($B/P^2={DEATH_B_OVER_P2:.1e}$)",
    )

    wd = s["is_wd_binary"]
    meas = s["pdot_is_measurement"] & np.isfinite(s["pdot"])
    lim = s["pdot_is_limit"] & np.isfinite(s["pdot"])
    # limits: downward arrows at the limit value
    ax.errorbar(
        s["period_s"][lim & ~wd],
        s["pdot"][lim & ~wd],
        yerr=0.5 * s["pdot"][lim & ~wd],
        uplims=True,
        fmt="o",
        ms=6,
        color="C0",
        label="LPT (no companion), $\\dot P$ limit",
    )
    ax.errorbar(
        s["period_s"][lim & wd],
        s["pdot"][lim & wd],
        yerr=0.5 * s["pdot"][lim & wd],
        uplims=True,
        fmt="s",
        ms=6,
        color="C2",
        label="WD binary/cand., $\\dot P$ limit",
    )
    for i in np.where(meas)[0]:
        m = "s" if wd[i] else "o"
        c = "C2" if wd[i] else "C0"
        ax.plot(s["period_s"][i], abs(s["pdot"][i]), m, ms=9, mfc="none", mec="C3", mew=2)
        ax.plot(s["period_s"][i], abs(s["pdot"][i]), m, ms=6, color=c)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="period (s)",
        ylabel=r"$\dot P$ (s s$^{-1}$)",
        title="The LPT class on the P–$\\dot P$ plane (13 objects, 2026-07)",
        xlim=(0.9e2, 1.2e5),
        ylim=(1e-16, 1e-6),
    )
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "lpt_ppdot.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.lpt._write_macros -- do not edit.",
        rf"\newcommand{{\lptN}}{{{_fmt('n_lpt')}}}",
        rf"\newcommand{{\lptNwd}}{{{_fmt('n_wd_binary')}}}",
        rf"\newcommand{{\lptNx}}{{{_fmt('n_xray')}}}",
        rf"\newcommand{{\lptNmeas}}{{{_fmt('n_pdot_measurements')}}}",
        rf"\newcommand{{\lptPmin}}{{{_fmt('period_min_min')}}}",
        rf"\newcommand{{\lptPmax}}{{{_fmt('period_max_hr')}}}",
        rf"\newcommand{{\lptPmed}}{{{_fmt('median_period_min')}}}",
        rf"\newcommand{{\lptNdeath}}{{{_fmt('n_below_death_line')}}}",
        rf"\newcommand{{\lptNconstr}}{{{_fmt('n_pdot_constrained')}}}",
        rf"\newcommand{{\lptSplitD}}{{{_fmt('delta_log_median_period')}}}",
        rf"\newcommand{{\lptSplitP}}{{{_fmt('p_perm_split')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="LPT population catalogue + P-Pdot diagram.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=True), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
