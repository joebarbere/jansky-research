"""The first long-period-transient (LPT) population catalogue and P--Pdot diagram (plan 35).

The confirmed LPT sample grew from 2 (2022) to 16 objects by mid-2026, and the class review (Rea,
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
        "pdot_type": ptype,
        "pdot_is_limit": ptype == "upper_limit",
        "pdot_is_measurement": np.isin(ptype, ("measurement", "measurement_disputed")),
        "is_wd_binary": np.array([r["binary_status"] in ("yes", "candidate") for r in rows]),
        # the raw status, because "unknown" is NOT "no": counting unknown-companion objects with
        # the no-companion class biased the split test toward its own null (round-10 referee)
        "binary_status": np.array([r["binary_status"] for r in rows]),
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

    At this sample size the test is genuinely EXACT: all C(n, n_a) partitions of the pooled
    log-periods are enumerated (11,440 at 16 choose 7) and the two-sided probability is the
    fraction with |median difference| >= observed. (An earlier version ran 20,000 random shuffles
    and called the result "exact" --- a Monte Carlo estimate whose third digit was seed noise.)
    Falls back to Monte Carlo only above ``_EXACT_LIMIT`` partitions.
    """
    from itertools import combinations
    from math import comb

    a = np.log10(period_s[is_wd])
    b = np.log10(period_s[~is_wd])
    if a.size < 2 or b.size < 2:
        return {"delta_log_median": float("nan"), "p_perm": float("nan"), "method": "n/a"}
    obs = np.median(a) - np.median(b)
    pool = np.concatenate([a, b])
    n = pool.size
    n_a = a.size
    if comb(n, n_a) <= _EXACT_LIMIT:
        count = 0
        total = 0
        idx = np.arange(n)
        for sel in combinations(idx, n_a):
            mask = np.zeros(n, bool)
            mask[list(sel)] = True
            d = np.median(pool[mask]) - np.median(pool[~mask])
            if abs(d) >= abs(obs) - 1e-12:
                count += 1
            total += 1
        return {
            "delta_log_median": round(float(obs), 3),
            "p_perm": round(count / total, 4),
            "method": f"exact ({total} partitions)",
        }
    rng = np.random.default_rng(0)
    shuffled = pool.copy()
    count = 0
    n_perm = 20000
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        d = np.median(shuffled[:n_a]) - np.median(shuffled[n_a:])
        if abs(d) >= abs(obs):
            count += 1
    return {
        "delta_log_median": round(float(obs), 3),
        "p_perm": round(count / n_perm, 4),
        "method": f"monte carlo ({n_perm})",
    }


_EXACT_LIMIT = 200_000  # enumerate exactly below this many partitions


def label_sensitivity(s: dict) -> list[dict]:
    """The split statistic under every defensible companion labelling.

    Three objects' companion status is UNKNOWN, and the published headline counted them with the
    no-companion class --- including two of the three longest periods in the sample, which biases
    the WD-vs-rest contrast toward zero. The primary comparison excludes unknowns; the variants
    bound the labelling systematic (the round-10 referee's table).
    """
    status = s["binary_status"]
    period = s["period_s"]
    wd = np.isin(status, ("yes", "candidate"))
    unknown = status == "unknown"
    out = []
    for label, mask_a, keep in (
        ("unknowns_excluded", wd, ~unknown),
        ("as_published_unknowns_with_rest", wd, np.ones(period.size, bool)),
        ("unknowns_as_binary", wd | unknown, np.ones(period.size, bool)),
        ("confirmed_only", status == "yes", np.isin(status, ("yes", "no"))),
    ):
        st = period_split_stat(period[keep], mask_a[keep])
        out.append(
            {
                "labelling": label,
                "n_wd": int(mask_a[keep].sum()),
                "n_rest": int((~mask_a[keep]).sum()),
                "delta_log_median": st["delta_log_median"],
                "p": st["p_perm"],
            }
        )
    return out


def death_line_margins(s: dict) -> dict:
    """How far below the death line each constrained object sits --- and why none could sit above.

    "9/9 below the death line" reads as nine independent confirmations; it is one structural fact.
    The death line scales as P^3, so at minute-to-hour periods the Pdot needed to sit ABOVE it
    implies fields orders of magnitude beyond any magnetar. This reports each object's margin
    (Pdot_needed / Pdot_constraint), the minimum, the count for the no-companion subset (the one
    that carries the theoretical claim), and a sweep of the death-valley constant showing the
    count is 9/9 across the literature range.
    """
    from typing import Any

    rows: list[dict[str, Any]] = []
    for i in range(s["period_s"].size):
        pd = s["pdot"][i]
        if not (np.isfinite(pd) and pd > 0):
            continue
        p = float(s["period_s"][i])
        needed = float(death_line(np.array([p]))[0])
        rows.append(
            {
                "name": str(s["name"][i]),
                "period_s": p,
                "pdot": float(pd),
                "pdot_needed": float(f"{needed:.2e}"),
                "margin": float(f"{needed / pd:.3g}"),
                "is_wd_binary": bool(s["is_wd_binary"][i]),
            }
        )
    margins = [float(r["margin"]) for r in rows]
    nonwd = [r for r in rows if not r["is_wd_binary"]]
    sweep = []
    for b_over_p2 in (5e10, 1e11, 2e11, 5e11, 1e12, 2e12):
        n_below = 0
        min_m = float("inf")
        for r in rows:
            needed = b_over_p2 / DEATH_B_OVER_P2 * float(r["pdot_needed"])
            pdv = float(r["pdot"])
            if pdv < needed:
                n_below += 1
            min_m = min(min_m, needed / pdv)
        sweep.append(
            {
                "b_over_p2": float(f"{b_over_p2:.0e}"),
                "n_below": n_below,
                "min_margin": float(f"{min_m:.3g}"),
            }
        )
    return {
        "objects": rows,
        "min_margin": min(margins) if margins else None,
        "max_margin": max(margins) if margins else None,
        "n_below_nonbinary": sum(1 for r in nonwd if float(r["pdot"]) < float(r["pdot_needed"])),
        "n_nonbinary_constrained": len(nonwd),
        "valley_sweep": sweep,
    }


def split_power(
    delta_log: float,
    n_a: int,
    n_b: int,
    *,
    width_dex: float = 1.75,
    n_sims: int = 400,
    alpha: float = 0.05,
    seed: int = 0,
) -> dict:
    """Power of the exact split test against OVERLAPPING classes offset by ``delta_log`` dex.

    The earlier "the test has power" demonstration injected DISJOINT classes ~1 dex apart ---
    roughly six times the observed offset --- so it could not fail. This measures the detection
    probability at the offset actually observed: two log-uniform classes each spanning the
    sample's own log-period width, shifted by ``delta_log``, scored by the same exact test.
    """
    from itertools import combinations

    rng = np.random.default_rng(seed)
    n = n_a + n_b
    all_idx = frozenset(range(n))
    sel_list = list(combinations(range(n), n_a))
    sel_idx = np.array(sel_list, dtype=int)
    rest_idx = np.array([sorted(all_idx - set(sel)) for sel in sel_list], dtype=int)
    hits = 0
    for _ in range(n_sims):
        a = rng.uniform(0.0, width_dex, n_a) + delta_log
        b = rng.uniform(0.0, width_dex, n_b)
        pool = np.concatenate([a, b])
        obs = np.median(a) - np.median(b)
        med_a = np.median(pool[sel_idx], axis=1)
        med_b = np.median(pool[rest_idx], axis=1)
        p = float(np.mean(np.abs(med_a - med_b) >= abs(obs) - 1e-12))
        if p < alpha:
            hits += 1
    return {
        "delta_log": delta_log,
        "n_a": n_a,
        "n_b": n_b,
        "width_dex": width_dex,
        "n_sims": n_sims,
        "power": round(hits / n_sims, 3),
    }


def synthetic_lpt_population(
    n: int = 16, *, n_wd: int = 7, split_min: float = 78.0, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Round-trip fixture at the REAL sample composition (16 objects, 7 binaries).

    NOTE the injected classes are disjoint at ``split_min`` --- an easy target that verifies the
    machinery, not a power statement at the observed offset; :func:`split_power` is that.
    """
    rng = np.random.default_rng(seed)
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
    """Build the population table, the split statistics, and the class P--Pdot diagram."""

    s = load_sample()
    pop = population_table(s)
    # primary: unknown-companion objects EXCLUDED (counting them with the no-companion class
    # biased the contrast toward zero); the full labelling table bounds the systematic
    labels = label_sensitivity(s)
    primary = labels[0]
    published = labels[1]
    # round-trip check at the real composition: the injected disjoint split must register
    p_syn, wd_syn = synthetic_lpt_population()
    split_syn = period_split_stat(p_syn, wd_syn)
    margins = death_line_margins(s)
    power = {
        "at_primary_offset": split_power(primary["delta_log_median"], 7, 6),
        "at_published_offset": split_power(published["delta_log_median"], 7, 9),
    }

    metrics = {
        "source": f"vendored cross-checked sample ({pop['n_lpt']} LPTs, per-value provenance)",
        **pop,
        "delta_log_median_period": primary["delta_log_median"],
        "p_split": primary["p"],
        "split_primary_n": [primary["n_wd"], primary["n_rest"]],
        "label_sensitivity": labels,
        "p_perm_synthetic_split": split_syn["p_perm"],
        "split_method": period_split_stat(s["period_s"], s["is_wd_binary"])["method"],
        "death_line_margins": margins,
        "split_power": power,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "lpt_metrics.json")
    _figure(s, op / "papers" / "lpt" / "figures")
    _write_macros(metrics, op / "papers" / "lpt" / "generated" / "macros.tex")
    _write_table(s, op / "papers" / "lpt" / "generated" / "table.tex")
    return metrics


ARXIV_TO_BIBKEY = {
    "astro-ph/0503052": "hyman2005",
    "2503.08033": "hurleywalker2022",
    "2503.08036": "hurleywalker2023",
    "2407.12266": "caleb2024",
    "2407.07480": "dong2025a",
    "2411.16606": "wang2025",
    "2408.11536": "deruiter2025",
    "2408.15757": "hurleywalker2024",
    "2501.09133": "lee2025",
    "2507.13453": "anumarlapudi2025",
    "2507.05139": "dong2025b",
    "2507.14448": "mcsweeney2025",
    "2606.04232": "rose2026",
    "2603.07857": "pritchard2026",
    "2606.20067": "wang2026vaster",
}


def _write_table(s: dict, path) -> None:
    """Emit the catalogue table the paper's contribution claim promises: all 16 objects, cited.

    Unknown companion status is shown as "unknown", never folded into "no" --- the labelling
    distinction the split statistic turns on.
    """
    from pathlib import Path as _P

    order = np.argsort(s["period_s"])
    lines = [
        "% Auto-generated by jansky_research.lpt._write_table -- do not edit by hand.",
        r"\begin{deluxetable*}{lrrlclc}",
        r"\tablecaption{The 16 confirmed long-period transients, per-value provenance."
        r" $\dot P$ types: m = measurement, l = upper limit, c = consistent with zero,"
        r" (d) = disputed. \label{tab:lpts}}",
        r"\tablehead{\colhead{Name} & \colhead{$P$ (min)} & \colhead{$\dot P$ (s\,s$^{-1}$)} &"
        r" \colhead{type} & \colhead{Companion} & \colhead{X-ray} & \colhead{Discovery}}",
        r"\startdata",
    ]
    abbr = {
        "upper_limit": "l",
        "measurement": "m",
        "measurement_disputed": "m(d)",
        "consistent_zero": "c",
    }
    for i in order:
        pd = s["pdot"][i]
        if np.isfinite(pd):
            mant, exp = f"{pd:.2e}".split("e")
            pd_s = rf"${float(mant):g}\times10^{{{int(exp)}}}$"
        else:
            pd_s = r"\nodata"
        ptype = abbr.get(str(s["pdot_type"][i]), r"\nodata")
        key = ARXIV_TO_BIBKEY.get(str(s["arxiv"][i]), "")
        cite = rf"\citet{{{key}}}" if key else r"\nodata"
        name = str(s["name"][i]).replace("_", r"\_")
        lines.append(
            f"{name} & {s['period_s'][i] / 60.0:.1f} & {pd_s} & {ptype} & "
            f"{s['binary_status'][i]} & {'yes' if s['xray'][i] else 'no'} & {cite} \\\\"
        )
    lines += [r"\enddata", r"\end{deluxetable*}"]
    p = _P(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
        title=f"The LPT class on the P–$\\dot P$ plane ({s['period_s'].size} objects)",
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
        rf"\newcommand{{\lptSplitP}}{{{_fmt('p_split')}}}",
        rf"\newcommand{{\lptSynSplitP}}{{{_fmt('p_perm_synthetic_split')}}}",
    ]
    labels = {row["labelling"]: row for row in m.get("label_sensitivity") or []}
    pub = labels.get("as_published_unknowns_with_rest", {})
    asbin = labels.get("unknowns_as_binary", {})
    conf = labels.get("confirmed_only", {})
    dm = m.get("death_line_margins") or {}
    sweep = dm.get("valley_sweep") or []
    sw_min = min((r["min_margin"] for r in sweep), default=None)
    pw = m.get("split_power") or {}

    def _d(dic: dict, key: str) -> str:
        v = dic.get(key)
        return "--" if v is None else str(v)

    lines += [
        rf"\newcommand{{\lptSplitDPub}}{{{_d(pub, 'delta_log_median')}}}",
        rf"\newcommand{{\lptSplitPPub}}{{{_d(pub, 'p')}}}",
        rf"\newcommand{{\lptSplitDAsBin}}{{{_d(asbin, 'delta_log_median')}}}",
        rf"\newcommand{{\lptSplitPAsBin}}{{{_d(asbin, 'p')}}}",
        rf"\newcommand{{\lptSplitDConf}}{{{_d(conf, 'delta_log_median')}}}",
        rf"\newcommand{{\lptSplitPConf}}{{{_d(conf, 'p')}}}",
        rf"\newcommand{{\lptMinMargin}}{{{_d(dm, 'min_margin')}}}",
        rf"\newcommand{{\lptMaxMargin}}{{{_d(dm, 'max_margin')}}}",
        rf"\newcommand{{\lptNdeathNonbin}}{{{_d(dm, 'n_below_nonbinary')}}}",
        rf"\newcommand{{\lptNnonbinConstr}}{{{_d(dm, 'n_nonbinary_constrained')}}}",
        rf"\newcommand{{\lptValleyMinMargin}}{{{'--' if sw_min is None else sw_min}}}",
        rf"\newcommand{{\lptPowerPrimary}}{{{_d(pw.get('at_primary_offset', {}), 'power')}}}",
        rf"\newcommand{{\lptPowerPub}}{{{_d(pw.get('at_published_offset', {}), 'power')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


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
