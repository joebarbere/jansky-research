"""Inner-heliosphere type III bursts: tracking the electron beam beyond the corona (Wind/WAVES).

A type III electron beam does not stop at the top of the corona --- it streams out along the open
field, exciting radio emission at the falling plasma frequency as it goes. Space-based receivers see
this as a slow drift from a few MHz down to tens of kHz over tens of minutes; the high-frequency
(RAD2, ~1--14 MHz) part traces the beam through the **inner heliosphere** (a few to ~10 R_sun, near
the Alfven surface), and RAD1 (down to ~20 kHz) would follow it on toward 1 AU. This slice fits the
drift in a Wind/WAVES dynamic spectrum and inverts it, via a **heliospheric** density model (Leblanc,
Dulk & Bougeret 1998), to the beam's outward radial speed and the heliocentric distance it reaches ---
the wider-distance companion to the coronal ``solarbursts`` slice (which used the Newkirk corona, valid
only to a few solar radii).

Reuses ``solarbursts``' dynamic-spectrum tools (background subtraction, burst windowing, the
per-channel ridge detector, the robust fit) and ``jansky.solar.density_from_plasma_frequency``; adds
the Leblanc model and the Wind/WAVES Level-2 CDF fetch. Pure NumPy with a synthetic offline fixture;
the real fetch needs the ``windwaves`` extra (``cdflib``) and is network-gated.
"""

from __future__ import annotations

import numpy as np

from . import solarbursts

__all__ = [
    "C_KMS",
    "R_AU_RSUN",
    "SPEED_GRID",
    "beam_speed",
    "speed_grid",
    "emission_radius",
    "fetch_windwaves",
    "leblanc_density",
    "leblanc_radius",
    "run",
    "synthetic_ip_burst",
]

C_KMS = 299792.458
R_AU_RSUN = 215.0  # 1 AU in solar radii


def leblanc_density(r_rsun: np.ndarray) -> np.ndarray:
    """Heliospheric electron density (cm⁻³) at heliocentric radius ``r_rsun`` (Leblanc et al. 1998).

    :math:`n_e(r) = 3.3\\times10^{5} r^{-2} + 4.1\\times10^{6} r^{-4} + 8.0\\times10^{7} r^{-6}`
    (r in solar radii), normalised to ~7.2 cm⁻³ at 1 AU. The :math:`r^{-2}` term dominates far from
    the Sun (a constant-speed solar wind); the steeper terms matter near the corona.
    """
    r = np.asarray(r_rsun, float)
    return 3.3e5 * r**-2.0 + 4.1e6 * r**-4.0 + 8.0e7 * r**-6.0


def leblanc_radius(n_e_cm3: np.ndarray) -> np.ndarray:
    """Invert :func:`leblanc_density`: heliocentric radius (R⊙) for a density, numerically.

    The model is monotonic in :math:`r`, so we interpolate on a fine grid from 1.3 to 250 R⊙.
    Densities outside the grid clamp to its ends.
    """
    rg = np.logspace(np.log10(1.3), np.log10(250.0), 4000)
    ng = leblanc_density(rg)  # strictly decreasing with r
    # np.interp needs increasing xp: reverse so log10(density) increases
    return np.interp(np.log10(np.asarray(n_e_cm3, float)), np.log10(ng[::-1]), rg[::-1])


def emission_radius(
    freq_mhz: np.ndarray, *, harmonic: int = 2, density_scale: float = 1.0
) -> np.ndarray:
    """Heliocentric radius (R⊙) of emission at ``freq_mhz``, via plasma frequency and the Leblanc model.

    The observed frequency is the ``harmonic`` of the local plasma frequency, so
    :math:`f_p = f/\\mathrm{harmonic}` gives the density
    (``jansky.solar.density_from_plasma_frequency``) and hence the radius (:func:`leblanc_radius`).
    ``density_scale`` multiplies the Leblanc profile: an event driving a factor-$k$ enhancement
    places the same plasma frequency at the radius where $k\\,n_{\\rm Leblanc}(r)$ matches, i.e.
    ``leblanc_radius(n_e / k)``. Note the exact degeneracy this exposes: ``harmonic=1`` with
    ``density_scale=4`` equals ``harmonic=2`` with ``density_scale=1``, because $f_p^2 \\propto n$
    — the emission-mode and density-enhancement systematics are one axis, not two.
    """
    from jansky import solar

    fp = np.asarray(freq_mhz, float) / harmonic
    return leblanc_radius(solar.density_from_plasma_frequency(fp) / density_scale)


def beam_speed(
    ridge_freqs_mhz: np.ndarray,
    ridge_times_s: np.ndarray,
    *,
    harmonic: int = 2,
    density_scale: float = 1.0,
) -> dict:
    """Outward beam speed from the drift ridge, via the Leblanc heliospheric density model.

    Maps each ridge frequency to a heliocentric radius (:func:`emission_radius`) and fits ONE
    point per distinct time sample (the column median) — the independent unit under the coarse
    cadence — with a leave-one-column-out jackknife error. A per-point fit is
    leverage-dominated (hundreds of channels stack into a few samples), and iterating a
    residual clip on points collapses onto the crowded near-Sun columns; both failure modes
    were referee-caught on the committed evidence. The estimator bracket (inverse regression
    on columns; the naive all-points OLS) is returned alongside.
    """
    f = np.asarray(ridge_freqs_mhz, float)
    t = np.asarray(ridge_times_s, float)
    r = emission_radius(f, harmonic=harmonic, density_scale=density_scale)
    nan = float("nan")
    cols = np.unique(t)
    if r.size < 3 or cols.size < 3:
        return {
            "speed_kms": nan,
            "speed_c": nan,
            "speed_c_se": nan,
            "speed_c_inverse": nan,
            "speed_c_points": nan,
            "r_lo": nan,
            "r_hi": nan,
            "r_hi_au": nan,
            "r2": nan,
            "n_used": 0,
            "n_time_cols": int(cols.size),
        }
    # One point per DISTINCT TIME SAMPLE (column median), no cross-column clipping. The
    # column is the independent unit: the coarse cadence stacks hundreds of channels into a
    # few samples, so a per-point fit is leverage-dominated, and iterating a residual clip on
    # points collapses onto the crowded near-Sun columns and rejects the sparse low-frequency
    # columns that carry the drift (measured live on the STEREO ridge: the converged per-point
    # clip returned 0.031 c on 2 surviving columns against 0.150 c from the column fit).
    col_r = np.asarray([float(np.median(r[t == c])) for c in cols])
    slope, icpt = np.polyfit(cols, col_r, 1)
    speed_kms = abs(slope) * 695700.0  # R_sun -> km
    model = slope * cols + icpt
    ss_res = float(np.sum((col_r - model) ** 2))
    ss_tot = float(np.sum((col_r - col_r.mean()) ** 2))
    # leave-one-column-out jackknife
    jk = []
    if cols.size >= 4:
        for i in range(cols.size):
            m = np.arange(cols.size) != i
            s, _ = np.polyfit(cols[m], col_r[m], 1)
            jk.append(abs(s))
        k = len(jk)
        jarr = np.asarray(jk)
        se = float(np.sqrt((k - 1) / k * np.sum((jarr - jarr.mean()) ** 2)))
    else:
        se = nan
    # estimator bracket: inverse regression on columns, and the naive all-points OLS
    inv = nan
    if np.ptp(col_r) > 0:
        s_tr, _ = np.polyfit(col_r, cols, 1)
        inv = abs(1.0 / s_tr) if s_tr != 0 else nan
    s_pts, _ = np.polyfit(t, r, 1)
    to_c = 695700.0 / C_KMS
    return {
        "speed_kms": float(speed_kms),
        "speed_c": float(speed_kms / C_KMS),
        "fit_slope": float(slope),
        "fit_intercept": float(icpt),
        "speed_c_se": float(se * to_c) if np.isfinite(se) else nan,
        "speed_c_inverse": float(inv * to_c) if np.isfinite(inv) else nan,
        "speed_c_points": float(abs(s_pts) * to_c),
        "r_lo": float(np.min(r)),
        "r_hi": float(np.max(r)),
        "r_hi_au": float(np.max(r) / R_AU_RSUN),
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else nan,
        "n_used": int(r.size),
        "n_time_cols": int(cols.size),
    }


#: The (harmonic x density-scale) systematics grid: emission mode is a factor of two in plasma
#: frequency and an event-driven density enhancement is "a factor of several"; f_p^2 ~ n makes
#: the two exactly degenerate (harmonic 1 at scale 4 == harmonic 2 at scale 1), so the grid's
#: distinct rows are what the data can actually distinguish: nothing. The bracket is the result.
SPEED_GRID = ((1, 1.0), (2, 1.0), (2, 2.0), (2, 4.0))


def speed_grid(ridge_freqs_mhz: np.ndarray, ridge_times_s: np.ndarray) -> list[dict]:
    """Beam speed over the (harmonic, density-scale) grid, from one fitted ridge."""
    out = []
    for harmonic, scale in SPEED_GRID:
        spd = beam_speed(ridge_freqs_mhz, ridge_times_s, harmonic=harmonic, density_scale=scale)
        out.append(
            {
                "harmonic": harmonic,
                "density_scale": scale,
                "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
                "r_lo_rsun": round(spd["r_lo"], 2) if np.isfinite(spd["r_lo"]) else None,
                "r_hi_rsun": round(spd["r_hi"], 2) if np.isfinite(spd["r_hi"]) else None,
                "r_hi_au": round(spd["r_hi_au"], 3) if np.isfinite(spd["r_hi_au"]) else None,
                "r2": round(spd["r2"], 3) if np.isfinite(spd["r2"]) else None,
                "n_used": spd["n_used"],
            }
        )
    return out


def synthetic_ip_burst(
    *,
    speed_c: float = 0.15,
    r0_rsun: float = 2.0,
    harmonic: int = 2,
    f_lo_mhz: float = 0.2,
    f_hi_mhz: float = 14.0,
    n_freq: int = 256,
    duration_s: float = 1800.0,
    n_time: int = 600,
    width_dex: float = 0.04,
    amp: float = 12.0,
    noise: float = 1.0,
    seed: int = 0,
) -> dict:
    """Synthetic interplanetary type III with an injected beam speed, via the Leblanc forward model.

    A beam climbs from ``r0_rsun`` at ``speed_c`` × c; at each instant the Leblanc density sets the
    plasma frequency and hence the (harmonic) emission frequency, tracing a slow high-to-low drift over
    a logarithmic frequency grid (Wind/WAVES spans ~0.02–14 MHz). Built from the same Leblanc mapping
    the analysis inverts, so a clean burst round-trips. Returns ``data`` (n_freq × n_time), ``freqs``
    (MHz, descending), ``times`` (s), and the injected ``truth_speed_c``.
    """
    from jansky import solar

    rng = np.random.default_rng(seed)
    freqs = np.logspace(np.log10(f_hi_mhz), np.log10(f_lo_mhz), n_freq)  # descending
    times = np.linspace(0.0, duration_s, n_time)
    v_rsun_per_s = speed_c * C_KMS / 695700.0
    r_t = r0_rsun + v_rsun_per_s * times
    f_ridge = harmonic * solar.plasma_frequency(leblanc_density(r_t))  # MHz, decreasing
    logf = np.log10(freqs)
    data = rng.normal(0.0, noise, (n_freq, n_time))
    for j, fr in enumerate(f_ridge):
        if fr <= 0:
            continue
        data[:, j] += amp * np.exp(-0.5 * ((logf - np.log10(fr)) / width_dex) ** 2)
    return {"data": data, "freqs": freqs, "times": times, "truth_speed_c": speed_c}


def fetch_windwaves(
    date_yyyymmdd: str, *, receiver: str = "rad2"
) -> dict:  # pragma: no cover - network
    """Fetch a Wind/WAVES Level-2 radio dynamic spectrum (RAD1 or RAD2) from SPDF for one day.

    Downloads ``wi_l2_wav_{receiver}_{date}_v01.cdf`` from the public SPDF archive and returns
    ``data`` (n_freq × n_time, PSD), ``freqs`` (MHz, descending), and ``times`` (seconds from the file
    start). Needs the ``windwaves`` extra (``cdflib``).
    """
    import re

    import cdflib
    import requests

    yyyy = date_yyyymmdd[:4]
    base = f"https://spdf.gsfc.nasa.gov/pub/data/wind/waves/{receiver}_l2/{yyyy}/"
    idx = requests.get(base, timeout=60).text
    pat = rf"wi_l2_wav_{receiver}_{date_yyyymmdd}_v[0-9]+\.cdf"
    m = re.findall(pat, idx)
    if not m:
        raise RuntimeError(f"no Wind/WAVES {receiver} file for {date_yyyymmdd}")
    raw = requests.get(base + m[0], timeout=120).content
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".cdf") as fh:
        fh.write(raw)
        fh.flush()
        c = cdflib.CDF(fh.name)
        freq = np.asarray(c.varget("FREQUENCY"), float) / 1e6  # Hz -> MHz
        psd = np.asarray(c.varget("PSD_V2_SP"), float)  # (time, freq)
        ep = cdflib.cdfepoch.to_datetime(c.varget("Epoch"))
    times = (ep - ep[0]) / np.timedelta64(1, "s")
    # return as (freq, time) with freq descending, to match the synthetic / solarbursts convention
    order = np.argsort(freq)[::-1]
    return {
        "data": psd.T[order],
        "freqs": freq[order],
        "times": np.asarray(times, float),
        "t0_utc": str(ep[0]),
    }


def run(
    out: str = ".",
    *,
    offline: bool = True,
    date: str | None = None,
    receiver: str = "rad2",
    harmonic: int = 2,
    pad_s: float = 1200.0,
) -> dict:
    """Full slice: fit an interplanetary type III drift and report the beam speed and reach."""
    from pathlib import Path

    if offline or date is None:
        burst = synthetic_ip_burst(harmonic=harmonic)
        source = "synthetic"
        truth: float | None = burst["truth_speed_c"]
    else:  # pragma: no cover - network
        burst = fetch_windwaves(date, receiver=receiver)
        source = f"Wind/WAVES {receiver.upper()} {date}"
        truth = None

    window = solarbursts.find_burst_window(burst["data"], burst["times"], pad_s=pad_s)
    rf, rt = solarbursts.detect_burst_ridge(
        burst["data"], burst["freqs"], burst["times"], window=window
    )
    spd = beam_speed(rf, rt, harmonic=harmonic)
    metrics: dict = _speed_metrics(source, rf, rt, spd, harmonic, pad_s)
    # burst peak epoch, so the flare association is committed evidence rather than prose
    if "t0_utc" in burst:  # pragma: no cover - real data only
        clean = solarbursts.background_subtract(burst["data"])
        t_pk = float(burst["times"][int(np.argmax(clean.sum(axis=0)))])
        metrics["burst_peak_utc"] = str(
            np.datetime64(burst["t0_utc"].split(".")[0]) + np.timedelta64(int(round(t_pk)), "s")
        )
    if truth is not None:
        metrics["truth_speed_c"] = truth
        if np.isfinite(spd["speed_c"]):
            metrics["recovery_ratio"] = round(spd["speed_c"] / truth, 3) if truth else None

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "windwaves_metrics.json")
    if not offline:  # pragma: no cover - the real ridge is committed evidence
        _write_ridge(rf, rt, source, pad_s, op / "results" / "windwaves_ridge.csv")
    _figure(burst, rf, rt, harmonic, op / "papers" / "windwaves" / "figures")
    _write_macros(metrics, op / "papers" / "windwaves" / "generated" / "macros.tex")
    return metrics


def _speed_metrics(source, rf, rt, spd, harmonic, pad_s) -> dict:
    """The shared metrics block for the heliospheric drift slices (windwaves + swaves)."""

    def _r(key, nd=4):
        v = spd.get(key)
        return round(v, nd) if v is not None and np.isfinite(v) else None

    return {
        "source": source,
        "n_ridge": int(rf.size),
        "n_used": spd["n_used"],
        "n_time_cols": spd["n_time_cols"],
        "r2": _r("r2", 3),
        "f_lo_mhz": round(float(np.min(rf)), 4) if rf.size else None,
        "f_hi_mhz": round(float(np.max(rf)), 3) if rf.size else None,
        "harmonic": harmonic,
        "pad_s": float(pad_s),
        "snr_threshold": 5.0,
        "fit_estimator": "column-median OLS with leave-one-column-out jackknife",
        "r_lo_rsun": _r("r_lo", 2),
        "r_hi_rsun": _r("r_hi", 2),
        "r_hi_au": _r("r_hi_au", 3),
        "speed_kms": int(round(spd["speed_kms"], -2)) if np.isfinite(spd["speed_kms"]) else None,
        "speed_c": _r("speed_c", 4),
        "speed_c_se": _r("speed_c_se", 4),
        "speed_c_inverse": _r("speed_c_inverse", 4),
        "speed_c_points": _r("speed_c_points", 4),
        "speed_grid": speed_grid(rf, rt),
    }


def _write_ridge(rf, rt, source, pad_s, path) -> None:
    import csv as _csv
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        fh.write(f"# {source}; pad_s={pad_s:g}; snr_threshold=5\n")
        w = _csv.writer(fh)
        w.writerow(["freq_mhz", "time_s"])
        w.writerows(zip(np.round(rf, 4), np.round(rt, 2), strict=True))


def _figure(burst, rf, rt, harmonic, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    clean = solarbursts.background_subtract(burst["data"])
    freqs, times = burst["freqs"], burst["times"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    ax1.pcolormesh(times, freqs, clean, cmap="inferno", shading="auto")
    ax1.plot(rt, rf, ".", color="cyan", ms=2, label="ridge")
    ax1.set(
        xlabel="time (s)", ylabel="frequency (MHz)", yscale="log", title="Interplanetary type III"
    )
    ax1.legend(loc="upper right", fontsize=8)
    if rf.size >= 2:
        # draw the SAME fit the paper quotes: one point per time sample (column medians)
        r = emission_radius(rf, harmonic=harmonic)
        spd = beam_speed(rf, rt, harmonic=harmonic)
        cols = np.unique(rt)
        col_r = np.asarray([float(np.median(r[rt == c])) for c in cols])
        ax2.plot(rt, r, ".", color="0.75", ms=2, label="per-channel")
        ax2.plot(cols, col_r, "o", color="C0", ms=4, label="per-sample median")
        if np.isfinite(spd.get("fit_slope", float("nan"))):
            ax2.plot(
                cols,
                spd["fit_slope"] * cols + spd["fit_intercept"],
                "-",
                color="C3",
                lw=1,
                label="column fit",
            )
        ax2.legend(fontsize=7)
        ax2.set(xlabel="time (s)", ylabel=r"heliocentric radius ($R_\odot$)", title="Beam track")
    fig.tight_layout()
    fig.savefig(out / "ipburst.pdf")
    plt.close(fig)


def _write_macros(m: dict, path, *, prefix: str = "ww", paper: str = "windwaves") -> None:
    """Namespaced macros (``<prefix>Syn*``/``<prefix>Real*``), shared by windwaves and swaves."""
    from pathlib import Path

    real = not str(m.get("source", "")).lower().startswith("synthetic")
    ns, other = (f"{prefix}Real", f"{prefix}Syn") if real else (f"{prefix}Syn", f"{prefix}Real")

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    grid = {(g["harmonic"], g["density_scale"]): g for g in m.get("speed_grid", [])}
    fund = grid.get((1, 1.0), {})
    speeds = [g["speed_c"] for g in m.get("speed_grid", []) if g.get("speed_c") is not None]
    reaches = [g["r_hi_rsun"] for g in m.get("speed_grid", []) if g.get("r_hi_rsun") is not None]

    def _g(v, nd=4):
        return "--" if v is None else str(round(v, nd))

    values = (
        ("Nridge", _fmt("n_ridge")),
        ("Nused", _fmt("n_used")),
        ("NtimeCols", _fmt("n_time_cols")),
        ("Rsq", _fmt("r2")),
        ("Flo", _fmt("f_lo_mhz")),
        ("Fhi", _fmt("f_hi_mhz")),
        ("Rlo", _fmt("r_lo_rsun")),
        ("Rhi", _fmt("r_hi_rsun")),
        ("RhiAU", _fmt("r_hi_au")),
        ("SpeedKms", _fmt("speed_kms")),
        ("SpeedC", _fmt("speed_c")),
        ("SpeedCSe", _fmt("speed_c_se")),
        ("SpeedCInverse", _fmt("speed_c_inverse")),
        ("SpeedCPoints", _fmt("speed_c_points")),
        ("SpeedCFund", _g(fund.get("speed_c"))),
        ("RhiFund", _g(fund.get("r_hi_rsun"), 2)),
        ("RhiAUFund", _g(fund.get("r_hi_au"), 3)),
        ("GridSpeedLo", _g(min(speeds)) if speeds else "--"),
        ("GridSpeedHi", _g(max(speeds)) if speeds else "--"),
        ("GridReachLo", _g(min(reaches), 1) if reaches else "--"),
        ("GridReachHi", _g(max(reaches), 1) if reaches else "--"),
        ("PeakUTC", _fmt("burst_peak_utc")),
        ("Truth", _fmt("truth_speed_c")),
        ("Ratio", _fmt("recovery_ratio")),
    )
    lines = [
        f"% Auto-generated by jansky_research.{paper}._write_macros -- do not edit by hand.",
        "% Mode-dependent values are namespaced (Syn*/Real*); the inactive namespace holds",
        "% placeholders so offline CI and real runs never collide.",
        rf"\newcommand{{\{prefix}Source}}{{{m['source']}}}",
        rf"\newcommand{{\{prefix}Harmonic}}{{{m['harmonic']}}}",
    ]
    for suffix, value in values:
        lines.append(rf"\newcommand{{\{ns}{suffix}}}{{{value}}}")
        lines.append(rf"\newcommand{{\{other}{suffix}}}{{--}}")
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

    p = argparse.ArgumentParser(description="Interplanetary type III beam speed (Wind/WAVES).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--date", help="YYYYMMDD")
    p.add_argument("--receiver", default="rad2", choices=["rad1", "rad2"])
    p.add_argument("--harmonic", type=int, default=2)
    p.add_argument("--pad", type=float, default=1200.0)
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline or not args.date,
        date=args.date,
        receiver=args.receiver,
        harmonic=args.harmonic,
        pad_s=args.pad,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
