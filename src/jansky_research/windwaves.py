"""Interplanetary type III bursts: tracking the electron beam into the heliosphere (Wind/WAVES).

A type III electron beam does not stop at the top of the corona --- it streams out along the open
field into interplanetary space, exciting radio emission at the falling plasma frequency all the way
toward 1 AU. Space-based receivers see this as a slow drift from a few MHz down to tens of kHz over
tens of minutes. This slice fits that drift in a Wind/WAVES dynamic spectrum and inverts it, via a
**heliospheric** density model (Leblanc, Dulk & Bougeret 1998), to the beam's outward speed and the
heliocentric distance it reaches --- the interplanetary companion to the coronal ``solarbursts`` slice
(which used the Newkirk corona, valid only to a few solar radii).

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
    "beam_speed",
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


def emission_radius(freq_mhz: np.ndarray, *, harmonic: int = 2) -> np.ndarray:
    """Heliocentric radius (R⊙) of emission at ``freq_mhz``, via plasma frequency and the Leblanc model.

    The observed frequency is the ``harmonic`` of the local plasma frequency, so
    :math:`f_p = f/\\mathrm{harmonic}` gives the density
    (``jansky.solar.density_from_plasma_frequency``) and hence the radius (:func:`leblanc_radius`).
    """
    from jansky import solar

    fp = np.asarray(freq_mhz, float) / harmonic
    return leblanc_radius(solar.density_from_plasma_frequency(fp))


def beam_speed(
    ridge_freqs_mhz: np.ndarray, ridge_times_s: np.ndarray, *, harmonic: int = 2
) -> dict:
    """Outward beam speed from the drift ridge, via the Leblanc heliospheric density model.

    Maps each ridge frequency to a heliocentric radius (:func:`emission_radius`), robustly fits radius
    versus time (reusing ``solarbursts._robust_linfit`` for outlier rejection), and returns the radial
    speed in km/s and in units of :math:`c`, the radius range (R⊙ and AU), the fit :math:`R^2`, and the
    number of ridge points used.
    """
    f = np.asarray(ridge_freqs_mhz, float)
    t = np.asarray(ridge_times_s, float)
    r = emission_radius(f, harmonic=harmonic)
    nan = float("nan")
    if r.size < 3 or np.ptp(t) == 0:
        return {
            "speed_kms": nan,
            "speed_c": nan,
            "r_lo": nan,
            "r_hi": nan,
            "r_hi_au": nan,
            "r2": nan,
            "n_used": 0,
        }
    dr_dt, icpt, keep = solarbursts._robust_linfit(t, r)
    speed_kms = abs(dr_dt) * 695700.0  # R_sun -> km
    rk, tk = r[keep], t[keep]
    ss_res = float(np.sum((rk - (dr_dt * tk + icpt)) ** 2))
    ss_tot = float(np.sum((rk - np.mean(rk)) ** 2))
    return {
        "speed_kms": float(speed_kms),
        "speed_c": float(speed_kms / C_KMS),
        "r_lo": float(np.min(rk)),
        "r_hi": float(np.max(rk)),
        "r_hi_au": float(np.max(rk) / R_AU_RSUN),
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else nan,
        "n_used": int(keep.sum()),
    }


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
    return {"data": psd.T[order], "freqs": freq[order], "times": np.asarray(times, float)}


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
    import json
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
    metrics: dict = {
        "source": source,
        "n_ridge": int(rf.size),
        "n_used": spd["n_used"],
        "r2": round(spd["r2"], 3) if np.isfinite(spd["r2"]) else None,
        "f_lo_mhz": round(float(np.min(rf)), 3) if rf.size else None,
        "f_hi_mhz": round(float(np.max(rf)), 3) if rf.size else None,
        "harmonic": harmonic,
        "r_lo_rsun": round(spd["r_lo"], 2) if np.isfinite(spd["r_lo"]) else None,
        "r_hi_rsun": round(spd["r_hi"], 2) if np.isfinite(spd["r_hi"]) else None,
        "r_hi_au": round(spd["r_hi_au"], 3) if np.isfinite(spd["r_hi_au"]) else None,
        "speed_kms": round(spd["speed_kms"], 1) if np.isfinite(spd["speed_kms"]) else None,
        "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
    }
    if truth is not None:
        metrics["truth_speed_c"] = truth
        if np.isfinite(spd["speed_c"]):
            metrics["recovery_ratio"] = round(spd["speed_c"] / truth, 3) if truth else None

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "windwaves_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(burst, rf, rt, harmonic, op / "papers" / "windwaves" / "figures")
    _write_macros(metrics, op / "papers" / "windwaves" / "generated" / "macros.tex")
    return metrics


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
        r = emission_radius(rf, harmonic=harmonic)
        ax2.plot(rt, r, "o", color="C0", ms=3)
        slope, icpt, keep = solarbursts._robust_linfit(rt, r)
        ax2.plot(rt, slope * rt + icpt, "-", color="C3", lw=1)
        ax2.set(xlabel="time (s)", ylabel=r"heliocentric radius ($R_\odot$)", title="Beam track")
    fig.tight_layout()
    fig.savefig(out / "ipburst.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.windwaves._write_macros -- do not edit by hand.",
        rf"\newcommand{{\wwSource}}{{{m['source']}}}",
        rf"\newcommand{{\wwNridge}}{{{m['n_ridge']}}}",
        rf"\newcommand{{\wwNused}}{{{_fmt('n_used')}}}",
        rf"\newcommand{{\wwRsq}}{{{_fmt('r2')}}}",
        rf"\newcommand{{\wwFlo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\wwFhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\wwHarmonic}}{{{m['harmonic']}}}",
        rf"\newcommand{{\wwRlo}}{{{_fmt('r_lo_rsun')}}}",
        rf"\newcommand{{\wwRhi}}{{{_fmt('r_hi_rsun')}}}",
        rf"\newcommand{{\wwRhiAU}}{{{_fmt('r_hi_au')}}}",
        rf"\newcommand{{\wwSpeedKms}}{{{_fmt('speed_kms')}}}",
        rf"\newcommand{{\wwSpeedC}}{{{_fmt('speed_c')}}}",
        rf"\newcommand{{\wwTruth}}{{{_fmt('truth_speed_c')}}}",
        rf"\newcommand{{\wwRatio}}{{{_fmt('recovery_ratio')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
