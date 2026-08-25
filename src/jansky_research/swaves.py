"""Type III to true interplanetary distances with STEREO/WAVES (HFR).

The Wind/WAVES slice tracked a type III electron beam only to the Alfven surface (~10 R_sun), because
RAD2 stops at 1 MHz. STEREO/WAVES (the SWAVES instrument; Bougeret et al. 2008) reaches down to
**0.125 MHz** with its High Frequency Receiver, where the emission comes from genuinely
**interplanetary** (super-Alfvenic, >~20 R_sun) plasma --- out to ~0.4 AU. This slice fits a
STEREO/WAVES type III drift across the full HFR band and inverts it, via the Leblanc heliospheric
density model, to the beam's outward speed and the true interplanetary distance it reaches.

It reuses the Wind/WAVES Leblanc tooling (``windwaves.emission_radius`` / ``beam_speed``) and the
``solarbursts`` dynamic-spectrum pipeline almost wholesale; the new code is the STEREO/WAVES one-minute
HFR ASCII parser (no auth, no CDF). Pure NumPy with a synthetic offline fixture; the real fetch is
network-gated.
"""

from __future__ import annotations

import numpy as np

from . import solarbursts, windwaves

__all__ = [
    "fetch_swaves",
    "parse_swaves_ascii",
    "run",
    "synthetic_ip_burst",
]


def parse_swaves_ascii(text: str) -> dict:
    """Parse a STEREO/WAVES one-minute HFR ASCII file → ``data`` (freq × time), ``freqs``, ``times``.

    The file's first line is the frequency axis (kHz); a background row follows; each subsequent data
    row is a minute index followed by one intensity per frequency. Returns frequencies in MHz
    (descending), times in seconds from the file start, and the intensity matrix as ``(n_freq, n_time)``.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    freqs_khz = np.array(lines[0].split(), float)
    nf = freqs_khz.size
    minutes: list[float] = []
    rows: list[list[float]] = []
    for ln in lines[1:]:
        parts = ln.split()
        if (
            len(parts) == nf + 1
        ):  # minute index + one intensity per frequency (skip the background row)
            minutes.append(float(parts[0]))
            rows.append([float(x) for x in parts[1:]])
    data = np.asarray(rows, float).T  # (freq, time)
    freqs_mhz = freqs_khz / 1000.0
    order = np.argsort(freqs_mhz)[::-1]  # descending, to match the solarbursts/windwaves convention
    return {
        "data": data[order],
        "freqs": freqs_mhz[order],
        "times": np.asarray(minutes, float) * 60.0,
    }


def synthetic_ip_burst(*, speed_c: float = 0.15, harmonic: int = 2, seed: int = 0) -> dict:
    """Synthetic interplanetary type III over the STEREO/WAVES HFR band (0.125--16 MHz).

    Thin wrapper over :func:`windwaves.synthetic_ip_burst` with HFR defaults (lower frequencies, a
    longer drift reaching interplanetary distances), so a clean burst round-trips the Leblanc inversion.
    """
    return windwaves.synthetic_ip_burst(
        speed_c=speed_c,
        harmonic=harmonic,
        seed=seed,
        f_lo_mhz=0.125,
        f_hi_mhz=16.0,
        r0_rsun=1.6,
        duration_s=2400.0,
        n_time=320,
        n_freq=319,
        width_dex=0.05,
    )


def fetch_swaves(
    date_yyyymmdd: str, *, spacecraft: str = "a"
) -> dict:  # pragma: no cover - network
    """Fetch a STEREO/WAVES one-minute HFR ASCII dynamic spectrum for one day from SPDF (no auth).

    ``spacecraft`` is ``"a"`` (ahead) or ``"b"`` (behind). Returns the parsed dynamic spectrum
    (see :func:`parse_swaves_ascii`).
    """
    import re

    import requests

    side = "ahead" if spacecraft.lower() == "a" else "behind"
    yyyy = date_yyyymmdd[:4]
    base = f"https://spdf.gsfc.nasa.gov/pub/data/stereo/{side}/swaves/one-minute/ascii/hfr/{yyyy}/"
    idx = requests.get(base, timeout=60).text
    pat = rf"stereo-{spacecraft.lower()}_swaves_hfr_average_{date_yyyymmdd}_v[0-9]+\.txt"
    m = re.findall(pat, idx)
    if not m:
        raise RuntimeError(f"no STEREO-{spacecraft.upper()} SWAVES HFR file for {date_yyyymmdd}")
    out = parse_swaves_ascii(requests.get(base + m[0], timeout=120).text)
    # minutes count from file start = 00:00 UT of the file's day
    d = date_yyyymmdd
    out["t0_utc"] = f"{d[:4]}-{d[4:6]}-{d[6:]}T00:00:00"
    return out


def run(
    out: str = ".",
    *,
    offline: bool = True,
    date: str | None = None,
    spacecraft: str = "a",
    harmonic: int = 2,
    pad_s: float = 2400.0,
) -> dict:
    """Full slice: fit a STEREO/WAVES type III drift and report the beam speed and interplanetary reach."""
    from pathlib import Path

    if offline or date is None:
        burst = synthetic_ip_burst(harmonic=harmonic)
        source = "synthetic"
        truth: float | None = burst["truth_speed_c"]
    else:  # pragma: no cover - network
        burst = fetch_swaves(date, spacecraft=spacecraft)
        source = f"STEREO-{spacecraft.upper()}/WAVES HFR {date}"
        truth = None

    window = solarbursts.find_burst_window(burst["data"], burst["times"], pad_s=pad_s)
    rf, rt = solarbursts.detect_burst_ridge(
        burst["data"], burst["freqs"], burst["times"], window=window
    )
    spd = windwaves.beam_speed(rf, rt, harmonic=harmonic)
    # the one-minute cadence makes many high-frequency channels share a timestamp, so the
    # *effective* number of independent time samples is the time-column count in the fit
    metrics: dict = windwaves._speed_metrics(source, rf, rt, spd, harmonic, pad_s)
    metrics["n_time_bins"] = metrics["n_time_cols"]
    if "t0_utc" in burst:  # pragma: no cover - real data only
        clean = solarbursts.background_subtract(burst["data"])
        t_pk = float(burst["times"][int(np.argmax(clean.sum(axis=0)))])
        metrics["burst_peak_utc"] = str(
            np.datetime64(burst["t0_utc"]) + np.timedelta64(int(round(t_pk)), "s")
        )
    if truth is not None:
        metrics["truth_speed_c"] = truth
        if np.isfinite(spd["speed_c"]):
            metrics["recovery_ratio"] = round(spd["speed_c"] / truth, 3) if truth else None

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "swaves_metrics.json")
    if not offline:  # pragma: no cover - the real ridge is committed evidence
        windwaves._write_ridge(rf, rt, source, pad_s, op / "results" / "swaves_ridge.csv")
    _figure(burst, rf, rt, harmonic, op / "papers" / "swaves" / "figures")
    _write_macros(metrics, op / "papers" / "swaves" / "generated" / "macros.tex")
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
    ax1.pcolormesh(times / 60.0, freqs, clean, cmap="inferno", shading="auto")
    ax1.plot(rt / 60.0, rf, ".", color="cyan", ms=2, label="ridge")
    ax1.set(
        xlabel="time (min)", ylabel="frequency (MHz)", yscale="log", title="STEREO/WAVES type III"
    )
    ax1.legend(loc="upper right", fontsize=8)
    if rf.size >= 2:
        # draw the SAME fit the paper quotes: one point per one-minute sample
        r = windwaves.emission_radius(rf, harmonic=harmonic)
        spd = windwaves.beam_speed(rf, rt, harmonic=harmonic)
        cols = np.unique(rt)
        col_r = np.asarray([float(np.median(r[rt == c])) for c in cols])
        au = windwaves.R_AU_RSUN
        ax2.plot(rt / 60.0, r / au, ".", color="0.75", ms=2, label="per-channel")
        ax2.plot(cols / 60.0, col_r / au, "o", color="C0", ms=4, label="per-sample median")
        if np.isfinite(spd.get("fit_slope", float("nan"))):
            ax2.plot(
                cols / 60.0,
                (spd["fit_slope"] * cols + spd["fit_intercept"]) / au,
                "-",
                color="C3",
                lw=1,
                label="sample fit",
            )
        ax2.legend(fontsize=7)
        ax2.set(xlabel="time (min)", ylabel="heliocentric distance (AU)", title="Beam track")
    fig.tight_layout()
    fig.savefig(out / "swaves.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    """Namespaced swSyn*/swReal* macros via the shared heliospheric writer."""
    windwaves._write_macros(m, path, prefix="sw", paper="swaves")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Type III to interplanetary distances (STEREO/WAVES HFR)."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--date", help="YYYYMMDD")
    p.add_argument("--spacecraft", default="a", choices=["a", "b"])
    p.add_argument("--harmonic", type=int, default=2)
    p.add_argument("--pad", type=float, default=2400.0)
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline or not args.date,
        date=args.date,
        spacecraft=args.spacecraft,
        harmonic=args.harmonic,
        pad_s=args.pad,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
