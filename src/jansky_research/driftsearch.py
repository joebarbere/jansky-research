"""CPU-only SETI Doppler-drift injection-recovery benchmark.

Technosignature searches look for narrowband signals that drift in frequency (the transmitter's
relative acceleration). Many surveys quote an injection-recovery efficiency, but there is no shared,
reproducible, CPU-only benchmark to compare detectors on the *same* reference set
(see ``survey/literature.md``). This module provides one, built on the pure-NumPy drift search in
``jansky.seti``: inject synthetic drifting tones over a grid of signal-to-noise ratios and drift
rates, run the brute-force de-drift search, and measure the recovered fraction
$P_\\mathrm{detect}(\\mathrm{SNR}, \\dot f)$ — plus the noise-only false-positive rate that
calibrates the detection threshold.

Everything is offline and seedable, so the benchmark is fully reproducible. Pointing the same
detector at real data (the Voyager-1 file, :func:`validate_voyager`) is an honest *negative* check:
it bounds where this synthetic-tuned teaching detector works (injected tones) and where it does not
(the real, drifting Voyager carrier amid a band-centre DC spike).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from jansky import seti

__all__ = [
    "RecoveryResult",
    "completeness_snr",
    "false_positive_rate",
    "injection_recovery",
    "locate_carrier",
    "measure_drift",
    "noise_peak_stats",
    "run",
    "run_voyager",
    "validate_voyager",
]


@dataclass(frozen=True)
class RecoveryResult:
    """Injection-recovery efficiency over a grid of injected SNR and drift rate."""

    inj_snrs: np.ndarray
    drift_rates: np.ndarray
    p_detect: np.ndarray  # shape (n_snr, n_drift): recovered fraction per cell
    threshold: float
    n_trials: int
    false_positive_rate: float
    completeness_snr_50: float
    completeness_snr_90: float


def _detect(snr, drift, *, search_drifts, threshold, n_time, n_freq, noise, present, rng) -> bool:
    """One trial: inject (or not) a tone, run the drift search, return whether it clears threshold."""
    wf = seti.drifting_tone(
        n_time,
        n_freq,
        drift_rate=drift,
        snr=snr,
        noise=noise,
        present=present,
        seed=int(rng.integers(2**31)),
    )
    return bool(seti.drift_search(wf, search_drifts).best_snr > threshold)


def false_positive_rate(
    *,
    n_trials=400,
    threshold=10.0,
    n_time=64,
    n_freq=512,
    noise=1.0,
    search_drifts=None,
    seed=0,
) -> float:
    """Fraction of *noise-only* waterfalls whose best drift-search S/N exceeds ``threshold``."""
    if search_drifts is None:
        search_drifts = np.linspace(-1.0, 1.0, 41)
    rng = np.random.default_rng(seed)
    hits = sum(
        _detect(
            10.0,
            0.0,
            search_drifts=search_drifts,
            threshold=threshold,
            n_time=n_time,
            n_freq=n_freq,
            noise=noise,
            present=False,
            rng=rng,
        )
        for _ in range(n_trials)
    )
    return hits / n_trials


def completeness_snr(inj_snrs: np.ndarray, p_mean: np.ndarray, level: float = 0.5) -> float:
    """Injected SNR at which the recovered fraction first crosses ``level`` (linear interpolation)."""
    inj_snrs = np.asarray(inj_snrs, float)
    p_mean = np.asarray(p_mean, float)
    if p_mean.max() < level:
        return float("nan")
    # np.interp needs a monotonic-increasing xp; sort by p_mean so finite-trial
    # non-monotonicity at the low-SNR end can't silently return a wrong crossing.
    order = np.argsort(p_mean)
    return float(np.interp(level, p_mean[order], inj_snrs[order]))


def injection_recovery(
    inj_snrs: np.ndarray,
    *,
    drift_rates: np.ndarray | None = None,
    n_trials: int = 20,
    threshold: float = 10.0,
    n_time: int = 64,
    n_freq: int = 512,
    noise: float = 1.0,
    search_drifts: np.ndarray | None = None,
    fpr_trials: int = 400,
    seed: int = 0,
) -> RecoveryResult:
    """Measure detector recovery efficiency over the injected (SNR × drift) grid.

    For each cell, ``n_trials`` drifting tones are injected and searched; the recovered fraction is
    the detection probability. Returns the full matrix plus drift-averaged 50% and 90% completeness
    SNRs and the noise-only false-positive rate at this ``threshold``.
    """
    if drift_rates is None:
        drift_rates = np.array([0.0, 0.3, 0.6])
    if search_drifts is None:
        search_drifts = np.linspace(-1.0, 1.0, 41)
    inj_snrs = np.asarray(inj_snrs, float)
    rng = np.random.default_rng(seed)
    p = np.zeros((inj_snrs.size, drift_rates.size))
    for i, snr in enumerate(inj_snrs):
        for j, dr in enumerate(drift_rates):
            hits = sum(
                _detect(
                    snr,
                    dr,
                    search_drifts=search_drifts,
                    threshold=threshold,
                    n_time=n_time,
                    n_freq=n_freq,
                    noise=noise,
                    present=True,
                    rng=rng,
                )
                for _ in range(n_trials)
            )
            p[i, j] = hits / n_trials
    p_mean = p.mean(axis=1)
    fpr = false_positive_rate(
        n_trials=fpr_trials,
        threshold=threshold,
        n_time=n_time,
        n_freq=n_freq,
        noise=noise,
        search_drifts=search_drifts,
        seed=seed + 1,
    )
    return RecoveryResult(
        inj_snrs=inj_snrs,
        drift_rates=drift_rates,
        p_detect=p,
        threshold=threshold,
        n_trials=n_trials,
        false_positive_rate=fpr,
        completeness_snr_50=completeness_snr(inj_snrs, p_mean, 0.5),
        completeness_snr_90=completeness_snr(inj_snrs, p_mean, 0.9),
    )


def noise_peak_stats(
    *,
    n_draws: int = 400,
    n_time: int = 64,
    n_freq: int = 512,
    noise: float = 1.0,
    search_drifts: np.ndarray | None = None,
    seed: int = 1,
) -> dict:
    """Distribution of the noise-only best drift-search S/N — the honest FPR statement.

    "0 false positives in 400 trials at threshold 10" could not have failed: measured, the
    noise-only best S/N sits ~18 sigma below that threshold. The informative quantities are the
    distribution's mean/sd/p99/max and the one-sided Clopper–Pearson bound at the stated trial
    count, all committed here rather than asserted.
    """
    if search_drifts is None:
        search_drifts = np.linspace(-1.0, 1.0, 41)
    rng = np.random.default_rng(seed)
    best = []
    for _ in range(n_draws):
        wf = seti.drifting_tone(
            n_time,
            n_freq,
            drift_rate=0.0,
            snr=0.0,
            noise=noise,
            present=False,
            seed=int(rng.integers(2**31)),
        )
        best.append(float(seti.drift_search(wf, search_drifts).best_snr))
    arr = np.asarray(best)
    return {
        "n_draws": n_draws,
        "mean": round(float(arr.mean()), 3),
        "sd": round(float(arr.std(ddof=1)), 3),
        "p99": round(float(np.percentile(arr, 99)), 3),
        "max": round(float(arr.max()), 3),
        # one-sided 95% CP upper bound on the FPR given zero exceedances in n_draws
        "fpr_upper_95_one_sided": round(1.0 - 0.05 ** (1.0 / n_draws), 5),
    }


def run(out: str = ".", *, n_trials: int = 30, threshold: float = 10.0, seed: int = 0) -> dict:
    """Compute the injection-recovery benchmark; write metrics + a recovery heatmap. Returns metrics.

    ``n_trials`` defaults to 30, matching the CLI and the committed evidence (an earlier
    signature default of 100 disagreed with both).
    """
    from pathlib import Path

    inj_snrs = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0])
    n_time, n_freq, fpr_trials = 64, 512, 400
    search_drifts = np.linspace(-1.0, 1.0, 41)
    res = injection_recovery(
        inj_snrs,
        n_trials=n_trials,
        threshold=threshold,
        n_time=n_time,
        n_freq=n_freq,
        search_drifts=search_drifts,
        fpr_trials=fpr_trials,
        seed=seed,
    )
    # the off-grid check: injecting half a grid step off the searched drifts must not move
    # the completeness materially (dedrift rounds to integer channel shifts, so it does not)
    off = injection_recovery(
        inj_snrs,
        drift_rates=np.array([0.325, 0.475]),
        n_trials=n_trials,
        threshold=threshold,
        n_time=n_time,
        n_freq=n_freq,
        search_drifts=search_drifts,
        fpr_trials=4,  # not the statistic of interest here
        seed=seed,
    )
    metrics = {
        "source": "synthetic injection-recovery benchmark (offline, seeded)",
        "config": {
            "n_time": n_time,
            "n_freq": n_freq,
            "n_search_drifts": int(search_drifts.size),
            "search_drift_min": float(search_drifts.min()),
            "search_drift_max": float(search_drifts.max()),
            "noise": 1.0,
            "fpr_trials": fpr_trials,
            "seed": seed,
        },
        "inj_snrs": res.inj_snrs.tolist(),
        "drift_rates": res.drift_rates.tolist(),
        # the FULL matrix, not just the drift average: "flat across drift" must be auditable
        "p_detect": res.p_detect.tolist(),
        "p_detect_mean": res.p_detect.mean(axis=1).tolist(),
        "completeness_snr_50_per_drift": [
            completeness_snr(res.inj_snrs, res.p_detect[:, j], 0.5)
            for j in range(res.drift_rates.size)
        ],
        "threshold": res.threshold,
        "n_trials": res.n_trials,
        "false_positive_rate": res.false_positive_rate,
        "noise_peak_stats": noise_peak_stats(
            n_draws=fpr_trials,
            n_time=n_time,
            n_freq=n_freq,
            search_drifts=search_drifts,
            seed=seed + 1,
        ),
        "off_grid_check": {
            "drift_rates": off.drift_rates.tolist(),
            "completeness_snr_50": off.completeness_snr_50,
        },
        "completeness_snr_50": res.completeness_snr_50,
        "completeness_snr_90": res.completeness_snr_90,
    }
    op = Path(out)
    paper = op / "papers" / "driftsearch"
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "drift_metrics.json")
    _heatmap(res, paper / "figures")
    _write_macros(metrics, paper / "generated" / "macros.tex")
    return metrics


def _write_macros(m: dict | None, path, voyager: dict | None = None) -> None:
    """Emit the FULL macro union (benchmark + Voyager namespaces) with placeholders.

    Both legs write the same file, so each emits every name — its own values plus ``--`` for
    the other leg's — and ``preserve_live_macros`` restores the real values on merge (the
    repo's two-namespace pattern; a writer that emitted only its own names would DELETE the
    other leg's numbers).
    """
    from pathlib import Path

    def _b(fmt: str, key: str, scale=None) -> str:
        if m is None or key not in m:
            return "--"
        v = m[key]
        return format(scale(v) if scale else v, fmt)

    lines = [
        "% Auto-generated by jansky_research.driftsearch._write_macros — do not edit by hand.",
        rf"\newcommand{{\dsBenchSource}}{{{m['source'] if m else '--'}}}",
        rf"\newcommand{{\dsThreshold}}{{{_b('.0f', 'threshold')}}}",
        rf"\newcommand{{\dsNtrials}}{{{_b('d', 'n_trials')}}}",
        rf"\newcommand{{\dsCfifty}}{{{_b('.1f', 'completeness_snr_50')}}}",
        rf"\newcommand{{\dsCninety}}{{{_b('.1f', 'completeness_snr_90')}}}",
        rf"\newcommand{{\dsNdrift}}{{{len(m['drift_rates']) if m else '--'}}}",
        rf"\newcommand{{\dsDriftMax}}{{{max(m['drift_rates']) if m else '--'}}}",
    ]
    if m:
        cfg = m["config"]
        nps = m["noise_peak_stats"]
        lines += [
            rf"\newcommand{{\dsNtime}}{{{cfg['n_time']}}}",
            rf"\newcommand{{\dsNfreq}}{{{cfg['n_freq']}}}",
            rf"\newcommand{{\dsNsearchDrifts}}{{{cfg['n_search_drifts']}}}",
            rf"\newcommand{{\dsFprTrials}}{{{cfg['fpr_trials']}}}",
            rf"\newcommand{{\dsNoiseMean}}{{{nps['mean']:.1f}}}",
            rf"\newcommand{{\dsNoiseSd}}{{{nps['sd']:.2f}}}",
            rf"\newcommand{{\dsNoiseMax}}{{{nps['max']:.1f}}}",
            rf"\newcommand{{\dsFprBound}}{{{100 * nps['fpr_upper_95_one_sided']:.2f}}}",
            rf"\newcommand{{\dsOffGridCfifty}}{{{m['off_grid_check']['completeness_snr_50']:.2f}}}",
        ]
    else:
        lines += [
            rf"\newcommand{{\ds{k}}}{{--}}"
            for k in (
                "Ntime",
                "Nfreq",
                "NsearchDrifts",
                "FprTrials",
                "NoiseMean",
                "NoiseSd",
                "NoiseMax",
                "FprBound",
                "OffGridCfifty",
            )
        ]
    if voyager:
        car = voyager["carrier"]
        lines += [
            rf"\newcommand{{\dsVoySource}}{{{voyager['source']}}}",
            rf"\newcommand{{\dsVoySnr}}{{{car['snr']:.0f}}}",
            rf"\newcommand{{\dsVoyFreq}}{{{car['freq_mhz']:.5f}}}",
            rf"\newcommand{{\dsVoyDriftHz}}{{{car['measured_drift_hz_s']:.3f}}}",
            rf"\newcommand{{\dsVoyBestDrift}}{{{car['search_best_drift_chan_per_sample']:.2f}}}",
            rf"\newcommand{{\dsVoyDriftChan}}{{{car['measured_drift_chan_per_sample']:.3f}}}",
            rf"\newcommand{{\dsVoyDcSnr}}{{{voyager['dc_spike']['search_snr']:.2g}}}",
            rf"\newcommand{{\dsVoyLegacySnr}}{{{voyager['legacy_asserted']['snr']:.2f}}}",
            rf"\newcommand{{\dsVoyLegacyFreq}}{{{voyager['legacy_asserted']['freq_mhz']}}}",
            rf"\newcommand{{\dsVoyBlankSnr}}{{{voyager['blank_snr']:.2f}}}",
        ]
    else:
        lines += [
            rf"\newcommand{{\dsVoy{k}}}{{--}}"
            for k in (
                "Source",
                "Snr",
                "Freq",
                "DriftHz",
                "BestDrift",
                "DriftChan",
                "DcSnr",
                "LegacySnr",
                "LegacyFreq",
                "BlankSnr",
            )
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a run may only ADD information, so an
    # offline rebuild can never blank a real value (report.preserve_live_macros).
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _heatmap(res: RecoveryResult, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 3.5))

    # pcolormesh with explicit cell edges: the injected S/N rows are NOT uniformly spaced
    # (0.25..1.5, 2.0, 3.0), and imshow's linear extent drew the 50% crossing ~0.5 in S/N
    # above where the data put it (round-11 referee)
    def _edges(v: np.ndarray) -> np.ndarray:
        mid = 0.5 * (v[1:] + v[:-1])
        return np.concatenate([[v[0] - (mid[0] - v[0])], mid, [v[-1] + (v[-1] - mid[-1])]])

    im = ax.pcolormesh(
        _edges(res.drift_rates),
        _edges(res.inj_snrs),
        res.p_detect,
        vmin=0,
        vmax=1,
        cmap="viridis",
        shading="flat",
    )
    ax.set(
        xlabel="injected drift rate (chan/sample)",
        ylabel="injected S/N",
        title=f"Drift-search recovery $P_{{det}}$ (FPR={res.false_positive_rate:.2g})",
    )
    fig.colorbar(im, ax=ax, label=r"$P_\mathrm{detect}$")
    fig.tight_layout()
    fig.savefig(out / "drift_recovery.pdf")
    plt.close(fig)


# The frequency an earlier version of this module ASSERTED for the Voyager-1 carrier and
# searched at. In this file that frequency maps to blank sky: the actual carrier sits ~0.92 MHz
# away (v/c ~ 33 km/s at X band — the scale of a topocentric/barycentric frame difference), a
# lesson kept on the record. The carrier is now LOCATED in the data, never asserted.
LEGACY_ASSERTED_CARRIER_MHZ = 8420.216


def locate_carrier(wf: np.ndarray, *, dc_halfwidth: int = 2048) -> dict:
    """Locate the brightest narrowband non-DC feature by MAD-z of the time-averaged spectrum.

    The band-centre DC spike (channel N/2) is orders of magnitude brighter than any real tone,
    so a brightest-channel search reports the artifact; masking the DC region and normalising by
    the median absolute deviation finds the brightest genuine narrowband feature. Returns the
    channel, its MAD-z, and the DC spike's own channel and z for the cautionary comparison.
    """
    col = np.asarray(wf, float).mean(axis=0)
    n = col.size
    med = float(np.median(col))
    mad = 1.4826 * float(np.median(np.abs(col - med)))
    z = (col - med) / mad if mad > 0 else np.zeros_like(col)
    dc_chan = int(np.argmax(z))
    zm = z.copy()
    lo = max(0, n // 2 - dc_halfwidth)
    zm[lo : n // 2 + dc_halfwidth + 1] = -np.inf
    chan = int(np.argmax(zm))
    return {
        "channel": chan,
        "mad_z": round(float(z[chan]), 1),
        "dc_channel": dc_chan,
        "dc_mad_z": round(float(z[dc_chan]), 1),
        "dc_is_band_centre": bool(abs(dc_chan - n // 2) <= dc_halfwidth),
    }


def measure_drift(wf: np.ndarray, channel: int, *, halfwidth: int = 200) -> dict:
    """Drift of a narrowband feature: linear fit to the per-sample peak channel near ``channel``."""
    sub = np.asarray(wf, float)[:, channel - halfwidth : channel + halfwidth + 1]
    peaks = np.argmax(sub, axis=1) + channel - halfwidth
    t = np.arange(peaks.size, dtype=float)
    slope, intercept = np.polyfit(t, peaks.astype(float), 1)
    resid = peaks - (slope * t + intercept)
    return {
        "chan_per_sample": round(float(slope), 4),
        "fit_rms_chan": round(float(np.std(resid)), 3),
        "n_samples": int(peaks.size),
    }


def validate_voyager(
    path=None, *, window: int = 4096
) -> dict:  # pragma: no cover - net + optional deps
    """Real-data check of the detector on the Breakthrough Listen Voyager-1 file.

    Two results, both measured from the file rather than asserted. (1) The recovery: the
    carrier is LOCATED in the data (brightest non-DC narrowband feature by MAD-z), its drift is
    measured from the per-sample peak walk, and the ``jansky.seti`` drift search evaluated in a
    window around it recovers it at S/N ~10^3 with a best-fit drift matching the measured walk.
    (2) The caution this slice has always carried: the band-centre DC-spike artifact is ~10^3x
    brighter than the carrier, so a brightest-channel search reports the artifact at S/N ~10^5
    and is wrong. An earlier version asserted the carrier frequency
    (:data:`LEGACY_ASSERTED_CARRIER_MHZ`) instead of locating it; in this file that frequency is
    blank sky, and the resulting "null" was a targeting error — the legacy numbers are kept in
    the output as the record of that lesson. Requires the optional ``voyager`` extra
    (``h5py`` + ``hdf5plugin``).
    """
    import h5py
    import hdf5plugin  # noqa: F401 - registers the bitshuffle filter

    from . import data as _data

    if path is None:
        path = _data.fetch("voyager1-h5")
    with h5py.File(path, "r") as f:
        wf = np.asarray(f["data"][:]).squeeze().astype(float)
        fch1 = float(f["data"].attrs["fch1"])
        foff = float(f["data"].attrs["foff"])
        tsamp = float(f["data"].attrs["tsamp"])
    n = wf.shape[1]
    drifts = np.linspace(-8.0, 8.0, 321)

    def _snr(center: int) -> tuple[float, float]:
        lo = max(0, center - window // 2)
        sub = wf[:, lo : lo + window]
        r = seti.drift_search(sub - np.median(sub), drifts)
        return float(r.best_snr), float(r.best_drift)

    loc = locate_carrier(wf)
    drift = measure_drift(wf, loc["channel"])
    carrier_snr, carrier_best_drift = _snr(loc["channel"])
    dc_snr, _ = _snr(loc["dc_channel"])
    legacy_chan = int(round((LEGACY_ASSERTED_CARRIER_MHZ - fch1) / foff))
    legacy_snr, _ = _snr(legacy_chan)
    blank_chan = 500_000 if abs(500_000 - n // 2) > window else n // 4
    blank_snr, _ = _snr(blank_chan)
    return {
        "source": "Breakthrough Listen Voyager-1 GBT file (real data)",
        "is_real": True,
        "window_channels": window,
        "n_channels": n,
        "carrier": {
            "channel": loc["channel"],
            "freq_mhz": round(fch1 + foff * loc["channel"], 5),
            "mad_z": loc["mad_z"],
            "measured_drift_chan_per_sample": drift["chan_per_sample"],
            "measured_drift_hz_s": round(drift["chan_per_sample"] * foff * 1e6 / tsamp, 4),
            "drift_fit_rms_chan": drift["fit_rms_chan"],
            "snr": round(carrier_snr, 1),
            "search_best_drift_chan_per_sample": carrier_best_drift,
        },
        "dc_spike": {
            "channel": loc["dc_channel"],
            "is_band_centre": loc["dc_is_band_centre"],
            "mad_z": loc["dc_mad_z"],
            "search_snr": round(dc_snr, 1),
        },
        "legacy_asserted": {
            "freq_mhz": LEGACY_ASSERTED_CARRIER_MHZ,
            "channel": legacy_chan,
            "snr": round(legacy_snr, 2),
            "note": "the earlier 'Voyager null' searched here: blank sky",
        },
        "blank_snr": round(blank_snr, 2),
        "recovered": bool(carrier_snr > blank_snr + 3.0),
    }


def run_voyager(out: str = ".", path=None) -> dict:  # pragma: no cover - net + optional deps
    """Run the real-data Voyager leg and COMMIT its evidence (the leg previously had none)."""
    from pathlib import Path

    from .report import write_results

    m = validate_voyager(path)
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    write_results(m, op / "results" / "drift_voyager.json")
    _write_macros(None, op / "papers" / "driftsearch" / "generated" / "macros.tex", voyager=m)
    return m


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="SETI drift-search injection-recovery benchmark.")
    p.add_argument("--out", default=".")
    p.add_argument("--n-trials", type=int, default=30)
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, n_trials=args.n_trials), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
