"""CPU-only SETI Doppler-drift injection-recovery benchmark.

Technosignature searches look for narrowband signals that drift in frequency (the transmitter's
relative acceleration). Many surveys quote an injection-recovery efficiency, but there is no shared,
reproducible, CPU-only benchmark to compare detectors on the *same* reference set
(see ``survey/literature.md``). This module provides one, built on the pure-NumPy drift search in
``jansky.seti``: inject synthetic drifting tones over a grid of signal-to-noise ratios and drift
rates, run the brute-force de-drift search, and measure the recovered fraction
$P_\\mathrm{detect}(\\mathrm{SNR}, \\dot f)$ — plus the noise-only false-positive rate that
calibrates the detection threshold.

Everything is offline and seedable, so the benchmark is fully reproducible. The same `jansky.seti`
detector can then be pointed at real data (e.g. the Voyager-1 carrier) as an external validation.
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
    "run",
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
    return float(np.interp(level, p_mean, inj_snrs))


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


def run(out: str = ".", *, n_trials: int = 30, threshold: float = 10.0, seed: int = 0) -> dict:
    """Compute the injection-recovery benchmark; write metrics + a recovery heatmap. Returns metrics."""
    import json
    from pathlib import Path

    inj_snrs = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0])
    res = injection_recovery(inj_snrs, n_trials=n_trials, threshold=threshold, seed=seed)
    metrics = {
        "inj_snrs": res.inj_snrs.tolist(),
        "drift_rates": res.drift_rates.tolist(),
        "p_detect_mean": res.p_detect.mean(axis=1).tolist(),
        "threshold": res.threshold,
        "n_trials": res.n_trials,
        "false_positive_rate": res.false_positive_rate,
        "completeness_snr_50": res.completeness_snr_50,
        "completeness_snr_90": res.completeness_snr_90,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "drift_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _heatmap(res, op / "paper" / "figures")
    return metrics


def _heatmap(res: RecoveryResult, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    im = ax.imshow(
        res.p_detect,
        origin="lower",
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap="viridis",
        extent=[
            res.drift_rates.min(),
            res.drift_rates.max(),
            res.inj_snrs.min(),
            res.inj_snrs.max(),
        ],
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


def validate_voyager(
    path=None, *, window: int = 4096
) -> dict:  # pragma: no cover - net + optional deps
    """Recover the Voyager-1 carrier from the Breakthrough Listen open-data file (external validation).

    Runs the *same* ``jansky.seti`` drift detector used by the benchmark on a real, known narrowband
    signal — the Voyager-1 spacecraft downlink near 8420 MHz. Requires the optional ``voyager`` extra
    (``h5py`` + ``hdf5plugin``; the file is bitshuffle-compressed). Returns the detected tone
    frequency, drift, and S/N versus a blank-window control.
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
    drifts = np.linspace(-3.0, 3.0, 121)
    k = int(np.argmax(wf.mean(axis=0)))
    lo = max(0, k - window // 2)
    carrier = wf[:, lo : lo + window]
    blank = wf[:, 500_000 : 500_000 + window]
    res = seti.drift_search(carrier - np.median(carrier), drifts)
    blank_snr = float(seti.drift_search(blank - np.median(blank), drifts).best_snr)
    return {
        "freq_mhz": fch1 + k * foff,
        "channel": k,
        "best_snr": float(res.best_snr),
        "best_drift": float(res.best_drift),
        "blank_snr": blank_snr,
    }


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
