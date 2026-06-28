"""Sub-threshold radio stacking of a population in VLASS, with injection-recovery calibration.

Most members of an optically/IR-selected population are fainter than a radio survey's single-source
detection limit, but their *average* flux is measurable by **image-plane stacking**: at N known
positions thermal noise averages down as $N^{-1/2}$ while a coherent sub-threshold signal adds, so the
stacked image reveals the population mean (White et al. 2007; Karim et al. 2011). A stacked flux is
only believable once the **bias** is calibrated --- snapshot/CLEAN and residual-vs-restored effects
corrupt raw stacks --- so this module pairs a robust median stack with an **injection-recovery** step
that measures the recovered/injected ratio and de-biases the result.

Reuses the project's verified VLASS CADC-SODA cutout path (the ``radio-cutout`` skill / ``vlass``) and
the ``vlass.measure_image_flux`` peak+annulus pattern. Pure NumPy + a synthetic offline fixture.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "gaussian_psf",
    "injection_recovery",
    "measure_stacked_flux",
    "median_stack",
    "run",
    "synthetic_population",
]


def gaussian_psf(size: int, fwhm_pix: float, amp: float = 1.0) -> np.ndarray:
    """A centred 2-D Gaussian PSF stamp of given FWHM (pixels) and peak amplitude."""
    sigma = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    c = (size - 1) / 2.0
    yy, xx = np.mgrid[0:size, 0:size]
    return amp * np.exp(-((xx - c) ** 2 + (yy - c) ** 2) / (2.0 * sigma**2))


def median_stack(cutouts: np.ndarray, *, sigma: float = 3.0, maxiters: int = 3) -> np.ndarray:
    """Pixel-wise sigma-clipped **median** of N centred cutout stamps (robust to bright interlopers).

    ``cutouts`` is shape ``(N, H, W)``. The sigma-clip down each pixel column rejects the rare bright
    neighbour or artefact (White et al. preferred the median over the mean for exactly this reason),
    and the median averages the thermal noise down by $\\sqrt{N}$.
    """
    from astropy.stats import sigma_clip

    arr = np.asarray(cutouts, float)
    clipped = sigma_clip(arr, sigma=sigma, maxiters=maxiters, axis=0, masked=True)
    return np.ma.median(clipped, axis=0).filled(np.nan)


def measure_stacked_flux(
    stack: np.ndarray, *, search_pix: float = 3.0, annulus_pix: tuple[float, float] = (10.0, 22.0)
) -> dict[str, float]:
    """Central peak, annulus RMS, and SNR of a stacked stamp (mirrors ``vlass.measure_image_flux``)."""
    a = np.asarray(stack, float)
    ny, nx = a.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - cx, yy - cy)
    near = (rr <= search_pix) & np.isfinite(a)
    nan = float("nan")
    peak = float(np.nanmax(a[near])) if near.any() else nan
    ann = a[(rr > annulus_pix[0]) & (rr < annulus_pix[1]) & np.isfinite(a)]
    rms = float(np.std(ann)) if ann.size > 20 else nan
    return {"peak": peak, "rms": rms, "snr": peak / rms if (rms and rms > 0) else nan}


def injection_recovery(
    background: np.ndarray, inject_amp: float, *, fwhm_pix: float = 2.5, sigma: float = 3.0
) -> dict[str, float]:
    """Calibrate the stacking bias: inject a known PSF into each background cutout, stack, measure.

    Injects a ``gaussian_psf`` of peak ``inject_amp`` at the centre of every ``(N, H, W)`` background
    cutout, median-stacks, and measures the recovered central peak. Returns the injected amplitude, the
    recovered peak, and the **ratio** ``recovered/injected`` --- the multiplicative bias to divide a
    measured stacked flux by. (On real VLASS cutouts this absorbs the CLEAN/residual flux bias.)
    """
    bg = np.asarray(background, float)
    psf = gaussian_psf(bg.shape[1], fwhm_pix, inject_amp)
    stack = median_stack(bg + psf[None, :, :], sigma=sigma)
    rec = measure_stacked_flux(stack)["peak"]
    return {
        "injected": float(inject_amp),
        "recovered": float(rec),
        "ratio": float(rec / inject_amp) if inject_amp else float("nan"),
    }


def synthetic_population(
    n_sources: int = 600,
    *,
    source_flux: float = 0.05,
    noise: float = 0.12,
    size: int = 51,
    fwhm_pix: float = 2.5,
    seed: int = 0,
) -> np.ndarray:
    """Synthetic stack of a sub-threshold population: a faint central source + noise per cutout.

    Each of ``n_sources`` cutouts is a centred Gaussian of peak ``source_flux`` (well below the
    per-cutout ``noise``, so individually undetected) plus Gaussian noise. The stack of all N recovers
    ``source_flux`` at high SNR. Returns the ``(N, size, size)`` cube.
    """
    rng = np.random.default_rng(seed)
    psf = gaussian_psf(size, fwhm_pix, source_flux)
    return psf[None, :, :] + rng.normal(0.0, noise, (n_sources, size, size))


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: stack a (synthetic or fetched) population, calibrate with injection-recovery, write."""
    import json
    from pathlib import Path

    if offline:
        cutouts = synthetic_population()
        source = "synthetic"
        injected_truth = 0.05
    else:  # pragma: no cover - network
        raise NotImplementedError("real VLASS-SE fetch + target list is wired in the next step")

    stack = median_stack(cutouts)
    meas = measure_stacked_flux(stack)
    # injection-recovery: inject the measured-level flux into background (noise-only) cutouts
    bg = (
        cutouts - np.nanmedian(cutouts, axis=0)[None, :, :]
    )  # remove the stacked signal -> background
    cal = injection_recovery(bg, meas["peak"])
    debiased = meas["peak"] / cal["ratio"] if cal["ratio"] else float("nan")
    metrics = {
        "source": source,
        "n_stacked": int(cutouts.shape[0]),
        "stacked_peak": round(meas["peak"], 4),
        "stacked_rms": round(meas["rms"], 4),
        "stacked_snr": round(meas["snr"], 1),
        "recovery_ratio": round(cal["ratio"], 3),
        "debiased_flux": round(debiased, 4),
    }
    if offline:
        metrics["injected_truth"] = injected_truth

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "stacking_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(stack, op / "papers" / "stacking" / "figures")
    _write_macros(metrics, op / "papers" / "stacking" / "generated" / "macros.tex")
    return metrics


def _figure(stack: np.ndarray, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4, 3.6))
    im = ax.imshow(np.asarray(stack, float), origin="lower", cmap="inferno")
    fig.colorbar(im, ax=ax, label="stacked brightness")
    ax.set(title="Median-stacked image", xlabel="pixel", ylabel="pixel")
    fig.tight_layout()
    fig.savefig(out / "stack.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.stacking._write_macros — do not edit by hand.",
        rf"\newcommand{{\stSource}}{{{m['source']}}}",
        rf"\newcommand{{\stN}}{{{m['n_stacked']}}}",
        rf"\newcommand{{\stPeak}}{{{m['stacked_peak']}}}",
        rf"\newcommand{{\stSNR}}{{{m['stacked_snr']}}}",
        rf"\newcommand{{\stRatio}}{{{m['recovery_ratio']}}}",
        rf"\newcommand{{\stDebiased}}{{{m['debiased_flux']}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Sub-threshold radio stacking with injection-recovery.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
