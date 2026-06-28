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
    "fetch_population",
    "fetch_se_cutout",
    "gaussian_psf",
    "injection_recovery",
    "measure_stacked_flux",
    "median_stack",
    "run",
    "stack_in_bins",
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
    cutout, median-stacks, and measures the recovered central peak **above the no-injection baseline**
    (so it works on real cutouts that already hold the faint population signal). Returns the injected
    amplitude, the recovered excess, and the **ratio** ``recovered/injected`` --- the multiplicative
    bias to divide a measured stacked flux by. (On real VLASS cutouts this absorbs the flux-scale bias.)
    """
    bg = np.asarray(background, float)

    # measure at the exact centre (the injected source is centred), so the common noise cancels
    def _centre(a: np.ndarray) -> float:
        ny, nx = a.shape
        return float(a[ny // 2, nx // 2])

    psf = gaussian_psf(bg.shape[1], fwhm_pix, inject_amp)
    base = _centre(median_stack(bg, sigma=sigma))
    rec = _centre(median_stack(bg + psf[None, :, :], sigma=sigma)) - base
    return {
        "injected": float(inject_amp),
        "recovered": float(rec),
        "ratio": float(rec / inject_amp) if inject_amp else float("nan"),
    }


def stack_in_bins(
    cutouts: np.ndarray, values: np.ndarray, *, n_bins: int = 3, min_per_bin: int = 10
) -> list[dict]:
    """Stack the cutouts in ``n_bins`` quantile bins of ``values``, injection-recovering each bin.

    Turns one stacked number into a population *trend*: split the cube into equal-count bins of the
    binning property (e.g. optical magnitude), median-stack and injection-recover each, and return a
    per-bin dict with ``n``, the value range/median, the stacked peak/SNR, the recovery ratio, and the
    de-biased flux. Bins with fewer than ``min_per_bin`` sources are skipped.
    """
    arr = np.asarray(cutouts, float)
    vals = np.asarray(values, float)
    good = np.isfinite(vals)
    arr, vals = arr[good], vals[good]
    edges = np.quantile(vals, np.linspace(0.0, 1.0, n_bins + 1))
    out: list[dict] = []
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        mask = (vals >= lo) & (vals <= hi) if b == n_bins - 1 else (vals >= lo) & (vals < hi)
        if int(mask.sum()) < min_per_bin:
            continue
        sub = arr[mask]
        meas = measure_stacked_flux(median_stack(sub))
        amp = (
            5.0 * meas["rms"] if (meas["rms"] and meas["rms"] > 0) else 5.0 * float(np.nanstd(sub))
        )
        cal = injection_recovery(sub, amp)
        out.append(
            {
                "n": int(mask.sum()),
                "value_lo": float(lo),
                "value_hi": float(hi),
                "value_med": float(np.median(vals[mask])),
                "peak": meas["peak"],
                "snr": meas["snr"],
                "ratio": cal["ratio"],
                "debiased": meas["peak"] / cal["ratio"] if cal["ratio"] else float("nan"),
            }
        )
    return out


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


def fetch_se_cutout(
    ra: float, dec: float, *, size_pix: int = 51, search_deg: float = 0.006
) -> np.ndarray | None:  # pragma: no cover - network
    """One VLASS Single-Epoch Stokes-I cutout (mJy/beam) at ``(ra, dec)`` via CADC SODA, or None.

    ``get_image_list`` returns server-side **cutout** URLs; we pick the SE Stokes-I ``tt0`` product and
    download it (a small stamp), then trim to a fixed ``size_pix`` square centred on the source so all
    cutouts stack on a common grid. None if there is no SE image or the download fails.
    """
    import io

    import numpy as _np
    import requests
    from astropy import units as _u
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS
    from astroquery.cadc import Cadc

    pos = SkyCoord(ra, dec, unit="deg")
    rad = search_deg * _u.deg
    try:
        cadc = Cadc()
        urls = cadc.get_image_list(cadc.query_region(pos, radius=rad, collection="VLASS"), pos, rad)
        se = [u for u in urls if ".se." in u and ".I." in u and "tt0" in u]
        if not se:
            return None
        data = requests.get(se[0], timeout=120).content
        if b"SIMPLE" not in data[:80]:
            return None
        with fits.open(io.BytesIO(data)) as hd:
            img = _np.squeeze(_np.asarray(hd[0].data, float)) * 1e3  # Jy/beam -> mJy/beam
            w = WCS(hd[0].header).celestial
        cut = Cutout2D(img, pos, (size_pix, size_pix), wcs=w, mode="partial", fill_value=_np.nan)
        return cut.data if cut.data.shape == (size_pix, size_pix) else None
    except Exception:
        return None


def fetch_population(
    center, radius_deg: float, *, max_sources: int = 300
) -> tuple:  # pragma: no cover - network
    """Cone-search SDSS DR16 quasars (VizieR ``VII/289``); returns ra, dec, and i-band magnitude."""
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "imag"])
    v.ROW_LIMIT = max_sources
    res = v.query_region(center, radius=radius_deg * _u.deg, catalog="VII/289/dr16q")
    t = res[0]
    return (
        _np.asarray(t["RAJ2000"], float),
        _np.asarray(t["DEJ2000"], float),
        _np.asarray(t["imag"], float),
    )


def run(
    center=None,
    radius_deg: float = 3.0,
    out: str = ".",
    *,
    offline: bool = True,
    max_sources: int = 300,
) -> dict:
    """Full slice: stack a (synthetic or fetched) population, calibrate with injection-recovery, write."""
    import json
    from pathlib import Path

    values: np.ndarray
    if offline or center is None:
        cutouts = synthetic_population()
        values = np.asarray(np.random.default_rng(0).uniform(18.0, 21.0, cutouts.shape[0]))  # i-mag
        source = "synthetic"
        injected_truth: float | None = 0.05
    else:  # pragma: no cover - network
        ra, dec, imag = fetch_population(center, radius_deg, max_sources=max_sources)
        pairs = [
            (c, m)
            for c, m in (
                (fetch_se_cutout(float(r), float(d)), float(m))
                for r, d, m in zip(ra, dec, imag, strict=True)
            )
            if c is not None
        ]
        if len(pairs) < 20:
            raise RuntimeError(f"only {len(pairs)} VLASS-SE cutouts fetched; need more for a stack")
        cutouts = np.asarray([c for c, _ in pairs])
        values = np.asarray([m for _, m in pairs])
        source = f"SDSS DR16Q x VLASS-SE @ ({center.ra.deg:.1f},{center.dec.deg:.1f})"
        injected_truth = None

    stack = median_stack(cutouts)
    meas = measure_stacked_flux(stack)
    # injection-recovery on the actual cutouts (baseline-subtracted): inject at a clean detectable level
    inject_amp = (
        5.0 * meas["rms"] if (meas["rms"] and meas["rms"] > 0) else 5.0 * float(np.nanstd(cutouts))
    )
    cal = injection_recovery(cutouts, inject_amp)
    debiased = meas["peak"] / cal["ratio"] if cal["ratio"] else float("nan")
    # magnitude-binned trend: turn one number into the radio-optical luminosity relation
    bins = sorted(stack_in_bins(cutouts, values, n_bins=3), key=lambda b: b["value_med"])
    binned: list[dict] = [
        {
            "imag_med": round(b["value_med"], 2),
            "n": b["n"],
            "debiased_uJy": round(1e3 * b["debiased"], 1),
            "snr": round(b["snr"], 1),
        }
        for b in bins
    ]
    metrics = {
        "source": source,
        "n_stacked": int(cutouts.shape[0]),
        "stacked_peak": round(meas["peak"], 4),
        "stacked_rms": round(meas["rms"], 4),
        "stacked_snr": round(meas["snr"], 1),
        "recovery_ratio": round(cal["ratio"], 3),
        "debiased_flux": round(debiased, 4),
        "n_bins": len(bins),
        "bins": binned,
    }
    if injected_truth is not None:
        metrics["injected_truth"] = injected_truth

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "stacking_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(stack, binned, op / "papers" / "stacking" / "figures")
    _write_macros(metrics, op / "papers" / "stacking" / "generated" / "macros.tex")
    return metrics


def _figure(stack: np.ndarray, bins: list[dict], out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.6))
    im = ax1.imshow(np.asarray(stack, float), origin="lower", cmap="inferno")
    fig.colorbar(im, ax=ax1, label="mJy/beam")
    ax1.set(title="Median-stacked image", xlabel="pixel", ylabel="pixel")
    if bins:
        mag = [b["imag_med"] for b in bins]
        flux = [b["debiased_uJy"] for b in bins]
        ax2.plot(mag, flux, "o-", color="C0")
        ax2.set(
            xlabel=r"median $i$ magnitude",
            ylabel=r"mean radio flux ($\mu$Jy/beam)",
            title="Radio--optical trend",
        )
        ax2.invert_xaxis()  # brighter (smaller mag) to the right
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
        rf"\newcommand{{\stNbins}}{{{m.get('n_bins', 0)}}}",
    ]
    bins = m.get("bins", [])
    if bins:
        bright, faint = bins[0], bins[-1]  # bins sorted by median i-mag (brightest first)
        lines += [
            rf"\newcommand{{\stBrightMag}}{{{bright['imag_med']}}}",
            rf"\newcommand{{\stBrightFlux}}{{{bright['debiased_uJy']}}}",
            rf"\newcommand{{\stFaintMag}}{{{faint['imag_med']}}}",
            rf"\newcommand{{\stFaintFlux}}{{{faint['debiased_uJy']}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    from astropy.coordinates import SkyCoord

    p = argparse.ArgumentParser(description="Sub-threshold radio stacking with injection-recovery.")
    p.add_argument("--ra", type=float, help="field-centre RA (deg)")
    p.add_argument("--dec", type=float, help="field-centre Dec (deg)")
    p.add_argument("--radius", type=float, default=3.0, help="cone radius (deg)")
    p.add_argument("--max-sources", type=int, default=300)
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    metrics = run(center, args.radius, args.out, offline=args.offline, max_sources=args.max_sources)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
