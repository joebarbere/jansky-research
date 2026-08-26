"""Sub-threshold radio stacking of a population in VLASS, with forced photometry and controls.

Most members of an optically/IR-selected population are fainter than a radio survey's single-source
detection limit, but their *median* flux is measurable by **image-plane stacking**: at N known
positions thermal noise averages down as $N^{-1/2}$ while a coherent sub-threshold signal adds, so
the stacked image reveals the population's central flux (White et al. 2007; Karim et al. 2011).

Two lessons this module now encodes (round-9 referee):

- **The stacked flux is a FORCED measurement at the known position.** A searched peak within even a
  3-pixel radius on beam-correlated noise reads +1.57 sigma and is positive ~99% of the time --- the
  stokesv lesson. ``measure_stacked_flux`` reads the central pixel; the searched maximum is kept
  only as a labelled diagnostic.
- **A shift-equivariant estimator makes same-plane injection an identity.** Adding one PSF plane to
  every cutout leaves the sigma-clip mask unchanged, so recovered == injected for ANY input and the
  old "ratio" measured nothing. ``injection_recovery`` now injects at random sub-pixel offsets
  (astrometric scatter) at the measured amplitude, so the test can fail; ``run`` additionally
  stacks off-source control cutouts (30 arcsec offsets) through the identical pipeline, which a
  centred annulus cannot substitute for.

Reuses the project's verified VLASS CADC-SODA cutout path (the ``radio-cutout`` skill / ``vlass``).
Pure NumPy + a synthetic offline fixture.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "fetch_population",
    "fetch_se_cutout",
    "gaussian_psf",
    "individually_detected",
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
    """FORCED central-pixel flux, annulus RMS, and SNR of a stacked stamp.

    The stack is made at known positions, so the measurement is the central pixel --- a genuinely
    forced photometry that goes negative about half the time on pure noise. The maximum within
    ``search_pix`` is returned only as ``peak_searched``, a labelled diagnostic: on beam-correlated
    noise it reads +1.57 sigma and is positive ~99% of the time, so it must never be the headline.
    """
    a = np.asarray(stack, float)
    ny, nx = a.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - cx, yy - cy)
    nan = float("nan")
    flux = float(a[int(round(cy)), int(round(cx))])
    if not np.isfinite(flux):
        flux = nan
    near = (rr <= search_pix) & np.isfinite(a)
    peak = float(np.nanmax(a[near])) if near.any() else nan
    ann = a[(rr > annulus_pix[0]) & (rr < annulus_pix[1]) & np.isfinite(a)]
    rms = float(np.std(ann)) if ann.size > 20 else nan
    return {
        "flux": flux,
        "peak_searched": peak,
        "rms": rms,
        "snr": flux / rms if (rms and rms > 0 and np.isfinite(flux)) else nan,
    }


def individually_detected(
    cutouts: np.ndarray,
    *,
    thresh_sigma: float = 5.0,
    search_pix: float = 3.0,
    annulus_pix: tuple[float, float] = (10.0, 22.0),
) -> np.ndarray:
    """Which cutouts hold an individually-detected source at the target position.

    A stack of "individually undetected" sources must actually exclude the detected ones --- a
    sentence in a paper is not a flux cut. Per cutout: the searched maximum within ``search_pix``
    (searching is correct HERE, it is a detection test and conservative) against that cutout's own
    annulus RMS. Returns a boolean mask, True = detected.
    """
    arr = np.asarray(cutouts, float)
    n, ny, nx = arr.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - cx, yy - cy)
    near = rr <= search_pix
    ann_m = (rr > annulus_pix[0]) & (rr < annulus_pix[1])
    out = np.zeros(n, dtype=bool)
    for i in range(n):
        a = arr[i]
        fin_near = near & np.isfinite(a)
        fin_ann = ann_m & np.isfinite(a)
        if not fin_near.any() or fin_ann.sum() < 20:
            continue
        rms = float(np.std(a[fin_ann]))
        if rms > 0 and float(np.nanmax(a[fin_near])) > thresh_sigma * rms:
            out[i] = True
    return out


def injection_recovery(
    background: np.ndarray,
    inject_amp: float,
    *,
    fwhm_pix: float = 2.5,
    sigma: float = 3.0,
    jitter_pix: float = 0.3,
    n_trials: int = 8,
    seed: int = 0,
) -> dict[str, float]:
    """An injection test that can fail: per-cutout sub-pixel offsets, measured amplitude, forced read.

    The earlier version added the SAME PSF plane to every cutout and differenced two clipped
    medians: sigma-clip is shift-equivariant, so the mask never changes and recovered == injected
    identically --- a ratio of 1.0 for any input, including pure noise. This version injects each
    cutout's PSF at an independent random sub-pixel offset (``jitter_pix``, the assumed astrometric
    scatter between catalogue and image in pixels --- a documented assumption, not a fit), stacks,
    and reads the FORCED central pixel above the no-injection baseline. Pixelization and centring
    losses now show up, the sigma-clip can interact, and the result carries a spread over
    ``n_trials`` independent jitter draws. Inject at the measured amplitude, not a comfortable 5x
    the RMS.
    """
    bg = np.asarray(background, float)
    n, size, _ = bg.shape
    c = (size - 1) / 2.0
    grid = np.mgrid[0:size, 0:size]
    yy, xx = grid[0], grid[1]
    s = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    rng = np.random.default_rng(seed)

    def _centre(a: np.ndarray) -> float:
        return float(a[int(round(c)), int(round(c))])

    base = _centre(median_stack(bg, sigma=sigma))
    ratios = []
    for _ in range(n_trials):
        dx = np.asarray(rng.normal(0.0, jitter_pix, n))
        dy = np.asarray(rng.normal(0.0, jitter_pix, n))
        psf = inject_amp * np.exp(
            -(
                (xx[None, :, :] - c - dx[:, None, None]) ** 2
                + (yy[None, :, :] - c - dy[:, None, None]) ** 2
            )
            / (2.0 * s**2)
        )
        rec = _centre(median_stack(bg + psf, sigma=sigma)) - base
        ratios.append(rec / inject_amp if inject_amp else float("nan"))
    r = np.asarray(ratios, float)
    return {
        "injected": float(inject_amp),
        "jitter_pix": float(jitter_pix),
        "ratio": float(np.mean(r)),
        "ratio_sd": float(np.std(r, ddof=1)) if r.size > 1 else float("nan"),
        "n_trials": int(n_trials),
    }


def stack_in_bins(
    cutouts: np.ndarray, values: np.ndarray, *, n_bins: int = 3, min_per_bin: int = 10
) -> list[dict]:
    """Stack the cutouts in ``n_bins`` quantile bins of ``values``, with forced photometry per bin.

    Turns one stacked number into a population *trend*: split the cube into equal-count bins of the
    binning property (e.g. optical magnitude), median-stack each, and return a per-bin dict with
    ``n``, the value range/median, the FORCED central flux, the annulus RMS, and the SNR. (The old
    per-bin "recovery ratio" was an algebraic identity --- always 1.0 --- and is gone; the global
    jittered injection test in :func:`run` covers the estimator once.) Bins with fewer than
    ``min_per_bin`` sources are skipped.
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
        meas = measure_stacked_flux(median_stack(arr[mask]))
        out.append(
            {
                "n": int(mask.sum()),
                "value_lo": float(lo),
                "value_hi": float(hi),
                "value_med": float(np.median(vals[mask])),
                "flux": meas["flux"],
                "rms": meas["rms"],
                "snr": meas["snr"],
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
    flux_scatter_dex: float = 0.0,
    seed: int = 0,
) -> np.ndarray:
    """Synthetic stack of a sub-threshold population: a faint central source + noise per cutout.

    Each of ``n_sources`` cutouts is a centred Gaussian of peak ``source_flux`` (well below the
    per-cutout ``noise``, so individually undetected) plus Gaussian noise. ``flux_scatter_dex`` draws
    each source's flux from a log-normal about ``source_flux`` --- a skewed population, for which the
    clipped-MEDIAN stack recovers the median, not the (larger) mean; with the default 0 every source
    is identical and mean == median, which is exactly the blindness the round-9 referee flagged, so
    tests of the estimator's meaning must set it. Returns the ``(N, size, size)`` cube.
    """
    rng = np.random.default_rng(seed)
    psf = gaussian_psf(size, fwhm_pix, 1.0)
    if flux_scatter_dex > 0:
        fluxes = source_flux * 10.0 ** rng.normal(0.0, flux_scatter_dex, n_sources)
    else:
        fluxes = np.full(n_sources, source_flux)
    return fluxes[:, None, None] * psf[None, :, :] + rng.normal(0.0, noise, (n_sources, size, size))


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
    """Cone-search SDSS DR16 quasars (VizieR ``VII/289``); returns ra, dec, i-band mag, and redshift."""
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "imag", "z"])
    v.ROW_LIMIT = max_sources
    res = v.query_region(center, radius=radius_deg * _u.deg, catalog="VII/289/dr16q")
    t = res[0]
    return (
        _np.asarray(t["RAJ2000"], float),
        _np.asarray(t["DEJ2000"], float),
        _np.asarray(t["imag"], float),
        _np.asarray(t["z"], float),
    )


CONTROL_OFFSET_ARCSEC = 30.0  # off-source control positions: this far north of each target


def run(
    center=None,
    radius_deg: float = 3.0,
    out: str = ".",
    *,
    offline: bool = True,
    max_sources: int = 300,
) -> dict:
    """Full slice: stack a (synthetic or fetched) population with forced photometry and controls."""
    from pathlib import Path

    mags: np.ndarray
    redshifts: np.ndarray
    targets: list[dict] | None = None
    if offline or center is None:
        cutouts = synthetic_population()
        rng = np.random.default_rng(0)
        mags = np.asarray(rng.uniform(18.0, 21.0, cutouts.shape[0]))  # i-mag
        # a synthetic redshift trend (brighter radio at higher z) for the offline binning test
        redshifts = np.asarray(np.random.default_rng(1).uniform(0.5, 3.0, cutouts.shape[0]))
        source = "synthetic"
        injected_truth: float | None = 0.05
        n_queried = int(cutouts.shape[0])
        # the offline control: noise-only cutouts through the identical pipeline
        controls: np.ndarray | None = np.random.default_rng(2).normal(0.0, 0.12, cutouts.shape)
    else:  # pragma: no cover - network
        ra, dec, imag, zarr = fetch_population(center, radius_deg, max_sources=max_sources)
        n_queried = int(ra.size)
        rows = []
        for r, d, m, z in zip(ra, dec, imag, zarr, strict=True):
            c = fetch_se_cutout(float(r), float(d))
            ctl = (
                fetch_se_cutout(float(r), float(d) + CONTROL_OFFSET_ARCSEC / 3600.0)
                if c is not None
                else None
            )
            rows.append(
                {
                    "ra": float(r),
                    "dec": float(d),
                    "imag": float(m),
                    "z": float(z),
                    "cutout": c,
                    "control": ctl,
                }
            )
        got = [x for x in rows if x["cutout"] is not None]
        if len(got) < 20:
            raise RuntimeError(f"only {len(got)} VLASS-SE cutouts fetched; need more for a stack")
        cutouts = np.asarray([x["cutout"] for x in got])
        mags = np.asarray([x["imag"] for x in got])
        redshifts = np.asarray([x["z"] for x in got])
        ctl_list = [x["control"] for x in got if x["control"] is not None]
        controls = np.asarray(ctl_list) if len(ctl_list) >= 20 else None
        source = (
            f"SDSS DR16Q x VLASS-SE @ ({center.ra.deg:.1f},{center.dec.deg:.1f}) "
            f"r={radius_deg:g} deg"
        )
        injected_truth = None
        targets = rows

    # the stated sample cut, implemented: drop individually-detected sources before stacking
    det = individually_detected(cutouts)
    n_detected = int(det.sum())
    cutouts, mags, redshifts = cutouts[~det], mags[~det], redshifts[~det]

    stack = median_stack(cutouts)
    meas = measure_stacked_flux(stack)
    # the injection test that can fail: sub-pixel jitter, at the MEASURED amplitude
    inject_amp = meas["flux"] if (np.isfinite(meas["flux"]) and meas["flux"] > 0) else meas["rms"]
    cal = injection_recovery(cutouts, inject_amp)
    # the off-source control: identical pipeline at positions holding no source
    control: dict | None = None
    if controls is not None:
        cm = measure_stacked_flux(median_stack(controls))
        control = {
            "n": int(controls.shape[0]),
            "flux": round(cm["flux"], 4),
            "rms": round(cm["rms"], 4),
            "snr": round(cm["snr"], 1) if np.isfinite(cm["snr"]) else None,
        }
    mag_bins = sorted(stack_in_bins(cutouts, mags, n_bins=3), key=lambda b: b["value_med"])
    binned: list[dict] = [
        {
            "imag_med": round(b["value_med"], 2),
            "n": b["n"],
            "flux_uJy": round(1e3 * b["flux"], 1),
            "rms_uJy": round(1e3 * b["rms"], 1),
            "snr": round(b["snr"], 1),
        }
        for b in mag_bins
    ]
    z_bins = sorted(stack_in_bins(cutouts, redshifts, n_bins=3), key=lambda b: b["value_med"])
    binned_z: list[dict] = [
        {
            "z_med": round(b["value_med"], 3),
            "n": b["n"],
            "flux_uJy": round(1e3 * b["flux"], 1),
            "rms_uJy": round(1e3 * b["rms"], 1),
            "snr": round(b["snr"], 1),
        }
        for b in z_bins
    ]
    metrics = {
        "source": source,
        "n_queried": n_queried,
        "radius_deg": radius_deg,
        "max_sources": max_sources,
        "n_with_cutout": int(det.size),
        "n_detected_excluded": n_detected,
        "n_stacked": int(cutouts.shape[0]),
        "stacked_flux": round(meas["flux"], 4),
        "stacked_peak_searched": round(meas["peak_searched"], 4),
        "stacked_rms": round(meas["rms"], 4),
        "stacked_snr": round(meas["snr"], 1),
        "injection": {
            "amp": round(float(cal["injected"]), 4),
            "jitter_pix": cal["jitter_pix"],
            "ratio": round(cal["ratio"], 3),
            "ratio_sd": round(cal["ratio_sd"], 3),
            "n_trials": cal["n_trials"],
        },
        "control": control,
        "n_bins": len(mag_bins),
        "bins": binned,
        "n_zbins": len(z_bins),
        "zbins": binned_z,
    }
    if center is not None:  # pragma: no cover - network
        metrics["field_ra"] = round(float(center.ra.deg), 2)
        metrics["field_dec"] = round(float(center.dec.deg), 2)
    if injected_truth is not None:
        metrics["injected_truth"] = injected_truth

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import _results_are_real, write_results

    json_path = op / "results" / "stacking_metrics.json"
    write_artifacts = True
    if source == "synthetic":
        try:
            import json as _json

            write_artifacts = not (
                json_path.is_file() and _results_are_real(_json.loads(json_path.read_text()))
            )
        except Exception:
            write_artifacts = True

    write_results(metrics, json_path)
    if targets is not None:  # pragma: no cover - network
        _write_targets(op / "results" / "stacking_targets.csv", targets, det)
    if write_artifacts:
        _figure(stack, binned, binned_z, op / "papers" / "stacking" / "figures")
    _write_macros(metrics, op / "papers" / "stacking" / "generated" / "macros.tex")
    return metrics


def _write_targets(path, rows: list[dict], det: np.ndarray) -> None:  # pragma: no cover - network
    """Commit the full queried target list: the denominator the stack is drawn from."""
    import csv
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # det is aligned to the rows that HAVE a cutout, in order
    det_iter = iter(det)
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ra", "dec", "imag", "z", "has_cutout", "has_control", "detected_excluded"])
        for x in rows:
            has = x["cutout"] is not None
            d = bool(next(det_iter)) if has else ""
            w.writerow(
                [
                    f"{x['ra']:.5f}",
                    f"{x['dec']:.5f}",
                    f"{x['imag']:.2f}",
                    f"{x['z']:.3f}",
                    int(has),
                    int(x["control"] is not None),
                    int(d) if d != "" else "",
                ]
            )


def _figure(stack: np.ndarray, bins: list[dict], zbins: list[dict], out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 3.5))
    im = ax1.imshow(np.asarray(stack, float), origin="lower", cmap="inferno")
    fig.colorbar(im, ax=ax1, label="mJy/beam")
    ax1.set(title="Median-stacked image", xlabel="pixel", ylabel="pixel")
    if bins:
        mag = [b["imag_med"] for b in bins]
        flux = [b["flux_uJy"] for b in bins]
        err = [b["rms_uJy"] for b in bins]
        ax2.errorbar(mag, flux, yerr=err, fmt="o-", color="C0")
        ax2.set(
            xlabel=r"median $i$ magnitude",
            ylabel=r"median radio flux ($\mu$Jy/beam)",
            title="Flux vs. apparent magnitude",
        )
        ax2.invert_xaxis()  # brighter (smaller mag) to the right
    if zbins:
        z = [b["z_med"] for b in zbins]
        fluxz = [b["flux_uJy"] for b in zbins]
        errz = [b["rms_uJy"] for b in zbins]
        ax3.errorbar(z, fluxz, yerr=errz, fmt="s-", color="C1")
        ax3.set(
            xlabel=r"median redshift $z$",
            ylabel=r"median radio flux ($\mu$Jy/beam)",
            title="Flux vs. redshift",
        )
    fig.tight_layout()
    fig.savefig(out / "stack.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    inj = m.get("injection") or {}
    ctl = m.get("control") or {}

    def _g(dic: dict, key: str) -> str:
        val = dic.get(key)
        return "--" if val is None else str(val)

    def _mm(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    flux_ujy = round(1e3 * m["stacked_flux"], 1)
    rms_ujy = round(1e3 * m["stacked_rms"], 1)
    srch_ujy = round(1e3 * m["stacked_peak_searched"], 1)
    ctl_flux_ujy = "--" if ctl.get("flux") is None else str(round(1e3 * ctl["flux"], 1))
    ctl_rms_ujy = "--" if ctl.get("rms") is None else str(round(1e3 * ctl["rms"], 1))
    lines = [
        "% Auto-generated by jansky_research.stacking._write_macros — do not edit by hand.",
        rf"\newcommand{{\stSource}}{{{m['source']}}}",
        rf"\newcommand{{\stN}}{{{m['n_stacked']}}}",
        rf"\newcommand{{\stNqueried}}{{{_mm('n_queried')}}}",
        rf"\newcommand{{\stNwithCutout}}{{{_mm('n_with_cutout')}}}",
        rf"\newcommand{{\stNdetExcl}}{{{_mm('n_detected_excluded')}}}",
        rf"\newcommand{{\stRadius}}{{{_mm('radius_deg')}}}",
        rf"\newcommand{{\stMaxSources}}{{{_mm('max_sources')}}}",
        rf"\newcommand{{\stFieldRa}}{{{_mm('field_ra')}}}",
        rf"\newcommand{{\stFieldDec}}{{{_mm('field_dec')}}}",
        rf"\newcommand{{\stFlux}}{{{flux_ujy}}}",
        rf"\newcommand{{\stRms}}{{{rms_ujy}}}",
        rf"\newcommand{{\stSNR}}{{{m['stacked_snr']}}}",
        rf"\newcommand{{\stPeakSearched}}{{{srch_ujy}}}",
        rf"\newcommand{{\stInjRatio}}{{{_g(inj, 'ratio')}}}",
        rf"\newcommand{{\stInjRatioSD}}{{{_g(inj, 'ratio_sd')}}}",
        rf"\newcommand{{\stInjJitter}}{{{_g(inj, 'jitter_pix')}}}",
        rf"\newcommand{{\stCtlN}}{{{_g(ctl, 'n')}}}",
        rf"\newcommand{{\stCtlFlux}}{{{ctl_flux_ujy}}}",
        rf"\newcommand{{\stCtlRms}}{{{ctl_rms_ujy}}}",
        rf"\newcommand{{\stCtlSNR}}{{{_g(ctl, 'snr')}}}",
        rf"\newcommand{{\stNbins}}{{{m.get('n_bins', 0)}}}",
    ]
    bins = m.get("bins", [])
    if bins:
        bright, faint = bins[0], bins[-1]  # bins sorted by median i-mag (brightest first)
        lines += [
            rf"\newcommand{{\stBrightMag}}{{{bright['imag_med']}}}",
            rf"\newcommand{{\stBrightFlux}}{{{bright['flux_uJy']}}}",
            rf"\newcommand{{\stBrightSNR}}{{{bright['snr']}}}",
            rf"\newcommand{{\stFaintMag}}{{{faint['imag_med']}}}",
            rf"\newcommand{{\stFaintFlux}}{{{faint['flux_uJy']}}}",
            rf"\newcommand{{\stFaintSNR}}{{{faint['snr']}}}",
        ]
    zbins = m.get("zbins", [])
    if zbins:
        lowz, highz = zbins[0], zbins[-1]  # zbins sorted by median z (lowest first)
        # the brightest z-bin and the flux range, so the paper can describe a non-monotonic trend honestly
        peakz = max(zbins, key=lambda b: b["flux_uJy"])
        lines += [
            rf"\newcommand{{\stNzbins}}{{{m.get('n_zbins', 0)}}}",
            rf"\newcommand{{\stLowzZ}}{{{lowz['z_med']}}}",
            rf"\newcommand{{\stLowzFlux}}{{{lowz['flux_uJy']}}}",
            rf"\newcommand{{\stHighzZ}}{{{highz['z_med']}}}",
            rf"\newcommand{{\stHighzFlux}}{{{highz['flux_uJy']}}}",
            rf"\newcommand{{\stPeakzZ}}{{{peakz['z_med']}}}",
            rf"\newcommand{{\stPeakzFlux}}{{{peakz['flux_uJy']}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a run may only ADD information, so an
    # offline rebuild can never blank a real value (report.preserve_live_macros).
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


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
