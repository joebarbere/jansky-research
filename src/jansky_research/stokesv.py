"""Stokes-V (circular-polarization) selection of coherent radio emitters in RACS.

Circularly polarized radio emission is a near-unambiguous flag of a *coherent* process
(electron-cyclotron maser, coherent plasma emission). Extragalactic AGN are at most a fraction of a
percent circularly polarized, whereas flaring M dwarfs, ultracool/brown dwarfs, magnetic
chemically-peculiar stars, RS\\,CVn binaries and pulsars reach tens of percent. So a high
$|V|/I$ in a wide-field survey is a clean coherent-emitter finder (Pritchard et al. 2021).

The dominant false positive is instrumental: off-axis **Stokes-I$\\to$V leakage** puts spurious
circular polarization at the $\\sim$1\\% level on bright sources, rising to $\\sim$10\\% near beam
edges. We defeat it with a per-region **leakage floor** estimated from the (assumed unpolarised) field
population --- a candidate must clear $n_\\sigma\\times\\mathrm{median}(|V/I|)$, not merely the image
noise (the RACS convention, $n_\\sigma=7$). Genuine stellar association is then confirmed by **proper
motion**: a real radio star sits at the catalogue position *propagated to the radio epoch*, not at the
optical-epoch position --- a discriminant a chance background match fails.

This composes :mod:`jansky_research.spectra` (cross-match, two-point index) and the
:mod:`jansky_research.vlass` forced-photometry/vetting pattern rather than reimplementing them. Pure
NumPy + a synthetic offline fixture for tests; the real run does forced photometry of a curated
late-type-star / ultracool-dwarf target list in CASDA RACS Stokes-V cutouts.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "LEAKAGE_NSIGMA",
    "RACS_NU_GHZ",
    "classify_emitter",
    "fetch_racs_i",
    "fetch_radio_star_measurements",
    "fetch_radio_stars",
    "fractional_circular_pol",
    "handedness",
    "leakage_floor",
    "match_targets_to_radio",
    "proper_motion_confirm",
    "run",
    "select_circular_pol",
    "synthetic_field",
    "validate_srsc",
]

# RACS band reference frequencies (GHz): RACS-low, RACS-mid, RACS-high.
RACS_NU_GHZ = {"low": 0.8875, "mid": 1.3675, "high": 1.6555}
# Candidate threshold in units of the median field |V/I| (the RACS-low2 / Pritchard convention).
LEAKAGE_NSIGMA = 7.0
# VizieR catalogue IDs. RACS-low DR1 Stokes-I (Hale+2021, split into the extragalactic |b| cut and the
# galactic region) and the Sydney Radio Star Catalogue (Driessen+2024) for the recover-a-known.
VIZIER_RACS_LOW_I = ("J/other/PASA/38.58/galcut", "J/other/PASA/38.58/galreg")
VIZIER_SRSC = "J/other/PASA/41.84/stars"


def fractional_circular_pol(
    v_flux: np.ndarray,
    i_flux: np.ndarray,
    e_v: np.ndarray | None = None,
    e_i: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Fractional circular polarization $|V|/I$ and its 1-$\sigma$ error.

    $f=|V|/I$. The error propagates both flux errors without dividing by the (possibly zero or
    sign-flipping) $V$: $\sigma_f=\sqrt{(e_V/I)^2+(f\,e_I/I)^2}$.
    """
    v = np.asarray(v_flux, float)
    i = np.asarray(i_flux, float)
    frac = np.abs(v) / i
    if e_v is None or e_i is None:
        return frac, np.full_like(frac, np.nan)
    e_v = np.asarray(e_v, float)
    e_i = np.asarray(e_i, float)
    sigma = np.sqrt((e_v / i) ** 2 + (frac * e_i / i) ** 2)
    return frac, sigma


def handedness(v_flux: float) -> str:
    """Sense of circular polarization from the sign of $V$ (``RCP`` if $V>0$, else ``LCP``).

    Note: the absolute $V$ sign convention varies between ASKAP/RACS pipeline versions and epochs, so
    handedness should be cross-checked against the relevant survey paper before physical interpretation.
    """
    return "RCP" if v_flux > 0 else "LCP"


def leakage_floor(frac_pol: np.ndarray, *, n_sigma: float = LEAKAGE_NSIGMA) -> float:
    r"""Per-region leakage threshold from the field $|V/I|$ distribution.

    The field is dominated by unpolarised sources whose apparent circular polarization is
    instrumental Stokes-I$\to$V leakage, so $\mathrm{median}(|V/I|)$ estimates the local leakage
    level (robust to the rare genuine emitters). The credible candidate threshold is
    ``n_sigma`` $\times$ that median --- a candidate must clear *leakage*, not just image noise. Feed
    it the $|V/I|$ of **bright** ($I/\sigma_I\gg1$) sources only: on faint sources $|V/I|\approx
    \sigma_V/I$ is noise-, not leakage-, dominated and would inflate the floor (handle that noise with
    a separate $V$-SNR cut, see :func:`select_circular_pol`).
    """
    f = np.asarray(frac_pol, float)
    f = f[np.isfinite(f)]
    if f.size == 0:
        return float("nan")
    return float(n_sigma * np.median(f))


def select_circular_pol(
    i_flux: np.ndarray,
    v_flux: np.ndarray,
    e_i: np.ndarray,
    e_v: np.ndarray,
    *,
    leakage_threshold: float,
    v_snr_min: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Select circular-polarization candidates above the leakage floor *and* the $V$ SNR cut.

    A source is a candidate when $|V|/I>$ ``leakage_threshold`` (clears instrumental leakage) **and**
    $|V|/\sigma_V\ge$ ``v_snr_min`` (a real $V$ detection, not noise). Forced photometry of a known
    target can use a looser ``v_snr_min`` (3--4) than a blind 5$\sigma$ extraction. Returns
    ``(mask, frac_pol)``.
    """
    frac, _ = fractional_circular_pol(v_flux, i_flux, e_v, e_i)
    v_snr = np.abs(np.asarray(v_flux, float)) / np.asarray(e_v, float)
    mask = (frac > leakage_threshold) & (v_snr >= v_snr_min)
    return mask, frac


def proper_motion_confirm(
    ra_radio: np.ndarray,
    dec_radio: np.ndarray,
    ra_cat: np.ndarray,
    dec_cat: np.ndarray,
    pmra: np.ndarray,
    pmdec: np.ndarray,
    dt_yr: float,
    *,
    match_arcsec: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Confirm a stellar association by proper motion.

    Propagates each star's catalogue position to the radio epoch --- ``pmra`` is $\mu_{\alpha*}=
    \mu_\alpha\cos\delta$ (mas/yr), so the RA *coordinate* shift is $\mu_{\alpha*}\,\Delta t/\cos\delta$
    --- and compares it to the measured radio centroid. A genuine radio star matches the *propagated*
    position; a chance background match (or a static source) does not, for an appreciable proper
    motion. Returns ``(confirmed_mask, sep_arcsec)`` where ``sep`` is to the propagated position.
    """
    ra_cat = np.asarray(ra_cat, float)
    dec_cat = np.asarray(dec_cat, float)
    cosd = np.cos(np.radians(dec_cat))
    ra_prop = ra_cat + (np.asarray(pmra, float) * dt_yr / 1000.0 / 3600.0) / cosd
    dec_prop = dec_cat + np.asarray(pmdec, float) * dt_yr / 1000.0 / 3600.0
    dra = (np.asarray(ra_radio, float) - ra_prop) * cosd
    ddec = np.asarray(dec_radio, float) - dec_prop
    sep = np.sqrt(dra**2 + ddec**2) * 3600.0
    return sep <= match_arcsec, sep


def match_targets_to_radio(
    target_ra: np.ndarray,
    target_dec: np.ndarray,
    radio_ra: np.ndarray,
    radio_dec: np.ndarray,
    radio_i: np.ndarray,
    radio_ei: np.ndarray,
    *,
    radius_arcsec: float = 15.0,
) -> dict[str, np.ndarray]:
    """Attach each target star's nearest radio Stokes-I component (reuses ``spectra.crossmatch``).

    Returns per-target arrays: ``matched`` (a radio component within ``radius_arcsec``), the matched
    ``i_flux``/``e_i`` (NaN where unmatched), and ``sep_arcsec``. A target with no match is a radio
    non-detection in this catalogue --- where forced photometry on the images then sets an upper limit.
    """
    from .spectra import crossmatch

    n = np.asarray(target_ra).size
    out = {
        "matched": np.zeros(n, dtype=bool),
        "i_flux": np.full(n, np.nan),
        "e_i": np.full(n, np.nan),
        "sep_arcsec": np.full(n, np.nan),
    }
    if np.asarray(radio_ra).size == 0:
        return out
    it, ir, sep = crossmatch(target_ra, target_dec, radio_ra, radio_dec, radius_arcsec)
    out["matched"][it] = True
    out["i_flux"][it] = np.asarray(radio_i, float)[ir]
    out["e_i"][it] = np.asarray(radio_ei, float)[ir]
    out["sep_arcsec"][it] = sep
    return out


def fetch_radio_stars() -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Fetch the Sydney Radio Star Catalogue (Driessen+2024) from VizieR (the recover-a-known target list).

    Returns positions, Gaia proper motions, the detection ``method`` (which flags Stokes-V detections),
    and the SIMBAD identifier per star.
    """
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "pmRA", "pmDE", "Method", "Survey", "Simbad"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs(VIZIER_SRSC)[0]
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "pmra": _np.asarray(t["pmRA"], float),
        "pmdec": _np.asarray(t["pmDE"], float),
        "method": _np.asarray([str(x) for x in t["Method"]], dtype=object),
        "survey": _np.asarray([str(x) for x in t["Survey"]], dtype=object),
        "simbad": _np.asarray([str(x) for x in t["Simbad"]], dtype=object),
    }


def fetch_racs_i(center, radius_deg: float) -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Cone-search RACS-low DR1 Stokes-I (Hale+2021) on VizieR; returns ra/dec/i_flux/e_i/noise (mJy)."""
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "Fpk", "e_Fpk", "Noise"])
    v.ROW_LIMIT = -1
    res = v.query_region(center, radius=radius_deg * _u.deg, catalog=list(VIZIER_RACS_LOW_I))
    if not res:
        return {k: _np.array([]) for k in ("ra", "dec", "i_flux", "e_i", "noise")}
    ra, dec, fpk, efpk, noise = [], [], [], [], []
    for t in res:
        ra.append(_np.asarray(t["RAJ2000"], float))
        dec.append(_np.asarray(t["DEJ2000"], float))
        fpk.append(_np.asarray(t["Fpk"], float))
        efpk.append(_np.asarray(t["e_Fpk"], float))
        noise.append(_np.asarray(t["Noise"], float))
    return {
        "ra": _np.concatenate(ra),
        "dec": _np.concatenate(dec),
        "i_flux": _np.concatenate(fpk),
        "e_i": _np.concatenate(efpk),
        "noise": _np.concatenate(noise),
    }


def fetch_radio_star_measurements() -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Fetch the SRSC per-detection radio table (Driessen+2024): Stokes I & V peak fluxes per survey.

    Returns ra/dec, ``freq`` (MHz), ``survey``, ``i_flux``/``e_i`` (SpeakI), ``v_flux``/``e_v``
    (SpeakV), and ``rms_v`` (localrmsV) --- the real coherent-emitter measurements used by
    :func:`validate_srsc`.
    """
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    cats = v.get_catalogs("J/other/PASA/41.84")
    t = [c for c in cats if c.meta["name"].endswith("/radio")][0]
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "freq": _np.asarray(t["Freq"], float),
        "survey": _np.asarray([str(x) for x in t["Survey"]], dtype=object),
        "i_flux": _np.asarray(t["SpeakI"], float),
        "e_i": _np.asarray(t["e_SpeakI"], float),
        "v_flux": _np.asarray(t["SpeakV"], float),
        "e_v": _np.asarray(t["e_SpeakV"], float),
        "rms_v": _np.asarray(t["localrmsV"], float),
    }


def validate_srsc(*, survey_prefix: str = "RACS") -> dict:  # pragma: no cover - network
    """Recover-a-known: do the SRSC's known V-detected radio stars classify as circular emitters?

    Pulls the SRSC radio measurements (:func:`fetch_radio_star_measurements`), keeps detections from
    surveys whose name starts with ``survey_prefix`` (default the three RACS bands) that have both a
    Stokes I and a Stokes V peak flux, and runs :func:`fractional_circular_pol` /
    :func:`classify_emitter` on them. Coherent emitters should be strongly circularly polarized, so
    the fraction classified ``circular``/``highly_circular`` is the recovery. (This validates the
    selection logic on real data; the leakage-floor + forced-photometry path is exercised separately.)
    """
    import numpy as _np

    m = fetch_radio_star_measurements()
    sel = (
        _np.asarray([str(s).startswith(survey_prefix) for s in m["survey"]])
        & _np.isfinite(m["i_flux"])
        & (m["i_flux"] > 0)
        & _np.isfinite(m["v_flux"])
    )
    frac, _ = fractional_circular_pol(m["v_flux"][sel], m["i_flux"][sel])
    cls = [
        classify_emitter(vv, ii) for vv, ii in zip(m["v_flux"][sel], m["i_flux"][sel], strict=True)
    ]
    n = int(sel.sum())
    circ = sum(c in ("circular", "highly_circular") for c in cls)
    return {
        "n_detections": n,
        "median_frac_pol": float(_np.median(frac)) if n else float("nan"),
        "frac_circular": circ / n if n else 0.0,
        "n_highly_circular": sum(c == "highly_circular" for c in cls),
    }


def classify_emitter(
    v_flux: float, i_flux: float, *, strong: float = 0.3, weak: float = 0.06
) -> str:
    """Coarse circular-polarization class from $|V|/I$.

    ``highly_circular`` ($\\ge$ ``strong``: deep ECME-like polarization, a strong coherent-emitter
    candidate); ``circular`` (``weak``--``strong``); ``weak`` (below ``weak``, leakage-dominated
    regime); ``nan`` for non-finite input.
    """
    if not (np.isfinite(v_flux) and np.isfinite(i_flux)) or i_flux <= 0:
        return "nan"
    frac = abs(v_flux) / i_flux
    if frac >= strong:
        return "highly_circular"
    if frac >= weak:
        return "circular"
    return "weak"


def synthetic_field(
    n_stars: int = 500,
    *,
    pol_fraction: float = 0.05,
    leakage_scale: float = 0.008,
    dt_yr: float = 20.0,
    seed: int = 0,
) -> tuple[dict[str, np.ndarray], float]:
    """Synthetic forced-photometry Stokes-V target sample (offline fixture).

    Models a curated late-type-star target list: each star has a forced $(I,V)$ measurement at the
    radio epoch. A ``pol_fraction`` minority are genuine coherent emitters with deep circular
    polarization ($|V|/I\\sim0.2$--$0.8$); the rest show only instrumental leakage
    ($|V|/I\\sim$ ``leakage_scale``). Emitters have a real radio centroid at the *proper-motion-
    propagated* position; non-emitters' forced centroid scatters about the catalogue position.
    Returns ``(stars, dt_yr)`` where ``stars`` has ra/dec (catalogue epoch), pmra/pmdec, the radio
    centroid ra_radio/dec_radio, fluxes, and the ``is_emitter`` truth mask.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(40.0, 60.0, n_stars)
    dec = rng.uniform(-40.0, -20.0, n_stars)
    # nearby late-type stars have appreciable proper motions (tens to hundreds of mas/yr)
    pmra = rng.normal(0.0, 150.0, n_stars)
    pmdec = rng.normal(0.0, 150.0, n_stars)
    i_flux = 10.0 ** rng.uniform(0.7, 2.5, n_stars)  # ~5-300 mJy Stokes I
    e_i = 0.05 * i_flux + 0.25  # mJy (~5% + a noise floor)
    e_v = 0.25 * np.ones(n_stars)  # ~0.25 mJy V noise per RACS image

    is_emitter = rng.random(n_stars) < pol_fraction
    sign = rng.choice([-1.0, 1.0], n_stars)
    # leakage component: a half-normal |V/I| centred near leakage_scale for every star
    leak_frac = np.abs(rng.normal(0.0, leakage_scale, n_stars))
    v_flux = sign * leak_frac * i_flux + rng.normal(0.0, e_v)
    # genuine emitters: deep circular polarization replaces the leakage value
    frac_em = rng.uniform(0.2, 0.8, n_stars)
    v_flux[is_emitter] = sign[is_emitter] * frac_em[is_emitter] * i_flux[is_emitter]

    cosd = np.cos(np.radians(dec))
    ra_prop = ra + (pmra * dt_yr / 1000.0 / 3600.0) / cosd
    dec_prop = dec + pmdec * dt_yr / 1000.0 / 3600.0
    jit = lambda s: rng.normal(0.0, s / 3600.0, n_stars)  # noqa: E731
    # emitters sit at the PM-propagated position; non-emitter forced centroids scatter near catalogue
    ra_radio = np.where(is_emitter, ra_prop + jit(0.5), ra + jit(1.0))
    dec_radio = np.where(is_emitter, dec_prop + jit(0.5), dec + jit(1.0))

    stars = {
        "ra": ra,
        "dec": dec,
        "pmra": pmra,
        "pmdec": pmdec,
        "ra_radio": ra_radio,
        "dec_radio": dec_radio,
        "i_flux": i_flux,
        "v_flux": v_flux,
        "e_i": e_i,
        "e_v": e_v,
        "is_emitter": is_emitter,
    }
    return stars, dt_yr


def run(
    out: str = ".", *, offline: bool = True, v_snr_min: float = 5.0, i_snr_ref: float = 10.0
) -> dict:
    """Full slice: forced-photometry V over a target star list, select coherent emitters, write artifacts.

    Offline uses :func:`synthetic_field`; the real path (forced photometry of a curated target list in
    CASDA RACS Stokes-V cutouts) is wired separately. The leakage floor is estimated from the
    **bright** sources only ($I/\\sigma_I>$ ``i_snr_ref``), where $|V/I|$ reflects instrumental leakage
    rather than $V$ noise; the noise on faint targets is handled separately by the ``v_snr_min`` cut.
    Candidates clear the leakage floor *and* the $V$-SNR cut, and are confirmed by proper motion.
    """
    import json
    from pathlib import Path

    if offline:
        stars, dt_yr = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        raise NotImplementedError("real CASDA forced-photometry path is wired in the next step")

    frac, _ = fractional_circular_pol(stars["v_flux"], stars["i_flux"], stars["e_v"], stars["e_i"])
    bright = (
        stars["i_flux"] / stars["e_i"]
    ) > i_snr_ref  # leakage ~ |V/I| only where V noise is small
    threshold = leakage_floor(frac[bright] if bright.any() else frac)
    mask, _ = select_circular_pol(
        stars["i_flux"],
        stars["v_flux"],
        stars["e_i"],
        stars["e_v"],
        leakage_threshold=threshold,
        v_snr_min=v_snr_min,
    )
    pm_ok, sep = proper_motion_confirm(
        stars["ra_radio"],
        stars["dec_radio"],
        stars["ra"],
        stars["dec"],
        stars["pmra"],
        stars["pmdec"],
        dt_yr,
    )
    candidates = mask & pm_ok
    truth = stars["is_emitter"]
    n_cand = int(candidates.sum())
    metrics = {
        "source": source,
        "n_targets": int(truth.size),
        "n_bright_ref": int(bright.sum()),
        "leakage_floor_pct": round(100.0 * threshold, 3),
        "n_above_floor": int(mask.sum()),
        "n_candidates": n_cand,
        "n_injected": int(truth.sum()),
        "n_recovered": int((candidates & truth).sum()),
        "n_pm_rejected": int((mask & ~pm_ok).sum()),
        "purity": round(float((candidates & truth).sum()) / n_cand, 3) if n_cand else 0.0,
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "stokesv_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(stars, frac, threshold, candidates, op / "papers" / "stokesv" / "figures")
    _write_macros(metrics, op / "papers" / "stokesv" / "generated" / "macros.tex")
    return metrics


def _figure(
    stars: dict, frac: np.ndarray, threshold: float, candidates: np.ndarray, out_dir
) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    i = stars["i_flux"]
    other = ~candidates
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(i[other], frac[other], s=6, color="0.6", label="targets (leakage)")
    ax.scatter(i[candidates], frac[candidates], s=40, color="r", marker="*", label="V candidates")
    ax.axhline(
        threshold, color="b", ls="--", lw=0.8, label=rf"leakage floor ({threshold * 100:.1f}%)"
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel=r"Stokes $I$ (mJy)",
        ylabel=r"$|V|/I$",
        title="Circular-polarization selection",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "circular_pol.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.stokesv._write_macros — do not edit by hand.",
        rf"\newcommand{{\svSource}}{{{m['source']}}}",
        rf"\newcommand{{\svNtargets}}{{{m['n_targets']}}}",
        rf"\newcommand{{\svLeakFloorPct}}{{{m['leakage_floor_pct']}}}",
        rf"\newcommand{{\svNcandidates}}{{{m['n_candidates']}}}",
        rf"\newcommand{{\svNinjected}}{{{m.get('n_injected', 0)}}}",
        rf"\newcommand{{\svNrecovered}}{{{m.get('n_recovered', 0)}}}",
        rf"\newcommand{{\svNpmrejected}}{{{m.get('n_pm_rejected', 0)}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Select Stokes-V coherent radio emitters (RACS).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--v-snr-min", type=float, default=5.0)
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=True, v_snr_min=args.v_snr_min), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
