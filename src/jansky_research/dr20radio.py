"""First radio-counterpart census of the SDSS-V DR20 Black Hole Mapper quasars.

DR20 (Almeida et al., arXiv:2607.26149) delivers ~500k BHM spectroscopic objects including
the first optical SDSS spectra from the southern hemisphere. No radio cross-match of any
SDSS-V catalog exists (plans/88): this module matches the DR20 quasar table against VLASS
(north, local CIRADA epoch catalogs) and RACS (south), with the false-match rate *measured*
by position-shift trials and the radio-targeted open-fiber cartons excluded from fractions
(they were selected BECAUSE they are radio sources — counting them would be circular). Their
match rate against the SELECTING survey (RACS, southern leg) is the ~100% pipeline validation;
against VLASS at 3 GHz it is a cross-frequency detection fraction (steep-spectrum sources
selected at 144/888 MHz routinely fade below VLASS depth) — measured, not assumed.

Committed-real-results pattern: real legs write force-tracked ``results/dr20radio_*.json``;
paper macros come only from committed evidence; synthetic fixtures feed tests alone.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "RADIO_CARTON_PREFIXES",
    "SPALL_URL",
    "VLASS_DEC_LIMIT_DEG",
    "crossmatch",
    "detection_fraction",
    "false_match_rate",
    "log_luminosity_whz",
    "luminosity_matched_fractions",
    "parse_racs_csv",
    "run_north",
    "run_south",
    "select_quasars",
    "synthetic_survey",
    "synthetic_two_surveys",
    "wilson_interval",
]

SPALL_URL = (
    "https://data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/summary/allepoch/"
    "spAll-lite-v6_2_1-allepoch.fits.gz"
)
# Open-fiber cartons that TARGETED radio sources — excluded from detection fractions
# (circular) and used as a positive-control validation set instead.
RADIO_CARTON_PREFIXES = (
    "openfibertargets_bhm_racsradio",
    "openfibertargets_bhm_lofarradio",
)
VLASS_DEC_LIMIT_DEG = -40.0
# Common-luminosity-limit inputs (stated, not hidden): VLASS QL per-epoch reliability
# threshold 1.0 mJy (CIRADA catalog user guide); RACS-low DR1 95% point-source completeness
# 3.0 mJy (Hale et al. 2021). Canonical alpha = -0.7 for the K-correction.
VLASS_FREQ_GHZ = 3.0
VLASS_S_LIM_MJY = 1.0
RACS_FREQ_GHZ = 0.8875
RACS_S_LIM_MJY = 3.0
# Hale et al. 2021 quote ~3 mJy (source-count based) and ~5 mJy (simulation based) for the
# 95% completeness; both variants are computed so the paper can show the sensitivity.
RACS_S_LIM_CONSERVATIVE_MJY = 5.0
# The MIRROR variant, added 2026-08-12. Raising the RACS limit provably cannot move the
# north/south ratio (the common limit is the RACS one in both legs, so it rescales them
# identically) -- so the only "conservative" check in this module was one that could not fail.
# The VLASS side is the one with leverage: VLASS's 1 mJy is a per-epoch RELIABILITY threshold
# while RACS's 3 mJy is a 95% COMPLETENESS limit, which are not the same kind of number, and
# the northern fraction is steeply sensitive to its own cut.
VLASS_S_LIM_CONSERVATIVE_MJY: tuple[float, float] = (2.0, 3.0)
# Synthesised-beam FWHMs, quoted in the papers to justify the two match radii.
VLASS_BEAM_ARCSEC = 2.5
RACS_BEAM_ARCSEC = 25.0


def select_quasars(
    cls: np.ndarray,
    zwarning: np.ndarray,
    firstcarton: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(quasar_mask, radio_carton_mask)`` from spAll columns.

    A quasar is ``CLASS == 'QSO'`` with a clean redshift (``ZWARNING == 0``). The radio-carton
    mask flags objects whose ``FIRSTCARTON`` starts with any radio-selected open-fiber prefix;
    it is reported separately, never inside the census fractions.
    """
    cls = np.char.strip(np.asarray(cls, dtype=str))
    carton = np.char.strip(np.asarray(firstcarton, dtype=str))
    quasar = (cls == "QSO") & (np.asarray(zwarning) == 0)
    radio_carton = np.zeros(carton.size, dtype=bool)
    for prefix in RADIO_CARTON_PREFIXES:
        radio_carton |= np.char.startswith(carton, prefix)
    return quasar, radio_carton


def crossmatch(
    ra1_deg: np.ndarray,
    dec1_deg: np.ndarray,
    ra2_deg: np.ndarray,
    dec2_deg: np.ndarray,
    *,
    radius_arcsec: float = 2.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nearest-neighbour sky match of catalog 1 onto catalog 2.

    Returns ``(matched_mask, sep_arcsec, idx)`` for catalog-1 rows — ``idx`` is the index of
    the nearest catalog-2 source (for flux lookup), valid where ``matched_mask`` is True.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c1 = SkyCoord(np.asarray(ra1_deg, float) * u.deg, np.asarray(dec1_deg, float) * u.deg)
    c2 = SkyCoord(np.asarray(ra2_deg, float) * u.deg, np.asarray(dec2_deg, float) * u.deg)
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    sep_as = sep.arcsec
    return sep_as <= radius_arcsec, sep_as, np.asarray(idx)


def false_match_rate(
    ra_q: np.ndarray,
    dec_q: np.ndarray,
    ra_r: np.ndarray,
    dec_r: np.ndarray,
    *,
    radius_arcsec: float = 2.5,
    n_trials: int = 10,
    shift_arcmin: tuple[float, float] = (5.0, 30.0),
    seed: int = 0,
) -> dict:
    """Chance-coincidence match rate, MEASURED by rigid position-shift trials.

    Each trial shifts every quasar position by a random offset (uniform in ``shift_arcmin``,
    random direction — large enough to decorrelate, small enough to stay in the same source
    density) and re-runs the match. The mean shifted match fraction is the false-match rate
    to subtract from raw detection fractions; the std across trials is its uncertainty.
    """
    rng = np.random.default_rng(seed)
    rates = []
    dec_q = np.asarray(dec_q, float)
    ra_q = np.asarray(ra_q, float)
    for _ in range(n_trials):
        amp = rng.uniform(*shift_arcmin, size=ra_q.size) / 60.0
        ang = rng.uniform(0, 2 * np.pi, size=ra_q.size)
        dec_s = np.clip(dec_q + amp * np.sin(ang), -89.9, 89.9)
        ra_s = (ra_q + amp * np.cos(ang) / np.cos(np.deg2rad(dec_s))) % 360.0
        m, _, _ = crossmatch(ra_s, dec_s, ra_r, dec_r, radius_arcsec=radius_arcsec)
        rates.append(float(np.mean(m)))
    return {
        "rate": float(np.mean(rates)),
        "std": float(np.std(rates)),
        "n_trials": n_trials,
        "radius_arcsec": radius_arcsec,
    }


def wilson_interval(k: int, n: int, *, z: float = 1.0) -> tuple[float, float, float]:
    """Wilson score interval for a binomial fraction: ``(fraction, lo, hi)``."""
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return p, max(center - half, 0.0), min(center + half, 1.0)


def detection_fraction(values: np.ndarray, matched: np.ndarray, *, bins: np.ndarray) -> dict:
    """Matched fraction binned over ``values`` (e.g. redshift), with Wilson errors."""
    values = np.asarray(values, float)
    matched = np.asarray(matched, bool)
    out: dict[str, list] = {
        "bin_lo": [],
        "bin_hi": [],
        "n": [],
        "k": [],
        "frac": [],
        "lo": [],
        "hi": [],
    }
    for lo, hi in zip(bins[:-1], bins[1:], strict=False):
        sel = (values >= lo) & (values < hi)
        n, k = int(sel.sum()), int((matched & sel).sum())
        p, plo, phi = wilson_interval(k, n)
        out["bin_lo"].append(float(lo))
        out["bin_hi"].append(float(hi))
        out["n"].append(n)
        out["k"].append(k)
        out["frac"].append(None if n == 0 else round(p, 5))
        out["lo"].append(None if n == 0 else round(plo, 5))
        out["hi"].append(None if n == 0 else round(phi, 5))
    return out


def synthetic_survey(
    *,
    n_quasars: int = 2000,
    n_radio: int = 5000,
    counterpart_fraction: float = 0.12,
    n_carton: int = 50,
    area_deg: float = 10.0,
    radius_arcsec: float = 2.5,
    seed: int = 0,
) -> dict:
    """Synthetic quasar + radio catalogs with known truth, for the offline round-trip.

    A ``counterpart_fraction`` of ordinary quasars gets a radio source planted at its
    position (within 1"); ``n_carton`` extra quasars are radio-TARGETED (counterpart by
    construction) and carry a radio-carton ``FIRSTCARTON`` string — the circularity the
    census must exclude. Everything lives in an ``area_deg``-sized box at Dec ~ +10.
    """
    rng = np.random.default_rng(seed)
    ra_q = rng.uniform(180.0, 180.0 + area_deg, n_quasars)
    dec_q = rng.uniform(10.0, 10.0 + area_deg, n_quasars)
    z_q = rng.uniform(0.3, 4.0, n_quasars)
    ra_r = rng.uniform(180.0, 180.0 + area_deg, n_radio)
    dec_r = rng.uniform(10.0, 10.0 + area_deg, n_radio)
    is_cp = rng.random(n_quasars) < counterpart_fraction
    jitter = 1.0 / 3600.0
    ra_r = np.concatenate([ra_r, ra_q[is_cp] + rng.normal(0, jitter, is_cp.sum())])
    dec_r = np.concatenate([dec_r, dec_q[is_cp] + rng.normal(0, jitter, is_cp.sum())])
    # radio-carton quasars: targeted because radio-detected -> counterpart by construction
    ra_c = rng.uniform(180.0, 180.0 + area_deg, n_carton)
    dec_c = rng.uniform(10.0, 10.0 + area_deg, n_carton)
    ra_r = np.concatenate([ra_r, ra_c + rng.normal(0, jitter, n_carton)])
    dec_r = np.concatenate([dec_r, dec_c + rng.normal(0, jitter, n_carton)])
    cls = np.array(["QSO"] * (n_quasars + n_carton))
    zwarning = np.zeros(n_quasars + n_carton, dtype=int)
    carton = np.array(
        ["bhm_spiders_agn"] * n_quasars + ["openfibertargets_bhm_racsradio_boss"] * n_carton
    )
    return {
        "ra_q": np.concatenate([ra_q, ra_c]),
        "dec_q": np.concatenate([dec_q, dec_c]),
        "z_q": np.concatenate([z_q, rng.uniform(0.3, 4.0, n_carton)]),
        "cls": cls,
        "zwarning": zwarning,
        "firstcarton": carton,
        "ra_r": ra_r,
        "dec_r": dec_r,
        "true_fraction": counterpart_fraction,
        "radius_arcsec": radius_arcsec,
    }


RACS_TAP_SYNC = "https://casda.csiro.au/casda_vo_tools/tap/sync"
RACS_TABLE = "AS110.racs_dr1_sources_galacticcut_v2021_08_v02"  # 2,123,638 sources (Hale+ 2021)


def synthetic_two_surveys(
    *,
    fade_fraction: float = 0.35,
    seed: int = 0,
    **kwargs,
) -> dict:
    """Two-survey variant of :func:`synthetic_survey` — the increment-1 blind spot, fixed.

    The base survey acts as the SELECTING survey: every radio-carton quasar has a counterpart
    there by construction. The second ("matching") survey keeps only ``fade_fraction`` of the
    carton counterparts (spectral fading between observing frequencies) while ordinary
    counterparts carry over unchanged — so a census pipeline must see carton-match ~100%
    against the selecting survey but only ~``fade_fraction`` against the other one.
    """
    s = synthetic_survey(seed=seed, **kwargs)
    rng = np.random.default_rng(seed + 1)
    n_carton = int((np.char.startswith(s["firstcarton"], "openfibertargets")).sum())
    # the last n_carton radio sources are the carton counterparts (see synthetic_survey)
    keep = np.ones(s["ra_r"].size, dtype=bool)
    carton_slice = np.arange(s["ra_r"].size - n_carton, s["ra_r"].size)
    keep[carton_slice] = rng.random(n_carton) < fade_fraction
    s["ra_r2"] = s["ra_r"][keep]
    s["dec_r2"] = s["dec_r"][keep]
    s["fade_fraction"] = fade_fraction
    return s


def parse_racs_csv(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse a RACS TAP CSV chunk into ``(ra, dec, peak_flux_mjy)`` arrays."""
    ra, dec, flux = [], [], []
    for line in text.splitlines()[1:]:
        parts = line.split(",")
        if len(parts) != 3:
            continue
        try:
            ra.append(float(parts[0]))
            dec.append(float(parts[1]))
            flux.append(float(parts[2]))
        except ValueError:
            continue
    return np.array(ra), np.array(dec), np.array(flux)


def log_luminosity_whz(
    z: np.ndarray,
    s_mjy: np.ndarray,
    *,
    freq_ghz: float,
    alpha: float | np.ndarray = -0.7,
    ref_freq_ghz: float = 1.4,
) -> np.ndarray:
    """log10 rest-frame 1.4 GHz spectral luminosity (W/Hz) from observed flux density.

    K-corrected assuming a power law ``S ~ nu^alpha`` (alpha = -0.7, the canonical synchrotron
    slope — an assumption the paper states, not hides): ``L = 4 pi d_L^2 S (1+z)^-(1+alpha)
    (ref/freq)^alpha``.
    """
    from astropy.cosmology import Planck18

    dl_m = Planck18.luminosity_distance(np.asarray(z, float)).to("m").value
    s_si = np.asarray(s_mjy, float) * 1e-29  # mJy -> W m^-2 Hz^-1
    with np.errstate(divide="ignore", invalid="ignore"):
        lum = (
            4.0
            * np.pi
            * dl_m**2
            * s_si
            * (1.0 + np.asarray(z, float)) ** (-(1.0 + alpha))
            * (ref_freq_ghz / freq_ghz) ** alpha
        )
        return np.log10(lum)


ALPHA_SWEEP: tuple[float, ...] = (0.0, -0.35, -0.7, -1.0)
# The MEASURED spectral index and its +/-1 bootstrap SE, from the flux-complete joint-detection
# sample in results/dr20radio_alpha.json (run_alpha, 2026-08-12: 4190 quasars detected by both
# VLASS and RACS in the overlap band, above the flux where the joint-detection requirement
# cannot truncate the distribution). These replace the sweep as the paper's uncertainty on the
# K-correction: the sweep answered "how bad could it be", this answers "what is it".
# test_dr20radio.py asserts these stay equal to the committed measurement.
# Kaplan-Meier over ALL 6,626 RACS-detected band quasars, treating the 1,055 VLASS
# non-detections as the left-censored observations they are. This replaced the
# flux-complete-cut value (-0.7218) as the primary estimator on 2026-08-12: the cut is
# unbiased only over the alpha range it assumes and discards a quarter of the sample, while
# KM uses every object and is unbiased over the whole range. It lands steeper, as it must --
# the objects a cut throws away are the steep ones.
ALPHA_MEASURED: float = -0.7546
ALPHA_MEASURED_SE: float = 0.0118
# The same estimator restricted to the faintest flux bin. The index depends on flux (a real
# trend: it survives the censoring correction), and the census's detection decisions happen
# at the faint end, so this is the value appropriate to converting a flux LIMIT -- while
# -0.7546 is appropriate to the sample as a whole. The contrast is reported at both.
ALPHA_THRESHOLD_REGIME: float = -0.6095
ALPHA_MEASURED_SWEEP: tuple[float, ...] = (
    ALPHA_MEASURED - ALPHA_MEASURED_SE,
    ALPHA_MEASURED,
    ALPHA_MEASURED + ALPHA_MEASURED_SE,
    ALPHA_THRESHOLD_REGIME,
)
"""Spectral indices for the luminosity-matched contrast.

The K-correction is the only thing making the 3 GHz (VLASS) and 888 MHz (RACS) legs
comparable, and it is not a measured quantity -- so the contrast has to be reported as a
range across it, not at one value. At alpha = -0.7 the north's effective flux cut is
1.28 mJy against the south's 3.0 mJy; at alpha = 0 both become 3.0 mJy.
"""


def spectral_index(
    s_hi_mjy: np.ndarray, s_lo_mjy: np.ndarray, *, freq_hi_ghz: float, freq_lo_ghz: float
) -> np.ndarray:
    """Two-point spectral index alpha, in the convention S ~ nu^alpha.

    Negative alpha is a steep (falling) spectrum. Non-positive fluxes give NaN rather than an
    exception: a catalog peak flux can be <= 0 after deconvolution and such rows must drop out
    of a distribution, not poison it.
    """
    s_hi = np.asarray(s_hi_mjy, float)
    s_lo = np.asarray(s_lo_mjy, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.log(s_hi / s_lo) / np.log(freq_hi_ghz / freq_lo_ghz)
    return np.where((s_hi > 0) & (s_lo > 0), a, np.nan)


def kaplan_meier_median(
    values: np.ndarray, censored: np.ndarray, *, right_censored: bool = True
) -> dict:
    """Kaplan-Meier median of a censored sample, with no distributional assumption.

    ``censored[i]`` marks an observation known only to lie beyond ``values[i]`` (above it for
    ``right_censored``, below it otherwise). Left-censored data are handled by negating, which
    turns them into the right-censored case and back again.

    Why this instead of a completeness cut: throwing away everything below a flux where the
    censoring cannot bite is unbiased only over the alpha range the cut assumes, and it
    discards most of the sample. Kaplan-Meier keeps every object by using the non-detections
    as the information they actually are -- "this source is steeper than X" -- so the estimate
    is unbiased over the whole range.

    Returns the median, the number of events and censored points, and the fraction of the
    sample the estimator can resolve (the KM curve is undefined below the last event if the
    tail is entirely censored, in which case ``median`` is None rather than a guess).
    """
    v = np.asarray(values, float)
    c = np.asarray(censored, bool)
    ok = np.isfinite(v)
    v, c = v[ok], c[ok]
    if v.size == 0:
        return {"median": None, "n_events": 0, "n_censored": 0, "n": 0}
    t = v if right_censored else -v
    order = np.argsort(t, kind="stable")
    t, c = t[order], c[order]
    surv, n_at_risk, median_t = 1.0, t.size, None
    for i, ti in enumerate(t):
        if c[i]:
            n_at_risk -= 1
            continue
        # events are processed one at a time; ties fall out of the product identically
        surv *= 1.0 - 1.0 / n_at_risk
        n_at_risk -= 1
        if median_t is None and surv <= 0.5:
            median_t = float(ti)
    return {
        "median": None if median_t is None else float(median_t if right_censored else -median_t),
        "n_events": int((~c).sum()),
        "n_censored": int(c.sum()),
        "n": int(t.size),
        "survival_at_last_event": float(surv),
    }


def censored_median_bounds(value: np.ndarray, kind: np.ndarray) -> dict:
    """Distribution-free bounds on a median when data are censored from BOTH sides.

    ``kind`` is 0 for a measured value, -1 for left-censored (the truth is below ``value``)
    and +1 for right-censored (the truth is above it). Kaplan-Meier handles one side only;
    with both, no point estimate is identified without extra assumptions, but the median is
    bounded, and the bound is honest where a one-sided estimate is merely convenient.

    This exists because a one-sided estimate here is biased in a nameable direction: RACS
    detections without a VLASS counterpart are the steep sources, VLASS detections without a
    RACS counterpart are the flat ones, and dropping either tail tilts the median toward the
    other. Bounding uses both.
    """
    v = np.asarray(value, float)
    k = np.asarray(kind, int)
    ok = np.isfinite(v)
    v, k = v[ok], k[ok]
    n = v.size
    if n == 0:
        return {"lo": None, "hi": None, "n": 0}
    grid = np.unique(v)
    meas, left, right = v[k == 0], v[k == -1], v[k == 1]
    # F_lo: only what must lie below x. F_hi: everything that could.
    f_lo = np.array([((meas < x).sum() + (left <= x).sum()) / n for x in grid])
    f_hi = np.array([((meas < x).sum() + left.size + (right < x).sum()) / n for x in grid])
    # Step OUTWARD to the grid point before each crossing. The bound is over a discrete grid
    # of observed values, so the true median can sit between two of them; rounding inward
    # would report an interval that excludes the quantity it bounds, which is worse than a
    # slightly wide one.
    i_lo = int(np.argmax(f_hi >= 0.5)) if (f_hi >= 0.5).any() else None
    i_hi = int(np.argmax(f_lo >= 0.5)) if (f_lo >= 0.5).any() else None
    lo = None if i_lo is None else grid[max(i_lo - 1, 0)]
    hi = None if i_hi is None else grid[min(i_hi + 1, grid.size - 1)]
    return {
        "lo": None if lo is None else float(lo),
        "hi": None if hi is None else float(hi),
        "n": int(n),
        "n_measured": int(meas.size),
        "n_left_censored": int(left.size),
        "n_right_censored": int(right.size),
    }


def alpha_complete_limit_mjy(
    *,
    s_lim_hi_mjy: float,
    freq_hi_ghz: float,
    freq_lo_ghz: float,
    alpha_min: float = -1.5,
) -> float:
    """Low-frequency flux above which a source of any ``alpha >= alpha_min`` stays detectable.

    Requiring a detection in BOTH surveys truncates the alpha distribution: at the shallower
    survey's limit a steep source has already fallen below the deeper survey's limit at the
    higher frequency, so it is missing from the joint sample and the measured median comes out
    too flat. Above the flux returned here that truncation cannot operate for any
    ``alpha >= alpha_min``, so the sub-sample is unbiased over that range -- at the cost of
    sample size. Quote both, and never quote the joint-sample median alone.
    """
    return float(s_lim_hi_mjy * (freq_hi_ghz / freq_lo_ghz) ** (-alpha_min))


def luminosity_matched_per_source_alpha(
    z: np.ndarray,
    matched: np.ndarray,
    s_matched_mjy: np.ndarray,
    *,
    freq_ghz: float,
    s_lim_this_mjy: float,
    s_lim_other_mjy: float,
    other_freq_ghz: float,
    bins: np.ndarray,
    alpha_samples: np.ndarray,
    n_real: int = 20,
    seed: int = 0,
) -> dict:
    """The comparison with each quasar given its OWN spectral index, drawn from the measured
    distribution, instead of one median index for the whole population.

    This measures a bias, not a variance. Drawing per-source indices and looking at the spread
    across realizations would be self-deceiving: the fraction averages over ~10^5 sources, so
    realization scatter falls as 1/sqrt(N) and comes out negligible no matter how broad the
    index distribution is. That is the same trap as quoting a bootstrap SE for something the
    bootstrap cannot see. What the population scatter actually does is shift the answer,
    because the detection criterion is a threshold and a threshold is not linear in alpha: the
    steep and flat halves of the distribution do not cancel. So the quantity to report is the
    DIFFERENCE between this and the single-median-index result, with the (small) realization
    spread quoted alongside only to show it is small.
    """
    rng = np.random.default_rng(seed)
    z = np.asarray(z, float)
    matched = np.asarray(matched, bool)
    samples = np.asarray(alpha_samples, float)
    samples = samples[np.isfinite(samples)]
    totals = []
    for _ in range(n_real):
        a_i = rng.choice(samples, size=z.size, replace=True)
        lim_this = log_luminosity_whz(
            z, np.full(z.size, s_lim_this_mjy), freq_ghz=freq_ghz, alpha=a_i
        )
        lim_other = log_luminosity_whz(
            z, np.full(z.size, s_lim_other_mjy), freq_ghz=other_freq_ghz, alpha=a_i
        )
        lum = np.full(z.size, -np.inf)
        lum[matched] = log_luminosity_whz(
            z[matched], s_matched_mjy[matched], freq_ghz=freq_ghz, alpha=a_i[matched]
        )
        above = matched & (lum >= np.maximum(lim_this, lim_other))
        totals.append(float(np.mean(above[np.isfinite(z)])))
    return {
        "fraction_mean": float(np.mean(totals)),
        "fraction_realization_sd": float(np.std(totals)),
        "n_real": n_real,
        "n_alpha_samples": int(samples.size),
        "seed": seed,
    }


def luminosity_matched_fractions(
    z: np.ndarray,
    matched: np.ndarray,
    s_matched_mjy: np.ndarray,
    *,
    freq_ghz: float,
    s_lim_this_mjy: float,
    s_lim_other_mjy: float,
    other_freq_ghz: float,
    bins: np.ndarray,
    alpha: float = -0.7,
) -> dict:
    """Detection fractions above the COMMON luminosity limit of two surveys.

    At each quasar's z, both surveys' flux limits convert to rest-1.4 GHz luminosity limits;
    the common limit is the max. A quasar counts as detected only if matched AND its
    counterpart luminosity clears the common limit — making north (VLASS, 3 GHz) and south
    (RACS, 0.888 GHz) fractions comparable despite different depths and frequencies.
    """
    z = np.asarray(z, float)
    lim_this = log_luminosity_whz(
        z, np.full(z.size, s_lim_this_mjy), freq_ghz=freq_ghz, alpha=alpha
    )
    lim_other = log_luminosity_whz(
        z, np.full(z.size, s_lim_other_mjy), freq_ghz=other_freq_ghz, alpha=alpha
    )
    common = np.maximum(lim_this, lim_other)
    lum = np.full(z.size, -np.inf)
    lum[matched] = log_luminosity_whz(
        z[matched], s_matched_mjy[matched], freq_ghz=freq_ghz, alpha=alpha
    )
    above = np.asarray(matched, bool) & (lum >= common)
    out = detection_fraction(z, above, bins=bins)
    out["s_lim_this_mjy"] = s_lim_this_mjy
    out["s_lim_other_mjy"] = s_lim_other_mjy
    out["alpha"] = alpha
    return out


# ------------------------------------------------------------------------- real data legs


def fetch_spall(dest_dir: str = "data") -> str:  # pragma: no cover - network
    """Download the DR20 spAll-lite summary file (177 MiB, resumable); returns the path."""
    from pathlib import Path
    from urllib.request import Request, urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / SPALL_URL.rsplit("/", 1)[1]
    if path.exists():
        return str(path)
    part = path.with_suffix(path.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    req = Request(SPALL_URL)
    if offset:
        req.add_header("Range", f"bytes={offset}-")
    with urlopen(req, timeout=60) as r, open(part, "ab" if offset else "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    part.rename(path)
    return str(path)


def read_spall_quasars(path: str) -> dict:
    """Read the spAll-lite columns the census needs and apply :func:`select_quasars`."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        cls = d["CLASS"]
        zw = np.asarray(d["ZWARNING"])
        carton = d["FIRSTCARTON"]
        quasar, radio_carton = select_quasars(cls, zw, carton)
        carton_s = np.char.strip(np.asarray(carton, dtype=str))
        return {
            "ra": np.asarray(d["RACAT"], float)[quasar],
            "dec": np.asarray(d["DECCAT"], float)[quasar],
            "z": np.asarray(d["Z"], float)[quasar],
            "obs": np.char.strip(np.asarray(d["OBS"], dtype=str))[quasar],
            "radio_carton": radio_carton[quasar],
            "carton_racs": np.char.startswith(carton_s, "openfibertargets_bhm_racsradio")[quasar],
            "carton_lofar": np.char.startswith(carton_s, "openfibertargets_bhm_lofarradio")[quasar],
            "n_total_rows": int(len(d)),
        }


def load_vlass_positions(*, total_flux: bool = False) -> dict:  # pragma: no cover - local bulk
    """Thin wrapper preserved for callers that want integrated rather than peak flux."""
    return _load_vlass(total_flux=total_flux)


def _load_vlass(*, total_flux: bool = False) -> dict:  # pragma: no cover - local bulk files
    """Positions of quality-cut VLASS components from the local epoch catalogs.

    Applies the same cuts as the merged `vlass` slice: ``Duplicate_flag < 2``,
    ``Quality_flag in (0, 4)``, ``S_Code != 'E'``. Returns ``{"E2": (ra, dec, peak_mjy), ...}``.
    """
    import csv
    import gzip

    from astropy.io import fits

    out = {}
    ra, dec, fx = [], [], []
    with gzip.open("data/CIRADA_VLASS2QLv2_table1_components.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            try:
                if int(float(r["Duplicate_flag"])) >= 2:
                    continue
                if int(float(r["Quality_flag"])) not in (0, 4):
                    continue
                if r["S_Code"].strip() == "E":
                    continue
                ra.append(float(r["RA"]))
                dec.append(float(r["DEC"]))
                fx.append(float(r["Total_flux" if total_flux else "Peak_flux"]))
            except (KeyError, ValueError):
                continue
    out["E2"] = (np.array(ra), np.array(dec), np.array(fx))
    ra3, dec3, fx3 = [], [], []
    for name in ("data/QL3.1_components.fits", "data/QL3.2_components.fits"):
        with fits.open(name, memmap=True) as hdul:
            d = hdul[1].data
            # The E3 interim lists (VLASS Memo 22) carry a simplified schema: a binary
            # quality `Flag` (0 = good) and no Duplicate_flag/Quality_flag columns; empty
            # islands (S_Code 'E') are already absent.
            ok = np.asarray(d["Flag"]) == 0
            ra3.append(np.asarray(d["RA"], float)[ok])
            dec3.append(np.asarray(d["DEC"], float)[ok])
            fx3.append(np.asarray(d["Total_flux" if total_flux else "Peak_flux"], float)[ok])
    out["E3"] = (np.concatenate(ra3), np.concatenate(dec3), np.concatenate(fx3))
    return out


def run_north(
    out: str = ".", *, radius_arcsec: float = 2.5, n_shift_trials: int = 10
) -> dict:  # pragma: no cover - network + bulk local data (pure pieces tested offline)
    """Real leg A: DR20 quasars (Dec > -40) vs the local VLASS epoch catalogs.

    Writes the committed evidence file ``results/dr20radio_north.json``.
    """
    import json
    from pathlib import Path

    spall = fetch_spall()
    q = read_spall_quasars(spall)
    north = q["dec"] > VLASS_DEC_LIMIT_DEG
    census = north & ~q["radio_carton"]
    carton = north & q["radio_carton"]
    vlass = load_vlass_positions()
    _asamp = _alpha_samples()
    zbins = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0])
    epochs = {}
    matched_any = np.zeros(int(census.sum()), dtype=bool)
    s_any = np.zeros(int(census.sum()))
    for name, (ra_r, dec_r, fx_r) in vlass.items():
        m, _, idx = crossmatch(
            q["ra"][census], q["dec"][census], ra_r, dec_r, radius_arcsec=radius_arcsec
        )
        s_any = np.where(m, np.maximum(s_any, fx_r[idx]), s_any)
        fm = false_match_rate(
            q["ra"][census],
            q["dec"][census],
            ra_r,
            dec_r,
            radius_arcsec=radius_arcsec,
            n_trials=n_shift_trials,
        )

        def _cv(mask: np.ndarray, rr: np.ndarray = ra_r, dd: np.ndarray = dec_r) -> dict:
            mm, _, _ = crossmatch(
                q["ra"][mask], q["dec"][mask], rr, dd, radius_arcsec=radius_arcsec
            )
            return {
                "n": int(mm.size),
                "matched": int(mm.sum()),
                "fraction": round(float(np.mean(mm)), 4) if mm.size else None,
            }

        matched_any |= m
        p, lo, hi = wilson_interval(int(m.sum()), int(m.size))
        epochs[name] = {
            "n_census": int(m.size),
            "n_matched": int(m.sum()),
            "raw_fraction": round(p, 5),
            "wilson_lo": round(lo, 5),
            "wilson_hi": round(hi, 5),
            "false_match": fm,
            "corrected_fraction": round(p - fm["rate"], 5),
            "carton_validation": {
                "racsradio_cross_frequency": _cv(north & q["carton_racs"]),
                "lofarradio_cross_frequency": _cv(north & q["carton_lofar"]),
            },
            "n_radio_sources": int(ra_r.size),
        }
    p_any, lo_any, hi_any = wilson_interval(int(matched_any.sum()), int(matched_any.size))
    metrics = {
        "source": "SDSS-V DR20 spAll-lite v6_2_1 allepoch x VLASS QL E2+E3 (local CIRADA)",
        "n_spall_rows": q["n_total_rows"],
        "n_quasars_clean": int((q["z"] > -1).sum()),
        "n_north_census": int(census.sum()),
        "n_north_radio_carton_excluded": int(carton.sum()),
        "obs_breakdown_north_census": {
            o: int((q["obs"][census] == o).sum()) for o in np.unique(q["obs"][census])
        },
        "radius_arcsec": radius_arcsec,
        "epochs": epochs,
        "any_epoch": {
            "n_matched": int(matched_any.sum()),
            "raw_fraction": round(p_any, 5),
            "wilson_lo": round(lo_any, 5),
            "wilson_hi": round(hi_any, 5),
        },
        "fraction_vs_z_any_epoch": detection_fraction(q["z"][census], matched_any, bins=zbins),
        "luminosity_matched": luminosity_matched_fractions(
            q["z"][census],
            matched_any,
            s_any,
            freq_ghz=VLASS_FREQ_GHZ,
            s_lim_this_mjy=VLASS_S_LIM_MJY,
            s_lim_other_mjy=RACS_S_LIM_MJY,
            other_freq_ghz=RACS_FREQ_GHZ,
            bins=zbins,
        ),
        "luminosity_matched_alpha": {
            f"{a:g}": luminosity_matched_fractions(
                q["z"][census],
                matched_any,
                s_any,
                freq_ghz=VLASS_FREQ_GHZ,
                s_lim_this_mjy=VLASS_S_LIM_MJY,
                s_lim_other_mjy=RACS_S_LIM_MJY,
                other_freq_ghz=RACS_FREQ_GHZ,
                bins=zbins,
                alpha=a,
            )
            for a in (*ALPHA_SWEEP, *ALPHA_MEASURED_SWEEP)
        },
        "luminosity_matched_conservative": luminosity_matched_fractions(
            q["z"][census],
            matched_any,
            s_any,
            freq_ghz=VLASS_FREQ_GHZ,
            s_lim_this_mjy=VLASS_S_LIM_MJY,
            s_lim_other_mjy=RACS_S_LIM_CONSERVATIVE_MJY,
            other_freq_ghz=RACS_FREQ_GHZ,
            bins=zbins,
        ),
        "luminosity_matched_vlass_conservative": {
            f"{v:g}": luminosity_matched_fractions(
                q["z"][census],
                matched_any,
                s_any,
                freq_ghz=VLASS_FREQ_GHZ,
                s_lim_this_mjy=v,
                s_lim_other_mjy=RACS_S_LIM_MJY,
                other_freq_ghz=RACS_FREQ_GHZ,
                bins=zbins,
                alpha=ALPHA_MEASURED,
            )
            for v in VLASS_S_LIM_CONSERVATIVE_MJY
        },
        "luminosity_matched_per_source_alpha": (
            luminosity_matched_per_source_alpha(
                q["z"][census],
                matched_any,
                s_any,
                freq_ghz=VLASS_FREQ_GHZ,
                s_lim_this_mjy=VLASS_S_LIM_MJY,
                s_lim_other_mjy=RACS_S_LIM_MJY,
                other_freq_ghz=RACS_FREQ_GHZ,
                bins=zbins,
                alpha_samples=_asamp,
            )
            if _asamp.size
            else None
        ),
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "dr20radio_north.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


# RACS-low DR1 sky coverage, for the footprint diagnostic in run_south.
RACS_PLANE_CUT_DEG: float = 5.0
RACS_DEC_FLOOR_DEG: float = -85.0


def fetch_racs_positions(
    dest_dir: str = "data/racs_dr1",
    *,
    dec_min: float = -85.0,
    dec_max: float = 30.0,
    strip_deg: float = 1.0,
) -> dict:  # pragma: no cover - network
    """Fetch the RACS-low DR1 source positions by resumable 1-degree Dec strips.

    Each strip is cached as ``<dest>/strip_<dec>.csv`` (a completed strip is never re-fetched);
    the consolidated arrays are returned and cached as ``<dest>/racs_positions.npz``.
    """
    import time
    from pathlib import Path
    from urllib.parse import urlencode
    from urllib.request import urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    npz = dest / "racs_positions.npz"
    if npz.exists():
        d = np.load(npz)
        return {"ra": d["ra"], "dec": d["dec"], "flux": d["flux"]}
    ras, decs, fluxes = [], [], []
    lo = dec_min
    while lo < dec_max:
        hi = min(lo + strip_deg, dec_max)
        cache = dest / f"strip_{lo:+.0f}.csv"
        if not cache.exists():
            q = f"SELECT ra, dec, peak_flux FROM {RACS_TABLE} WHERE dec >= {lo} AND dec < {hi}"
            params = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
            for attempt in range(6):
                try:
                    with urlopen(f"{RACS_TAP_SYNC}?{params}", timeout=180) as r:
                        text = r.read().decode()
                    if not text.startswith("ra"):
                        raise OSError("unexpected TAP response")
                    cache.write_text(text)
                    break
                except (TimeoutError, OSError):
                    time.sleep(min(2**attempt, 60))
            else:
                raise OSError(f"RACS strip {lo} failed after retries")
            print(f"[dr20radio] RACS strip {lo:+.0f}: cached", flush=True)
        ra, dec, fx = parse_racs_csv(cache.read_text())
        ras.append(ra)
        decs.append(dec)
        fluxes.append(fx)
        lo = hi
    out = {
        "ra": np.concatenate(ras),
        "dec": np.concatenate(decs),
        "flux": np.concatenate(fluxes),
    }
    np.savez_compressed(npz, ra=out["ra"], dec=out["dec"], flux=out["flux"])
    return out


def fetch_racs_total_flux(
    dest_dir: str = "data/racs_dr1_total",
    *,
    dec_min: float = -40.0,
    dec_max: float = 30.0,
    strip_deg: float = 1.0,
) -> dict:  # pragma: no cover - network
    """RACS-low DR1 positions with BOTH peak and integrated flux, over the overlap band.

    The main position cache carries ``peak_flux`` only, which is all a detection census needs.
    A spectral index does not: alpha from peak fluxes across a 2.5" beam and a 25" one is
    biased steep, because the finer beam resolves out flux the coarser one keeps. Measuring
    that bias needs the integrated flux on both sides -- VLASS's is in the local component
    catalogs, RACS's is this column.

    Cached per 1-degree Dec strip exactly like ``fetch_racs_positions``; a completed strip is
    never re-fetched, so an interrupted run resumes.
    """
    import csv
    import time
    from pathlib import Path
    from urllib.parse import urlencode
    from urllib.request import urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    npz = dest / "racs_total.npz"
    if npz.exists():
        d = np.load(npz)
        return {k: d[k] for k in d}
    ra, dec, peak, total = [], [], [], []
    lo = dec_min
    while lo < dec_max:
        hi = min(lo + strip_deg, dec_max)
        path = dest / f"strip_{lo:+06.1f}.csv"
        if not path.exists():
            q = (
                f"SELECT ra, dec, peak_flux, total_flux_source FROM {RACS_TABLE} "
                f"WHERE dec >= {lo} AND dec < {hi}"
            )
            params = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
            for attempt in range(6):
                try:
                    with urlopen(f"{RACS_TAP_SYNC}?{params}", timeout=600) as fh:
                        path.write_bytes(fh.read())
                    break
                except Exception:  # noqa: BLE001 - retry transient TAP failures
                    if attempt == 5:
                        raise
                    time.sleep(5 * (attempt + 1))
        with path.open() as fh:
            for row in csv.DictReader(fh):
                try:
                    p_, t_ = float(row["peak_flux"]), float(row["total_flux_source"])
                except (KeyError, ValueError):
                    continue
                ra.append(float(row["ra"]))
                dec.append(float(row["dec"]))
                peak.append(p_)
                total.append(t_)
        lo = hi
    out = {
        "ra": np.array(ra),
        "dec": np.array(dec),
        "peak": np.array(peak),
        "total": np.array(total),
    }
    np.savez(npz, ra=out["ra"], dec=out["dec"], peak=out["peak"], total=out["total"])
    return out


def _alpha_samples(results_dir: str = "results") -> np.ndarray:  # pragma: no cover - real leg
    """The committed flux-complete alpha samples, or an empty array if run_alpha has not run.

    The census legs use these for the per-source-index comparison. Absence is not an error:
    run_alpha is a separate leg and a census run must still work before it has been executed.
    """
    import json
    from pathlib import Path

    path = Path(results_dir) / "dr20radio_alpha.json"
    if not path.exists():
        return np.array([])
    return np.asarray(json.loads(path.read_text()).get("alpha_samples_flux_complete", []), float)


def run_alpha(
    out: str = ".",
    *,
    vlass_radius_arcsec: float = 2.5,
    racs_radius_arcsec: float = 5.0,
    alpha_min_complete: float = -1.5,
    n_boot: int = 2000,
    seed: int = 0,
) -> dict:  # pragma: no cover - network + bulk local data (pure pieces tested offline)
    """Real leg C: MEASURE the spectral index instead of assuming it.

    The K-correction is the only thing making a 3 GHz and an 888 MHz survey comparable, and
    the north/south contrast in Section 4.3 turned out to be dominated by it: sweeping alpha
    over 0..-1 moved the gap from 0.23 to 1.66 percentage points. alpha need not be assumed.
    The overlap band (-40 < dec <= +30) is covered by BOTH surveys, so every quasar detected
    twice yields a two-point alpha, and the census's own median can replace the canonical
    -0.7.

    Two numbers come out, and the paper must quote both. The joint-detection sample is
    truncated -- a steep source at the RACS limit has already dropped below the VLASS limit at
    3 GHz -- so its median is biased flat. Above ``alpha_complete_limit_mjy`` that truncation
    cannot operate for any ``alpha >= alpha_min_complete``, giving a smaller unbiased sample.

    Writes ``results/dr20radio_alpha.json``.
    """
    import json
    from pathlib import Path

    rng = np.random.default_rng(seed)
    q = read_spall_quasars(fetch_spall())
    band = (q["dec"] > VLASS_DEC_LIMIT_DEG) & (q["dec"] <= 30.0) & ~q["radio_carton"]
    ra, dec = q["ra"][band], q["dec"][band]

    # VLASS: same construction as run_north -- brightest peak across the epoch catalogs. That
    # max-of-epochs is right for a detection census and WRONG for a flux ratio, where it is a
    # positively biased estimator (noise and real variability both push it up, and 3 GHz
    # variability is largest for compact flat-spectrum sources, so the bias is flatward).
    # Each epoch is therefore also kept separately, and the spread between them is reported.
    s_vlass = np.zeros(ra.size)
    m_vlass = np.zeros(ra.size, dtype=bool)
    per_epoch: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for _name, (ra_r, dec_r, fx_r) in load_vlass_positions().items():
        m, _, idx = crossmatch(ra, dec, ra_r, dec_r, radius_arcsec=vlass_radius_arcsec)
        s_vlass = np.where(m, np.maximum(s_vlass, fx_r[idx]), s_vlass)
        m_vlass |= m
        per_epoch[_name] = (m, np.where(m, fx_r[idx], 0.0))

    racs = fetch_racs_positions()
    m_racs, _, idx_r = crossmatch(
        ra, dec, racs["ra"], racs["dec"], radius_arcsec=racs_radius_arcsec
    )
    s_racs = np.where(m_racs, racs["flux"][idx_r], 0.0)

    both = m_vlass & m_racs & (s_vlass > 0) & (s_racs > 0)
    a = spectral_index(
        s_vlass[both], s_racs[both], freq_hi_ghz=VLASS_FREQ_GHZ, freq_lo_ghz=RACS_FREQ_GHZ
    )
    a = a[np.isfinite(a)]

    s_complete = alpha_complete_limit_mjy(
        s_lim_hi_mjy=VLASS_S_LIM_MJY,
        freq_hi_ghz=VLASS_FREQ_GHZ,
        freq_lo_ghz=RACS_FREQ_GHZ,
        alpha_min=alpha_min_complete,
    )
    unb = both & (s_racs > s_complete)
    a_unb = spectral_index(
        s_vlass[unb], s_racs[unb], freq_hi_ghz=VLASS_FREQ_GHZ, freq_lo_ghz=RACS_FREQ_GHZ
    )
    a_unb = a_unb[np.isfinite(a_unb)]

    def _block(x: np.ndarray, label: str) -> dict:
        if x.size == 0:
            return {"label": label, "n": 0}
        boot = np.array(
            [np.median(rng.choice(x, size=x.size, replace=True)) for _ in range(n_boot)]
        )
        return {
            "label": label,
            "n": int(x.size),
            "median": float(np.median(x)),
            "median_boot_se": float(np.std(boot)),
            "mean": float(np.mean(x)),
            "std": float(np.std(x)),
            "p16": float(np.percentile(x, 16)),
            "p84": float(np.percentile(x, 84)),
            "frac_steeper_than_canonical": float(np.mean(x < -0.7)),
        }

    def _alpha_for(mask: np.ndarray, s_hi: np.ndarray) -> np.ndarray:
        sel = mask & m_racs & (s_hi > 0) & (s_racs > 0)
        x = spectral_index(
            s_hi[sel], s_racs[sel], freq_hi_ghz=VLASS_FREQ_GHZ, freq_lo_ghz=RACS_FREQ_GHZ
        )
        return x[np.isfinite(x)]

    # (d) The estimator that needs no completeness cut at all. Every RACS-detected quasar
    # carries information about its index: if VLASS also saw it, alpha is measured; if not,
    # S_VLASS < the VLASS limit, so alpha is LEFT-CENSORED at
    # log(S_lim_VLASS / S_RACS) / log(nu_V / nu_R) -- "this source is steeper than X", which is
    # exactly the objects a completeness cut throws away. Kaplan-Meier uses all of them and is
    # unbiased over the whole alpha range instead of over the range the cut assumes.
    _lnr = np.log(VLASS_FREQ_GHZ / RACS_FREQ_GHZ)
    racs_sel = m_racs & (s_racs > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        a_obs = np.where(
            m_vlass & (s_vlass > 0),
            np.log(s_vlass / np.where(s_racs > 0, s_racs, np.nan)) / _lnr,
            np.nan,
        )
        a_lim = np.log(VLASS_S_LIM_MJY / np.where(s_racs > 0, s_racs, np.nan)) / _lnr
    is_cens = ~(m_vlass & (s_vlass > 0))
    a_km = np.where(is_cens, a_lim, a_obs)[racs_sel]
    km = kaplan_meier_median(a_km, is_cens[racs_sel], right_censored=False)
    _kmv, _kmc = a_km, is_cens[racs_sel]
    _kb = []
    for _ in range(n_boot // 4):  # KM is O(n log n) per draw; a quarter of the SE budget is ample
        j = rng.integers(0, _kmv.size, _kmv.size)
        r = kaplan_meier_median(_kmv[j], _kmc[j], right_censored=False)
        if r["median"] is not None:
            _kb.append(r["median"])
    km["median_boot_se"] = float(np.std(_kb)) if _kb else None
    km["n_boot"] = len(_kb)

    # (a) Does the completeness floor matter? alpha >= -1.5 is an ASSUMPTION, and the sample's
    # own p16 is steeper than it, so for the steepest sixth the truncation is still active.
    # Re-cutting at successively steeper floors is the version of this test that can fail: if
    # the median steepens monotonically, -0.72 is an upper bound on flatness, not a value.
    floors = {}
    for a_min in (-1.5, -2.0, -2.5):
        lim = alpha_complete_limit_mjy(
            s_lim_hi_mjy=VLASS_S_LIM_MJY,
            freq_hi_ghz=VLASS_FREQ_GHZ,
            freq_lo_ghz=RACS_FREQ_GHZ,
            alpha_min=a_min,
        )
        xa = _alpha_for(m_vlass & (s_racs > lim), s_vlass)
        floors[f"{a_min:g}"] = {
            "flux_cut_mjy": float(lim),
            "n": int(xa.size),
            "median": float(np.median(xa)) if xa.size else None,
        }

    # (b) Is the bright-sample index transferable to the threshold regime it is applied in?
    # Each bin is measured BOTH ways. The detections-only median is what a completeness cut
    # would report and is truncation-biased flat, worst in the faintest bin -- so a trend in
    # that column can be manufactured entirely by the censoring. The Kaplan-Meier column uses
    # the non-detections and is the one to read; the gap between the columns is the size of
    # the artefact.
    flux_bins = []
    for lo, hi in (
        (RACS_S_LIM_MJY, s_complete),
        (s_complete, 10.0),
        (10.0, 20.0),
        (20.0, float("inf")),
    ):
        binsel = racs_sel & (s_racs > lo) & (s_racs <= hi)
        xb = _alpha_for(m_vlass & (s_racs > lo) & (s_racs <= hi), s_vlass)
        kmb = kaplan_meier_median(
            np.where(is_cens, a_lim, a_obs)[binsel], is_cens[binsel], right_censored=False
        )
        flux_bins.append(
            {
                "s_racs_lo_mjy": float(lo),
                "s_racs_hi_mjy": None if not np.isfinite(hi) else float(hi),
                "n_detected": int(xb.size),
                "median_detected": float(np.median(xb)) if xb.size else None,
                "n_km": kmb["n"],
                "n_censored": kmb["n_censored"],
                "median_km": kmb["median"],
            }
        )

    # (c) How much of the index comes from taking the max over epochs?
    epochs_out = {}
    for name, (me, se) in per_epoch.items():
        xe = _alpha_for(me & (s_racs > s_complete), se)
        epochs_out[name] = {
            "n": int(xe.size),
            "median": float(np.median(xe)) if xe.size else None,
        }

    # (e) The resolution systematic, measured rather than estimated. alpha from PEAK fluxes
    # compares a 2.5" beam against a 25" one, and the finer beam resolves out flux the coarser
    # one keeps, biasing alpha steep. Repeating with integrated flux on both sides -- VLASS's
    # from the same component catalogs, RACS's from total_flux_source -- gives the size of it.
    integrated: dict = {"available": False}
    try:
        rt = fetch_racs_total_flux()
        sv_t = np.zeros(ra.size)
        mv_t = np.zeros(ra.size, dtype=bool)
        for _n, (rr, dd, ff) in load_vlass_positions(total_flux=True).items():
            mt, _, it = crossmatch(ra, dec, rr, dd, radius_arcsec=vlass_radius_arcsec)
            sv_t = np.where(mt, np.maximum(sv_t, ff[it]), sv_t)
            mv_t |= mt
        mr_t, _, ir_t = crossmatch(ra, dec, rt["ra"], rt["dec"], radius_arcsec=racs_radius_arcsec)
        sr_t = np.where(mr_t, rt["total"][ir_t], 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ao_t = np.where(
                mv_t & (sv_t > 0),
                np.log(sv_t / np.where(sr_t > 0, sr_t, np.nan)) / _lnr,
                np.nan,
            )
            al_t = np.log(VLASS_S_LIM_MJY / np.where(sr_t > 0, sr_t, np.nan)) / _lnr
        cens_t = ~(mv_t & (sv_t > 0))
        sel_t = mr_t & (sr_t > 0)
        km_t = kaplan_meier_median(
            np.where(cens_t, al_t, ao_t)[sel_t], cens_t[sel_t], right_censored=False
        )
        # Per-bin as well as global. The global shift is a median over a mostly-compact
        # sample and is guaranteed small; the question the paper actually leans on is whether
        # the steep BRIGHT bin is a spectrum or a beam, and only the per-bin version answers it.
        bins_t = []
        for lo, hi in (
            (RACS_S_LIM_MJY, s_complete),
            (s_complete, 10.0),
            (10.0, 20.0),
            (20.0, float("inf")),
        ):
            bs = sel_t & (sr_t > lo) & (sr_t <= hi)
            kb = kaplan_meier_median(
                np.where(cens_t, al_t, ao_t)[bs], cens_t[bs], right_censored=False
            )
            bins_t.append(
                {
                    "s_racs_lo_mjy": float(lo),
                    "s_racs_hi_mjy": None if not np.isfinite(hi) else float(hi),
                    "n_km": kb["n"],
                    "median_km": kb["median"],
                }
            )
        integrated = {
            "available": True,
            "kaplan_meier_median": km_t["median"],
            "n": km_t["n"],
            "n_censored": km_t["n_censored"],
            "shift_from_peak": float(km_t["median"] - km["median"]),
            "flux_bins": bins_t,
        }
    except Exception as e:  # noqa: BLE001 - the integrated cache is optional, absence is a result
        integrated = {"available": False, "error": str(e)}

    # (f) The mirror censoring. The Kaplan-Meier above conditions on a RACS detection, so it
    # keeps the steep sources a completeness cut loses -- but it drops the quasars VLASS saw
    # and RACS did not, which are the FLAT ones, right-censored at
    # log(S_VLASS / S_lim_RACS)/log(nu_V/nu_R). Dropping that tail tilts the median steep by
    # exactly the mechanism that tilts the joint-detection median flat. With both tails
    # censored no point estimate is identified, but the median is bounded, and the bound is
    # the honest object.
    with np.errstate(divide="ignore", invalid="ignore"):
        a_rlim = np.log(np.where(s_vlass > 0, s_vlass, np.nan) / RACS_S_LIM_MJY) / _lnr
    either = (m_racs & (s_racs > 0)) | (m_vlass & (s_vlass > 0))
    both_m = m_racs & (s_racs > 0) & m_vlass & (s_vlass > 0)
    kind = np.where(both_m, 0, np.where(m_racs & (s_racs > 0), -1, 1))
    val = np.where(both_m, a_obs, np.where(m_racs & (s_racs > 0), a_lim, a_rlim))
    bounds = censored_median_bounds(val[either], kind[either])

    metrics = {
        "source": f"SDSS-V DR20 spAll-lite x VLASS x RACS-low DR1, overlap band ({RACS_TABLE})",
        "kaplan_meier": km,
        "double_censored_bounds": bounds,
        "integrated_flux": integrated,
        "completeness_floor_sensitivity": floors,
        "flux_bins": flux_bins,
        "per_epoch": epochs_out,
        "band": "-40 < dec <= +30 (both surveys cover it)",
        "n_band_census": int(band.sum()),
        "n_vlass_matched": int(m_vlass.sum()),
        "n_racs_matched": int(m_racs.sum()),
        "n_both": int(both.sum()),
        "vlass_radius_arcsec": vlass_radius_arcsec,
        "racs_radius_arcsec": racs_radius_arcsec,
        "freq_hi_ghz": VLASS_FREQ_GHZ,
        "freq_lo_ghz": RACS_FREQ_GHZ,
        "joint_detection": _block(a, "all joint detections (truncation-biased flat)"),
        "flux_complete": _block(
            a_unb, f"S_RACS > {s_complete:.2f} mJy (unbiased for alpha >= {alpha_min_complete:g})"
        ),
        # The samples themselves, so the per-source-alpha comparison in the census legs draws
        # from committed evidence rather than from a re-run.
        "alpha_samples_flux_complete": [round(float(x), 4) for x in a_unb],
        "completeness_flux_mjy": s_complete,
        "alpha_min_complete": alpha_min_complete,
        "canonical_alpha": -0.7,
        "n_boot": n_boot,
        "seed": seed,
    }
    Path(out, "results").mkdir(parents=True, exist_ok=True)
    Path(out, "results", "dr20radio_alpha.json").write_text(json.dumps(metrics, indent=1))
    return metrics


def run_south(
    out: str = ".", *, radius_arcsec: float = 5.0, n_shift_trials: int = 10
) -> dict:  # pragma: no cover - network + bulk data (pure pieces tested offline)
    """Real leg B: the categorical first — DR20 quasars south of -40 deg vs RACS-low DR1.

    Also computes the overlap band (-40..+30) for the VLASS cross-check, and validates the
    racsradio carton against its SELECTING survey (expected ~100%). RACS-low's 25" beam
    motivates the wider default match radius (5"), with the false-match rate measured as
    always. Writes ``results/dr20radio_south.json``.
    """
    import json
    from pathlib import Path

    spall = fetch_spall()
    q = read_spall_quasars(spall)
    racs = fetch_racs_positions()
    _asamp = _alpha_samples()
    zbins = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0])

    def census_block(sel: np.ndarray, label: str) -> dict:
        cen = sel & ~q["radio_carton"]
        m, _, idx = crossmatch(
            q["ra"][cen], q["dec"][cen], racs["ra"], racs["dec"], radius_arcsec=radius_arcsec
        )
        s_m = np.where(m, racs["flux"][idx], 0.0)
        fm = false_match_rate(
            q["ra"][cen],
            q["dec"][cen],
            racs["ra"],
            racs["dec"],
            radius_arcsec=radius_arcsec,
            n_trials=n_shift_trials,
        )
        p, lo_w, hi_w = wilson_interval(int(m.sum()), int(m.size))
        return {
            "label": label,
            "n_census": int(m.size),
            "n_matched": int(m.sum()),
            "raw_fraction": round(p, 5),
            "wilson_lo": round(lo_w, 5),
            "wilson_hi": round(hi_w, 5),
            "false_match": fm,
            "corrected_fraction": round(p - fm["rate"], 5),
            "fraction_vs_z": detection_fraction(q["z"][cen], m, bins=zbins),
            "luminosity_matched": luminosity_matched_fractions(
                q["z"][cen],
                m,
                s_m,
                freq_ghz=RACS_FREQ_GHZ,
                s_lim_this_mjy=RACS_S_LIM_MJY,
                s_lim_other_mjy=VLASS_S_LIM_MJY,
                other_freq_ghz=VLASS_FREQ_GHZ,
                bins=zbins,
            ),
            "luminosity_matched_vlass_conservative": {
                f"{v:g}": luminosity_matched_fractions(
                    q["z"][cen],
                    m,
                    s_m,
                    freq_ghz=RACS_FREQ_GHZ,
                    s_lim_this_mjy=RACS_S_LIM_MJY,
                    s_lim_other_mjy=v,
                    other_freq_ghz=VLASS_FREQ_GHZ,
                    bins=zbins,
                    alpha=ALPHA_MEASURED,
                )
                for v in VLASS_S_LIM_CONSERVATIVE_MJY
            },
            "luminosity_matched_per_source_alpha": (
                luminosity_matched_per_source_alpha(
                    q["z"][cen],
                    m,
                    s_m,
                    freq_ghz=RACS_FREQ_GHZ,
                    s_lim_this_mjy=RACS_S_LIM_MJY,
                    s_lim_other_mjy=VLASS_S_LIM_MJY,
                    other_freq_ghz=VLASS_FREQ_GHZ,
                    bins=zbins,
                    alpha_samples=_asamp,
                )
                if _asamp.size
                else None
            ),
            "luminosity_matched_alpha": {
                f"{a:g}": luminosity_matched_fractions(
                    q["z"][cen],
                    m,
                    s_m,
                    freq_ghz=RACS_FREQ_GHZ,
                    s_lim_this_mjy=RACS_S_LIM_MJY,
                    s_lim_other_mjy=VLASS_S_LIM_MJY,
                    other_freq_ghz=VLASS_FREQ_GHZ,
                    bins=zbins,
                    alpha=a,
                )
                for a in (*ALPHA_SWEEP, *ALPHA_MEASURED_SWEEP)
            },
            "luminosity_matched_conservative": luminosity_matched_fractions(
                q["z"][cen],
                m,
                s_m,
                freq_ghz=RACS_FREQ_GHZ,
                s_lim_this_mjy=RACS_S_LIM_CONSERVATIVE_MJY,
                s_lim_other_mjy=VLASS_S_LIM_MJY,
                other_freq_ghz=VLASS_FREQ_GHZ,
                bins=zbins,
            ),
            "obs_breakdown": {o: int((q["obs"][cen] == o).sum()) for o in np.unique(q["obs"][cen])},
        }

    deep_south = q["dec"] <= VLASS_DEC_LIMIT_DEG
    overlap = (q["dec"] > VLASS_DEC_LIMIT_DEG) & (q["dec"] <= 30.0)

    # Carton validation split by SELECTING survey: racsradio cartons vs RACS is the true
    # ~100% pipeline validation; lofarradio (144 MHz-selected) vs RACS is a cross-frequency
    # fraction, exactly like the VLASS case in increment 1.
    def carton_block(mask: np.ndarray) -> dict:
        mm, _, _ = crossmatch(
            q["ra"][mask], q["dec"][mask], racs["ra"], racs["dec"], radius_arcsec=radius_arcsec
        )
        return {
            "n": int(mm.size),
            "matched": int(mm.sum()),
            "fraction": round(float(np.mean(mm)), 4) if mm.size else None,
        }

    # Footprint diagnostic (added 2026-08-12 after a referee asked what the denominator
    # includes). `deep_south` is a pure declination cut, so quasars in RACS-low DR1's
    # Galactic-plane hole (|b| <~ 5) or below its dec floor sit in the denominator of the
    # headline fraction as guaranteed non-detections. Reporting the count converts a silent
    # one-sided bias into a stated — and, as it turns out, negligible — bound.
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    # Same mask census_block applies, so this describes the census's own denominator rather
    # than the raw declination cut (they differ by the 10 excluded radio-carton objects).
    _cen = deep_south & ~q["radio_carton"]
    _b = SkyCoord(q["ra"][_cen] * u.deg, q["dec"][_cen] * u.deg).galactic.b.deg
    _uncovered = (np.abs(_b) < RACS_PLANE_CUT_DEG) | (q["dec"][_cen] < RACS_DEC_FLOOR_DEG)
    footprint = {
        "n_deep_south": int(_cen.sum()),
        "n_in_plane_hole": int((np.abs(_b) < RACS_PLANE_CUT_DEG).sum()),
        "n_below_dec_floor": int((q["dec"][_cen] < RACS_DEC_FLOOR_DEG).sum()),
        "n_uncovered": int(_uncovered.sum()),
        "uncovered_fraction": round(float(_uncovered.mean()), 6),
        "plane_cut_deg": RACS_PLANE_CUT_DEG,
        "dec_floor_deg": RACS_DEC_FLOOR_DEG,
    }

    in_racs_sky = q["dec"] <= 30.0
    carton_racs_v = carton_block(q["carton_racs"] & in_racs_sky)
    carton_lofar_v = carton_block(q["carton_lofar"] & in_racs_sky)
    metrics = {
        "source": f"SDSS-V DR20 spAll-lite x RACS-low DR1 ({RACS_TABLE})",
        "n_racs_sources": int(racs["ra"].size),
        "radius_arcsec": radius_arcsec,
        "deep_south": census_block(deep_south, "dec <= -40 (SDSS x RACS: categorical first)"),
        "racs_footprint": footprint,
        "overlap_band": census_block(overlap, "-40 < dec <= +30 (VLASS cross-check band)"),
        "carton_validation": {
            "racsradio_vs_racs_selecting_survey": carton_racs_v,
            "lofarradio_vs_racs_cross_frequency": carton_lofar_v,
        },
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "dr20radio_south.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def paper_assets(out: str = ".", *, results_dir: str = "results") -> None:  # pragma: no cover
    """Figures + macros for papers/dr20radio from the COMMITTED evidence JSONs only.

    The committed-real-results rule: no census runs here — absent evidence fails loudly.
    """
    import json
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = json.loads((Path(results_dir) / "dr20radio_north.json").read_text())
    s = json.loads((Path(results_dir) / "dr20radio_south.json").read_text())
    fig_dir = Path(out) / "papers" / "dr20radio" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    def _binned(block: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        mid = 0.5 * (np.asarray(block["bin_lo"]) + np.asarray(block["bin_hi"]))
        kk, nn = np.asarray(block["k"], float), np.asarray(block["n"], float)
        lo = np.array([wilson_interval(int(k), int(m))[1] for k, m in zip(kk, nn, strict=True)])
        hi = np.array([wilson_interval(int(k), int(m))[2] for k, m in zip(kk, nn, strict=True)])
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(nn > 0, kk / nn, np.nan)
        return mid, frac, lo, hi

    def _plot(ax: plt.Axes, block: dict, label: str, color: str, marker: str) -> None:
        mid, frac, lo, hi = _binned(block)
        ax.errorbar(
            mid,
            frac,
            yerr=[frac - lo, hi - frac],
            fmt=marker,
            color=color,
            label=label,
            capsize=2,
            ms=5,
        )

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    _plot(axes[0], n["fraction_vs_z_any_epoch"], "North: VLASS any-epoch (3 GHz)", "C0", "o")
    _plot(
        axes[0],
        s["deep_south"]["fraction_vs_z"],
        "South: RACS, $\\delta \\leq -40^\\circ$",
        "C3",
        "s",
    )
    _plot(axes[0], s["overlap_band"]["fraction_vs_z"], "Overlap band: RACS", "0.6", "^")
    axes[0].set_title("Raw detection fraction (survey-native depths)")
    _plot(axes[1], n["luminosity_matched"], "North above common $L$ limit", "C0", "o")
    _plot(axes[1], s["deep_south"]["luminosity_matched"], "South above common $L$ limit", "C3", "s")
    axes[1].set_title("Above the common 1.4 GHz luminosity limit")
    for ax in axes:
        ax.set_xlabel("redshift $z$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("radio-detected fraction")
    fig.tight_layout()
    fig.savefig(fig_dir / "dr20radio_fractions.pdf")
    plt.close(fig)

    ds, ov, ae = s["deep_south"], s["overlap_band"], n["any_epoch"]
    cv = s["carton_validation"]
    fp = s["racs_footprint"]
    al = json.loads((Path(results_dir) / "dr20radio_alpha.json").read_text())
    # Sweep keys for the measured alpha and its +/-1 SE, formatted exactly as the runs wrote
    # them. Derived from the same constants the runs used, so a drift in either shows up as a
    # KeyError here rather than as a silently wrong number in the paper.
    _am, _am_lo, _am_hi = (f"{a:g}" for a in (ALPHA_MEASURED, *sorted(ALPHA_MEASURED_SWEEP)[::2]))
    e2, e3 = n["epochs"]["E2"], n["epochs"]["E3"]
    lum_n = n["luminosity_matched"]
    lum_s = ds["luminosity_matched"]

    def _tot(b: dict) -> float:
        return float(sum(b["k"])) / float(sum(b["n"]))

    _gap_meas = 100 * (
        _tot(n["luminosity_matched_alpha"][_am]) - _tot(ds["luminosity_matched_alpha"][_am])
    )
    _gap_vc = [
        100
        * (
            _tot(n["luminosity_matched_vlass_conservative"][f"{v:g}"])
            - _tot(ds["luminosity_matched_vlass_conservative"][f"{v:g}"])
        )
        for v in VLASS_S_LIM_CONSERVATIVE_MJY
    ]
    # Span of the median across every systematic actually tested: completeness floor, epoch
    # choice, and flux bin. This is the honest scale of the uncertainty on alpha.
    _rat = lambda a, b: _tot(a) / _tot(b)  # noqa: E731
    _ratio_meas = _rat(n["luminosity_matched_alpha"][_am], ds["luminosity_matched_alpha"][_am])
    _ratios_lim = [
        _ratio_meas,
        _rat(n["luminosity_matched_conservative"], ds["luminosity_matched_conservative"]),
        *(
            _rat(
                n["luminosity_matched_vlass_conservative"][f"{v:g}"],
                ds["luminosity_matched_vlass_conservative"][f"{v:g}"],
            )
            for v in VLASS_S_LIM_CONSERVATIVE_MJY
        ),
    ]
    _ratios_alpha = [
        _rat(n["luminosity_matched_alpha"][k], ds["luminosity_matched_alpha"][k])
        for k in (f"{ALPHA_THRESHOLD_REGIME:g}", _am, "-1")
    ]
    _gap_faint = 100 * (
        _tot(n["luminosity_matched_alpha"][f"{ALPHA_THRESHOLD_REGIME:g}"])
        - _tot(ds["luminosity_matched_alpha"][f"{ALPHA_THRESHOLD_REGIME:g}"])
    )
    # How much of the raw flux trend was censoring rather than physics: the detections-only
    # column minus the Kaplan-Meier column, worst in the faintest bin.
    _bin_art = max(b["median_detected"] - b["median_km"] for b in al["flux_bins"])
    _alpha_medians = [
        al["kaplan_meier"]["median"],
        *(v["median"] for v in al["completeness_floor_sensitivity"].values()),
        *(v["median"] for v in al["per_epoch"].values()),
        *(b["median_km"] for b in al["flux_bins"]),
        al["integrated_flux"]["kaplan_meier_median"],
    ]
    _alpha_sys = max(_alpha_medians) - min(_alpha_medians)
    _ps_n = n["luminosity_matched_per_source_alpha"]
    _ps_s = ds["luminosity_matched_per_source_alpha"]
    _gap_se = [
        100 * (_tot(n["luminosity_matched_alpha"][k]) - _tot(ds["luminosity_matched_alpha"][k]))
        for k in (_am_lo, _am_hi)
    ]
    lines = [
        "% Auto-generated by jansky_research.dr20radio.paper_assets from the committed",
        "% results/dr20radio_*.json evidence — do not edit by hand.",
        rf"\newcommand{{\drNorthCensus}}{{{n['n_north_census']:,}}}".replace(",", r"{,}"),
        rf"\newcommand{{\drSouthCensus}}{{{ds['n_census']:,}}}".replace(",", r"{,}"),
        rf"\newcommand{{\drOverlapCensus}}{{{ov['n_census']:,}}}".replace(",", r"{,}"),
        rf"\newcommand{{\drNorthAnyPct}}{{{100 * ae['raw_fraction']:.2f}}}",
        rf"\newcommand{{\drNorthEtwoPct}}{{{100 * e2['corrected_fraction']:.2f}}}",
        rf"\newcommand{{\drNorthEthreePct}}{{{100 * e3['corrected_fraction']:.2f}}}",
        rf"\newcommand{{\drSouthPct}}{{{100 * ds['corrected_fraction']:.2f}}}",
        rf"\newcommand{{\drSouthMatched}}{{{ds['n_matched']:,}}}".replace(",", r"{,}"),
        rf"\newcommand{{\drOverlapPct}}{{{100 * ov['corrected_fraction']:.2f}}}",
        rf"\newcommand{{\drSouthFmPct}}{{{100 * ds['false_match']['rate']:.2f}}}",
        rf"\newcommand{{\drNorthFmPct}}{{{100 * e2['false_match']['rate']:.3f}}}",
        rf"\newcommand{{\drLumNorthPct}}{{{100 * _tot(lum_n):.2f}}}",
        rf"\newcommand{{\drLumSouthPct}}{{{100 * _tot(lum_s):.2f}}}",
        rf"\newcommand{{\drOverlapLumPct}}{{{100 * _tot(ov['luminosity_matched']):.2f}}}",
        # The contrast as a RANGE over spectral index. alpha is not measured for this sample
        # and the gap is dominated by it: 0.23 pp at alpha=0 against 1.66 pp at alpha=-1. The
        # 5 mJy variant cannot show this because it rescales both legs identically.
        # Two figures that were hard-typed in the paper until 2026-08-12, both derivable from
        # this same committed evidence. Out-of-range quasars fall outside the 0 <= z < 6 binning
        # and so enter the raw fractions but not the binned/luminosity-matched ones. The northern
        # census above dec +30 lies outside the RACS footprint entirely, which is why the
        # north-vs-overlap comparison is not a matched-sky one.
        rf"\newcommand{{\drOutZNorthPct}}{{{100 * (1 - sum(n['fraction_vs_z_any_epoch']['n']) / n['n_north_census']):.2f}}}",
        rf"\newcommand{{\drOutZSouthPct}}{{{100 * (1 - sum(ds['fraction_vs_z']['n']) / ds['n_census']):.2f}}}",
        rf"\newcommand{{\drNorthOutsideRacsPct}}{{{100 * (1 - ov['n_census'] / n['n_north_census']):.1f}}}",
        rf"\newcommand{{\drAlphaLo}}{{{min(ALPHA_SWEEP):g}}}",
        rf"\newcommand{{\drAlphaHi}}{{{max(ALPHA_SWEEP):g}}}",
        rf"\newcommand{{\drLumNorthFlatPct}}{{{100 * _tot(n['luminosity_matched_alpha']['0']):.2f}}}",
        rf"\newcommand{{\drLumSouthFlatPct}}{{{100 * _tot(ds['luminosity_matched_alpha']['0']):.2f}}}",
        rf"\newcommand{{\drGapFlatPp}}{{{100 * (_tot(n['luminosity_matched_alpha']['0']) - _tot(ds['luminosity_matched_alpha']['0'])):.2f}}}",
        rf"\newcommand{{\drGapSteepPp}}{{{100 * (_tot(n['luminosity_matched_alpha']['-1']) - _tot(ds['luminosity_matched_alpha']['-1'])):.2f}}}",
        rf"\newcommand{{\drGapFidPp}}{{{100 * (_tot(lum_n) - _tot(lum_s)):.2f}}}",
        # The steep end of the sweep. Emitted because the first revision reached for
        # \drLumNorthConsPct (the 5 mJy variant, 3.45%) to close the range "3.06--...", which is a
        # different axis entirely and understated the sweep by a full percentage point. A referee
        # caught it. If a range needs an endpoint, the endpoint gets its own macro.
        # The MEASURED spectral index and the contrast evaluated at it. This is what replaces
        # the sweep as the paper's statement: the sweep bounded how bad the assumption could
        # be, these numbers say what the assumption actually is.
        # KM is the primary estimator; the flux-complete cut is kept for comparison because the
        # difference between them IS the bias a cut leaves behind.
        rf"\newcommand{{\drAlphaMeas}}{{{al['kaplan_meier']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaMeasSe}}{{{al['kaplan_meier']['median_boot_se']:.2f}}}",
        rf"\newcommand{{\drAlphaMeasN}}{{{al['kaplan_meier']['n']}}}",
        rf"\newcommand{{\drAlphaKmCens}}{{{al['kaplan_meier']['n_censored']}}}",
        rf"\newcommand{{\drAlphaKmEvents}}{{{al['kaplan_meier']['n_events']}}}",
        rf"\newcommand{{\drAlphaCutMed}}{{{al['flux_complete']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaCutN}}{{{al['flux_complete']['n']}}}",
        rf"\newcommand{{\drAlphaFaint}}{{{ALPHA_THRESHOLD_REGIME:.2f}}}",
        rf"\newcommand{{\drGapFaintPp}}{{{_gap_faint:.2f}}}",
        rf"\newcommand{{\drAlphaBinArt}}{{{_bin_art:.2f}}}",
        # The RATIO, which is the scale-free contrast. Raising either survey's limit rescales
        # both fractions and shrinks the pp gap without touching the ratio -- so a pp gap that
        # moves under a deeper cut is measuring normalisation, not contrast. Reporting the
        # ratio is what makes the limit variants interpretable instead of misleading.
        rf"\newcommand{{\drRatioMeas}}{{{_ratio_meas:.2f}}}",
        rf"\newcommand{{\drRatioLimLo}}{{{min(_ratios_lim):.2f}}}",
        rf"\newcommand{{\drRatioLimHi}}{{{max(_ratios_lim):.2f}}}",
        rf"\newcommand{{\drRatioLimSpreadPct}}{{{100 * (max(_ratios_lim) / min(_ratios_lim) - 1):.0f}}}",
        rf"\newcommand{{\drRatioAlphaFlat}}{{{_rat(n['luminosity_matched_alpha']['0'], ds['luminosity_matched_alpha']['0']):.2f}}}",
        rf"\newcommand{{\drRatioAlphaLo}}{{{min(_ratios_alpha):.2f}}}",
        rf"\newcommand{{\drRatioAlphaHi}}{{{max(_ratios_alpha):.2f}}}",
        # The selection bracket: alpha depends on WHICH survey selects the sample, and that
        # spread is larger than every other systematic term combined.
        rf"\newcommand{{\drAlphaBoundLo}}{{{al['double_censored_bounds']['lo']:.2f}}}",
        rf"\newcommand{{\drAlphaBoundHi}}{{{al['double_censored_bounds']['hi']:.2f}}}",
        rf"\newcommand{{\drAlphaBoundN}}{{{al['double_censored_bounds']['n']}}}",
        rf"\newcommand{{\drAlphaRightCens}}{{{al['double_censored_bounds']['n_right_censored']}}}",
        rf"\newcommand{{\drAlphaBinCInt}}{{{al['integrated_flux']['flux_bins'][3]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaBinCShift}}{{{al['integrated_flux']['flux_bins'][3]['median_km'] - al['flux_bins'][3]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaIntKm}}{{{al['integrated_flux']['kaplan_meier_median']:.2f}}}",
        rf"\newcommand{{\drAlphaIntShift}}{{{al['integrated_flux']['shift_from_peak']:.3f}}}",
        rf"\newcommand{{\drAlphaMeasSd}}{{{al['flux_complete']['std']:.2f}}}",
        rf"\newcommand{{\drAlphaMeasPlo}}{{{al['flux_complete']['p16']:.2f}}}",
        rf"\newcommand{{\drAlphaMeasPhi}}{{{al['flux_complete']['p84']:.2f}}}",
        rf"\newcommand{{\drAlphaJointMed}}{{{al['joint_detection']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaJointN}}{{{al['joint_detection']['n']}}}",
        rf"\newcommand{{\drAlphaComplMjy}}{{{al['completeness_flux_mjy']:.1f}}}",
        rf"\newcommand{{\drAlphaBandN}}{{{al['n_band_census']}}}",
        rf"\newcommand{{\drAlphaVlassRad}}{{{al['vlass_radius_arcsec']:g}}}",
        rf"\newcommand{{\drAlphaRacsRad}}{{{al['racs_radius_arcsec']:g}}}",
        rf"\newcommand{{\drAlphaVlassBeam}}{{{VLASS_BEAM_ARCSEC:g}}}",
        rf"\newcommand{{\drAlphaRacsBeam}}{{{RACS_BEAM_ARCSEC:g}}}",
        rf"\newcommand{{\drAlphaTruncPct}}{{{100 * (1 - al['flux_complete']['n'] / al['joint_detection']['n']):.0f}}}",
        rf"\newcommand{{\drLumNorthMeasPct}}{{{100 * _tot(n['luminosity_matched_alpha'][_am]):.2f}}}",
        rf"\newcommand{{\drLumSouthMeasPct}}{{{100 * _tot(ds['luminosity_matched_alpha'][_am]):.2f}}}",
        rf"\newcommand{{\drGapMeasPp}}{{{100 * (_tot(n['luminosity_matched_alpha'][_am]) - _tot(ds['luminosity_matched_alpha'][_am])):.2f}}}",
        # min/max, not lo-alpha/hi-alpha: a steeper alpha widens the gap, so keying these to
        # the alpha endpoints put the larger number in the macro named "Lo".
        rf"\newcommand{{\drGapMeasLoPp}}{{{min(_gap_se):.2f}}}",
        rf"\newcommand{{\drGapMeasHiPp}}{{{max(_gap_se):.2f}}}",
        # The population-scatter check: each quasar given its OWN index drawn from the measured
        # distribution, instead of one median index for everyone. It shifts both fractions by
        # ~0.3 pp and the GAP by only 0.03 -- so the contrast survives the scatter but the
        # absolute fractions do not, and the paper has to say both.
        # The VLASS-limit mirror sweep, and the systematic budget on the measured index. The
        # bootstrap SE (0.015) is the SMALLEST of these terms by an order of magnitude; quoting
        # it alone as the uncertainty was the same mistake in a new place.
        rf"\newcommand{{\drNorthEffCut}}{{{RACS_S_LIM_MJY * (VLASS_FREQ_GHZ / RACS_FREQ_GHZ) ** ALPHA_MEASURED:.2f}}}",
        rf"\newcommand{{\drSlimVlassConsA}}{{{VLASS_S_LIM_CONSERVATIVE_MJY[0]:g}}}",
        rf"\newcommand{{\drSlimVlassConsB}}{{{VLASS_S_LIM_CONSERVATIVE_MJY[1]:g}}}",
        rf"\newcommand{{\drGapVlassConsA}}{{{_gap_vc[0]:.2f}}}",
        rf"\newcommand{{\drGapVlassConsB}}{{{_gap_vc[1]:.2f}}}",
        rf"\newcommand{{\drGapVlassDropPct}}{{{100 * (1 - _gap_vc[1] / _gap_meas):.0f}}}",
        rf"\newcommand{{\drAlphaFloorB}}{{{al['completeness_floor_sensitivity']['-2']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaFloorC}}{{{al['completeness_floor_sensitivity']['-2.5']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaFloorBCut}}{{{al['completeness_floor_sensitivity']['-2']['flux_cut_mjy']:.1f}}}",
        rf"\newcommand{{\drAlphaFloorCCut}}{{{al['completeness_floor_sensitivity']['-2.5']['flux_cut_mjy']:.1f}}}",
        rf"\newcommand{{\drAlphaEpochA}}{{{al['per_epoch']['E2']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaEpochB}}{{{al['per_epoch']['E3']['median']:.2f}}}",
        rf"\newcommand{{\drAlphaBinFaint}}{{{al['flux_bins'][0]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaBinA}}{{{al['flux_bins'][1]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaBinB}}{{{al['flux_bins'][2]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaBinC}}{{{al['flux_bins'][3]['median_km']:.2f}}}",
        rf"\newcommand{{\drAlphaSysSpan}}{{{_alpha_sys:.2f}}}",
        # Endpoints, not a half-width. Writing the span as "+/- 0.35" asserts twice the range
        # the checks actually cover, and contradicted the paper's own "-0.6 to -0.9".
        rf"\newcommand{{\drAlphaSysLo}}{{{min(_alpha_medians):.2f}}}",
        rf"\newcommand{{\drAlphaSysHi}}{{{max(_alpha_medians):.2f}}}",
        rf"\newcommand{{\drAlphaMeasSeFull}}{{{al['kaplan_meier']['median_boot_se']:.3f}}}",
        rf"\newcommand{{\drAlphaMeasFull}}{{{al['kaplan_meier']['median']:.3f}}}",
        rf"\newcommand{{\drAlphaComplFracPct}}{{{100 * al['flux_complete']['n'] / al['n_band_census']:.1f}}}",
        rf"\newcommand{{\drLumNorthPsPct}}{{{100 * _ps_n['fraction_mean']:.2f}}}",
        rf"\newcommand{{\drLumSouthPsPct}}{{{100 * _ps_s['fraction_mean']:.2f}}}",
        rf"\newcommand{{\drGapPsPp}}{{{100 * (_ps_n['fraction_mean'] - _ps_s['fraction_mean']):.2f}}}",
        rf"\newcommand{{\drGapPsShiftPp}}{{{100 * ((_ps_n['fraction_mean'] - _ps_s['fraction_mean']) - (_tot(n['luminosity_matched_alpha'][_am]) - _tot(ds['luminosity_matched_alpha'][_am]))):.2f}}}",
        rf"\newcommand{{\drPsRealSd}}{{{100 * max(_ps_n['fraction_realization_sd'], _ps_s['fraction_realization_sd']):.3f}}}",
        rf"\newcommand{{\drLumNorthSteepPct}}{{{100 * _tot(n['luminosity_matched_alpha']['-1']):.2f}}}",
        rf"\newcommand{{\drLumSouthSteepPct}}{{{100 * _tot(ds['luminosity_matched_alpha']['-1']):.2f}}}",
        # The two ratios the retired 5 mJy check moves between. Quoted so the "it cannot test the
        # asymmetry" argument is checkable rather than asserted.
        rf"\newcommand{{\drRatioFid}}{{{_tot(lum_n) / _tot(lum_s):.4f}}}",
        # What fraction of the southern census sits in sky RACS-low DR1 never covered. A referee
        # asked; the answer is small, but "small" had to be measured rather than assumed.
        rf"\newcommand{{\drSouthUncovPct}}{{{100 * fp['uncovered_fraction']:.3f}}}",
        rf"\newcommand{{\drSouthUncovN}}{{{fp['n_uncovered']}}}",
        rf"\newcommand{{\drRacsPlaneCutDeg}}{{{fp['plane_cut_deg']:g}}}",
        rf"\newcommand{{\drRatioCons}}{{{_tot(n['luminosity_matched_conservative']) / _tot(ds['luminosity_matched_conservative']):.4f}}}",
        rf"\newcommand{{\drLumNorthConsPct}}{{{100 * _tot(n['luminosity_matched_conservative']):.2f}}}",
        rf"\newcommand{{\drLumSouthConsPct}}{{{100 * _tot(ds['luminosity_matched_conservative']):.2f}}}",
        rf"\newcommand{{\drCartonRacsPct}}{{{100 * cv['racsradio_vs_racs_selecting_survey']['fraction']:.0f}}}",
        rf"\newcommand{{\drCartonRacsN}}{{{cv['racsradio_vs_racs_selecting_survey']['n']}}}",
        rf"\newcommand{{\drCartonLofarPct}}{{{100 * cv['lofarradio_vs_racs_cross_frequency']['fraction']:.0f}}}",
        rf"\newcommand{{\drCartonLofarN}}{{{cv['lofarradio_vs_racs_cross_frequency']['n']}}}",
        rf"\newcommand{{\drCartonVlassRacsPct}}{{{100 * e2['carton_validation']['racsradio_cross_frequency']['fraction']:.0f}}}",
        rf"\newcommand{{\drCartonVlassLofarPct}}{{{100 * e2['carton_validation']['lofarradio_cross_frequency']['fraction']:.0f}}}",
        rf"\newcommand{{\drSlimRacsCons}}{{{RACS_S_LIM_CONSERVATIVE_MJY:.1f}}}",
        rf"\newcommand{{\drRacsSrc}}{{{s['n_racs_sources']:,}}}".replace(",", r"{,}"),
        rf"\newcommand{{\drSlimVlass}}{{{VLASS_S_LIM_MJY:.1f}}}",
        rf"\newcommand{{\drSlimRacs}}{{{RACS_S_LIM_MJY:.1f}}}",
    ]
    gen = Path(out) / "papers" / "dr20radio" / "generated"
    gen.mkdir(parents=True, exist_ok=True)
    (gen / "macros.tex").write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="DR20 BHM x VLASS/RACS radio census (plan 88).")
    p.add_argument("--out", default=".")
    p.add_argument("--north", action="store_true", help="run the VLASS northern leg")
    p.add_argument("--south", action="store_true", help="run the RACS southern leg")
    p.add_argument("--paper", action="store_true", help="figures + macros from committed evidence")
    args = p.parse_args(argv)
    if args.paper:
        paper_assets(args.out)
        return 0
    if args.south:
        m = run_south(args.out)
        slim = {
            k: (
                {kk: vv for kk, vv in v.items() if kk != "fraction_vs_z"}
                if isinstance(v, dict) and "fraction_vs_z" in v
                else v
            )
            for k, v in m.items()
        }
        print(json.dumps(slim, indent=2))
        return 0
    if args.north:
        m = run_north(args.out)
        slim = {k: v for k, v in m.items() if k != "fraction_vs_z_any_epoch"}
        print(json.dumps(slim, indent=2))
        return 0
    p.error("choose a mode: --north or --south")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
