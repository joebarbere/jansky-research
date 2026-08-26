"""Galactic RM structure functions from SPICE-RACS --- the `rmsky` slice at 10x the data (plan 36).

SPICE-RACS DR2 (arXiv:2605.16917; ~2.5--3.4x10^5 RMs over 87.5% of sky, 6.7 deg^-2) was released
without a systematic Galactic structure-function analysis by latitude. This module supplies the
tooling: the second-order **RM structure function** :math:`\\mathrm{SF}(\\delta\\theta) =
\\langle[\\mathrm{RM}(\\hat n_1)-\\mathrm{RM}(\\hat n_2)]^2\\rangle` in angular-separation bins,
**noise-debiased** by subtracting :math:`\\langle\\sigma_1^2+\\sigma_2^2\\rangle` per bin (the
measurement-error term that otherwise floors the SF), with bootstrap errors --- computed per
Galactic-latitude bin so the turbulence amplitude and any coherence-scale break can be compared
across the disc--halo transition. The latitude profile / quadrant machinery reuses `rmsky`.

GATE 0 (2026-07-02): DR2's catalogue is public on the CSIRO DAP (collection csiro:64891,
`spice-racs.dr2.fits.gz`, 4.97 GB, no auth); the verified bounded first leg is DR1 on CASDA TAP ---
`AS110.spice_racs_dr1_corrected_cut_v02`, **24,758 rows** (live count) with `l, b, rm, rm_err,
snr_polint` columns. Offline, a synthetic Gaussian RM screen with a KNOWN coherence scale and a
known latitude enhancement drives the recover-a-known.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .rmsky import _ratio_bootstrap_se, enhancement_ratio

__all__ = [
    "latitude_ladder",
    "load_spice_racs_dr2",
    "matched_flag_sf",
    "spatial_block_jackknife",
    "structure_function",
    "synthetic_rm_screen",
    "fetch_spice_racs_dr1",
    "run",
]

SPICE_DR1_TABLE = "AS110.spice_racs_dr1_corrected_cut_v02"
CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap"
DR2_DAP = "https://data.csiro.au/collection/csiro:64891"

#: DR1 comparison counts the paper narrates: raw snr_polint>=8 rows from the CASDA TAP leg,
#: and the DR1 release paper's fully-flagged RM count (Thomson et al. 2023). Committed here so
#: the prose numbers trace to evidence rather than living only in a sentence.
DR1_SNR_CUT_ROWS = 7707
DR1_PAPER_FLAGGED_ROWS = 5818


def structure_function(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    rm: np.ndarray,
    rm_err: np.ndarray,
    *,
    bins_deg: np.ndarray | None = None,
    n_boot: int = 100,
    max_pairs: int = 2_000_000,
    seed: int = 0,
    block_ids: np.ndarray | None = None,
) -> dict:
    r"""Noise-debiased second-order RM structure function with bootstrap errors.

    For every source pair within the separation bins:
    :math:`\mathrm{SF}(\delta\theta) = \langle(\mathrm{RM}_1-\mathrm{RM}_2)^2\rangle -
    \langle\sigma_1^2+\sigma_2^2\rangle` --- the second term removes the measurement-noise floor
    (Haverkorn et al. 2004 convention). Pairs are randomly subsampled to ``max_pairs`` (recorded)
    so the cost stays quadratic-safe; bootstrap resamples *sources* (not pairs) to respect the
    correlated pair structure. Returns bin centres, SF, its bootstrap SE, pair counts, and the
    subsample fraction.

    ``block_ids`` (one label per source) additionally returns a leave-one-block-out jackknife SE
    for the plateau (median of the full sample's three largest-separation finite bins). The i.i.d.
    source bootstrap resamples within correlated patches, so on a real sky it is a shot-noise-only
    lower bound; the block jackknife is the field-level uncertainty (the slice's recorded lesson,
    previously applied to the headline ratio and to nothing else).
    """
    rng = np.random.default_rng(seed)
    ra = np.radians(np.asarray(ra_deg, float))
    dec = np.radians(np.asarray(dec_deg, float))
    rm = np.asarray(rm, float)
    var = np.asarray(rm_err, float) ** 2
    n = rm.size
    if bins_deg is None:
        bins_deg = np.logspace(-1, 1.3, 12)  # 0.1 -- 20 deg
    bins = np.radians(np.asarray(bins_deg, float))

    # all pairs (i<j) when tractable; RANDOM pair draws for large n (triu at n~2.5e5 would
    # need ~3e10 index entries -- hundreds of GB). Random pairs are an unbiased SF estimator.
    n_pairs_all = n * (n - 1) // 2
    if n <= 3000:
        i_idx, j_idx = np.triu_indices(n, k=1)
        if n_pairs_all > max_pairs:
            keep = rng.choice(n_pairs_all, max_pairs, replace=False)
            i_idx, j_idx = i_idx[keep], j_idx[keep]
    else:
        n_draw = int(min(max_pairs, n_pairs_all))
        i_idx = rng.integers(0, n, n_draw)
        j_idx = rng.integers(0, n, n_draw)
        good = i_idx != j_idx
        i_idx, j_idx = i_idx[good], j_idx[good]

    # angular separation via the haversine formula (stable at small angles)
    sdlat = np.sin((dec[j_idx] - dec[i_idx]) / 2.0)
    sdlon = np.sin((ra[j_idx] - ra[i_idx]) / 2.0)
    h = sdlat**2 + np.cos(dec[i_idx]) * np.cos(dec[j_idx]) * sdlon**2
    sep = 2.0 * np.arcsin(np.sqrt(np.clip(h, 0.0, 1.0)))

    d2 = (rm[i_idx] - rm[j_idx]) ** 2
    nvar = var[i_idx] + var[j_idx]
    which = np.digitize(sep, bins) - 1
    nb = len(bins) - 1

    def sf_of(mask_sources: np.ndarray) -> np.ndarray:
        m = mask_sources[i_idx] & mask_sources[j_idx]
        out = np.full(nb, np.nan)
        for b in range(nb):
            sel = m & (which == b)
            if sel.sum() >= 20:
                out[b] = d2[sel].mean() - nvar[sel].mean()
        return out

    all_mask = np.ones(n, bool)
    sf = sf_of(all_mask)
    boots = np.full((n_boot, nb), np.nan)
    for k in range(n_boot):
        pick = rng.integers(0, n, n)
        mask = np.zeros(n, bool)
        mask[np.unique(pick)] = True  # source-level resample (unique-set approximation)
        boots[k] = sf_of(mask)
    se = np.nanstd(boots, axis=0)
    counts = np.array([(which == b).sum() for b in range(nb)])
    centres = np.degrees(np.sqrt(bins[:-1] * bins[1:]))
    out = {
        "sep_deg": centres,
        "sf": sf,
        "sf_err": se,
        "n_pairs": counts,
        # ordered random draws double-count unique pairs; fraction quoted per UNIQUE pair
        "pair_fraction": min(1.0, 0.5 * max_pairs / max(1, n_pairs_all)),
    }
    if block_ids is not None:
        ids = np.asarray(block_ids)
        fin = np.where(np.isfinite(sf))[0]
        plats = []
        if fin.size >= 3:
            idx3 = fin[-3:]  # the plateau bins are FIXED by the full sample, then re-evaluated
            for u in np.unique(ids):
                p = sf_of(ids != u)[idx3]
                if np.isfinite(p).any():
                    plats.append(float(np.nanmedian(p)))
        vals = np.asarray(plats, float)
        k = vals.size
        out["plateau_jk_blocks"] = int(k)
        out["plateau_jk_se"] = (
            float(np.sqrt((k - 1) / k * np.sum((vals - vals.mean()) ** 2)))
            if k >= 5
            else float("nan")
        )
    return out


def matched_flag_sf(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    rm: np.ndarray,
    rm_err: np.ndarray,
    flag_ok: np.ndarray,
    *,
    bins_deg: np.ndarray | None = None,
    n_boot: int = 200,
    max_pairs: int = 2_000_000,
    seed: int = 0,
) -> dict:
    r"""Flagged-vs-unflagged SF comparison on ONE shared random pair set, with uncertainties.

    The paper's quality-flag claim compared half-plateau breaks from two *independent* random
    pair draws, neither committed, neither uncertainty-bearing --- an 11-valued statistic read
    off noisy curves (a referee blocker). Here the pair draw happens once over ALL sources; the
    flagged curve keeps only pairs whose both endpoints pass ``flag_ok``, so every difference
    between the two curves is attributable to the flagged-out sources, not to sampling. A shared
    source bootstrap re-evaluates both curves and both breaks per replicate, giving percentile
    intervals for each break and for their difference, plus the fraction of replicates in which
    the unflagged break is the smaller --- the honest strength of "the flag removes small-scale
    contamination". Both full curves (sf, sf_err, n_pairs per bin) are returned for committing.
    """
    rng = np.random.default_rng(seed)
    ra = np.radians(np.asarray(ra_deg, float))
    dec = np.radians(np.asarray(dec_deg, float))
    rm = np.asarray(rm, float)
    var = np.asarray(rm_err, float) ** 2
    ok = np.asarray(flag_ok, bool)
    n = rm.size
    if bins_deg is None:
        bins_deg = np.logspace(-1, 1.3, 12)
    bins = np.radians(np.asarray(bins_deg, float))
    n_pairs_all = n * (n - 1) // 2
    if n <= 3000:
        i_idx, j_idx = np.triu_indices(n, k=1)
        if n_pairs_all > max_pairs:
            keep = rng.choice(n_pairs_all, max_pairs, replace=False)
            i_idx, j_idx = i_idx[keep], j_idx[keep]
    else:  # pragma: no cover - big data only
        n_draw = int(min(max_pairs, n_pairs_all))
        i_idx = rng.integers(0, n, n_draw)
        j_idx = rng.integers(0, n, n_draw)
        good = i_idx != j_idx
        i_idx, j_idx = i_idx[good], j_idx[good]
    sdlat = np.sin((dec[j_idx] - dec[i_idx]) / 2.0)
    sdlon = np.sin((ra[j_idx] - ra[i_idx]) / 2.0)
    h = sdlat**2 + np.cos(dec[i_idx]) * np.cos(dec[j_idx]) * sdlon**2
    sep = 2.0 * np.arcsin(np.sqrt(np.clip(h, 0.0, 1.0)))
    d2 = (rm[i_idx] - rm[j_idx]) ** 2
    nvar = var[i_idx] + var[j_idx]
    which = np.digitize(sep, bins) - 1
    nb = len(bins) - 1
    centres = np.degrees(np.sqrt(bins[:-1] * bins[1:]))

    def sf_of(mask_sources: np.ndarray) -> np.ndarray:
        m = mask_sources[i_idx] & mask_sources[j_idx]
        out = np.full(nb, np.nan)
        for b in range(nb):
            sel = m & (which == b)
            if sel.sum() >= 20:
                out[b] = d2[sel].mean() - nvar[sel].mean()
        return out

    all_mask = np.ones(n, bool)
    curves = {"unflagged": sf_of(all_mask), "flagged": sf_of(ok)}
    boots = {k: np.full((n_boot, nb), np.nan) for k in curves}
    breaks = {k: np.full(n_boot, np.nan) for k in curves}
    for r in range(n_boot):
        pick = rng.integers(0, n, n)
        mask = np.zeros(n, bool)
        mask[np.unique(pick)] = True  # shared source resample: both curves see the SAME sky draw
        for k, base in (("unflagged", mask), ("flagged", mask & ok)):
            boots[k][r] = sf_of(base)
            breaks[k][r] = _sf_break(centres, boots[k][r])

    def _pairs_of(mask_sources: np.ndarray) -> np.ndarray:
        m = mask_sources[i_idx] & mask_sources[j_idx]
        return np.array([(m & (which == b)).sum() for b in range(nb)])

    def _pcts(x: np.ndarray) -> tuple[float, float]:
        x = x[np.isfinite(x)]
        if x.size < 10:
            return float("nan"), float("nan")
        return float(np.percentile(x, 16)), float(np.percentile(x, 84))

    out: dict[str, Any] = {"sep_deg": centres, "n_total": n, "n_flagged": int(ok.sum())}
    for k, base in (("unflagged", all_mask), ("flagged", ok)):
        lo16, hi84 = _pcts(breaks[k])
        out[k] = {
            "sf": curves[k],
            "sf_err": np.nanstd(boots[k], axis=0),
            "n_pairs": _pairs_of(base),
            "break_deg": _sf_break(centres, curves[k]),
            "break_p16": lo16,
            "break_p84": hi84,
        }
    both = np.isfinite(breaks["unflagged"]) & np.isfinite(breaks["flagged"])
    out["n_boot_finite"] = int(both.sum())
    out["frac_boot_unflagged_break_smaller"] = (
        float(np.mean(breaks["unflagged"][both] < breaks["flagged"][both]))
        if both.any()
        else float("nan")
    )
    return out


#: Field realizations averaged for the offline recover-a-known. One realization's bootstrap
#: measures sampling noise inside that field, not the field-to-field scatter, and on this
#: screen the two differ by ~3x.
N_SYNTHETIC_REALIZATIONS = 30


def synthetic_rm_screen(
    n_sources: int = 1500,
    *,
    coherence_deg: float = 2.0,
    amp_high_b: float = 15.0,
    plane_boost: float = 5.0,
    noise: float = 2.0,
    seed: int = 0,
) -> dict:
    r"""A synthetic RM sky with a KNOWN coherence scale and a known plane enhancement.

    Sources are scattered over a patch; the RM field is a sum of Gaussian blobs of angular size
    ``coherence_deg`` (so the SF rises up to ~ the coherence scale and saturates at
    :math:`2\sigma_\mathrm{RM}^2` beyond it), with the RM amplitude boosted by ``plane_boost``
    at low latitude. Gaussian measurement noise ``noise`` (rad/m^2) is added and recorded in
    ``rm_err`` --- the debiasing target the SF must remove.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(0.0, 40.0, n_sources)
    dec = rng.uniform(-20.0, 20.0, n_sources)
    gb = dec  # patch geometry: treat dec as latitude for the fixture
    # Gaussian-blob random field with the requested coherence scale
    n_blob = 220
    bra = rng.uniform(-5.0, 45.0, n_blob)
    bdec = rng.uniform(-25.0, 25.0, n_blob)
    bamp = rng.normal(0.0, amp_high_b, n_blob)
    rm_true = np.zeros(n_sources)
    for k in range(n_blob):
        d2 = (ra - bra[k]) ** 2 * np.cos(np.radians(dec)) ** 2 + (dec - bdec[k]) ** 2
        rm_true += bamp[k] * np.exp(-0.5 * d2 / coherence_deg**2)
    rm_true *= 1.0 + (plane_boost - 1.0) * np.exp(-0.5 * (gb / 5.0) ** 2)
    rm_err = np.full(n_sources, noise)
    rm_obs = rm_true + rng.normal(0.0, noise, n_sources)
    return {
        "ra": ra,
        "dec": dec,
        "gal_b": gb,
        "gal_l": ra,
        "rm": rm_obs,
        "rm_err": rm_err,
        "coherence_deg": coherence_deg,
    }


def fetch_spice_racs_dr1(
    *, snr_min: float = 8.0, max_rows: int = 30000
) -> dict:  # pragma: no cover - network
    """Fetch the SPICE-RACS DR1 RM catalogue from CASDA TAP (verified 24,758 rows, no auth)."""
    from astroquery.utils.tap.core import TapPlus

    tap = TapPlus(url=CASDA_TAP)
    q = (
        f"SELECT TOP {int(max_rows)} ra, dec, l, b, rm, rm_err, snr_polint "
        f"FROM {SPICE_DR1_TABLE} WHERE snr_polint >= {snr_min} AND rm_err > 0"
    )
    t = tap.launch_job(q).get_results()
    return {
        "ra": np.asarray(t["ra"], float),
        "dec": np.asarray(t["dec"], float),
        "gal_l": np.asarray(t["l"], float),
        "gal_b": np.asarray(t["b"], float),
        "rm": np.asarray(t["rm"], float),
        "rm_err": np.asarray(t["rm_err"], float),
    }


DR2_LOCAL = Path("data/spice-racs.dr2.fits")  # gunzipped from the DAP .gz so astropy can memmap


def load_spice_racs_dr2(
    path: str | Path = DR2_LOCAL, *, snr_min: float = 8.0, goodrm: bool = True, dedup: bool = True
) -> dict:  # pragma: no cover - needs the 5 GB DAP file
    """Load the public SPICE-RACS DR2 catalogue FITS (CSIRO DAP csiro:64891, no auth).

    Column names follow the DR1 convention (rm, rm_err, snr_polint, l, b); any variant casing is
    resolved by lookup, and the column each cut actually used is **recorded** in the returned
    ``meta`` dict together with the full cut cascade (raw rows, S/N cut, finite-error cut,
    ``goodRM_flag``, dedup) — a referee round found the sample definition unauditable without it.

    ``dedup=True`` keeps one row per ``cat_id``: RACS tiles overlap, so ~26% of goodRM rows are
    repeat observations of the same component (333,173 rows but 246,508 unique ``cat_id``, which
    is the release's own published post-dedup 8-sigma count). The kept row is the one observed
    closest to its tile centre (best PSF/leakage behaviour), with highest ``snr_polint`` as the
    tie-break. Duplicates are not merely redundant for pair statistics: a pair of rows for the
    SAME source at zero separation with independent noise puts pure noise power into the smallest
    SF bins, and duplicated sources are double-weighted in every median.
    """
    from astropy.io import fits
    from astropy.table import Table

    with fits.open(path, memmap=True) as hdul:
        t = Table(hdul[1].data)
    cols = {c.lower(): c for c in t.colnames}
    used: dict[str, str] = {}

    def col(field: str, *names):
        for nm in names:
            if nm in cols:
                used[field] = cols[nm]
                return np.asarray(t[cols[nm]], float)
        raise KeyError(f"none of {names} in DR2 table (has: {sorted(cols)[:40]}...)")

    rm = col("rm", "rm")
    rm_err = col("rm_err", "rm_err", "e_rm", "rm_err_obs")
    snr = col("snr", "snr_polint", "snr_pi", "snr")
    gl = col("gal_l", "l", "gal_l", "glon")
    gb = col("gal_b", "b", "gal_b", "glat")
    ra = col("ra", "ra", "ra_deg")
    dec = col("dec", "dec", "dec_deg")
    meta: dict = {"snr_column": used["snr"], "n_raw": int(rm.size), "snr_min": float(snr_min)}
    m = snr >= snr_min
    meta["n_after_snr"] = int(m.sum())
    m &= np.isfinite(rm) & (rm_err > 0)
    meta["n_after_finite"] = int(m.sum())
    # the survey's own quality selection: goodRM_flag folds in leakage/snr/StokesI-fit rejection
    # (GATE-2 delta: without it the sample is leakage-contaminated, worst near the plane)
    if goodrm and "goodrm_flag" in cols:
        m &= np.asarray(t[cols["goodrm_flag"]], bool)
    meta["goodrm_applied"] = bool(goodrm and "goodrm_flag" in cols)
    meta["n_after_goodrm"] = int(m.sum())
    idx = np.where(m)[0]
    if dedup and "cat_id" in cols:
        cat = np.asarray(t[cols["cat_id"]])[idx]
        septc = (
            np.asarray(t[cols["separation_tile_centre"]], float)[idx]
            if "separation_tile_centre" in cols
            else np.zeros(idx.size)
        )
        # smallest tile-centre separation first, then highest S/N; first occurrence per cat_id wins
        order = np.lexsort((-snr[idx], septc))
        _, first = np.unique(cat[order], return_index=True)
        idx = idx[order[np.sort(first)]]
        meta["dedup"] = "one row per cat_id: min separation_tile_centre, then max snr_polint"
    meta["n_final"] = int(idx.size)
    out = {
        "ra": ra[idx],
        "dec": dec[idx],
        "gal_l": gl[idx],
        "gal_b": gb[idx],
        "rm": rm[idx],
        "rm_err": rm_err[idx],
        "meta": meta,
    }
    if "goodrm_flag" in cols:
        # per-kept-row goodRM status, so an unflagged load can compare flagged vs unflagged
        # sources on ONE shared pair set (matched_flag_sf) instead of two independent draws
        out["goodrm_ok"] = np.asarray(t[cols["goodrm_flag"]], bool)[idx]
    return out


def spatial_block_jackknife(
    gal_l: np.ndarray,
    gal_b: np.ndarray,
    stat_fn,
    *,
    block_deg: float = 10.0,
    min_blocks: int = 5,
) -> dict:
    """Leave-one-sky-block-out jackknife for a statistic of a spatially correlated field.

    An i.i.d. source bootstrap on an RM sky with a ~2 deg coherence length resamples
    *within* correlated patches, so it measures shot noise and understates the field-level
    uncertainty — the slice's own recorded lesson, which its real-data headline then repeated.
    Blocks of ``block_deg`` (several coherence lengths) are the exchangeable unit. ``stat_fn``
    receives a boolean keep-mask over sources and returns the statistic. Returns the full-sample
    statistic, the jackknife SE, and the block count.
    """
    gal_l = np.asarray(gal_l, float)
    gal_b = np.asarray(gal_b, float)
    ids = np.floor(gal_l / block_deg) * 1000 + np.floor((gal_b + 90.0) / block_deg)
    uniq = np.unique(ids)
    full = float(stat_fn(np.ones(gal_l.size, bool)))
    vals = np.asarray([stat_fn(ids != u) for u in uniq], float)
    vals = vals[np.isfinite(vals)]
    k = vals.size
    # blocks outside the statistic's support (e.g. mid-latitude blocks for a plane/pole ratio)
    # leave it exactly unchanged; report how many actually move it, so "601 blocks" is honest
    n_eff = int(np.sum(vals != full))
    if k < min_blocks:
        return {"stat": full, "se": float("nan"), "n_blocks": int(k), "n_blocks_effective": n_eff}
    se = float(np.sqrt((k - 1) / k * np.sum((vals - vals.mean()) ** 2)))
    return {"stat": full, "se": se, "n_blocks": int(k), "n_blocks_effective": n_eff}


def latitude_ladder(
    s: dict,
    *,
    b_edges: tuple = (0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 90.0),
    max_pairs: int = 500_000,
    n_boot: int = 40,
    jk_block_deg: float | None = 10.0,
) -> dict:
    """SF plateau (and sqrt(plateau/2) = RM dispersion) per |b| bin --- the fluctuation-power profile.

    The two-bin disc--halo split showed a factor-~23 contrast; the ladder resolves HOW the RM
    fluctuation power falls with latitude. Each bin's plateau is the median of the SF's three
    largest-separation finite bins; per-bin source counts are reported (thin bins are honest
    NaNs). The intrinsic+extragalactic floor is latitude-independent, so the ladder's SHAPE is
    Galactic even though each absolute value is an upper bound.

    Each plateau carries TWO errors: ``plateau_err`` (i.i.d. source bootstrap --- shot noise
    only, a lower bound on a correlated sky) and ``plateau_jk_err`` (leave-one-sky-block-out
    jackknife over ``jk_block_deg`` blocks --- the field-level error, the same estimator the
    headline ratio uses). A referee round found every ladder error was the bootstrap the paper
    itself condemns as ~11x too small on the headline; ``sigma_rm_jk_err`` propagates the
    jackknife instead.
    """
    ab = np.abs(np.asarray(s["gal_b"], float))
    out: dict[str, list] = {
        "b_lo": [],
        "b_hi": [],
        "n": [],
        "plateau": [],
        "plateau_err": [],
        "plateau_jk_err": [],
        "jk_blocks": [],
        "sigma_rm": [],
        "sigma_rm_err": [],
        "sigma_rm_jk_err": [],
        "n_pairs": [],
        "pair_fraction": [],
    }
    for lo, hi in zip(b_edges[:-1], b_edges[1:], strict=True):
        m = (ab >= lo) & (ab < hi)
        out["b_lo"].append(lo)
        out["b_hi"].append(hi)
        out["n"].append(int(m.sum()))
        if m.sum() < 200:
            for k in (
                "plateau",
                "plateau_err",
                "plateau_jk_err",
                "sigma_rm",
                "sigma_rm_err",
                "sigma_rm_jk_err",
                "pair_fraction",
            ):
                out[k].append(float("nan"))
            out["n_pairs"].append(0)
            out["jk_blocks"].append(0)
            continue
        blocks = None
        if jk_block_deg is not None:
            gl = np.asarray(s["gal_l"], float)[m]
            gb = np.asarray(s["gal_b"], float)[m]
            blocks = np.floor(gl / jk_block_deg) * 1000 + np.floor((gb + 90.0) / jk_block_deg)
        sf = structure_function(
            s["ra"][m],
            s["dec"][m],
            s["rm"][m],
            s["rm_err"][m],
            max_pairs=max_pairs,
            n_boot=n_boot,
            block_ids=blocks,
        )
        good = np.isfinite(sf["sf"])
        plat = float(np.nanmedian(sf["sf"][good][-3:])) if good.sum() >= 3 else float("nan")
        perr = float(np.nanmedian(sf["sf_err"][good][-3:])) if good.sum() >= 3 else float("nan")
        pjk = float(sf.get("plateau_jk_se", float("nan")))
        out["plateau"].append(plat)
        out["plateau_err"].append(perr)
        out["plateau_jk_err"].append(pjk)
        out["jk_blocks"].append(int(sf.get("plateau_jk_blocks", 0)))
        sig = float(np.sqrt(plat / 2.0)) if plat > 0 else float("nan")
        out["sigma_rm"].append(sig)
        # error propagation through sigma = sqrt(plateau/2): dsigma = dplateau / (4 sigma) * ...
        out["sigma_rm_err"].append(
            float(perr / (2.0 * np.sqrt(2.0 * plat))) if plat > 0 else float("nan")
        )
        out["sigma_rm_jk_err"].append(
            float(pjk / (2.0 * np.sqrt(2.0 * plat)))
            if plat > 0 and np.isfinite(pjk)
            else float("nan")
        )
        out["n_pairs"].append(int(np.sum(sf["n_pairs"])))
        out["pair_fraction"].append(float(sf["pair_fraction"]))
    lad: dict[str, Any] = {k: np.asarray(v) for k, v in out.items()}
    # first-order floor subtraction: the highest-|b| bin estimates the latitude-independent
    # intrinsic+extragalactic floor; sigma_gal = sqrt(sigma^2 - floor^2) (NaN where <= floor)
    fin = np.isfinite(lad["sigma_rm"])
    floor = lad["sigma_rm"][fin][-1]
    floor_err = lad["sigma_rm_err"][fin][-1]
    with np.errstate(invalid="ignore"):
        lad["sigma_gal"] = np.sqrt(np.clip(lad["sigma_rm"] ** 2 - floor**2, 0.0, None))
    lad["floor_sigma"] = float(floor)
    lad["floor_sigma_err"] = float(floor_err)
    lad["floor_sigma_jk_err"] = float(lad["sigma_rm_jk_err"][fin][-1])
    # Floor sensitivity: the subtraction is licensed by ONE bin's plateau, so quote how the
    # plane-bin Galactic dispersion moves when the floor is taken instead from the neighbouring
    # high-|b| bin (the plausible alternative) and across the floor's own +/-1 sigma.
    plane_sig = lad["sigma_rm"][fin][0]
    alt_floor = lad["sigma_rm"][fin][-2] if fin.sum() >= 2 else float("nan")

    def _gal(f: float) -> float:
        return float(np.sqrt(plane_sig**2 - f**2)) if plane_sig > f else float("nan")

    lad["sigma_gal_plane"] = _gal(floor)
    lad["sigma_gal_plane_floor_alt"] = _gal(alt_floor)
    lad["sigma_gal_plane_floor_lo"] = _gal(floor + floor_err)
    lad["sigma_gal_plane_floor_hi"] = _gal(max(floor - floor_err, 0.0))
    lad["alt_floor_sigma"] = float(alt_floor)
    return lad


def _sf_break(sep_deg: np.ndarray, sf: np.ndarray) -> float:
    """Crude coherence-scale estimate: separation where the SF first reaches half its plateau."""
    good = np.isfinite(sf) & (sf > 0)
    if good.sum() < 4:
        return float("nan")
    plateau = np.nanmedian(sf[good][-3:])
    rising = np.where(good & (sf >= 0.5 * plateau))[0]
    return float(sep_deg[rising[0]]) if rising.size else float("nan")


def run(out: str = ".", *, offline: bool = True, dr2: bool = False) -> dict:
    """Offline: SF + latitude recover-a-known on the synthetic screen; real: SPICE-RACS DR1."""
    from pathlib import Path

    # the synthetic recover-a-known ALWAYS runs (its macros must never be overwritten by
    # real values -- GATE-2 delta caught exactly that collision)
    if offline:
        s = synthetic_rm_screen()
        source = "synthetic RM screen"
    elif dr2:  # pragma: no cover - needs the local DAP file
        s = load_spice_racs_dr2()
        s["coherence_deg"] = float("nan")
        source = "SPICE-RACS DR2 goodRM (CSIRO DAP csiro:64891)"
    else:  # pragma: no cover - network
        s = fetch_spice_racs_dr1()
        s["coherence_deg"] = float("nan")
        source = f"SPICE-RACS DR1 ({SPICE_DR1_TABLE})"

    lo = np.abs(s["gal_b"]) < 10.0
    hi = np.abs(s["gal_b"]) > 10.0
    sf_lo = structure_function(s["ra"][lo], s["dec"][lo], s["rm"][lo], s["rm_err"][lo])
    sf_hi = structure_function(s["ra"][hi], s["dec"][hi], s["rm"][hi], s["rm_err"][hi])
    pole = 15.0 if offline else 60.0  # the synthetic patch spans only ±20 deg
    ratio = enhancement_ratio(s["rm"], s["gal_b"], pole_deg=pole)
    ratio_se = _ratio_bootstrap_se(s["rm"], s["gal_b"], pole_deg=pole)
    # The i.i.d. source bootstrap resamples within correlated patches (coherence ~2 deg), so on
    # the real sky it understates the field-level error; the sky-block jackknife is the honest
    # uncertainty on the real leg, exactly as the seed ensemble is on the synthetic one.
    jk = spatial_block_jackknife(
        s["gal_l"],
        s["gal_b"],
        lambda m: enhancement_ratio(s["rm"][m], s["gal_b"][m], pole_deg=pole),
    )
    # The jackknife's own tuning knob, varied: 10 deg was the only block size ever run while
    # the paper's Honest Limits says plane gradients span "tens of degrees". Quote the largest.
    jk_sweep: dict[str, float | None] = {}
    for bd in (5.0, 10.0, 15.0, 20.0, 30.0):
        j = spatial_block_jackknife(
            s["gal_l"],
            s["gal_b"],
            lambda m: enhancement_ratio(s["rm"][m], s["gal_b"][m], pole_deg=pole),
            block_deg=bd,
        )
        jk_sweep[f"{bd:.0f}"] = round(float(j["se"]), 2) if np.isfinite(j["se"]) else None
    sweep_vals = [v for v in jk_sweep.values() if v is not None]
    jk_se_max = max(sweep_vals) if sweep_vals else None
    ratio_ens = ratio_ens_se = None
    ens: list[float] = []
    if offline:
        # The bootstrap resamples sources WITHIN one fixed field realization, so it measures
        # sampling noise and not the realization variance of a correlated random field. On
        # this screen it understates the true scatter by ~3x, and the default seed happens to
        # sit high -- which turned a validation claim into a lucky draw. Quote the ensemble.
        ens = [
            enhancement_ratio((t := synthetic_rm_screen(seed=k))["rm"], t["gal_b"], pole_deg=pole)
            for k in range(N_SYNTHETIC_REALIZATIONS)
        ]
        ratio_ens = float(np.mean(ens))
        ratio_ens_se = float(np.std(ens, ddof=1))
    break_lo = _sf_break(sf_lo["sep_deg"], sf_lo["sf"])
    break_hi = _sf_break(sf_hi["sep_deg"], sf_hi["sf"])

    # Definitional variants, mirrored from the sibling rmsky slice where each was free to move
    # the answer (alt bins moved rmsky's ratio +39%): never previously run on DR2, where the
    # quoted error is 10%.
    variants: dict[str, float | None] = {}
    if not offline:  # pragma: no cover - real legs only

        def _r2(x: float) -> float | None:
            return round(float(x), 2) if np.isfinite(x) else None

        variants["enhancement_ratio_alt_bins_5_70"] = _r2(
            enhancement_ratio(s["rm"], s["gal_b"], plane_deg=5.0, pole_deg=70.0)
        )
        c300 = np.abs(s["rm"]) < 300.0
        variants["n_cut300"] = int(c300.sum())
        variants["enhancement_ratio_cut300"] = _r2(
            enhancement_ratio(s["rm"][c300], s["gal_b"][c300], pole_deg=pole)
        )
        med_e = float(np.median(s["rm_err"]))
        ecut = s["rm_err"] < med_e
        variants["rm_err_median"] = round(med_e, 2)
        variants["enhancement_ratio_erm_cut"] = _r2(
            enhancement_ratio(s["rm"][ecut], s["gal_b"][ecut], pole_deg=pole)
        )

    ladder = latitude_ladder(s) if not offline else None  # pragma: no cover - big data only
    plat_lo = float(np.nanmedian(sf_lo["sf"][-3:]))
    plat_hi = float(np.nanmedian(sf_hi["sf"][-3:]))
    metrics = {
        "source": source,
        "is_real": not offline,
        "n_sources": int(s["rm"].size),
        "enhancement_ratio": round(float(ratio), 2),
        "enhancement_ratio_se": round(float(ratio_se), 2),
        "enhancement_ratio_jackknife_se": round(float(jk["se"]), 2),
        "n_jackknife_blocks": jk["n_blocks"],
        "n_jackknife_blocks_effective": jk["n_blocks_effective"],
        "jackknife_se_by_block_deg": jk_sweep,
        "enhancement_ratio_jackknife_se_max": jk_se_max,
        # Offline only: the honest uncertainty on a recover-a-known over a random field.
        "enhancement_ratio_ensemble": None if ratio_ens is None else round(ratio_ens, 2),
        "enhancement_ratio_ensemble_sd": None if ratio_ens_se is None else round(ratio_ens_se, 2),
        "enhancement_ratio_ensemble_sem": None
        if ratio_ens_se is None
        else round(ratio_ens_se / np.sqrt(N_SYNTHETIC_REALIZATIONS), 2),
        "n_realizations": N_SYNTHETIC_REALIZATIONS if offline else None,
        "sf_plateau_low_b": round(plat_lo, 1),
        "sf_plateau_high_b": round(plat_hi, 1),
        # the disc/halo fluctuation-POWER contrast the docstring's "factor ~23" refers to;
        # previously computed nowhere (the referee had to derive it from two other macros)
        "sf_plateau_ratio": round(plat_lo / plat_hi, 1) if plat_hi > 0 else None,
        "sf_break_low_b_deg": round(break_lo, 2) if np.isfinite(break_lo) else None,
        "sf_break_high_b_deg": round(break_hi, 2) if np.isfinite(break_hi) else None,
        # None (JSON null), not a bare NaN literal: RFC 8259 has no NaN
        "true_coherence_deg": (
            float(s["coherence_deg"]) if np.isfinite(s.get("coherence_deg", float("nan"))) else None
        ),
    }
    metrics.update(variants)
    if "meta" in s:  # pragma: no cover - big data only
        metrics["sample_cascade"] = s["meta"]
        # The DR1 comparison counts the paper narrates (raw S/N-cut rows vs the DR1 paper's
        # fully-flagged RM count) belong in committed evidence, not only in prose.
        metrics["dr1_snr_cut_rows"] = DR1_SNR_CUT_ROWS
        metrics["dr1_paper_flagged_rows"] = DR1_PAPER_FLAGGED_ROWS
    if ladder is not None:  # pragma: no cover - big data only
        fin = np.isfinite(ladder["sigma_rm"])
        fin_idx = np.where(fin)[0]
        floor_i = int(fin_idx[-1])  # the bin the floor is READ FROM: its sigma_gal is zero by
        # construction, a definition rather than a measurement -- record null, not 0.0

        def _jkv(i: int) -> float | None:
            v = float(ladder["sigma_rm_jk_err"][i])
            return round(v, 2) if np.isfinite(v) else None

        metrics.update(
            {
                "ladder_bins": [
                    {
                        "b": f"{ladder['b_lo'][i]:.0f}-{ladder['b_hi'][i]:.0f}",
                        "n": int(ladder["n"][i]),
                        "n_pairs": int(ladder["n_pairs"][i]),
                        "pair_fraction": round(float(ladder["pair_fraction"][i]), 4),
                        "plateau": round(float(ladder["plateau"][i]), 1),
                        "plateau_err": round(float(ladder["plateau_err"][i]), 1),
                        "plateau_jk_err": (
                            round(float(ladder["plateau_jk_err"][i]), 1)
                            if np.isfinite(ladder["plateau_jk_err"][i])
                            else None
                        ),
                        "jk_blocks": int(ladder["jk_blocks"][i]),
                        "sigma_rm": round(float(ladder["sigma_rm"][i]), 1),
                        "sigma_rm_err": round(float(ladder["sigma_rm_err"][i]), 2),
                        "sigma_rm_jk_err": _jkv(i),
                        "sigma_gal": (
                            None if i == floor_i else round(float(ladder["sigma_gal"][i]), 1)
                        ),
                    }
                    for i in range(len(ladder["n"]))
                    if fin[i]
                ],
                "sigma_rm_plane": round(float(ladder["sigma_rm"][fin][0]), 1),
                "sigma_rm_pole": round(float(ladder["sigma_rm"][fin][-1]), 1),
                "ladder_floor_sigma": round(ladder["floor_sigma"], 1),
                "ladder_floor_sigma_err": round(ladder["floor_sigma_err"], 2),
                "ladder_floor_sigma_jk_err": (
                    round(ladder["floor_sigma_jk_err"], 2)
                    if np.isfinite(ladder["floor_sigma_jk_err"])
                    else None
                ),
                "ladder_alt_floor_sigma": round(ladder["alt_floor_sigma"], 1),
                "sigma_gal_plane": round(ladder["sigma_gal_plane"], 1),
                "sigma_gal_plane_floor_alt": round(ladder["sigma_gal_plane_floor_alt"], 1),
                "sigma_gal_plane_floor_lo": round(ladder["sigma_gal_plane_floor_lo"], 1),
                "sigma_gal_plane_floor_hi": round(ladder["sigma_gal_plane_floor_hi"], 1),
            }
        )
        # Floor sensitivity where it bites (the second-highest-|b| bin), against the measured
        # floor and across the literature floor range. Previously derived inside _write_macros
        # only, so the two Honest-Limits numbers appeared in no committed JSON.
        if fin.sum() >= 2:
            mid_i = int(fin_idx[-2])
            sig_mid = float(ladder["sigma_rm"][mid_i])
            metrics["sigma_gal_mid_bin"] = (
                f"{ladder['b_lo'][mid_i]:.0f}-{ladder['b_hi'][mid_i]:.0f}"
            )
            metrics["sigma_gal_mid"] = round(float(ladder["sigma_gal"][mid_i]), 1)
            for name, f in zip(("lit_lo", "lit_hi"), LITERATURE_FLOOR_RANGE, strict=True):
                metrics[f"sigma_gal_mid_floor_{name}"] = round(
                    float(np.sqrt(max(sig_mid**2 - f**2, 0.0))), 1
                )
    flag_cmp = None
    if dr2 and not offline:  # pragma: no cover - big data only
        # The quality-flag comparison behind the paper's leakage claim, on a MATCHED pair set:
        # one pair draw over the unflagged sample, the flagged curve keeping only pairs whose
        # both endpoints pass goodRM. The previous version compared breaks from two independent
        # uncommitted draws with no uncertainty -- a referee blocker.
        su = load_spice_racs_dr2(goodrm=False)
        hiu = np.abs(su["gal_b"]) > 10.0
        flag_cmp = matched_flag_sf(
            su["ra"][hiu],
            su["dec"][hiu],
            su["rm"][hiu],
            su["rm_err"][hiu],
            su["goodrm_ok"][hiu],
        )
        metrics["unflagged_n_sources"] = int(su["rm"].size)

        def _leg(d: dict) -> dict:
            return {
                "sf": [None if not np.isfinite(x) else round(float(x), 1) for x in d["sf"]],
                "sf_err": [None if not np.isfinite(x) else round(float(x), 1) for x in d["sf_err"]],
                "n_pairs": [int(x) for x in d["n_pairs"]],
                "break_deg": round(d["break_deg"], 2) if np.isfinite(d["break_deg"]) else None,
                "break_p16": round(d["break_p16"], 2) if np.isfinite(d["break_p16"]) else None,
                "break_p84": round(d["break_p84"], 2) if np.isfinite(d["break_p84"]) else None,
            }

        metrics["flag_comparison_high_b"] = {
            "sep_deg": [round(float(x), 3) for x in flag_cmp["sep_deg"]],
            "n_total": flag_cmp["n_total"],
            "n_flagged": flag_cmp["n_flagged"],
            "n_boot_finite": flag_cmp["n_boot_finite"],
            "frac_boot_unflagged_break_smaller": round(
                flag_cmp["frac_boot_unflagged_break_smaller"], 3
            ),
            "flagged": _leg(flag_cmp["flagged"]),
            "unflagged": _leg(flag_cmp["unflagged"]),
        }
        # keep the headline keys, now sourced from the matched comparison
        metrics["unflagged_sf_break_high_b_deg"] = metrics["flag_comparison_high_b"]["unflagged"][
            "break_deg"
        ]
        metrics["flagged_sf_break_high_b_deg"] = metrics["flag_comparison_high_b"]["flagged"][
            "break_deg"
        ]
        metrics["unflagged_sf_plateau_high_b"] = round(
            float(np.nanmedian(flag_cmp["unflagged"]["sf"][-3:])), 1
        )
        # Does the goodRM selection itself vary with |b| (depolarization-selected disc
        # populations)? The bound is one line: the flag-removed fraction per ladder bin.
        abu = np.abs(su["gal_b"])
        flag_by_bin = []
        edges = (0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 90.0)
        for lo_b, hi_b in zip(edges[:-1], edges[1:], strict=True):
            mb = (abu >= lo_b) & (abu < hi_b)
            n_all = int(mb.sum())
            n_ok = int((mb & su["goodrm_ok"]).sum())
            flag_by_bin.append(
                {
                    "b": f"{lo_b:.0f}-{hi_b:.0f}",
                    "n_unflagged": n_all,
                    "n_goodrm": n_ok,
                    "pct_removed": round(100.0 * (n_all - n_ok) / n_all, 2) if n_all else None,
                }
            )
        metrics["goodrm_removed_by_bin"] = flag_by_bin
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "rmstructure_metrics.json")
    if offline:
        # The per-seed ratios behind the ensemble mean/SD: committed so "3.15 +/- 1.11 over 30
        # realizations" is auditable, not just asserted (a referee round flagged the omission).
        write_results(
            {
                "source": "synthetic RM screen seed ensemble (recover-a-known)",
                "pole_deg": pole,
                "injected_plane_boost": 5.0,
                "ratios": [round(float(x), 3) for x in ens],
                # The error-model recover-a-known the paper leads with, committed rather than
                # living only in macros: the seed-0 block jackknife against the true scatter.
                "seed0_jackknife_se": round(float(jk["se"]), 2),
                "seed0_jackknife_se_by_block_deg": jk_sweep,
                "ensemble_sd": None if ratio_ens_se is None else round(ratio_ens_se, 2),
            },
            op / "results" / "rmstructure_synthetic.json",
        )
    # Per-leg figure names: JSON and macros are merge-guarded but a shared figure name meant a
    # direct offline run in the repo root silently swapped the real figure for the synthetic one.
    _figure(
        s,
        sf_lo,
        sf_hi,
        op / "papers" / "rmstructure" / "figures",
        ladder=ladder,
        name="rmstructure_syn" if offline else "rmstructure_real",
    )
    _write_macros(metrics, op / "papers" / "rmstructure" / "generated" / "macros.tex")
    return metrics


def _figure(s, sf_lo, sf_hi, out_dir, *, ladder=None, name: str = "rmstructure_syn") -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ncols = 3 if ladder is not None else 2
    fig, axes = plt.subplots(1, ncols, figsize=(4.7 * ncols, 3.9))
    ax1, ax2 = axes[0], axes[1]
    # robust symmetric colour scale: fixed +/-40 saturated against a ~200 rad/m^2 plane sigma,
    # making the headline contrast invisible in the one figure the paper prints
    v = float(np.percentile(np.abs(s["rm"]), 98.0))
    sc = ax1.scatter(s["gal_l"], s["gal_b"], c=s["rm"], s=4, cmap="RdBu_r", vmin=-v, vmax=v)
    fig.colorbar(sc, ax=ax1, label="RM (rad m$^{-2}$)")
    ax1.set(xlabel="Galactic longitude (deg)", ylabel="Galactic latitude (deg)", title="RM sky")
    for sf, lab, c in ((sf_lo, "|b| < 10°", "C3"), (sf_hi, "|b| > 10°", "C0")):
        g = np.isfinite(sf["sf"])
        ax2.errorbar(
            sf["sep_deg"][g],
            sf["sf"][g],
            yerr=sf["sf_err"][g],
            fmt="o-",
            ms=3,
            color=c,
            lw=1,
            label=lab,
        )
    ax2.set(
        xscale="log",
        yscale="log",
        xlabel=r"separation $\delta\theta$ (deg)",
        ylabel=r"SF(RM) (rad$^2$ m$^{-4}$)",
        title="Noise-debiased structure function",
    )
    ax2.legend(fontsize=8)
    if ladder is not None:  # pragma: no cover - big data only
        ax3 = axes[2]
        fin = np.isfinite(ladder["sigma_rm"])
        bc = 0.5 * (ladder["b_lo"][fin] + ladder["b_hi"][fin])
        jk = ladder["sigma_rm_jk_err"][fin]
        boot = ladder["sigma_rm_err"][fin]
        yerr = np.where(np.isfinite(jk), jk, boot)
        ax3.errorbar(bc, ladder["sigma_rm"][fin], yerr=yerr, fmt="s-", color="k", ms=4, lw=1)
        ax3.set(
            yscale="log",
            xlabel=r"$|b|$ bin centre (deg)",
            ylabel=r"$\sigma_\mathrm{RM}$ (rad m$^{-2}$)",
            title="Latitude ladder (jackknife errors)",
        )
    fig.tight_layout()
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


#: Literature intrinsic+extragalactic RM floor range (rad/m^2) used for the floor-sensitivity
#: macros (Mao et al. 2010; the Taylor et al. 2009 discussion) — the plausible span a reader
#: might adopt instead of this paper's measured polar plateau.
LITERATURE_FLOOR_RANGE = (9.0, 15.0)

# Ratio-family macros are display-rounded to one decimal so the quoted value and its error
# carry matching precision (a referee flagged 11.17 +/- 0.1).
_ONE_DP = {
    "enhancement_ratio",
    "enhancement_ratio_se",
    "enhancement_ratio_jackknife_se",
    "enhancement_ratio_jackknife_se_max",
    "enhancement_ratio_ensemble",
    "enhancement_ratio_ensemble_sd",
    "enhancement_ratio_ensemble_sem",
    "enhancement_ratio_alt_bins_5_70",
    "enhancement_ratio_cut300",
    "enhancement_ratio_erm_cut",
    # the floor is quoted at one decimal, so its errors display at one decimal too
    "ladder_floor_sigma_err",
    "ladder_floor_sigma_jk_err",
}


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        if val is None:
            return "--"
        if key in _ONE_DP:
            return f"{float(val):.1f}"
        return str(val)

    def _casc(key: str) -> str:
        c = m.get("sample_cascade") or {}
        val = c.get(key)
        return "--" if val is None else str(val)

    def _flag(leg: str, key: str) -> str:
        d = (m.get("flag_comparison_high_b") or {}).get(leg) or {}
        val = d.get(key)
        return "--" if val is None else str(val)

    def _flag_top(key: str) -> str:
        val = (m.get("flag_comparison_high_b") or {}).get(key)
        return "--" if val is None else str(val)

    # Display-time summaries of the per-bin goodRM removal (all inputs live in the JSON): the
    # Honest Limits item needs only the span, not six numbers.
    fb = [
        b["pct_removed"]
        for b in (m.get("goodrm_removed_by_bin") or [])
        if b.get("pct_removed") is not None
    ]
    if fb:
        m = dict(m)
        m["goodrm_removed_pct_min"] = round(min(fb), 2)
        m["goodrm_removed_pct_max"] = round(max(fb), 2)

    pref = "rmsReal" if m.get("is_real") else "rmsSyn"
    lines = [
        "% Auto-generated by jansky_research.rmstructure._write_macros -- do not edit.",
        "% Synthetic (rmsSyn*) and real (rmsReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders so offline CI and real runs never collide.",
        # \rmsSource (read first by the merge guard) carries the LIVE run's provenance; the
        # per-namespace sources below it are for display, so the paper can cite each leg's
        # provenance without the shared marker inverting when the other leg rebuilds.
        rf"\newcommand{{\rmsSource}}{{{m['source']}}}",
    ]
    for ns in ("rmsSyn", "rmsReal"):
        live = ns == pref
        g = (lambda k: _fmt(k)) if live else (lambda k: "--")
        fl = _flag if live else (lambda leg, k: "--")
        ft = _flag_top if live else (lambda k: "--")
        lines += [
            rf"\newcommand{{\{ns}Source}}{{{m['source'] if live else '--'}}}",
            rf"\newcommand{{\{ns}N}}{{{g('n_sources')}}}",
            rf"\newcommand{{\{ns}Ratio}}{{{g('enhancement_ratio')}}}",
            rf"\newcommand{{\{ns}RatioSe}}{{{g('enhancement_ratio_se')}}}",
            rf"\newcommand{{\{ns}RatioJkSe}}{{{g('enhancement_ratio_jackknife_se')}}}",
            rf"\newcommand{{\{ns}JkBlocks}}{{{g('n_jackknife_blocks')}}}",
            rf"\newcommand{{\{ns}RatioEns}}{{{g('enhancement_ratio_ensemble')}}}",
            rf"\newcommand{{\{ns}RatioEnsSd}}{{{g('enhancement_ratio_ensemble_sd')}}}",
            rf"\newcommand{{\{ns}RatioEnsSem}}{{{g('enhancement_ratio_ensemble_sem')}}}",
            rf"\newcommand{{\{ns}NReal}}{{{g('n_realizations')}}}",
            rf"\newcommand{{\{ns}PlatLo}}{{{g('sf_plateau_low_b')}}}",
            rf"\newcommand{{\{ns}PlatHi}}{{{g('sf_plateau_high_b')}}}",
            rf"\newcommand{{\{ns}BreakLo}}{{{g('sf_break_low_b_deg')}}}",
            rf"\newcommand{{\{ns}BreakHi}}{{{g('sf_break_high_b_deg')}}}",
            rf"\newcommand{{\{ns}SigPlane}}{{{g('sigma_rm_plane')}}}",
            rf"\newcommand{{\{ns}SigPole}}{{{g('sigma_rm_pole')}}}",
            rf"\newcommand{{\{ns}FloorSig}}{{{g('ladder_floor_sigma')}}}",
            rf"\newcommand{{\{ns}FloorSigErr}}{{{g('ladder_floor_sigma_err')}}}",
            rf"\newcommand{{\{ns}AltFloorSig}}{{{g('ladder_alt_floor_sigma')}}}",
            rf"\newcommand{{\{ns}SigGalPlane}}{{{g('sigma_gal_plane')}}}",
            rf"\newcommand{{\{ns}SigGalPlaneAltFloor}}{{{g('sigma_gal_plane_floor_alt')}}}",
            rf"\newcommand{{\{ns}SigGalPlaneFloorLo}}{{{g('sigma_gal_plane_floor_lo')}}}",
            rf"\newcommand{{\{ns}SigGalPlaneFloorHi}}{{{g('sigma_gal_plane_floor_hi')}}}",
            rf"\newcommand{{\{ns}SigGalMidBin}}{{{g('sigma_gal_mid_bin')}}}",
            rf"\newcommand{{\{ns}SigGalMid}}{{{g('sigma_gal_mid')}}}",
            rf"\newcommand{{\{ns}SigGalMidFloorLitLo}}{{{g('sigma_gal_mid_floor_lit_lo')}}}",
            rf"\newcommand{{\{ns}SigGalMidFloorLitHi}}{{{g('sigma_gal_mid_floor_lit_hi')}}}",
            rf"\newcommand{{\{ns}BreakHiUnflagged}}{{{g('unflagged_sf_break_high_b_deg')}}}",
            rf"\newcommand{{\{ns}NUnflagged}}{{{g('unflagged_n_sources')}}}",
            rf"\newcommand{{\{ns}NRowsRaw}}{{{_casc('n_raw') if live else '--'}}}",
            rf"\newcommand{{\{ns}NAfterSnr}}{{{_casc('n_after_snr') if live else '--'}}}",
            rf"\newcommand{{\{ns}NAfterFinite}}{{{_casc('n_after_finite') if live else '--'}}}",
            rf"\newcommand{{\{ns}NAfterGoodrm}}{{{_casc('n_after_goodrm') if live else '--'}}}",
            rf"\newcommand{{\{ns}PlatRatio}}{{{g('sf_plateau_ratio')}}}",
            rf"\newcommand{{\{ns}RatioJkSeMax}}{{{g('enhancement_ratio_jackknife_se_max')}}}",
            rf"\newcommand{{\{ns}JkBlocksEff}}{{{g('n_jackknife_blocks_effective')}}}",
            rf"\newcommand{{\{ns}RatioAltBins}}{{{g('enhancement_ratio_alt_bins_5_70')}}}",
            rf"\newcommand{{\{ns}RatioCutThree}}{{{g('enhancement_ratio_cut300')}}}",
            rf"\newcommand{{\{ns}RatioErmCut}}{{{g('enhancement_ratio_erm_cut')}}}",
            rf"\newcommand{{\{ns}FloorSigJkErr}}{{{g('ladder_floor_sigma_jk_err')}}}",
            rf"\newcommand{{\{ns}BreakHiFlagM}}{{{fl('flagged', 'break_deg')}}}",
            rf"\newcommand{{\{ns}BreakHiFlagLo}}{{{fl('flagged', 'break_p16')}}}",
            rf"\newcommand{{\{ns}BreakHiFlagHi}}{{{fl('flagged', 'break_p84')}}}",
            rf"\newcommand{{\{ns}BreakHiUnflagLo}}{{{fl('unflagged', 'break_p16')}}}",
            rf"\newcommand{{\{ns}BreakHiUnflagHi}}{{{fl('unflagged', 'break_p84')}}}",
            rf"\newcommand{{\{ns}BreakFracSmall}}{{{ft('frac_boot_unflagged_break_smaller')}}}",
            rf"\newcommand{{\{ns}FlagRemPctMin}}{{{g('goodrm_removed_pct_min')}}}",
            rf"\newcommand{{\{ns}FlagRemPctMax}}{{{g('goodrm_removed_pct_max')}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and would
    # otherwise blank the other mode's macros with '--'. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="SPICE-RACS RM structure functions by latitude.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--dr2", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline, dr2=args.dr2), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
