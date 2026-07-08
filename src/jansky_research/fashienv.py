"""FASHI environment statistics: the HI mass function nobody has split by environment (plan 45).

The FAST All Sky HI Survey (FASHI) DR2 (arXiv:2606.31539; 156,411 sources) published only a
GLOBAL Schechter HI mass function (HIMF); the sole environment paper (arXiv:2510.22902) used
230 DR1-era group galaxies. GATE-0 (2026-07-08): **DR2's catalogue is not yet public** (release
targeted ~Aug 2026 at zcp521.github.io/fashi; the table links 404 today), so this is the DR1
FIRST LEG --- the same DR1-while-DR2-embargoed pattern the merged `rmstructure` slice used.
DR1 (VizieR J/other/SCPMA/67.19511/table2; **41,741 HI sources**, Dec>-14, z<0.09, precomputed
`logMass`/`Dist`) is still ~180x the only prior environment sample, and its environment-split HIMF
is equally unrun. The DR2 swap is a one-line source change once it publishes.

Scope (honest, self-contained; corrects plan 45 where FASHI lacks the needed inputs):

- **Environment-split HIMF** (the headline): the 1/Vmax HIMF of FASHI HI masses computed
  separately in voids vs walls (Douglass VoidFinder) and in group-member vs field environments
  (Tempel groups), each fit with a Schechter function. **The relative void-wall / group-field
  knee difference is the robust deliverable** (both bins share the same 1/Vmax method, so the
  comparison is insensitive to the absolute completeness); the ABSOLUTE faint-end slope from a
  simple 1/Vmax is steeper than the published FASHI global HIMF, which used the survey's full
  completeness function --- stated, not papered over. Anchor: Moorman+2014 found the void HIMF
  knee suppressed by ~0.1-0.2 dex vs walls; we recover the same sign and magnitude.
- The plan's "gas fraction at fixed M*" and "HI-deficiency vs clustercentric radius" are DROPPED:
  FASHI carries no stellar masses or optical diameters (both need optical counterparts), and the
  raw median-HI-vs-radius of DETECTED sources is selection-biased (deficient galaxies drop out of
  a flux-limited HI sample). Void-wall and group-field HIMFs are the cleanest self-contained
  statements the data support.

Cross-match catalogues (VizieR, all resolve): Tempel+2017 SDSS groups (J/A+A/602/A100: galaxy
GroupID/Ngal + group R200/M200); Douglass+2023 voids (J/ApJS/265/7: VoidFinder spheres +
V2/VIDE/REVOLVER galaxy membership). Footprint caveat (GATE-0): the group/void catalogues are
SDSS-based (Dec>~0, z<~0.05 volume-limited), so the cross-matched sample is the SDSS-cap subset,
reported honestly, not the full 41,741.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "comoving_xyz",
    "void_membership",
    "assign_groups",
    "clustercentric_radius",
    "vmax_1vmax",
    "himf",
    "fit_schechter",
    "schechter",
    "synthetic_environment_catalogue",
    "run",
]

# FASHI DR1 catalogue (VizieR) and cross-match catalogues
FASHI_DR1_VIZIER = "J/other/SCPMA/67.19511/table2"  # 41,741 extragalactic HI sources
TEMPEL_VIZIER = "J/A+A/602/A100"  # SDSS groups (table1 galaxies, table2 groups)
DOUGLASS_VIZIER = "J/ApJS/265/7"  # voids (table1 VoidFinder spheres, table5 V2 membership)
FASHI_DR2_PORTAL = "https://zcp521.github.io/fashi"  # DR2 catalogue, public ~Aug 2026

C_KM_S = 299792.458
H0 = 70.0  # km/s/Mpc, the FASHI DR1 distance convention (h70)
FASHI_FLUX_LIMIT = 0.30  # Jy km/s integrated-flux limit for the 1/Vmax weighting. A single-cut
# order-of-magnitude value (FAST is deeper than ALFALFA's ~0.7 Jy km/s); it is NOT FASHI's full
# completeness function, so the ABSOLUTE HIMF (esp. the faint-end slope) is only approximate ---
# the relative environment offsets, which share this cut, are what the slice reports.

# Planck2018 deceleration: q0 = Om/2 - OL = 0.315/2 - 0.685 = -0.527 (matches the Douglass frame)
_Q0 = -0.527


def _comoving_distance_mpc(z: np.ndarray, h0: float = H0) -> np.ndarray:
    r"""Low-z comoving distance, 2nd-order in z: :math:`(c/H_0)\,z\,[1 - (1+q_0)z/2]`.

    Uses the Planck2018 :math:`q_0=-0.527` (not the EdS 0.5) so the void-membership geometry
    matches the Douglass+2023 catalogue's own cosmology; the residual vs the full integral is
    <0.02% at z<0.09. With ``h0=100`` the result is Mpc/h (H0-independent) --- the Douglass
    frame; with ``h0=70`` it is physical Mpc (the FASHI DR1 convention).
    """
    z = np.asarray(z, float)
    return (C_KM_S / h0) * z * (1.0 - 0.5 * (1.0 + _Q0) * z)


def comoving_xyz(
    ra_deg: np.ndarray, dec_deg: np.ndarray, z: np.ndarray, *, h0: float = H0
) -> np.ndarray:
    """(RA, Dec, z) -> comoving Cartesian, the frame for void-sphere membership.

    Use ``h0=100`` to match the Douglass void spheres' Mpc/h frame; ``h0=70`` for physical Mpc.
    """
    d = _comoving_distance_mpc(z, h0)
    ra = np.radians(np.asarray(ra_deg, float))
    dec = np.radians(np.asarray(dec_deg, float))
    return np.stack(
        [d * np.cos(dec) * np.cos(ra), d * np.cos(dec) * np.sin(ra), d * np.sin(dec)], axis=1
    )


def void_membership(
    gal_xyz: np.ndarray, sphere_xyz: np.ndarray, sphere_radius: np.ndarray
) -> np.ndarray:
    """Boolean void membership: True where a galaxy lies inside ANY void sphere (VoidFinder).

    The standard VoidFinder test (Douglass+2023 table1 maximal spheres): a galaxy is a void
    member iff its comoving distance to some sphere centre is less than that sphere's radius.
    Chunked so a 10^4 x 10^4 distance matrix never materialises.
    """
    gal = np.asarray(gal_xyz, float)
    cen = np.asarray(sphere_xyz, float)
    rad = np.asarray(sphere_radius, float)
    inside = np.zeros(len(gal), bool)
    step = 2000
    for a in range(0, len(gal), step):
        g = gal[a : a + step]
        d2 = ((g[:, None, :] - cen[None, :, :]) ** 2).sum(axis=2)  # (chunk, n_sphere)
        inside[a : a + step] = (d2 < rad[None, :] ** 2).any(axis=1)
    return inside


def assign_groups(
    gal_ra: np.ndarray,
    gal_dec: np.ndarray,
    gal_cz: np.ndarray,
    grp_ra: np.ndarray,
    grp_dec: np.ndarray,
    grp_cz: np.ndarray,
    grp_r200_mpc: np.ndarray,
    *,
    dv_max: float = 1000.0,
) -> np.ndarray:
    """Nearest group within R200 (projected) and ``dv_max`` (km/s): index into groups, or -1.

    Projected separation uses the small-angle sky distance at the group's redshift; a galaxy is
    assigned to the group whose (projected sep / R200) is smallest among groups passing the
    velocity cut. Returns a per-galaxy group index (-1 = field).
    """
    ra = np.radians(np.asarray(gal_ra, float))
    dec = np.radians(np.asarray(gal_dec, float))
    cz = np.asarray(gal_cz, float)
    gra = np.radians(np.asarray(grp_ra, float))
    gdec = np.radians(np.asarray(grp_dec, float))
    gcz = np.asarray(grp_cz, float)
    r200 = np.asarray(grp_r200_mpc, float)
    out = np.full(len(cz), -1, int)
    for i in range(len(cz)):
        dv = np.abs(cz[i] - gcz)
        near = dv <= dv_max
        if not near.any():
            continue
        # angular separation (haversine), projected to Mpc at the group's distance
        sdlat = np.sin((gdec - dec[i]) / 2.0)
        sdlon = np.sin((gra - ra[i]) / 2.0)
        h = sdlat**2 + np.cos(dec[i]) * np.cos(gdec) * sdlon**2
        sep_rad = 2.0 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))
        sep_mpc = sep_rad * (gcz / H0)
        ratio = np.where(near, sep_mpc / np.maximum(r200, 1e-6), np.inf)
        j = int(np.argmin(ratio))
        if ratio[j] <= 1.0:
            out[i] = j
    return out


def clustercentric_radius(
    gal_ra: np.ndarray,
    gal_dec: np.ndarray,
    gal_cz: np.ndarray,
    grp_idx: np.ndarray,
    grp_ra: np.ndarray,
    grp_dec: np.ndarray,
    grp_cz: np.ndarray,
    grp_r200_mpc: np.ndarray,
) -> np.ndarray:
    """Projected R/R200 for each group member (NaN for field galaxies, grp_idx == -1)."""
    ra = np.radians(np.asarray(gal_ra, float))
    dec = np.radians(np.asarray(gal_dec, float))
    gra = np.radians(np.asarray(grp_ra, float))
    gdec = np.radians(np.asarray(grp_dec, float))
    gcz = np.asarray(grp_cz, float)
    r200 = np.asarray(grp_r200_mpc, float)
    out = np.full(len(ra), np.nan)
    m = np.asarray(grp_idx, int) >= 0
    for i in np.where(m)[0]:
        j = grp_idx[i]
        sdlat = np.sin((gdec[j] - dec[i]) / 2.0)
        sdlon = np.sin((gra[j] - ra[i]) / 2.0)
        h = sdlat**2 + np.cos(dec[i]) * np.cos(gdec[j]) * sdlon**2
        sep_mpc = 2.0 * np.arcsin(np.sqrt(np.clip(h, 0, 1))) * (gcz[j] / H0)
        out[i] = sep_mpc / max(r200[j], 1e-6)
    return out


def vmax_1vmax(
    log_mhi: np.ndarray,
    dist_mpc: np.ndarray,
    flux: np.ndarray,
    *,
    flux_limit: float = FASHI_FLUX_LIMIT,
    d_min: float = 5.0,
    d_max: float = 386.0,  # ~ z=0.09 at H0=70
) -> np.ndarray:
    """Per-source maximum comoving volume (Mpc^3) it could occupy and stay above ``flux_limit``.

    Flux scales as D^-2, so d_lim = d * sqrt(flux / flux_limit); V_max ∝ (min(d_lim, d_max)^3 -
    d_min^3). The survey solid angle cancels in the HIMF normalisation per bin, so V_max is
    quoted per unit steradian and the caller supplies the area (`himf`). Sources fainter than
    the limit get V_max via their own flux (they define the local completeness), never negative.
    """
    d = np.asarray(dist_mpc, float)
    f = np.asarray(flux, float)
    d_lim = d * np.sqrt(np.maximum(f, 1e-12) / flux_limit)
    d_hi = np.clip(d_lim, d_min, d_max)
    return np.maximum(d_hi**3 - d_min**3, 0.0) / 3.0  # per steradian


def himf(
    log_mhi: np.ndarray,
    vmax_per_sr: np.ndarray,
    *,
    area_sr: float,
    bins: np.ndarray | None = None,
) -> dict:
    """1/Vmax HI mass function: phi(logM) = (1/dlogM) Σ_i 1/V_i, Poisson errors per bin.

    ``vmax_per_sr`` from :func:`vmax_1vmax`; ``area_sr`` is the survey solid angle so the volume
    per source is ``area_sr * vmax_per_sr``. Returns bin centres, phi (Mpc^-3 dex^-1), its error,
    and per-bin counts.
    """
    lm = np.asarray(log_mhi, float)
    v = np.asarray(vmax_per_sr, float) * area_sr
    if bins is None:
        bins = np.arange(6.5, 11.01, 0.25)
    dlog = float(bins[1] - bins[0])
    which = np.digitize(lm, bins) - 1
    nb = len(bins) - 1
    phi = np.zeros(nb)
    err = np.zeros(nb)
    cnt = np.zeros(nb, int)
    for b in range(nb):
        sel = (which == b) & (v > 0)
        if sel.any():
            w = 1.0 / v[sel]
            phi[b] = w.sum() / dlog
            err[b] = np.sqrt((w**2).sum()) / dlog  # Poisson-weighted
            cnt[b] = int(sel.sum())
    centres = 0.5 * (bins[:-1] + bins[1:])
    return {"logm": centres, "phi": phi, "phi_err": err, "counts": cnt, "dlog": dlog}


def schechter(logm: np.ndarray, log_phi_star: float, log_m_star: float, alpha: float) -> np.ndarray:
    r"""Schechter HIMF in log-mass: :math:`\phi\,d\log M = \ln 10\,\phi^*\,x^{\alpha+1} e^{-x}`,
    with :math:`x = 10^{\log M - \log M^*}` (Mpc^-3 dex^-1)."""
    x = 10.0 ** (np.asarray(logm, float) - log_m_star)
    return np.log(10.0) * (10.0**log_phi_star) * x ** (alpha + 1.0) * np.exp(-x)


def fit_schechter(h: dict, *, p0: tuple = (-2.5, 9.9, -1.3)) -> dict:
    """Least-squares Schechter fit (log phi*, logM*, alpha) to a 1/Vmax HIMF (log-space)."""
    from scipy.optimize import curve_fit

    lm = h["logm"]
    phi = h["phi"]
    err = h["phi_err"]
    good = (phi > 0) & (h["counts"] >= 3) & np.isfinite(err) & (err > 0)
    if good.sum() < 4:
        return {
            "log_phi_star": np.nan,
            "log_m_star": np.nan,
            "alpha": np.nan,
            "n_bins": int(good.sum()),
        }

    def model(lmv, lps, lms, a):
        return np.log10(np.maximum(schechter(lmv, lps, lms, a), 1e-30))

    sigma = err[good] / (np.log(10.0) * phi[good])  # error on log10(phi)
    try:
        popt, pcov = curve_fit(
            model, lm[good], np.log10(phi[good]), p0=p0, sigma=sigma, maxfev=20000
        )
    except (RuntimeError, ValueError):
        return {
            "log_phi_star": np.nan,
            "log_m_star": np.nan,
            "alpha": np.nan,
            "n_bins": int(good.sum()),
        }
    perr = np.sqrt(np.diag(pcov))
    return {
        "log_phi_star": float(popt[0]),
        "log_m_star": float(popt[1]),
        "alpha": float(popt[2]),
        "log_m_star_err": float(perr[1]),
        "alpha_err": float(perr[2]),
        "n_bins": int(good.sum()),
    }


def synthetic_environment_catalogue(
    n: int = 150000,
    *,
    void_frac: float = 0.15,
    wall_logmstar: float = 9.95,
    void_logmstar: float = 9.70,  # voids: suppressed knee (Moorman+2014)
    wall_alpha: float = -1.25,
    void_alpha: float = -1.45,  # voids: steeper faint end
    area_sr: float = 2.0,
    seed: int = 0,
) -> dict:
    """A flux-limited mock with environment-dependent Schechter HIMFs --- the recover-a-known.

    Draws HI masses from a wall or void Schechter (the two differing in knee mass and faint
    slope, the injected signal), places them uniformly in a comoving shell, converts to an
    integrated flux via M_HI ∝ flux·D^2, and applies the survey flux limit. The pipeline must
    recover the two injected Schechter parameter sets separately by environment.
    """
    rng = np.random.default_rng(seed)
    is_void = rng.random(n) < void_frac

    def draw_logm(size, lms, a):
        # inverse-CDF-free rejection draw from the Schechter in logM over [6.5, 11.2]
        grid = np.linspace(6.5, 11.2, 400)
        pdf = schechter(grid, 0.0, lms, a)
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        return np.interp(rng.random(size), cdf, grid)

    logm = np.empty(n)
    logm[is_void] = draw_logm(int(is_void.sum()), void_logmstar, void_alpha)
    logm[~is_void] = draw_logm(int((~is_void).sum()), wall_logmstar, wall_alpha)
    # uniform in comoving volume out to d_max: d ∝ U^(1/3)
    d_max = 386.0
    dist = d_max * rng.random(n) ** (1.0 / 3.0)
    dist = np.clip(dist, 5.0, d_max)
    # integrated flux: log S = logM - log(2.356e5) - 2 log D  (M_HI = 2.356e5 D^2 S)
    flux = 10.0 ** (logm - np.log10(2.356e5) - 2.0 * np.log10(dist))
    z = dist * H0 / C_KM_S
    detected = flux >= FASHI_FLUX_LIMIT
    return {
        "log_mhi": logm[detected],
        "dist_mpc": dist[detected],
        "flux": flux[detected],
        "z": z[detected],
        "is_void": is_void[detected],
        "area_sr": area_sr,
        "truth": {
            "wall": (wall_logmstar, wall_alpha),
            "void": (void_logmstar, void_alpha),
        },
    }


def fetch_fashi_dr1() -> dict:  # pragma: no cover - network
    """Fetch FASHI DR1 (VizieR J/other/SCPMA/67.19511/table2): 41,741 HI sources."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "cz", "z", "W50", "Ssum", "Dist", "logMass"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs(FASHI_DR1_VIZIER)[0]
    return {
        "ra": np.asarray(t["RAJ2000"], float),
        "dec": np.asarray(t["DEJ2000"], float),
        "cz": np.asarray(t["cz"], float),
        "z": np.asarray(t["z"], float),
        "w50": np.asarray(t["W50"], float),
        "flux": np.asarray(t["Ssum"], float) / 1000.0,  # mJy km/s -> Jy km/s
        "dist_mpc": np.asarray(t["Dist"], float),
        "log_mhi": np.asarray(t["logMass"], float),
    }


def fetch_tempel_groups() -> dict:  # pragma: no cover - network
    """Fetch Tempel+2017 SDSS groups: member GroupID/Ngal (table1) + group R200/M200 (table2)."""
    from astroquery.vizier import Vizier

    vg = Vizier(columns=["GroupID", "Ngal", "RAJ2000", "DEJ2000", "Dist.c", "R200", "M200"])
    vg.ROW_LIMIT = -1
    grp = vg.get_catalogs(f"{TEMPEL_VIZIER}/table2")[0]
    vgal = Vizier(columns=["GroupID", "Ngal", "RAJ2000", "DEJ2000", "zobs"])
    vgal.ROW_LIMIT = -1
    gal = vgal.get_catalogs(f"{TEMPEL_VIZIER}/table1")[0]
    return {
        "grp_id": np.asarray(grp["GroupID"], int),
        "grp_ra": np.asarray(grp["RAJ2000"], float),
        "grp_dec": np.asarray(grp["DEJ2000"], float),
        "grp_cz": np.asarray(grp["Dist.c"], float) * H0,  # Mpc -> km/s (Hubble)
        "grp_r200": np.asarray(grp["R200"], float),
        "grp_ngal": np.asarray(grp["Ngal"], int),
        "gal_ra": np.asarray(gal["RAJ2000"], float),
        "gal_dec": np.asarray(gal["DEJ2000"], float),
        "gal_cz": np.asarray(gal["zobs"], float) * C_KM_S,
    }


def fetch_voidfinder_spheres() -> dict:  # pragma: no cover - network
    """Fetch Douglass+2023 VoidFinder maximal spheres (J/ApJS/265/7/table1)."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["x", "y", "z", "Rad", "Cosmo"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs(f"{DOUGLASS_VIZIER}/table1")[0]
    # one cosmology only (Planck2018); x,y,z,Rad are Mpc/h in the standard equatorial frame
    # (verified: RA/Dec reconstructed from x,y,z match the catalogue's own RA/Dec exactly)
    m = np.asarray([str(c).strip().startswith("Planck") for c in t["Cosmo"]])
    xyz = np.stack(
        [np.asarray(t["x"], float), np.asarray(t["y"], float), np.asarray(t["z"], float)], axis=1
    )
    return {"sphere_xyz": xyz[m], "sphere_radius": np.asarray(t["Rad"], float)[m]}


def _offset_stats(name: str, fit_a: dict, fit_b: dict) -> dict:
    """Knee offset (a - b), its combined error, and significance in sigma --- pipeline-made."""
    da = fit_a.get("log_m_star", np.nan)
    db = fit_b.get("log_m_star", np.nan)
    ea = fit_a.get("log_m_star_err", np.nan)
    eb = fit_b.get("log_m_star_err", np.nan)
    off = da - db
    err = float(np.sqrt(ea**2 + eb**2))
    return {
        f"{name}_knee_offset": round(off, 3),
        f"{name}_knee_offset_err": round(err, 3),
        f"{name}_knee_offset_sigma": round(abs(off) / err, 2) if err > 0 else None,
    }


def _himf_and_fit(cat_logm, cat_dist, cat_flux, area_sr, mask=None):
    lm = cat_logm if mask is None else cat_logm[mask]
    dd = cat_dist if mask is None else cat_dist[mask]
    ff = cat_flux if mask is None else cat_flux[mask]
    vmax = vmax_1vmax(lm, dd, ff)
    h = himf(lm, vmax, area_sr=area_sr)
    fit = fit_schechter(h)
    return h, fit


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: environment-split HIMF recover-a-known on the mock; real: FASHI DR1 x groups/voids."""
    import json

    if offline:
        cat = synthetic_environment_catalogue()
        area = cat["area_sr"]
        _h_all, fit_all = _himf_and_fit(cat["log_mhi"], cat["dist_mpc"], cat["flux"], area)
        _h_v, fit_void = _himf_and_fit(
            cat["log_mhi"], cat["dist_mpc"], cat["flux"], area, mask=cat["is_void"]
        )
        _h_w, fit_wall = _himf_and_fit(
            cat["log_mhi"], cat["dist_mpc"], cat["flux"], area, mask=~cat["is_void"]
        )
        metrics = {
            "source": "synthetic environment-split flux-limited mock",
            "is_real": False,
            "n_sources": int(cat["log_mhi"].size),
            "himf_global": {k: round(v, 3) for k, v in fit_all.items() if isinstance(v, float)},
            "himf_void": {k: round(v, 3) for k, v in fit_void.items() if isinstance(v, float)},
            "himf_wall": {k: round(v, 3) for k, v in fit_wall.items() if isinstance(v, float)},
            "void_knee_offset": round(fit_void["log_m_star"] - fit_wall["log_m_star"], 3),
            "true_void_logmstar": cat["truth"]["void"][0],
            "true_void_alpha": cat["truth"]["void"][1],
            "true_wall_logmstar": cat["truth"]["wall"][0],
            "true_wall_alpha": cat["truth"]["wall"][1],
            "true_knee_offset": round(cat["truth"]["void"][0] - cat["truth"]["wall"][0], 3),
        }
        figdata = (_h_v, _h_w, fit_void, fit_wall)
    else:  # pragma: no cover - needs the VizieR catalogues
        metrics, figdata = _real_leg()

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "fashienv_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(figdata, op / "papers" / "fashienv" / "figures")
    _write_macros(metrics, op / "papers" / "fashienv" / "generated" / "macros.tex")
    return metrics


def _real_leg():  # pragma: no cover - network + VizieR catalogues
    """FASHI DR1 x Tempel groups x Douglass voids: environment-split HIMF + R/R200 gradient."""
    fashi = fetch_fashi_dr1()
    area = 7600.0 * (np.pi / 180.0) ** 2  # FASHI DR1 sky area ~7600 deg^2 (Zhang+2024); note
    # this only sets the never-quoted phi* normalisation and cancels in every knee/slope offset
    finite = np.isfinite(fashi["log_mhi"]) & np.isfinite(fashi["dist_mpc"]) & (fashi["flux"] > 0)
    lm, dd, ff = fashi["log_mhi"][finite], fashi["dist_mpc"][finite], fashi["flux"][finite]
    ra, dec, z, cz = (fashi[k][finite] for k in ("ra", "dec", "z", "cz"))

    _h_all, fit_all = _himf_and_fit(lm, dd, ff, area)

    # voids (VoidFinder spheres, SDSS-cap overlap only). The spheres are Mpc/h, so place the
    # FASHI galaxies in the same H0-independent Mpc/h frame (h0=100) for the membership test.
    spheres = fetch_voidfinder_spheres()
    gal_xyz = comoving_xyz(ra, dec, z, h0=100.0)
    in_void = void_membership(gal_xyz, spheres["sphere_xyz"], spheres["sphere_radius"])
    # only galaxies within the void catalogue's comoving volume can be classified: a galaxy
    # outside every sphere but inside the SDSS-cap volume is a wall galaxy; a galaxy beyond that
    # volume (most of FASHI's southern/high-z sky) is unclassifiable and excluded -- the honest
    # footprint caveat, quantified by n_classifiable_void.
    classifiable = _within_void_footprint(gal_xyz, spheres["sphere_xyz"])
    _h_v, fit_void = _himf_and_fit(lm, dd, ff, area, mask=classifiable & in_void)
    _h_w, fit_wall = _himf_and_fit(lm, dd, ff, area, mask=classifiable & ~in_void)

    # group-member vs field HIMF (Tempel groups; a galaxy within R200+dv of any group is a
    # "group member", the rest are "field") -- the density-split HIMF, cleaner than the
    # selection-biased median-HI-vs-radius the plan proposed (dropped; see the module docstring)
    grp = fetch_tempel_groups()
    gidx = assign_groups(ra, dec, cz, grp["grp_ra"], grp["grp_dec"], grp["grp_cz"], grp["grp_r200"])
    in_group = gidx >= 0
    # restrict the field/group split to the same SDSS-cap footprint (else "field" is dominated
    # by FASHI's southern sky where no Tempel groups exist -- an apples-to-oranges comparison)
    cap = classifiable
    _h_g, fit_group = _himf_and_fit(lm, dd, ff, area, mask=cap & in_group)
    _h_f, fit_field = _himf_and_fit(lm, dd, ff, area, mask=cap & ~in_group)

    metrics = {
        "source": "FASHI DR1 (VizieR J/other/SCPMA/67.19511) x Tempel+2017 groups x Douglass+2023 voids",
        "is_real": True,
        "n_sources": int(lm.size),
        "n_classifiable_void": int(classifiable.sum()),
        "n_in_void": int((classifiable & in_void).sum()),
        "n_group_members": int((cap & in_group).sum()),
        "himf_global": {k: round(v, 3) for k, v in fit_all.items() if isinstance(v, float)},
        "himf_void": {k: round(v, 3) for k, v in fit_void.items() if isinstance(v, float)},
        "himf_wall": {k: round(v, 3) for k, v in fit_wall.items() if isinstance(v, float)},
        "himf_group": {k: round(v, 3) for k, v in fit_group.items() if isinstance(v, float)},
        "himf_field": {k: round(v, 3) for k, v in fit_field.items() if isinstance(v, float)},
        **_offset_stats("void", fit_void, fit_wall),
        **_offset_stats("group", fit_group, fit_field),
        "dr2_followon": "swap fetch_fashi_dr1 -> DR2 table when it publishes (~Aug 2026)",
    }
    return metrics, (_h_v, _h_w, fit_void, fit_wall)


def _within_void_footprint(gal_xyz, sphere_xyz, pad_mpc=20.0):  # pragma: no cover - network path
    """Galaxies inside the convex extent of the void survey (bounding box + pad); crude but honest."""
    lo = sphere_xyz.min(axis=0) - pad_mpc
    hi = sphere_xyz.max(axis=0) + pad_mpc
    return np.all((gal_xyz >= lo) & (gal_xyz <= hi), axis=1)


def _figure(figdata, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    h_v, h_w, fit_v, fit_w = figdata
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    for h, fit, c, lab in ((h_w, fit_w, "C0", "wall"), (h_v, fit_v, "C3", "void")):
        good = h["counts"] >= 3
        ax.errorbar(
            h["logm"][good],
            h["phi"][good],
            yerr=h["phi_err"][good],
            fmt="o",
            color=c,
            ms=4,
            label=f"{lab} ($\\log M^*$={fit.get('log_m_star', float('nan')):.2f}, "
            f"$\\alpha$={fit.get('alpha', float('nan')):.2f})",
        )
        if np.isfinite(fit.get("log_m_star", np.nan)):
            xs = np.linspace(7.0, 10.8, 100)
            ax.plot(
                xs,
                schechter(xs, fit["log_phi_star"], fit["log_m_star"], fit["alpha"]),
                color=c,
                lw=1,
            )
    ax.set(
        xlabel=r"$\log_{10}(M_{\rm HI}/M_\odot)$",
        ylabel=r"$\phi$ (Mpc$^{-3}$ dex$^{-1}$)",
        yscale="log",
        title="Environment-split HI mass function",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fashienv.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    def g(d: str, key: str | None) -> str:
        sub = m.get(d)
        v = sub.get(key) if isinstance(sub, dict) and key is not None else sub
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "feReal" if m.get("is_real") else "feSyn"
    lines = [
        "% Auto-generated by jansky_research.fashienv._write_macros -- do not edit.",
        "% Synthetic (feSyn*) and real (feReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under feReal* (an offline rebuild resets feReal* to placeholders by design).",
        rf"\newcommand{{\feSource}}{{{m['source']}}}",
        rf"\newcommand{{\feN}}{{{m['n_sources']}}}",
    ]
    for ns in ("feSyn", "feReal"):
        live = ns == pref
        for macro, d, key in (
            ("VoidLogMStar", "himf_void", "log_m_star"),
            ("VoidAlpha", "himf_void", "alpha"),
            ("WallLogMStar", "himf_wall", "log_m_star"),
            ("WallAlpha", "himf_wall", "alpha"),
            ("GlobalLogMStar", "himf_global", "log_m_star"),
            ("GlobalAlpha", "himf_global", "alpha"),
            ("GroupLogMStar", "himf_group", "log_m_star"),
            ("FieldLogMStar", "himf_field", "log_m_star"),
            ("VoidKneeOffset", "void_knee_offset", None),
            ("VoidKneeErr", "void_knee_offset_err", None),
            ("VoidKneeSigma", "void_knee_offset_sigma", None),
            ("GroupKneeOffset", "group_knee_offset", None),
            ("GroupKneeErr", "group_knee_offset_err", None),
            ("GroupKneeSigma", "group_knee_offset_sigma", None),
            ("NInVoid", "n_in_void", None),
            ("NGroupMembers", "n_group_members", None),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(d, key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="FASHI DR1 environment-split HI mass function.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
