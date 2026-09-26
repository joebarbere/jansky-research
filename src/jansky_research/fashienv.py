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
    "void_membership_holes",
    "cmb_to_helio_cz",
    "random_void_positions",
    "assign_groups",
    "clustercentric_radius",
    "vmax_1vmax",
    "vmax_from_catalogue",
    "load_fashi_dr2",
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
FASHI_DR2_PORTAL = "https://zcp521.github.io/fashi.html"  # links the CSTCloud share below
# DR2 (arXiv:2606.31539; Sci. China PMA 2026-09-08) is NOT on VizieR; the tables are on a public
# CSTCloud share. Download = POST shareGetInfo (list files) then shareDownloadRequest (a signed
# URL valid 24 h). Cached at data/fashi_dr2/.
FASHI_DR2_SHARE = "XfSZ82LQc4"
FASHI_DR2_TABLE = "Table2_FASHI_DR2_Extragalactic_Hi_Source_Catalog.csv"
FASHI_DR2_AREA_DEG2 = 19482.0  # the area DR2's own Vmax column is computed over (paper Sect. 5)
FASHI_DR2_C_MIN = 0.5  # "only galaxies above the 50% flux completeness limit" (DR2 HIMF sample)

C_KM_S = 299792.458
H0 = 70.0  # km/s/Mpc, the FASHI DR1 distance convention (h70)
FASHI_FLUX_LIMIT = 0.30  # Jy km/s integrated-flux limit for the 1/Vmax weighting. A single-cut
# order-of-magnitude value (FAST is deeper than ALFALFA's ~0.7 Jy km/s); it is NOT FASHI's full
# completeness function, so the ABSOLUTE HIMF (esp. the faint-end slope) is only approximate ---
# the relative environment offsets, which share this cut, are what the slice reports.

# Planck2018 deceleration: q0 = Om/2 - OL = 0.315/2 - 0.685 = -0.527 (matches the Douglass frame)
_Q0 = -0.527


def _comoving_distance_mpc(z: np.ndarray, h0: float = H0, q0: float = _Q0) -> np.ndarray:
    r"""Low-z comoving distance, 2nd-order in z: :math:`(c/H_0)\,z\,[1 - (1+q_0)z/2]`.

    Uses the Planck2018 :math:`q_0=-0.527` (not the EdS 0.5) so the void-membership geometry
    matches the Douglass+2023 catalogue's own cosmology; the residual vs the full integral is
    <0.02% at z<0.09. With ``h0=100`` the result is Mpc/h (H0-independent) --- the Douglass
    frame; with ``h0=70`` it is physical Mpc (the FASHI DR1 convention).
    """
    z = np.asarray(z, float)
    return (C_KM_S / h0) * z * (1.0 - 0.5 * (1.0 + q0) * z)


def comoving_xyz(
    ra_deg: np.ndarray, dec_deg: np.ndarray, z: np.ndarray, *, h0: float = H0, q0: float = _Q0
) -> np.ndarray:
    """(RA, Dec, z) -> comoving Cartesian, the frame for void-sphere membership.

    Use ``h0=100`` to match the Douglass void spheres' Mpc/h frame; ``h0=70`` for physical Mpc.
    """
    d = _comoving_distance_mpc(z, h0, q0)
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


def void_membership_holes(
    gal_xyz: np.ndarray, hole_xyz: np.ndarray, hole_radius: np.ndarray
) -> np.ndarray:
    """True where a galaxy lies inside ANY hole (sphere); KD-tree on galaxies, one query per hole.

    Same test as :func:`void_membership`, scaled for the full VoidFinder hole list (a void is the
    union of all its holes, Douglass+2023 table2) and for repeated calls in the random-void null.
    """
    from scipy.spatial import cKDTree

    gal = np.asarray(gal_xyz, float)
    inside = np.zeros(len(gal), bool)
    if len(gal) == 0 or len(hole_xyz) == 0:
        return inside
    tree = cKDTree(gal)
    for hits in tree.query_ball_point(
        np.asarray(hole_xyz, float), r=np.asarray(hole_radius, float)
    ):
        if hits:
            inside[hits] = True
    return inside


# Sun's velocity relative to the CMB (Planck 2018 dipole): 369.82 km/s toward (l, b) =
# (264.021, 48.253) deg. A galaxy's heliocentric cz = CMB-frame cz - v_sun . n_hat.
CMB_DIPOLE_KMS = 369.82
CMB_DIPOLE_LB = (264.021, 48.253)


def cmb_to_helio_cz(ra_deg: np.ndarray, dec_deg: np.ndarray, cz_cmb: np.ndarray) -> np.ndarray:
    """Convert CMB-frame recession velocities to heliocentric ones (km/s)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    apex = SkyCoord(l=CMB_DIPOLE_LB[0] * u.deg, b=CMB_DIPOLE_LB[1] * u.deg, frame="galactic")
    pos = SkyCoord(np.asarray(ra_deg, float) * u.deg, np.asarray(dec_deg, float) * u.deg)
    cos_t = np.cos(pos.separation(apex).radian)
    return np.asarray(cz_cmb, float) - CMB_DIPOLE_KMS * cos_t


def _unit(xyz: np.ndarray) -> np.ndarray:
    return xyz / np.linalg.norm(xyz, axis=-1, keepdims=True)


def _sky_cells(ra_deg: np.ndarray, dec_deg: np.ndarray, cell_deg: float) -> set:
    """Occupied equal-area cells (RA bins x sin(Dec) bins) -- a footprint map without healpy."""
    ds = np.radians(cell_deg)
    nra = int(round(360.0 / cell_deg))
    ira = (np.asarray(ra_deg, float) // cell_deg).astype(int) % nra
    idec = np.floor(np.sin(np.radians(np.asarray(dec_deg, float))) / ds).astype(int)
    return set(zip(ira.tolist(), idec.tolist(), strict=True))


def _rotation_to(a: np.ndarray, b: np.ndarray, roll: float) -> np.ndarray:
    """Rotation matrix taking unit vector ``a`` onto unit vector ``b``, then rolling about ``b``."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    r1 = np.eye(3) + vx + vx @ vx / (1.0 + c) if c > -1 + 1e-12 else -np.eye(3)
    k = b
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    r2 = np.eye(3) + np.sin(roll) * kx + (1 - np.cos(roll)) * kx @ kx
    return r2 @ r1


def random_void_positions(
    hole_xyz: np.ndarray,
    void_id: np.ndarray,
    footprint_ra: np.ndarray,
    footprint_dec: np.ndarray,
    rng: np.random.Generator,
    *,
    cell_deg: float = 2.0,
) -> np.ndarray:
    """Move every void RIGIDLY (all its holes together) to a random direction inside the footprint.

    Each void keeps its radial distance, size and internal shape; only its sky position (and a
    random roll about the line of sight) changes. The target direction is uniform on the sphere,
    accepted only where the footprint (cells occupied by ``footprint_ra/dec``) has coverage. This
    is the random-void null: it preserves the void population's radial profile and volume
    fraction, so any void-wall offset it produces is the estimator's and geometry's own bias.
    """
    cells = _sky_cells(footprint_ra, footprint_dec, cell_deg)
    hole_xyz = np.asarray(hole_xyz, float)
    out = np.empty_like(hole_xyz)
    for vid in np.unique(void_id):
        sel = np.flatnonzero(void_id == vid)
        centre = _unit(hole_xyz[sel].mean(axis=0))
        while True:
            z = rng.uniform(-1, 1)
            phi = rng.uniform(0, 2 * np.pi)
            tgt = np.array([np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z])
            ra = np.degrees(np.arctan2(tgt[1], tgt[0])) % 360.0
            dec = np.degrees(np.arcsin(tgt[2]))
            if next(iter(_sky_cells(np.array([ra]), np.array([dec]), cell_deg))) in cells:
                break
        rot = _rotation_to(centre, tgt, rng.uniform(0, 2 * np.pi))
        out[sel] = hole_xyz[sel] @ rot.T
    return out


def void_members(
    gal_xyz: np.ndarray, hole_xyz: np.ndarray, hole_radius: np.ndarray, void_id: np.ndarray
) -> dict:
    """Per-void member galaxy indices, and how many voids contain each galaxy."""
    from scipy.spatial import cKDTree

    tree = cKDTree(np.asarray(gal_xyz, float))
    members: dict[int, set] = {}
    for c, r, vid in zip(
        np.asarray(hole_xyz, float), np.asarray(hole_radius, float), void_id, strict=True
    ):
        hits = tree.query_ball_point(c, r)
        if hits:
            members.setdefault(int(vid), set()).update(hits)
    n_in = np.zeros(len(gal_xyz), int)
    for m in members.values():
        n_in[list(m)] += 1
    return {"members": {k: np.fromiter(v, int) for k, v in members.items()}, "n_voids": n_in}


def void_jackknife(
    cat: dict,
    area: float,
    in_void: np.ndarray,
    classifiable: np.ndarray,
    vm: dict,
    weights: dict,
) -> dict:
    """Delete-one-VOID jackknife of the void-wall knee offset under several weightings at once.

    ``vm`` is :func:`void_members` output; deleting void k turns every galaxy that ONLY void k
    contains into a wall galaxy. ``weights`` maps a label to ``(base_mask, vmax_or_None)``; the
    same deletions are applied to every weighting, so the jackknife error of their DIFFERENCE is
    a paired error (what the weighting comparison needs), not a quadrature guess.
    """
    lm, dd, ff = cat["log_mhi"], cat["dist_mpc"], cat["flux"]
    offs: dict[str, list] = {k: [] for k in weights}
    occupied = [k for k, m in vm["members"].items() if np.any(classifiable[m])]
    for k in occupied:
        m = vm["members"][k]
        iv = in_void.copy()
        iv[m[vm["n_voids"][m] == 1]] = False
        row = {}
        for lab, (base, vmax) in weights.items():
            _hv, fv = _himf_and_fit(lm, dd, ff, area, mask=base & classifiable & iv, vmax=vmax)
            _hw, fw = _himf_and_fit(lm, dd, ff, area, mask=base & classifiable & ~iv, vmax=vmax)
            row[lab] = fv.get("log_m_star", np.nan) - fw.get("log_m_star", np.nan)
        if all(np.isfinite(v) for v in row.values()):
            for lab, v in row.items():
                offs[lab].append(v)

    def jk(a: np.ndarray) -> float:
        n = a.size
        return float(np.sqrt((n - 1) / n * np.sum((a - a.mean()) ** 2))) if n > 2 else np.nan

    out: dict = {"n_voids_occupied": len(occupied), "n_ok": len(next(iter(offs.values())))}
    arrs = {k: np.asarray(v) for k, v in offs.items()}
    for lab, a in arrs.items():
        out[f"{lab}_jackknife_err"] = round(jk(a), 4)
        if a.size:
            out[f"{lab}_mean_offset"] = round(float(a.mean()), 4)
            out[f"{lab}_min_offset"] = round(float(a.min()), 4)
            out[f"{lab}_max_offset"] = round(float(a.max()), 4)
    labs = list(arrs)
    if len(labs) == 2:
        out[f"{labs[1]}_minus_{labs[0]}_jackknife_err"] = round(
            jk(arrs[labs[1]] - arrs[labs[0]]), 4
        )
    return out


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


def vmax_from_catalogue(
    vmax_mpc3: np.ndarray,
    completeness: np.ndarray,
    *,
    survey_area_deg2: float = FASHI_DR2_AREA_DEG2,
) -> np.ndarray:
    """Per-steradian EFFECTIVE volume from a catalogue's own Vmax and completeness: C * Vmax / Omega.

    FASHI DR2 publishes, per source, the comoving Vmax over its survey area and the completeness
    C at that source's flux and line width; the survey's HIMF weights each galaxy by
    1/(C * Vmax). Returned per steradian so it drops into :func:`himf` exactly like
    :func:`vmax_1vmax` (the caller's ``area_sr`` multiplies it back). Non-finite or non-positive
    inputs give 0, which :func:`himf` excludes.
    """
    v = np.asarray(vmax_mpc3, float)
    c = np.asarray(completeness, float)
    omega = survey_area_deg2 * (np.pi / 180.0) ** 2
    good = np.isfinite(v) & np.isfinite(c) & (v > 0) & (c > 0)
    return np.where(good, c * v, 0.0) / omega


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
    resid = (np.log10(phi[good]) - model(lm[good], *popt)) / sigma
    dof = max(int(good.sum()) - 3, 1)
    return {
        "log_phi_star": float(popt[0]),
        "log_m_star": float(popt[1]),
        "alpha": float(popt[2]),
        "log_m_star_err": float(perr[1]),
        "alpha_err": float(perr[2]),
        "n_bins": int(good.sum()),
        # Fit quality, committed so a missed bin is visible in the evidence, not only the figure.
        "red_chi2": float(np.sum(resid**2) / dof),
        "corr_mstar_alpha": float(pcov[1, 2] / (perr[1] * perr[2]))
        if perr[1] * perr[2] > 0
        else np.nan,
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


def load_fashi_dr2(path: str | Path) -> dict:
    """Read the FASHI DR2 Table 2 CSV into the same keys as :func:`fetch_fashi_dr1`, plus
    ``completeness`` and ``vmax_mpc3``.

    Column map (DR1 VizieR -> DR2 CSV): RAJ2000->ra, DEJ2000->dec, cz->v_opt, z->z_opt,
    W50->W_50, Ssum->S_sum (both mJy km/s), Dist->distance, logMass->mass. Note DR2's
    ``distance`` convention differs from DR1's (matched sources are ~2.5% farther in DR2), so
    DR1 and DR2 absolute masses are not interchangeable; relative offsets within one release are.
    """
    import csv

    cols = {
        "ra": "ra", "dec": "dec", "cz": "v_opt", "z": "z_opt", "w50": "W_50", "flux": "S_sum",
        "dist_mpc": "distance", "log_mhi": "mass", "completeness": "completeness",
        "vmax_mpc3": "Vmax",
    }  # fmt: skip
    vals: dict[str, list[float]] = {k: [] for k in cols}
    with Path(path).open(newline="") as fh:
        for r in csv.DictReader(fh):
            for k, c in cols.items():
                try:
                    vals[k].append(float(r[c]))
                except (TypeError, ValueError):
                    vals[k].append(np.nan)
    out = {k: np.asarray(v, float) for k, v in vals.items()}
    out["flux"] = out["flux"] / 1000.0  # mJy km/s -> Jy km/s, as for DR1
    return out


def fetch_fashi_dr2(cache_dir: str | Path | None = None) -> dict:  # pragma: no cover - network
    """FASHI DR2 Table 2 (156,411 sources) from the public CSTCloud share, cached on disk."""
    import json
    import urllib.request

    from . import data as _data

    d = Path(cache_dir) if cache_dir else _data.data_dir() / "fashi_dr2"
    target = d / FASHI_DR2_TABLE
    if not target.exists():
        d.mkdir(parents=True, exist_ok=True)

        def _post(api: str, body: dict) -> dict:
            req = urllib.request.Request(
                f"https://pan.cstcloud.cn/s/api/{api}",
                data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())

        info = _post("shareGetInfo", {"shareId": FASHI_DR2_SHARE})
        fid = next(
            (f["fid"] for f in _walk_share_files(info) if f.get("name") == FASHI_DR2_TABLE), None
        )
        if fid is None:
            raise RuntimeError("FASHI DR2 table not found in the CSTCloud share listing")
        dl = _post("shareDownloadRequest", {"shareId": FASHI_DR2_SHARE, "fid": fid})
        url = next(v for v in _walk_values(dl) if isinstance(v, str) and v.startswith("http"))
        urllib.request.urlretrieve(url, target)
    return load_fashi_dr2(target)


def _walk_share_files(obj):  # pragma: no cover - network helper
    if isinstance(obj, dict):
        if "fid" in obj and "name" in obj:
            yield obj
        for v in obj.values():
            yield from _walk_share_files(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk_share_files(v)


def _walk_values(obj):  # pragma: no cover - network helper
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _walk_values(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk_values(v)
    else:
        yield obj


def fetch_tempel_groups() -> dict:  # pragma: no cover - network
    """Fetch Tempel+2017 SDSS groups: member GroupID/Ngal (table1) + group R200/M200 (table2)."""
    from astroquery.vizier import Vizier

    vg = Vizier(columns=["GroupID", "Ngal", "RAJ2000", "DEJ2000", "zcmb", "R200", "M200"])
    vg.ROW_LIMIT = -1
    grp = vg.get_catalogs(f"{TEMPEL_VIZIER}/table2")[0]
    vgal = Vizier(columns=["GroupID", "Ngal", "RAJ2000", "DEJ2000", "zobs"])
    vgal.ROW_LIMIT = -1
    gal = vgal.get_catalogs(f"{TEMPEL_VIZIER}/table1")[0]
    return {
        "grp_id": np.asarray(grp["GroupID"], int),
        "grp_ra": np.asarray(grp["RAJ2000"], float),
        "grp_dec": np.asarray(grp["DEJ2000"], float),
        # FASHI v_opt is heliocentric; Tempel gives CMB-frame redshifts. The old
        # Dist.c * H0 used H0=70 against the catalogue's 67.8 and the wrong frame (referee
        # 2026-09-26, finding 10: a few hundred km/s against a +/-1000 km/s window).
        "grp_cz": cmb_to_helio_cz(
            np.asarray(grp["RAJ2000"], float),
            np.asarray(grp["DEJ2000"], float),
            np.asarray(grp["zcmb"], float) * C_KM_S,
        ),
        "grp_r200": np.asarray(grp["R200"], float),
        "grp_ngal": np.asarray(grp["Ngal"], int),
        "gal_ra": np.asarray(gal["RAJ2000"], float),
        "gal_dec": np.asarray(gal["DEJ2000"], float),
        "gal_cz": np.asarray(gal["zobs"], float) * C_KM_S,
    }


def fetch_voidfinder_spheres(*, all_holes: bool = True) -> dict:  # pragma: no cover - network
    """Douglass+2023 VoidFinder voids, Planck2018 cosmology, in Mpc/h equatorial Cartesian.

    ``all_holes=True`` (default since the DR2 referee round) returns every hole of every void
    (table2): a VoidFinder void is the UNION of its holes, so using only the maximal sphere
    (table1) classifies void-outskirt galaxies as wall. ``void_id`` groups holes into voids, for
    the delete-one-void jackknife and for moving whole voids in the random-void null.
    """
    from astroquery.vizier import Vizier

    table = "table2" if all_holes else "table1"
    v = Vizier(columns=["x", "y", "z", "Rad", "Cosmo", "void"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs(f"{DOUGLASS_VIZIER}/{table}")[0]
    # (verified 2026-07: RA/Dec reconstructed from x,y,z match the catalogue's own RA/Dec)
    m = np.asarray([str(c).strip().startswith("Planck") for c in t["Cosmo"]])
    xyz = np.stack(
        [np.asarray(t["x"], float), np.asarray(t["y"], float), np.asarray(t["z"], float)], axis=1
    )
    return {
        "sphere_xyz": xyz[m],
        "sphere_radius": np.asarray(t["Rad"], float)[m],
        "void_id": np.asarray(t["void"], int)[m],
        "table": table,
    }


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


def _himf_and_fit(cat_logm, cat_dist, cat_flux, area_sr, mask=None, vmax=None):
    """HIMF + Schechter fit. ``vmax`` (per sr, e.g. :func:`vmax_from_catalogue`) overrides the
    single-flux-cut :func:`vmax_1vmax` weighting when given."""
    lm = cat_logm if mask is None else cat_logm[mask]
    if vmax is not None:
        vmax = vmax if mask is None else vmax[mask]
    else:
        dd = cat_dist if mask is None else cat_dist[mask]
        ff = cat_flux if mask is None else cat_flux[mask]
        vmax = vmax_1vmax(lm, dd, ff)
    h = himf(lm, vmax, area_sr=area_sr)
    fit = fit_schechter(h)
    return h, fit


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: environment-split HIMF recover-a-known on the mock; real: FASHI DR1 x groups/voids."""

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
    from .report import write_results

    write_results(metrics, op / "results" / "fashienv_metrics.json")
    _figure(figdata, op / "papers" / "fashienv" / "figures")
    _write_macros(metrics, op / "papers" / "fashienv" / "generated" / "macros.tex")
    return metrics


def void_jackknife_offset(
    lm: np.ndarray,
    dd: np.ndarray,
    ff: np.ndarray,
    area: float,
    *,
    gal_xyz: np.ndarray,
    sphere_xyz: np.ndarray,
    sphere_radius: np.ndarray,
    classifiable: np.ndarray,
    vmax: np.ndarray | None = None,
) -> dict:  # pragma: no cover - real leg only
    """Delete-one-void jackknife on the void-wall knee offset.

    The fit errors this paper quotes are Poisson counting noise *within one realization* of
    large-scale structure. The void bin is a modest number of voids in a single SDSS volume,
    so the dominant uncertainty on a void HIMF knee is how much voids differ from each other,
    which counting noise cannot see -- structurally the same error as quoting a bootstrap SE
    for a quantity the bootstrap resamples inside (see CLAUDE.md).

    Dropping each void sphere in turn and refitting gives the sample variance directly. The
    jackknife error is sqrt((n-1)/n * sum (x_i - xbar)^2).
    """
    n_v = int(sphere_radius.size)
    # Membership computed once as a (galaxy, sphere) matrix: deleting a void is then a column
    # mask rather than a full recomputation, which turns ~37 minutes into ~1.
    d = np.linalg.norm(gal_xyz[:, None, :] - sphere_xyz[None, :, :], axis=2)
    member = d <= sphere_radius[None, :]
    del d
    # Only voids that actually contain a classifiable galaxy can change the answer; deleting an
    # empty one is a no-op, and including those would deflate the jackknife spread toward zero.
    occupied = np.flatnonzero(member[classifiable].any(axis=0))
    offsets = []
    for k in occupied:
        keep = np.ones(n_v, bool)
        keep[k] = False
        in_void_k = member[:, keep].any(axis=1)
        _hv, fv = _himf_and_fit(lm, dd, ff, area, mask=classifiable & in_void_k, vmax=vmax)
        _hw, fw = _himf_and_fit(lm, dd, ff, area, mask=classifiable & ~in_void_k, vmax=vmax)
        if np.isfinite(fv.get("log_m_star", np.nan)) and np.isfinite(fw.get("log_m_star", np.nan)):
            offsets.append(float(fv["log_m_star"] - fw["log_m_star"]))
    if len(offsets) < 3:
        return {
            "n_voids": n_v,
            "n_occupied": int(occupied.size),
            "n_ok": len(offsets),
            "jackknife_err": None,
        }
    a = np.asarray(offsets)
    mean = float(a.mean())
    err = float(np.sqrt((len(a) - 1) / len(a) * np.sum((a - mean) ** 2)))
    return {
        "n_voids": n_v,
        "n_occupied": int(occupied.size),
        "n_ok": int(a.size),
        "mean_offset": round(mean, 4),
        "jackknife_err": round(err, 4),
        "min_offset": round(float(a.min()), 4),
        "max_offset": round(float(a.max()), 4),
    }


def _environment_split(cat: dict, base: np.ndarray, vmax, area: float, env: dict) -> dict:
    """Void/wall and group/field HIMF fits + knee offsets for one catalogue and one weighting."""
    lm, dd, ff = cat["log_mhi"], cat["dist_mpc"], cat["flux"]
    cl, iv, ig = env["classifiable"], env["in_void"], env["in_group"]
    _h_all, f_all = _himf_and_fit(lm, dd, ff, area, mask=base, vmax=vmax)
    h_v, f_v = _himf_and_fit(lm, dd, ff, area, mask=base & cl & iv, vmax=vmax)
    h_w, f_w = _himf_and_fit(lm, dd, ff, area, mask=base & cl & ~iv, vmax=vmax)
    _h_g, f_g = _himf_and_fit(lm, dd, ff, area, mask=base & cl & ig, vmax=vmax)
    _h_f, f_f = _himf_and_fit(lm, dd, ff, area, mask=base & cl & ~ig, vmax=vmax)

    def rnd(f: dict) -> dict:
        return {
            k: (round(v, 3) if isinstance(v, float) else v)
            for k, v in f.items()
            if isinstance(v, (float, int))
        }

    n = int(base.sum())
    n_wall = int((base & cl & ~iv).sum())
    return {
        "n_sources": n,
        "n_classifiable_void": int((base & cl).sum()),
        "n_in_void": int((base & cl & iv).sum()),
        "n_wall": n_wall,
        # The comparison bin's share of the whole sample: the bounding-box classifiable region
        # makes "wall" close to "everything", which the paper must quantify, not assert.
        "wall_pct_of_sample": round(100.0 * n_wall / n, 1) if n else None,
        "wall_minus_global_logmstar": round(
            f_w.get("log_m_star", np.nan) - f_all.get("log_m_star", np.nan), 3
        ),
        "n_group_members": int((base & cl & ig).sum()),
        "n_field": int((base & cl & ~ig).sum()),
        "himf_global": rnd(f_all),
        "himf_void": rnd(f_v),
        "himf_wall": rnd(f_w),
        "himf_group": rnd(f_g),
        "himf_field": rnd(f_f),
        **_offset_stats("void", f_v, f_w),
        **_offset_stats("group", f_g, f_f),
        "_fig": (h_v, h_w, f_v, f_w),
    }


def _environments(cat: dict, voids: dict, grp: dict, *, q0: float = _Q0) -> dict:
    """Void membership over all holes (Douglass Mpc/h frame), classifiability, group membership."""
    xyz = comoving_xyz(cat["ra"], cat["dec"], cat["z"], h0=100.0, q0=q0)
    gidx = assign_groups(
        cat["ra"],
        cat["dec"],
        cat["cz"],
        grp["grp_ra"],
        grp["grp_dec"],
        grp["grp_cz"],
        grp["grp_r200"],
    )
    return {
        "xyz": xyz,
        "in_void": void_membership_holes(xyz, voids["sphere_xyz"], voids["sphere_radius"]),
        "classifiable": _within_void_footprint(xyz, voids["sphere_xyz"]),
        "in_group": gidx >= 0,
    }


def _match_releases(
    dr1: dict, dr2: dict, *, sep_arcmin: float = 1.5, dv_kms: float = 100.0
) -> dict:
    """Match DR1 to DR2 by position and velocity; report the matched fraction and distance ratio."""
    from scipy.spatial import cKDTree

    def xyz(c):
        r, d = np.radians(c["ra"]), np.radians(c["dec"])
        return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])

    chord = 2.0 * np.sin(np.radians(sep_arcmin / 60.0) / 2.0)
    hits = cKDTree(xyz(dr2)).query_ball_point(xyz(dr1), r=chord)
    ratio, n_match = [], 0
    for i, js in enumerate(hits):
        if not js:
            continue
        js = np.asarray(js)
        dv = np.abs(dr2["cz"][js] - dr1["cz"][i])
        if dv.min() <= dv_kms:
            n_match += 1
            j = js[int(np.argmin(dv))]
            ratio.append(dr2["dist_mpc"][j] / dr1["dist_mpc"][i])
    r = np.asarray(ratio)
    return {
        "n_dr1": int(dr1["ra"].size),
        "n_matched": n_match,
        "matched_pct": round(100.0 * n_match / max(dr1["ra"].size, 1), 1),
        "median_dist_ratio_dr2_over_dr1": round(float(np.median(r)), 4) if r.size else None,
    }


def _clean(cat: dict) -> dict:
    """Finite mass/distance, positive flux, z > 0 (no comoving position otherwise)."""
    ok = (
        np.isfinite(cat["log_mhi"])
        & np.isfinite(cat["dist_mpc"])
        & (cat["flux"] > 0)
        & (cat["z"] > 0)
    )
    return {k: np.asarray(v)[ok] for k, v in cat.items()}


def _real_leg(n_null: int = 100):  # pragma: no cover - network + VizieR catalogues
    """FASHI DR2 x Tempel groups x Douglass voids: the environment-split HIMF.

    Headline weighting: DR2's own 1/(C * Vmax) on the C >= 0.5 sample. Since the DR2 referee
    round (2026-09-26): voids are full VoidFinder unions (all holes, table2), groups are in the
    heliocentric frame, the single-flux-cut weighting is ALSO run on the same C >= 0.5 sample with
    a paired delete-one-void jackknife on the difference, the EdS check reports how many
    galaxies change void status (by distance) instead of being read as fragility, and a
    random-void null (voids moved rigidly within the footprint) measures the estimator's and
    geometry's own bias with its sign and size.
    """
    voids = fetch_voidfinder_spheres(all_holes=True)
    spheres1 = fetch_voidfinder_spheres(all_holes=False)
    grp = fetch_tempel_groups()
    area = FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2  # DR2's own Vmax area

    raw = fetch_fashi_dr2()
    n_dr2 = int(raw["ra"].size)
    n_zle0 = int(np.sum(raw["z"] <= 0))
    cat = _clean(raw)
    env = _environments(cat, voids, grp)
    comp, vcat = cat["completeness"], cat["vmax_mpc3"]
    samp_b = np.isfinite(comp) & (comp >= FASHI_DR2_C_MIN) & np.isfinite(vcat) & (vcat > 0)
    vmax_b = vmax_from_catalogue(vcat, comp)
    vmax_a = vmax_1vmax(cat["log_mhi"], cat["dist_mpc"], cat["flux"])
    b = _environment_split(cat, samp_b, vmax_b, area, env)
    figdata = b.pop("_fig")
    a_same = _environment_split(cat, samp_b, vmax_a, area, env)
    a_same.pop("_fig")
    a_all = _environment_split(cat, np.ones(cat["ra"].size, bool), vmax_a, area, env)
    a_all.pop("_fig")

    # Old classification (maximal spheres only), for continuity with the DR1 paper.
    env_max = {
        **env,
        "in_void": void_membership_holes(
            env["xyz"], spheres1["sphere_xyz"], spheres1["sphere_radius"]
        ),
    }
    b_max = _environment_split(cat, samp_b, vmax_b, area, env_max)
    b_max.pop("_fig")

    # EdS: how many galaxies change void status, by redshift -- what the test actually probes.
    env_eds = _environments(cat, voids, grp, q0=0.5)
    b_eds = _environment_split(cat, samp_b, vmax_b, area, env_eds)
    b_eds.pop("_fig")
    cl = env["classifiable"] & samp_b
    flip = cl & (env["in_void"] != env_eds["in_void"])
    zb = [0.0, 0.02, 0.04, 0.06, 0.09]
    eds_flips = {
        f"z{lo:.2f}-{hi:.2f}": {
            "n_void": int((cl & env["in_void"] & (cat["z"] >= lo) & (cat["z"] < hi)).sum()),
            "n_changed": int((flip & (cat["z"] >= lo) & (cat["z"] < hi)).sum()),
        }
        for lo, hi in zip(zb[:-1], zb[1:], strict=True)
    }

    vm = void_members(env["xyz"], voids["sphere_xyz"], voids["sphere_radius"], voids["void_id"])
    jk = void_jackknife(
        cat, area, env["in_void"], env["classifiable"], vm,
        {"optA_same": (samp_b, vmax_a), "B": (samp_b, vmax_b)},
    )  # fmt: skip

    # Random-void null: voids moved rigidly within the Tempel (SDSS) footprint.
    rng = np.random.default_rng(64)
    null = []
    for _ in range(n_null):
        moved = random_void_positions(
            voids["sphere_xyz"], voids["void_id"], grp["gal_ra"], grp["gal_dec"], rng
        )
        iv = void_membership_holes(env["xyz"], moved, voids["sphere_radius"])
        _hv, fv = _himf_and_fit(
            cat["log_mhi"],
            cat["dist_mpc"],
            cat["flux"],
            area,
            mask=samp_b & env["classifiable"] & iv,
            vmax=vmax_b,
        )
        _hw, fw = _himf_and_fit(
            cat["log_mhi"],
            cat["dist_mpc"],
            cat["flux"],
            area,
            mask=samp_b & env["classifiable"] & ~iv,
            vmax=vmax_b,
        )
        d = fv.get("log_m_star", np.nan) - fw.get("log_m_star", np.nan)
        if np.isfinite(d):
            null.append(float(d))
    na = np.asarray(null)
    null_out = {
        "n_reps": n_null,
        "n_ok": int(na.size),
        "mean": round(float(na.mean()), 4) if na.size else None,
        "std": round(float(na.std(ddof=1)), 4) if na.size > 1 else None,
        # fraction of null offsets at least as negative as the measured one (one-sided p)
        "p_le_measured": round(float((np.sum(na <= b["void_knee_offset"]) + 1) / (na.size + 1)), 4)
        if na.size
        else None,
        "offsets": [round(x, 4) for x in null],
    }
    null_out["measured_minus_null_mean"] = (
        round(b["void_knee_offset"] - null_out["mean"], 3) if null_out["mean"] is not None else None
    )

    dr1 = _clean(fetch_fashi_dr1())
    dr1_env = _environments(dr1, voids, grp)
    d1 = _environment_split(
        dr1, np.ones(dr1["ra"].size, bool), None, 7600.0 * (np.pi / 180.0) ** 2, dr1_env
    )
    d1.pop("_fig")
    xm = _match_releases(dr1, cat)

    def keep(d: dict, *keys: str) -> dict:
        return {k: d[k] for k in keys if k in d}

    offs = (
        "n_sources",
        "void_knee_offset",
        "void_knee_offset_err",
        "void_knee_offset_sigma",
        "group_knee_offset",
        "group_knee_offset_err",
        "group_knee_offset_sigma",
        "himf_global",
        "n_in_void",
    )
    metrics = {
        "source": (
            "FASHI DR2 (arXiv:2606.31539, CSTCloud share) x Tempel+2017 groups x "
            "Douglass+2023 voids (all holes); weights 1/(C*Vmax) from the DR2 catalogue (C >= 0.5)"
        ),
        "is_real": True,
        "release": "DR2",
        "void_classification": "VoidFinder full voids (Douglass+2023 table2, all holes)",
        "n_dr2_catalogue": n_dr2,
        "n_dr2_z_nonpositive": n_zle0,
        "n_dropped_nonfinite_or_nonpositive_flux": int(n_dr2 - n_zle0 - cat["ra"].size),
        "n_finite_z_pos": int(cat["ra"].size),
        "weighting": "DR2 per-source completeness x Vmax (option B)",
        "c_min": FASHI_DR2_C_MIN,
        "dr2_paper_himf": {"log_m_star": 9.89, "log_m_star_err": 0.02, "alpha": -1.31,
                           "alpha_err": 0.02, "ref": "arXiv:2606.31539 abstract"},
        **b,
        "void_jackknife": jk,
        "random_void_null": null_out,
        "eds_void_knee_offset": b_eds["void_knee_offset"],
        "eds_void_knee_offset_err": b_eds["void_knee_offset_err"],
        "eds_void_knee_offset_sigma": b_eds["void_knee_offset_sigma"],
        "eds_void_status_changes": eds_flips,
        "maxsphere_void_knee_offset": b_max["void_knee_offset"],
        "maxsphere_void_knee_offset_err": b_max["void_knee_offset_err"],
        "maxsphere_n_in_void": b_max["n_in_void"],
        "optA_same_sample": keep(a_same, *offs),
        "optA_single_flux_cut": keep(a_all, *offs),
        "dr1_optA_single_flux_cut": keep(d1, *offs),
        "dr1_dr2_match": xm,
        # Same sample, only the weighting differs: the effect of the weighting alone.
        "weighting_effect_same_sample": round(b["void_knee_offset"] - a_same["void_knee_offset"], 3),
    }  # fmt: skip
    return metrics, figdata


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
            # the per-fit errors, so the paper stops hard-typing them (it quoted +/-0.08
            # against a committed 0.071, which then failed to reproduce the combined 0.087)
            ("VoidLogMStarErr", "himf_void", "log_m_star_err"),
            ("WallLogMStarErr", "himf_wall", "log_m_star_err"),
            ("NWall", "n_wall", None),
            ("NField", "n_field", None),
            ("GroupKneeOffset", "group_knee_offset", None),
            ("GroupKneeErr", "group_knee_offset_err", None),
            ("GroupKneeSigma", "group_knee_offset_sigma", None),
            ("NInVoid", "n_in_void", None),
            ("NGroupMembers", "n_group_members", None),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(d, key) if live else '--'}}}")
    if m.get("is_real") and m.get("release") == "DR2":
        # DR2-leg numbers the prose needs, all pipeline-made (no hand-typed comparisons).
        def nested(path_: str):
            cur: object = m
            for k in path_.split("."):
                cur = cur.get(k) if isinstance(cur, dict) else None
            return "--" if cur is None else cur

        for macro, key in (
            ("NDRTwo", "n_dr2_catalogue"),
            ("NZNonpos", "n_dr2_z_nonpositive"),
            ("CMin", "c_min"),
            ("WallPct", "wall_pct_of_sample"),
            ("WallMinusGlobal", "wall_minus_global_logmstar"),
            ("GlobalLogMStarErr", "himf_global.log_m_star_err"),
            ("GlobalAlphaErr", "himf_global.alpha_err"),
            ("EdsVoidKneeOffset", "eds_void_knee_offset"),
            ("EdsVoidKneeSigma", "eds_void_knee_offset_sigma"),
            ("OptAVoidKneeOffset", "optA_single_flux_cut.void_knee_offset"),
            ("OptAVoidKneeErr", "optA_single_flux_cut.void_knee_offset_err"),
            ("OptAGroupKneeOffset", "optA_single_flux_cut.group_knee_offset"),
            ("OptAGlobalAlpha", "optA_single_flux_cut.himf_global.alpha"),
            ("DROneN", "dr1_optA_single_flux_cut.n_sources"),
            ("DROneVoidKneeOffset", "dr1_optA_single_flux_cut.void_knee_offset"),
            ("DROneVoidKneeErr", "dr1_optA_single_flux_cut.void_knee_offset_err"),
            ("DROneVoidKneeSigma", "dr1_optA_single_flux_cut.void_knee_offset_sigma"),
            ("DROneGroupKneeOffset", "dr1_optA_single_flux_cut.group_knee_offset"),
            ("DROneGlobalAlpha", "dr1_optA_single_flux_cut.himf_global.alpha"),
            ("WeightingShiftPct", "weighting_shift_pct"),
            ("MatchedPct", "dr1_dr2_match.matched_pct"),
            ("NDropped", "n_dropped_nonfinite_or_nonpositive_flux"),
            ("JkErr", "void_jackknife.B_jackknife_err"),
            ("JkNOcc", "void_jackknife.n_voids_occupied"),
            ("WeightEffect", "weighting_effect_same_sample"),
            ("WeightEffectErr", "void_jackknife.B_minus_optA_same_jackknife_err"),
            ("OptASameVoidKneeOffset", "optA_same_sample.void_knee_offset"),
            ("OptASameGroupKneeOffset", "optA_same_sample.group_knee_offset"),
            ("MaxSphereVoidKneeOffset", "maxsphere_void_knee_offset"),
            ("MaxSphereNInVoid", "maxsphere_n_in_void"),
            ("NullReps", "random_void_null.n_ok"),
            ("NullMean", "random_void_null.mean"),
            ("NullStd", "random_void_null.std"),
            ("NullP", "random_void_null.p_le_measured"),
            ("NullExcess", "random_void_null.measured_minus_null_mean"),
            ("GlobalRedChi", "himf_global.red_chi2"),
            ("VoidRedChi", "himf_void.red_chi2"),
            ("WallRedChi", "himf_wall.red_chi2"),
            ("DistRatio", "dr1_dr2_match.median_dist_ratio_dr2_over_dr1"),
        ):
            lines.append(rf"\newcommand{{\feReal{macro}}}{{{nested(key)}}}")
    if m.get("is_real") and m.get("release") == "DR2":
        rn = m.get("random_void_null") or {}
        jk2 = m.get("void_jackknife") or {}
        ex, sd, je = rn.get("measured_minus_null_mean"), rn.get("std"), jk2.get("B_jackknife_err")
        we, wee = m.get("weighting_effect_same_sample"), jk2.get("B_minus_optA_same_jackknife_err")
        flips = m.get("eds_void_status_changes") or {}
        lo, hi = flips.get("z0.00-0.02") or {}, flips.get("z0.06-0.09") or {}

        def ratio(a_, b_, fmt="{:.1f}"):
            return fmt.format(abs(a_) / b_) if a_ is not None and b_ else "--"

        lines += [
            rf"\newcommand{{\feRealNullExcessSigma}}{{{ratio(ex, sd)}}}",
            rf"\newcommand{{\feRealNullExcessJkSigma}}{{{ratio(ex, je)}}}",
            rf"\newcommand{{\feRealWeightEffectSigma}}{{{ratio(we, wee)}}}",
            rf"\newcommand{{\feRealNullBiasPct}}{{{ratio(rn.get('mean'), abs(m.get('void_knee_offset') or 0) / 100 if m.get('void_knee_offset') else None, '{:.0f}')}}}",
            rf"\newcommand{{\feRealEdsChangedLowZPct}}{{{ratio(lo.get('n_changed'), lo.get('n_void', 0) / 100 if lo.get('n_void') else None, '{:.0f}')}}}",
            rf"\newcommand{{\feRealEdsChangedHighZPct}}{{{ratio(hi.get('n_changed'), hi.get('n_void', 0) / 100 if hi.get('n_void') else None, '{:.0f}')}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


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
