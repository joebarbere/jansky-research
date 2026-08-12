"""Independent replication of the Sofue & Kohno (2025) inner Milky Way rotation curve.

Sofue & Kohno (PASJ 77, 1335; arXiv:2509.23581) derive the modern inner rotation curve (RC) by
the terminal-velocity method (TVM): Gaussian decomposition of HI4PI / CO longitude–velocity
spectra (terminal velocity = the highest-|v| component), a velocity-dispersion correction
iterated so the curve meets the solar circular speed, Gaussian-weighted radial binning, and a
BH + Plummer bulge/disc + NFW halo decomposition whose local dark-matter density
(0.107 GeV/cm^3, halo-only, author-framed as a lower limit) sits ~3x below the ~0.3 consensus.

This module implements every stage independently (plans/86): the published ASCII tables from
the arXiv source are vendored as the offline anchor (``parse_paper_tables``), the estimators
(`gaussian_tvm` vs the `hi` slice's threshold rule) run head-to-head on the same spectra, and
``decompose_rc``/``rho_dm_local``/``sensitivity_scan`` reproduce and stress the dark-matter
number. Committed-real-results pattern: real legs write force-tracked ``results/innerrc_*.json``
and paper macros come only from those; synthetic spectra exist for tests alone.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "G_KPC",
    "R0_PC",
    "V0_KMS",
    "decompose_rc",
    "ew_asymmetry_fit",
    "gaussian_tvm",
    "paper_macros",
    "parse_paper_tables",
    "rc_from_terminal",
    "rho_dm_local_gev",
    "rotation_curve_weighted",
    "run_anchor",
    "calibrate_sigma",
    "sensitivity_scan",
    "synthetic_spectrum",
    "threshold_tvm",
]

# Galactic constants used by the paper (their eq. 6-9): keep identical for the anchor.
R0_PC = 8178.0
V0_KMS = 235.1
G_KPC = 4.30091e-6  # G in kpc (km/s)^2 / Msun
GEV_PER_MSUN_PC3 = 38.0  # 1 Msun/pc^3 = 38.0 GeV/cm^3 (mass-energy conversion, E=mc^2)


# ----------------------------------------------------------------------- anchor tables


def parse_paper_tables(table_dir: str) -> dict:
    """Parse the vendored arXiv-source RC tables into arrays.

    ``tab50[abc].tex`` hold the inner RC at dR = 50 pc; ``urcgatex[AB].tex`` the unified RC
    (TVM + VERA + VLBA + Gaia). Rows are ``R & V & dV \\\\`` with R in pc; the trailing solar
    row in the LaTeX wrapper is not part of the data files.
    """
    from pathlib import Path

    def read(names: list[str]) -> dict:
        rows = []
        for name in names:
            for line in (Path(table_dir) / name).read_text().splitlines():
                parts = [p.strip() for p in line.replace(r"\\", "").split("&")]
                if len(parts) == 3 and parts[0]:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        continue
        arr = np.array(rows)
        order = np.argsort(arr[:, 0])
        return {"R_pc": arr[order, 0], "V_kms": arr[order, 1], "dV_kms": arr[order, 2]}

    return {
        "inner": read(["tab50a.tex", "tab50b.tex", "tab50c.tex"]),
        "unified": read(["urcgatexA.tex", "urcgatexB.tex"]),
    }


# ------------------------------------------------------------------ terminal velocities


def synthetic_spectrum(
    vel_kms: np.ndarray,
    vterm_kms: float,
    *,
    sigma_v: float = 8.0,
    n_clouds: int = 6,
    noise_k: float = 0.15,
    peak_k: float = 4.0,
    seed: int = 0,
) -> np.ndarray:
    """A crowded synthetic HI/CO spectrum whose true terminal velocity is known.

    Emission is a set of Gaussian cloud components at velocities up to ``vterm_kms`` (the
    highest-velocity component is centred exactly there, with dispersion ``sigma_v``), plus
    noise — mimicking the paper's Fig. 2 decomposition targets.
    """
    rng = np.random.default_rng(seed)
    spec = np.zeros_like(vel_kms, dtype=float)
    centers = np.concatenate([rng.uniform(0.15, 0.85, size=n_clouds - 1) * vterm_kms, [vterm_kms]])
    for c in centers:
        amp = rng.uniform(0.4, 1.0) * peak_k
        width = rng.uniform(0.8, 1.4) * sigma_v
        spec += amp * np.exp(-0.5 * ((vel_kms - c) / width) ** 2)
    return spec + rng.normal(0.0, noise_k, size=vel_kms.size)


def gaussian_tvm(
    vel_kms: np.ndarray,
    spectrum: np.ndarray,
    *,
    min_amp_k: float = 1.0,
    max_components: int = 15,
    sign: int = +1,
    return_components: bool = False,
) -> float | tuple[float, list[tuple[float, float, float]]]:
    """Terminal velocity by Gaussian decomposition (the paper's Sec. 3.2 method).

    Deterministic greedy fit: repeatedly locate the strongest residual peak, fit a single
    Gaussian there (linearised moment fit over a local window), subtract, and stop at
    ``max_components`` or when the residual peak drops below ``min_amp_k``. The terminal
    velocity is the centre of the component with the largest ``sign * centre`` whose fitted
    amplitude exceeds ``min_amp_k``. No initial-parameter sensitivity — the failure mode of
    ad-hoc solver fits the paper's method is exposed to.
    """
    vel_kms = np.asarray(vel_kms, float)
    resid = np.asarray(spectrum, float).copy()
    dv = abs(float(np.median(np.diff(vel_kms))))
    centers = []
    for _ in range(max_components):
        i = int(np.argmax(resid))
        amp = resid[i]
        if amp < min_amp_k:
            break
        # moment fit in a window around the peak (robust to neighbours after subtraction)
        half = max(3, int(round(3 * 8.0 / dv)))
        lo, hi = max(0, i - half), min(vel_kms.size, i + half + 1)
        w = np.clip(resid[lo:hi], 0.0, None)
        if w.sum() <= 0:
            break
        c = float(np.sum(vel_kms[lo:hi] * w) / w.sum())
        var = float(np.sum((vel_kms[lo:hi] - c) ** 2 * w) / w.sum())
        width = max(math.sqrt(max(var, 1e-6)), dv)
        centers.append((c, amp, width))
        resid = resid - amp * np.exp(-0.5 * ((vel_kms - c) / width) ** 2)
    if not centers:
        return float("nan")
    # joint multi-Gaussian refinement from the greedy seeds: deterministic (seeds are
    # data-derived), and it un-blends neighbouring clouds the moment fit merges — the
    # classic failure mode of one-shot solver fits on crowded spectra.
    from scipy.optimize import curve_fit

    def multi(v, *p):
        out = np.zeros_like(v)
        for j in range(0, len(p), 3):
            out = out + p[j] * np.exp(-0.5 * ((v - p[j + 1]) / p[j + 2]) ** 2)
        return out

    p0, lo_b, hi_b = [], [], []
    for c, amp, width in centers:
        p0 += [amp, c, width]
        lo_b += [0.0, c - 3 * width, dv / 2]
        hi_b += [np.inf, c + 3 * width, 10 * width]
    try:
        popt, _ = curve_fit(
            multi,
            vel_kms,
            np.asarray(spectrum, float),
            p0=p0,
            bounds=(lo_b, hi_b),
            maxfev=20000,
        )
        refined = [
            (popt[j + 1], popt[j], popt[j + 2])
            for j in range(0, len(popt), 3)
            if popt[j] >= min_amp_k
        ]
        if refined:
            centers = refined
    except RuntimeError:
        pass  # keep greedy centers if refinement fails to converge
    best = max(centers, key=lambda t: sign * t[0])
    if return_components:
        return float(best[0]), [(float(c), float(a), float(w)) for c, a, w in centers]
    return float(best[0])


def threshold_tvm(vel_kms: np.ndarray, spectrum: np.ndarray, *, threshold_k: float = 2.0) -> float:
    """The `hi` slice's estimator: most positive velocity above a fixed brightness threshold.

    Kept verbatim for the head-to-head; documented to read high vs spectral fitting
    (McClure-Griffiths & Dickey 2016).
    """
    above = np.asarray(spectrum, float) > threshold_k
    if not above.any():
        return float("nan")
    return float(np.max(np.asarray(vel_kms, float)[above]))


def rc_from_terminal(
    longitudes_deg: np.ndarray,
    vterm_kms: np.ndarray,
    *,
    sigma_v_kms: float = 0.0,
    r0_pc: float = R0_PC,
    v0_kms: float = V0_KMS,
) -> tuple[np.ndarray, np.ndarray]:
    """Tangent-point mapping with the paper's dispersion correction (their eq. 10-12).

    ``R = R0 sin l``; ``V = (|v_term| - sigma_v) + V0 |sin l|`` — the correction subtracts the
    ISM velocity dispersion that inflates the envelope (15 km/s HI, 5 km/s CO in the paper).
    Works for either quadrant via ``abs``.
    """
    ell = np.deg2rad(np.asarray(longitudes_deg, float))
    s = np.abs(np.sin(ell))  # |sin l|: fourth-quadrant longitudes (l>180 or l<0) map to +R
    r = r0_pc * s
    v = (np.abs(np.asarray(vterm_kms, float)) - sigma_v_kms) + v0_kms * s
    return r, v


def rotation_curve_weighted(
    r_pc: np.ndarray,
    v_kms: np.ndarray,
    *,
    grid_pc: np.ndarray | None = None,
    dr_pc: float = 50.0,
    half_width_pc: float = 2.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The paper's Gaussian-weighted running average (their eq. 13-15).

    Returns ``(R_grid, V, dV)`` where each grid point averages all measurements with weights
    ``exp[-((R_i - R)/half_width)^2]`` and ``dV`` is the weighted standard deviation.
    """
    r_pc = np.asarray(r_pc, float)
    v_kms = np.asarray(v_kms, float)
    if grid_pc is None:
        grid_pc = np.arange(dr_pc, np.nanmax(r_pc) + dr_pc, dr_pc)
    out_v = np.full(grid_pc.size, np.nan)
    out_dv = np.full(grid_pc.size, np.nan)
    for k, rg in enumerate(grid_pc):
        w = np.exp(-(((r_pc - rg) / half_width_pc) ** 2))
        if w.sum() < 1e-12:
            continue
        m = float(np.sum(v_kms * w) / w.sum())
        out_v[k] = m
        out_dv[k] = math.sqrt(max(float(np.sum((v_kms - m) ** 2 * w) / w.sum()), 0.0))
    return grid_pc, out_v, out_dv


def calibrate_sigma(
    longitudes_deg: np.ndarray,
    vterm_kms: np.ndarray,
    *,
    solar_window_pc: tuple[float, float] = (7000.0, 8150.0),
    v0_kms: float = V0_KMS,
) -> float:
    """The paper's Sec.-3.3 dispersion calibration: choose sigma_v so the curve meets V0 at R0.

    V is linear in sigma_v, so the iteration has a closed form: sigma* is the mean excess of
    the uncorrected curve over V0 in the near-solar window (clipped at 0). Each estimator gets
    its own calibrated sigma, exactly as the paper calibrates to its own measurement.
    """
    r, v = rc_from_terminal(longitudes_deg, vterm_kms, sigma_v_kms=0.0, v0_kms=v0_kms)
    grid, vb, _ = rotation_curve_weighted(r, v, dr_pc=50.0, half_width_pc=25.0)
    sel = (grid >= solar_window_pc[0]) & (grid <= solar_window_pc[1]) & np.isfinite(vb)
    if not sel.any():
        return 0.0
    return float(max(np.mean(vb[sel]) - v0_kms, 0.0))


# ------------------------------------------------------------------------ decomposition


def _v_bh(r_pc: np.ndarray) -> np.ndarray:
    # their eq. 19: 4e6 Msun point mass -> 131.5 km/s at 1 pc
    return 131.5 / np.sqrt(np.maximum(r_pc, 1e-3))


def _v_plummer(r_pc: np.ndarray, v_i: float, a_pc: float) -> np.ndarray:
    # circular speed of a Plummer sphere, peaking near r ~ a*sqrt(2) (their eq. 20-21)
    x = np.maximum(r_pc, 1e-6) / a_pc
    return v_i * x / (1.0 + x**2) ** 0.75


def _v_nfw(r_pc: np.ndarray, v_h: float, h_pc: float) -> np.ndarray:
    x = np.maximum(r_pc, 1e-6) / h_pc
    g = np.log(1.0 + x) - x / (1.0 + x)
    return v_h * np.sqrt(np.maximum(g, 0.0) / x)


def _v_burkert(r_pc: np.ndarray, v_h: float, h_pc: float) -> np.ndarray:
    x = np.maximum(r_pc, 1e-6) / h_pc
    g = np.log(1.0 + x**2) / 2.0 + np.log(1.0 + x) - np.arctan(x)
    return v_h * np.sqrt(np.maximum(g, 0.0) / x)


HALO_MODELS = {"nfw": _v_nfw, "burkert": _v_burkert}


def rc_model(
    r_pc: np.ndarray,
    v_bulge: float,
    a_bulge: float,
    v_disc: float,
    a_disc: float,
    v_halo: float,
    h_halo: float,
    *,
    halo: str = "nfw",
) -> np.ndarray:
    """Total model rotation speed: BH + Plummer bulge + Plummer disc + halo in quadrature."""
    return np.sqrt(
        _v_bh(r_pc) ** 2
        + _v_plummer(r_pc, v_bulge, a_bulge) ** 2
        + _v_plummer(r_pc, v_disc, a_disc) ** 2
        + HALO_MODELS[halo](r_pc, v_halo, h_halo) ** 2
    )


# Fit bounds, hoisted out of decompose_rc so `bound_contact` can report against them.
# v_bulge's upper bound is the one that matters: the data themselves cap it, because
# components add in quadrature and the observed curve peaks at ~255 km/s, while the Plummer
# term peaks at 0.620 * v_bulge -> v_bulge <= ~410. The historic 800 is far above that, and
# variants that excise the inner curve (R > 2 kpc) leave the bulge unconstrained and run
# straight to it. Kept at 800 so the published fits are unchanged; railed variants are now
# flagged and excluded from quoted ranges instead of being silently counted as converged.
FIT_PARAM_NAMES = ("v_bulge", "a_bulge", "v_disc", "a_disc", "v_halo", "h_halo")
FIT_LOWER = (50.0, 50.0, 50.0, 1000.0, 10.0, 3000.0)
FIT_UPPER = (800.0, 2000.0, 600.0, 20000.0, 500.0, 100000.0)
BOUND_CONTACT_TOL = 0.01  # within 1% of the span of a bound counts as railed


def bound_contact(fit: dict, *, tol: float = BOUND_CONTACT_TOL) -> list[str]:
    """Names of fitted parameters sitting at (or within ``tol`` of) a bound.

    A `curve_fit` that does not raise is not the same as a converged fit: a parameter glued to
    a wall means the data did not determine it, and any quantity derived from it is reporting
    the bound rather than the measurement. Added 2026-08-12 after a referee found the quoted
    dark-matter-density maximum came from a variant with v_bulge at exactly 800.0.
    """
    railed = []
    for name, lo, hi in zip(FIT_PARAM_NAMES, FIT_LOWER, FIT_UPPER, strict=True):
        if name not in fit:
            continue
        v, span = float(fit[name]), hi - lo
        if v <= lo + tol * span or v >= hi - tol * span:
            railed.append(name)
    return railed


def decompose_rc(
    r_pc: np.ndarray,
    v_kms: np.ndarray,
    dv_kms: np.ndarray | None = None,
    *,
    halo: str = "nfw",
    p0: tuple[float, ...] = (400.0, 300.0, 320.0, 5500.0, 150.0, 20000.0),
) -> dict:
    """Least-squares decomposition into BH + bulge + disc + halo (the paper's Sec. 5.3)."""
    from scipy.optimize import curve_fit

    r_pc = np.asarray(r_pc, float)
    v_kms = np.asarray(v_kms, float)
    ok = np.isfinite(r_pc) & np.isfinite(v_kms)
    sigma = None
    if dv_kms is not None:
        dv = np.asarray(dv_kms, float)
        ok &= np.isfinite(dv) & (dv > 0)
        sigma = dv[ok]

    def f(r, vb, ab, vd, ad, vh, hh):
        return rc_model(r, vb, ab, vd, ad, vh, hh, halo=halo)

    popt, pcov = curve_fit(
        f,
        r_pc[ok],
        v_kms[ok],
        p0=p0,
        sigma=sigma,
        bounds=(list(FIT_LOWER), list(FIT_UPPER)),
        maxfev=20000,
    )
    perr = np.sqrt(np.diag(pcov))
    names = list(FIT_PARAM_NAMES)
    out: dict = {n: float(v) for n, v in zip(names, popt, strict=True)}
    out.update({f"d{n}": float(e) for n, e in zip(names, perr, strict=True)})
    out["halo"] = halo
    resid = v_kms[ok] - f(r_pc[ok], *popt)
    out["rms_kms"] = float(np.sqrt(np.mean(resid**2)))
    out["n_points"] = int(ok.sum())
    # chi2/N alongside rms: rms is not comparable across variants that fit different point
    # sets (the R > 2 kpc variants drop the 39 highest-residual inner rows), and the referee's
    # point about the outer curve is only visible in chi2.
    if sigma is not None:
        out["chi2_per_n"] = float(np.mean((resid / sigma) ** 2))
    out["rho_dm_gev"] = rho_dm_local_gev(out["v_halo"], out["h_halo"], halo=halo)
    # 1-sigma corners of rho_DM from the halo parameters' own uncertainties. This is the
    # honest carrier of any "compatible with X" statement — unlike the variant scan, it cannot
    # be manufactured by a bound.
    # All four corners, not an assumed pairing: rho is not monotonic in h (at R0 = 8.2 kpc it
    # *rises* with the NFW scale radius), so guessing which corner is extremal gets it wrong —
    # the (v+dv, h-dh) pairing gives 0.244 where the true 1-sigma maximum is 0.309.
    _corners = [
        rho_dm_local_gev(
            out["v_halo"] + a * out["dv_halo"], out["h_halo"] + b * out["dh_halo"], halo=halo
        )
        for a in (-1, 1)
        for b in (-1, 1)
    ]
    out["rho_dm_gev_lo"], out["rho_dm_gev_hi"] = float(min(_corners)), float(max(_corners))
    out["railed_params"] = bound_contact(out)
    return out


def rho_dm_local_gev(
    v_halo_kms: float, h_pc: float, *, halo: str = "nfw", r_pc: float = R0_PC
) -> float:
    """Local halo dark-matter density (GeV/cm^3) at ``r_pc`` from the fitted halo params.

    For NFW: rho0 = V_h^2 / (4 pi G h^2), rho(R) = rho0 / [x (1+x)^2] (the paper's eq. 23-26
    convention); Burkert analogously with rho(R) = rho0 / [(1+x)(1+x^2)].
    """
    h_kpc = h_pc / 1000.0
    rho0 = v_halo_kms**2 / (4.0 * math.pi * G_KPC * h_kpc**2)  # Msun / kpc^3
    x = r_pc / h_pc
    if halo == "nfw":
        rho = rho0 / (x * (1.0 + x) ** 2)
    else:
        rho = rho0 / ((1.0 + x) * (1.0 + x**2))
    return float(rho / 1e9 * GEV_PER_MSUN_PC3)  # Msun/kpc^3 -> Msun/pc^3 -> GeV/cm^3


def sensitivity_scan(r_pc: np.ndarray, v_kms: np.ndarray, dv_kms: np.ndarray | None = None) -> dict:
    """Map rho_DM's sensitivity to the decomposition choices (the paper-invited analysis).

    Varies: halo profile (NFW vs Burkert), inner fit boundary (all R vs R > 2 kpc, excising
    the bar), and error weighting (on/off). Returns every variant's fit + the rho_DM spread.
    """
    variants = {}
    for halo in HALO_MODELS:
        for rmin, tag in ((0.0, "all"), (2000.0, "R>2kpc")):
            for weighted in (True, False):
                sel = np.asarray(r_pc, float) >= rmin
                key = f"{halo}:{tag}:{'w' if weighted else 'u'}"
                try:
                    variants[key] = decompose_rc(
                        np.asarray(r_pc)[sel],
                        np.asarray(v_kms)[sel],
                        (np.asarray(dv_kms)[sel] if (weighted and dv_kms is not None) else None),
                        halo=halo,
                    )
                except Exception as e:  # noqa: BLE001 - a non-converging variant is a result
                    variants[key] = {"error": str(e)}
    # A variant only counts toward the quoted range if every fitted parameter is interior.
    # Before 2026-08-12 the range was taken over all eight, and its maximum came from a fit
    # with v_bulge at exactly the 800 km/s bound: excising R < 2 kpc removes all the data that
    # constrains the bulge, so (v_bulge, a_bulge) slide along the Plummer degeneracy ridge
    # until they hit walls. Railed variants are kept in the output — a reader should see them
    # — but they are reported separately and never set a quoted bound.
    ok_v = {k: v for k, v in variants.items() if "rho_dm_gev" in v and not v["railed_params"]}
    railed = {k: v["railed_params"] for k, v in variants.items() if v.get("railed_params")}
    rhos = [v["rho_dm_gev"] for v in ok_v.values()]
    all_rhos = [v["rho_dm_gev"] for v in variants.values() if "rho_dm_gev" in v]
    return {
        "variants": variants,
        "rho_dm_min_gev": float(min(rhos)),
        "rho_dm_max_gev": float(max(rhos)),
        "n_converged": len(rhos),
        "n_fitted": len(all_rhos),
        "railed_variants": railed,
        "rho_dm_min_gev_incl_railed": float(min(all_rhos)),
        "rho_dm_max_gev_incl_railed": float(max(all_rhos)),
    }


# -------------------------------------------------------------------------- E/W asymmetry


def ew_asymmetry_fit(r_pc: np.ndarray, dv_kms: np.ndarray) -> dict:
    """Fit the paper's damped sinusoid (their eq. 29) to an E-W velocity-difference curve.

    dV(R) = A exp(-R/L) sin(2 pi (R - R0s)/P); returns the fitted parameters and RMS.
    """
    from scipy.optimize import curve_fit

    def f(r, a, ell, r0s, p):
        return a * np.exp(-r / ell) * np.sin(2.0 * math.pi * (r - r0s) / p)

    ok = np.isfinite(np.asarray(r_pc, float)) & np.isfinite(np.asarray(dv_kms, float))
    popt, _ = curve_fit(
        f,
        np.asarray(r_pc, float)[ok],
        np.asarray(dv_kms, float)[ok],
        p0=(45.0, 3500.0, 4000.0, 4400.0),
        bounds=([1, 500, 0, 1000], [200, 20000, 8000, 20000]),
        maxfev=20000,
    )
    resid = np.asarray(dv_kms, float)[ok] - f(np.asarray(r_pc, float)[ok], *popt)
    return {
        "amp_kms": float(popt[0]),
        "damping_pc": float(popt[1]),
        "phase_pc": float(popt[2]),
        "period_pc": float(popt[3]),
        "rms_kms": float(np.sqrt(np.mean(resid**2))),
    }


# ------------------------------------------------------------------------------ anchor run


def run_anchor(out: str = ".", *, table_dir: str = "tests/data/sofue2025") -> dict:
    """The offline anchor leg: decompose the paper's own published RC tables.

    Reproducing their Table 1 fit (and the 0.107 GeV/cm^3 halo-only local DMD) from their own
    published curve licenses the variant analysis; the sensitivity scan then maps how far
    rho_DM moves under defensible alternative choices. Writes results/innerrc_anchor.json.
    """
    import json
    from pathlib import Path

    tables = parse_paper_tables(table_dir)
    uni = tables["unified"]
    fit = decompose_rc(uni["R_pc"], uni["V_kms"], uni["dV_kms"])
    scan = sensitivity_scan(uni["R_pc"], uni["V_kms"], uni["dV_kms"])
    # Their Table 1 solution evaluated against their own table: if its rms is comparable to
    # our best fit's, the two rho_DM values are degenerate alternatives, not a discrepancy.
    # CONVENTION: the paper's V_h (their eq. 24-25) absorbs 4pi into g(x), so their
    # rho0 = V_h^2/(G h^2) — a factor 4pi off ours, and their V_h = ours/sqrt(4pi). Convert
    # their Table-1 V_h into our convention for the curve, and use their rho0 for their rho.
    paper = dict(
        v_bulge=406.2,
        a_bulge=332.8,
        v_disc=322.4,
        a_disc=5624.8,
        v_halo=64.4 * math.sqrt(4 * math.pi),
        h_halo=22379.1,
    )
    v_paper = rc_model(uni["R_pc"], *paper.values())
    paper_rms = float(np.sqrt(np.mean((uni["V_kms"] - v_paper) ** 2)))
    # chi2/N for their solution on the same points, so it is comparable with the refit's, and
    # the outer-curve residual that the rms comparison hides: beyond 8 kpc their published halo
    # sits systematically below their own unified curve, which is exactly where rho_DM(R0) is
    # set. Reporting rms alone made a structured bias look like a wash.
    _res = (uni["V_kms"] - v_paper) / uni["dV_kms"]
    paper_chi2 = float(np.mean(_res**2))
    _outer = uni["R_pc"] > 8000.0
    paper_outer = {
        "n": int(_outer.sum()),
        "mean_resid_sigma": float(np.mean(_res[_outer])),
        "rms_kms": float(np.sqrt(np.mean((uni["V_kms"][_outer] - v_paper[_outer]) ** 2))),
    }
    paper_rho = rho_dm_local_gev(paper["v_halo"], paper["h_halo"])  # == their 0.107 if consistent
    metrics = {
        "source": "Sofue & Kohno 2025 published RC tables (arXiv:2509.23581 source, vendored)",
        "paper_table1_params": paper,
        "paper_table1_chi2_per_n": paper_chi2,
        "paper_table1_outer": paper_outer,
        "n_inner_rows": int(tables["inner"]["R_pc"].size),
        "n_unified_rows": int(uni["R_pc"].size),
        "anchor_fit": fit,
        "paper_table1_rms_kms": paper_rms,
        "paper_table1_rho_dm_gev": paper_rho,
        "paper_rho_dm_gev": 0.107,
        "rho_dm_ratio_vs_paper": fit["rho_dm_gev"] / 0.107,
        "r_max_pc": float(uni["R_pc"].max()),
        "sensitivity": {k: v for k, v in scan.items() if k != "variants"},
        # The whole fit, not a curated subset. This used to drop v_halo and h_halo — the only
        # two numbers rho_dm_gev is computed from — so a reader could not check the quoted
        # maximum, nor see that the variant producing it had a parameter on a bound.
        "sensitivity_variants": {k: v for k, v in scan["variants"].items()},
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "innerrc_anchor.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


# ------------------------------------------------------------------------- HI4PI real leg

HI4PI_BASE = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/594/A116/CUBES/GAL/CAR/"
# The galactic-plane (b=0-centred) CAR tiles covering the inner Galaxy, from cubes_gal.dat.
HI4PI_PLANE_TILES = {
    "CAR_E01.fits": 10.0,
    "CAR_E02.fits": 30.0,
    "CAR_E03.fits": 50.0,
    "CAR_E04.fits": 70.0,
    "CAR_E05.fits": 90.0,
    "CAR_E14.fits": 270.0,
    "CAR_E15.fits": 290.0,
    "CAR_E16.fits": 310.0,
    "CAR_E17.fits": 330.0,
    "CAR_E18.fits": 350.0,
}


def lv_from_cube(
    cube: np.ndarray,
    lon0_deg: float,
    dlon_deg: float,
    lat0_deg: float,
    dlat_deg: float,
    vel0_kms: float,
    dvel_kms: float,
    *,
    b_max_deg: float = 3.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Collapse a (v, b, l) HI cube into a |b|-averaged longitude-velocity diagram.

    Axis descriptions are the FITS reference values at pixel 0 (``lon0 + i*dlon`` etc.).
    Returns ``(glon_deg, vel_kms, T[v, l])``. Pure array function — offline-testable.
    """
    nv, nb, nl = cube.shape
    lat = lat0_deg + dlat_deg * np.arange(nb)
    sel = np.abs(lat) <= b_max_deg
    t_lv = np.nanmean(cube[:, sel, :], axis=1)
    glon = lon0_deg + dlon_deg * np.arange(nl)
    vel = vel0_kms + dvel_kms * np.arange(nv)
    return glon % 360.0, vel, t_lv


def tvm_spectrum(
    vel_kms: np.ndarray,
    spectrum: np.ndarray,
    *,
    sign: int,
    method: str = "gaussian",
    edge_threshold_k: float = 1.0,
    window_in_kms: float = 80.0,
    window_out_kms: float = 30.0,
) -> float:
    """Terminal velocity of one real HI spectrum, by either estimator, on the envelope window.

    Real plane spectra are ~100 K forests; fitting the whole profile buries the faint terminal
    component under bright inner-Galaxy emission. Both estimators therefore operate on a
    window around the emission edge: the largest ``sign * v`` where ``T > edge_threshold_k``,
    extended ``window_in_kms`` inward and ``window_out_kms`` outward. Within it, ``gaussian``
    decomposes and returns the outermost component centre (the paper's method); ``threshold``
    returns the 2 K envelope crossing (the `hi` slice's rule).
    """
    vel_kms = np.asarray(vel_kms, float)
    spectrum = np.asarray(spectrum, float)
    v_signed = sign * vel_kms
    above = spectrum > edge_threshold_k
    if not above.any():
        return float("nan")
    edge = float(np.max(v_signed[above]))
    win = (v_signed >= edge - window_in_kms) & (v_signed <= edge + window_out_kms)
    if win.sum() < 8:
        return float("nan")
    if method == "threshold":
        t = threshold_tvm(v_signed[win], spectrum[win], threshold_k=2.0)
        return sign * t
    g = gaussian_tvm(v_signed[win], spectrum[win], min_amp_k=edge_threshold_k, max_components=8)
    assert isinstance(g, float)
    return sign * g


def fetch_hi4pi_tile(name: str, dest_dir: str) -> str:  # pragma: no cover - network
    """Download one HI4PI CAR tile from CDS (resumable); returns the local path."""
    import time
    from pathlib import Path
    from urllib.request import Request, urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / name
    if path.exists():
        return str(path)
    part = path.with_suffix(".part")
    # Stall-proof: 60 s socket timeout per read, resume from the .part offset, retry with
    # backoff (a CDS transfer hung mid-stream for hours on 2026-08-04 without this).
    for attempt in range(8):
        offset = part.stat().st_size if part.exists() else 0
        req = Request(HI4PI_BASE + name)
        if offset:
            req.add_header("Range", f"bytes={offset}-")
        try:
            with urlopen(req, timeout=60) as r, open(part, "ab" if offset else "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            part.rename(path)
            return str(path)
        except (TimeoutError, OSError):
            time.sleep(min(2**attempt, 60))
    raise OSError(f"failed to fetch {name} after 8 attempts")


def run_hi4pi(  # pragma: no cover - network + real data (core functions tested offline)
    out: str = ".",
    *,
    dest_dir: str = "data/hi4pi",
    b_max_deg: float = 3.0,
    lon_step: int = 2,
    sigma_v_kms: float = 15.0,
    delete_after: bool = False,
    table_dir: str = "tests/data/sofue2025",
) -> dict:
    """The real-data leg: HI4PI tiles -> LV diagram -> both TVM estimators -> RC vs their Table 2.

    Writes ``results/innerrc_hi4pi.json`` (committed evidence). ``lon_step`` subsamples the
    0.0833-degree longitude grid (2 -> ~10 arcmin sampling, ~1200 sightlines).
    """
    import json
    import os
    from pathlib import Path

    from astropy.io import fits

    rows = []  # (glon, vterm_gauss, vterm_thresh)
    fig_lv = None  # (glon, vel, T) of the l~30 tile, kept for the LV-diagram figure
    for name in HI4PI_PLANE_TILES:
        path = fetch_hi4pi_tile(name, dest_dir)
        with fits.open(path, memmap=True) as hdul:
            h = hdul[0].header
            cube = np.asarray(hdul[0].data, float)
            # FITS pixel 0 world values per axis (CAR projection: linear)
            lon0 = h["CRVAL1"] + (1 - h["CRPIX1"]) * h["CDELT1"]
            lat0 = h["CRVAL2"] + (1 - h["CRPIX2"]) * h["CDELT2"]
            vel0 = (h["CRVAL3"] + (1 - h["CRPIX3"]) * h["CDELT3"]) / 1000.0
            dvel = h["CDELT3"] / 1000.0
            glon, vel, t_lv = lv_from_cube(
                cube, lon0, h["CDELT1"], lat0, h["CDELT2"], vel0, dvel, b_max_deg=b_max_deg
            )
            if name == "CAR_E02.fits":
                fig_lv = (glon.copy(), vel.copy(), t_lv.copy())
        for i in range(0, glon.size, lon_step):
            ell = glon[i]
            in_q1 = 5.0 <= ell <= 89.5
            in_q4 = 270.5 <= ell <= 355.0
            if not (in_q1 or in_q4):
                continue
            sign = +1 if in_q1 else -1
            spec = t_lv[:, i]
            if not np.isfinite(spec).any():
                continue
            vg = tvm_spectrum(vel, spec, sign=sign, method="gaussian")
            vt = tvm_spectrum(vel, spec, sign=sign, method="threshold")
            rows.append((float(ell), float(vg), float(vt)))
        if delete_after:
            os.remove(path)
        print(f"[innerrc] {name}: {len(rows)} sightlines cumulative", flush=True)

    arr = np.array(rows)
    ells, vgs, vts = arr[:, 0], arr[:, 1], arr[:, 2]
    ok = np.isfinite(vgs) & np.isfinite(vts)
    bias = np.abs(vts[ok]) - np.abs(vgs[ok])  # threshold minus gaussian, in |v| (envelope sense)
    bias_kms = {
        "median": float(np.median(bias)),
        "mean": float(np.mean(bias)),
        "p16": float(np.percentile(bias, 16)),
        "p84": float(np.percentile(bias, 84)),
    }
    # Their Sec.-3.3 calibration, per estimator (each meets V0 at R0 by construction; the
    # difference of the calibrated sigmas IS the estimator bias in the paper's own terms).
    sigma_gauss = calibrate_sigma(ells[ok], vgs[ok])
    sigma_thresh = calibrate_sigma(ells[ok], vts[ok])
    east = ok & ((ells > 0) & (ells < 180))
    west = ok & (ells > 180)
    r_e, v_e = rc_from_terminal(ells[east], vgs[east], sigma_v_kms=sigma_gauss)
    r_w, v_w = rc_from_terminal(ells[west], vgs[west], sigma_v_kms=sigma_gauss)
    r_t, v_t = rc_from_terminal(ells[ok], vts[ok], sigma_v_kms=sigma_thresh)
    grid = np.arange(50.0, 8200.0, 50.0)
    _, ve_b, dve_b = rotation_curve_weighted(r_e, v_e, grid_pc=grid, half_width_pc=25.0)
    _, vw_b, dvw_b = rotation_curve_weighted(r_w, v_w, grid_pc=grid, half_width_pc=25.0)
    both_r = np.concatenate([r_e, r_w])
    both_v = np.concatenate([v_e, v_w])
    _, v_b, dv_b = rotation_curve_weighted(both_r, both_v, grid_pc=grid, half_width_pc=25.0)

    paper = parse_paper_tables(table_dir)["inner"]
    v_paper_i = np.interp(grid, paper["R_pc"], paper["V_kms"])
    cmp_ok = np.isfinite(v_b) & (grid > 2000)  # outside the bar, where TVM is a mass tracer
    table2 = {
        "mean_abs_dv_kms": float(np.mean(np.abs(v_b[cmp_ok] - v_paper_i[cmp_ok]))),
        "median_dv_kms": float(np.median(v_b[cmp_ok] - v_paper_i[cmp_ok])),
        "n_bins": int(cmp_ok.sum()),
    }
    # the same comparison WITHOUT calibration (fixed sigma=15 as first attempted), committed so
    # the findings' pre-calibration number is traceable to pipeline output
    r_u, v_u = rc_from_terminal(ells[ok], vgs[ok], sigma_v_kms=15.0)
    _, vu_b, _ = rotation_curve_weighted(r_u, v_u, grid_pc=grid, half_width_pc=25.0)
    u_ok = np.isfinite(vu_b) & (grid > 2000)
    table2["median_dv_kms_fixed_sigma15"] = float(np.median(vu_b[u_ok] - v_paper_i[u_ok]))
    dv_ew = ve_b - vw_b
    # fit the mid-disc window (their Sec. 5.5 region): inside ~2 kpc the bar's non-circular
    # chaos rails any smooth fit — the first attempt hit its amplitude bound doing exactly that
    ew_ok = np.isfinite(dv_ew) & (grid > 2000) & (grid < 8000)
    ew = ew_asymmetry_fit(grid[ew_ok], dv_ew[ew_ok]) if ew_ok.sum() > 20 else {}
    ew["mid_disc_rms_kms"] = float(np.sqrt(np.mean(dv_ew[ew_ok] ** 2))) if ew_ok.sum() > 0 else None

    _, vt_b, _ = rotation_curve_weighted(r_t, v_t, grid_pc=grid, half_width_pc=25.0)
    metrics = {
        "source": f"HI4PI (CDS J/A+A/594/A116, CAR plane tiles), |b|<{b_max_deg:.1f} deg",
        "n_sightlines": int(ok.sum()),
        "estimator_bias_threshold_minus_gaussian_kms": bias_kms,
        "sigma_v_calibrated_gaussian_kms": round(sigma_gauss, 2),
        "sigma_v_calibrated_threshold_kms": round(sigma_thresh, 2),
        "table2_comparison_R>2kpc": table2,
        "ew_asymmetry_fit": ew,
        "rc_grid_pc": grid.tolist(),
        "rc_v_kms": [round(float(x), 2) if np.isfinite(x) else None for x in v_b],
        "rc_dv_kms": [round(float(x), 2) if np.isfinite(x) else None for x in dv_b],
        "rc_east_v_kms": [round(float(x), 2) if np.isfinite(x) else None for x in ve_b],
        "rc_west_v_kms": [round(float(x), 2) if np.isfinite(x) else None for x in vw_b],
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "innerrc_hi4pi.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figures_hi4pi(
        op / "papers" / "innerrc" / "figures",
        fig_lv=fig_lv,
        grid=grid,
        v_b=v_b,
        dv_b=dv_b,
        vt_b=vt_b,
        ve_b=ve_b,
        vw_b=vw_b,
        paper_r=paper["R_pc"],
        paper_v=paper["V_kms"],
        ew=ew,
        table_dir=table_dir,
    )
    return metrics


def _figures_hi4pi(  # pragma: no cover - exercised by the real leg
    out_dir,
    *,
    fig_lv,
    grid,
    v_b,
    dv_b,
    vt_b,
    ve_b,
    vw_b,
    paper_r,
    paper_v,
    ew,
    table_dir,
) -> None:
    """The source-paper-style figure set: LV diagram, decomposition, RC, E/W, anchor."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Fig A (their Fig. 1): the |b|-averaged LV diagram with the fitted terminal envelope
    if fig_lv is not None:
        glon, vel, t_lv = fig_lv
        fig, ax = plt.subplots(figsize=(7, 3.4))
        vmax = np.nanpercentile(t_lv, 99)
        ax.imshow(
            t_lv,
            origin="lower",
            aspect="auto",
            cmap="magma",
            vmin=0,
            vmax=vmax,
            extent=[glon.min(), glon.max(), vel.min(), vel.max()],
        )
        env_l, env_v = [], []
        for i in range(0, glon.size, 4):
            vt = tvm_spectrum(vel, t_lv[:, i], sign=+1, method="gaussian")
            if np.isfinite(vt):
                env_l.append(glon[i])
                env_v.append(vt)
        ax.plot(env_l, env_v, ".", ms=2.5, color="cyan", label="Gaussian-TVM terminal velocity")
        ax.set(
            xlabel="Galactic longitude (deg)",
            ylabel="$V_{LSR}$ (km/s)",
            ylim=(-100, 200),
            title="HI4PI longitude–velocity diagram, $|b|<3^\\circ$",
        )
        ax.legend(loc="upper right", fontsize=7)
        fig.tight_layout()
        fig.savefig(out / "innerrc_lv.pdf")
        plt.close(fig)

        # Fig B (their Fig. 2): one spectrum decomposed into Gaussian components
        i30 = int(np.argmin(np.abs(glon - 30.0)))
        spec = t_lv[:, i30]
        res = tvm_spectrum_components(vel, spec, sign=+1)
        if res is not None:
            vterm, comps, win = res
            fig, ax = plt.subplots(figsize=(6, 3.2))
            ax.step(vel[win], spec[win], where="mid", lw=0.9, color="k", label="HI4PI spectrum")
            vv = vel[win]
            for c, a, w in comps:
                ax.plot(vv, a * np.exp(-0.5 * ((vv - c) / max(w, 1e-3)) ** 2), lw=0.8, color="m")
            ax.axvline(vterm, color="c", ls="--", lw=1.2, label=f"$v_t$ = {vterm:.1f} km/s")
            ax.set(
                xlabel="$V_{LSR}$ (km/s)",
                ylabel="$T_B$ (K)",
                title=f"Envelope decomposition at $\\ell$ = {glon[i30]:.1f}$^\\circ$",
            )
            ax.legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(out / "innerrc_decomp.pdf")
            plt.close(fig)

    # Fig C (their Fig. 8): our binned RC vs their published Table 2, plus the threshold curve
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    okb = np.isfinite(v_b)
    ax.fill_between(grid[okb], (v_b - dv_b)[okb], (v_b + dv_b)[okb], color="C0", alpha=0.25, lw=0)
    ax.plot(grid[okb], v_b[okb], color="C0", lw=1.2, label="this work (Gaussian TVM)")
    okt = np.isfinite(vt_b)
    ax.plot(grid[okt], vt_b[okt], color="C2", lw=0.9, ls=":", label="threshold estimator")
    ax.plot(paper_r, paper_v, color="C3", lw=1.0, ls="--", label="Sofue & Kohno Table 2")
    ax.plot([R0_PC], [V0_KMS], "o", color="k", ms=5, label="Sun")
    ax.set(xlabel="R (pc)", ylabel="V (km/s)", xlim=(0, 8400), ylim=(80, 320))
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "innerrc_rc.pdf")
    plt.close(fig)

    # Fig D (their Fig. 14): East vs West curves and the damped-sinusoid fit
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 4.6), sharex=True)
    oke, okw = np.isfinite(ve_b), np.isfinite(vw_b)
    ax1.plot(grid[oke], ve_b[oke], color="C3", lw=0.9, label="East ($\\ell>0$)")
    ax1.plot(grid[okw], vw_b[okw], color="C0", lw=0.9, label="West ($\\ell<0$)")
    ax1.set(ylabel="V (km/s)")
    ax1.legend(fontsize=7)
    both = oke & okw
    dv_ew_line = ve_b - vw_b
    ax2.axhline(0, color="0.6", lw=0.6)
    ax2.plot(grid[both], dv_ew_line[both], color="k", lw=0.8, label="E $-$ W")
    if ew:
        model = (
            ew["amp_kms"]
            * np.exp(-grid / ew["damping_pc"])
            * np.sin(2 * np.pi * (grid - ew["phase_pc"]) / ew["period_pc"])
        )
        ax2.plot(grid, model, color="C1", lw=1.0, ls="--", label="eq.-29 damped sinusoid")
    ax2.set(xlabel="R (pc)", ylabel="$\\delta V$ (km/s)")
    ax2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "innerrc_ew.pdf")
    plt.close(fig)

    # Fig E (their Fig. 11): anchor decomposition of their unified RC, ours vs theirs
    tables = parse_paper_tables(table_dir)
    uni = tables["unified"]
    fit = decompose_rc(uni["R_pc"], uni["V_kms"], uni["dV_kms"])
    rr = np.linspace(20, 26000, 800)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.errorbar(
        uni["R_pc"],
        uni["V_kms"],
        yerr=uni["dV_kms"],
        fmt=".",
        ms=3,
        color="0.4",
        lw=0.5,
        label="their published unified RC",
    )
    ours = [fit[k] for k in ("v_bulge", "a_bulge", "v_disc", "a_disc", "v_halo", "h_halo")]
    ax.plot(
        rr,
        rc_model(rr, *ours),
        color="C0",
        lw=1.3,
        label=f"our refit ($\\rho_{{DM}}$={fit['rho_dm_gev']:.2f} GeV/cm$^3$)",
    )
    theirs = [406.2, 332.8, 322.4, 5624.8, 64.4 * math.sqrt(4 * math.pi), 22379.1]
    ax.plot(
        rr,
        rc_model(rr, *theirs),
        color="C3",
        lw=1.1,
        ls="--",
        label="their Table 1 ($\\rho_{DM}$=0.107)",
    )
    ax.plot(rr, _v_plummer(rr, ours[0], ours[1]), lw=0.7, ls=":", color="C0")
    ax.plot(rr, _v_plummer(rr, ours[2], ours[3]), lw=0.7, ls="-.", color="C0")
    ax.plot(rr, _v_nfw(rr, ours[4], ours[5]), lw=0.7, ls="--", color="C0")
    ax.set(xlabel="R (pc)", ylabel="V (km/s)", xscale="log", xlim=(30, 26000), ylim=(0, 320))
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "innerrc_decomposition.pdf")
    plt.close(fig)


def tvm_spectrum_components(vel_kms, spectrum, *, sign):  # pragma: no cover - figure helper
    """Like tvm_spectrum(gaussian) but also returns the fitted components and window mask."""
    vel_kms = np.asarray(vel_kms, float)
    spectrum = np.asarray(spectrum, float)
    v_signed = sign * vel_kms
    above = spectrum > 1.0
    if not above.any():
        return None
    edge = float(np.max(v_signed[above]))
    win = (v_signed >= edge - 80.0) & (v_signed <= edge + 30.0)
    if win.sum() < 8:
        return None
    res = gaussian_tvm(
        v_signed[win], spectrum[win], min_amp_k=1.0, max_components=8, return_components=True
    )
    vterm, comps = res
    comps = [(sign * c, a, w) for c, a, w in comps]
    return sign * vterm, comps, win


def paper_macros(
    out: str = ".", *, results_dir: str = "results", table_dir: str = "tests/data/sofue2025"
) -> str:
    """Generate papers/innerrc/generated/macros.tex from the COMMITTED evidence JSONs only.

    The committed-real-results rule: no run happens here — if the evidence files are absent
    this fails loudly rather than substituting anything.
    """
    import json
    from pathlib import Path

    a = json.loads((Path(results_dir) / "innerrc_anchor.json").read_text())
    h = json.loads((Path(results_dir) / "innerrc_hi4pi.json").read_text())
    fit, sens = a["anchor_fit"], a["sensitivity"]
    bias, ew = h["estimator_bias_threshold_minus_gaussian_kms"], h["ew_asymmetry_fit"]
    # Inside 2 kpc the committed HI-only curve and their (CO-dominated) inner table disagree,
    # and the paper claimed the opposite ("the bar-region inner peak reproduce fully"). Both
    # sides are already committed evidence, so the disagreement is measured here rather than
    # described: their inner table is vendored, our curve is in the hi4pi results JSON.
    _inner = parse_paper_tables(table_dir)["inner"]
    _g = np.asarray(h["rc_grid_pc"], float)
    _v = np.asarray(h["rc_v_kms"], float)
    _sel = np.isfinite(_v) & (_g < 2000.0)
    _theirs = np.interp(_g[_sel], _inner["R_pc"], _inner["V_kms"])
    _dv_inner = float(np.mean(_v[_sel] - _theirs))
    _ipk = int(np.nanargmax(_inner["V_kms"]))
    _opk = int(np.argmax(_v[_sel]))
    lines = [
        "% Auto-generated by jansky_research.innerrc.paper_macros from the committed",
        "% results/innerrc_*.json evidence — do not edit by hand.",
        rf"\newcommand{{\irRhoOurs}}{{{fit['rho_dm_gev']:.2f}}}",
        rf"\newcommand{{\irRhoTheirs}}{{{a['paper_table1_rho_dm_gev']:.3f}}}",
        rf"\newcommand{{\irRmsOurs}}{{{fit['rms_kms']:.1f}}}",
        rf"\newcommand{{\irRmsTheirs}}{{{a['paper_table1_rms_kms']:.1f}}}",
        rf"\newcommand{{\irVBulgeOurs}}{{{fit['v_bulge']:.0f}}}",
        rf"\newcommand{{\irABulgeOurs}}{{{fit['a_bulge']:.0f}}}",
        rf"\newcommand{{\irRhoMin}}{{{sens['rho_dm_min_gev']:.2f}}}",
        rf"\newcommand{{\irRhoMax}}{{{sens['rho_dm_max_gev']:.2f}}}",
        rf"\newcommand{{\irScanN}}{{{sens['n_converged']}}}",
        rf"\newcommand{{\irScanNFitted}}{{{sens['n_fitted']}}}",
        rf"\newcommand{{\irScanNRailed}}{{{len(sens['railed_variants'])}}}",
        # The honest carrier of "compatible with the consensus density": the primary fit's own
        # halo-parameter uncertainties, which no bound can manufacture. The variant scan cannot
        # play that role -- six of its eight variants have a parameter on a wall.
        rf"\newcommand{{\irRhoLo}}{{{fit['rho_dm_gev_lo']:.2f}}}",
        rf"\newcommand{{\irRhoHi}}{{{fit['rho_dm_gev_hi']:.2f}}}",
        rf"\newcommand{{\irChiOurs}}{{{fit['chi2_per_n']:.2f}}}",
        rf"\newcommand{{\irChiTheirs}}{{{a['paper_table1_chi2_per_n']:.2f}}}",
        rf"\newcommand{{\irOuterN}}{{{a['paper_table1_outer']['n']}}}",
        rf"\newcommand{{\irOuterBiasSig}}{{{a['paper_table1_outer']['mean_resid_sigma']:.2f}}}",
        # How far the fitted E/W period sits from the one their eq. 29 quotes. The paper called
        # this "replicates in period and phase"; it is a 36% disagreement.
        rf"\newcommand{{\irInnerDv}}{{{_dv_inner:.0f}}}",
        rf"\newcommand{{\irInnerN}}{{{int(_sel.sum())}}}",
        rf"\newcommand{{\irPeakTheirs}}{{{_inner['V_kms'][_ipk]:.0f}}}",
        rf"\newcommand{{\irPeakTheirsR}}{{{_inner['R_pc'][_ipk]:.0f}}}",
        rf"\newcommand{{\irPeakOurs}}{{{_v[_sel][_opk]:.0f}}}",
        rf"\newcommand{{\irPeakOursR}}{{{_g[_sel][_opk]:.0f}}}",
        # Both already committed; neither had a macro, so the paper could not cite the
        # uncalibrated offset that exposes how much of the Table-2 agreement is pinned, nor
        # the fit residual that shows the sinusoid explains under half the E/W variance.
        rf"\newcommand{{\irTabTwoFixedSigma}}{{{h['table2_comparison_R>2kpc']['median_dv_kms_fixed_sigma15']:.1f}}}",
        rf"\newcommand{{\irEwFitRms}}{{{ew['rms_kms']:.1f}}}",
        rf"\newcommand{{\irEwPeriodOffPct}}{{{100 * abs(ew['period_pc'] - 4400.0) / 4400.0:.0f}}}",
        rf"\newcommand{{\irUniRows}}{{{a['n_unified_rows']}}}",
        rf"\newcommand{{\irInnerRows}}{{{a['n_inner_rows']}}}",
        rf"\newcommand{{\irRmax}}{{{a['r_max_pc'] / 1000:.1f}}}",
        rf"\newcommand{{\irNsight}}{{{h['n_sightlines']}}}",
        rf"\newcommand{{\irBiasMed}}{{{bias['median']:.1f}}}",
        rf"\newcommand{{\irBiasPlo}}{{{bias['p16']:.0f}}}",
        rf"\newcommand{{\irBiasPhi}}{{{bias['p84']:.0f}}}",
        rf"\newcommand{{\irSigGauss}}{{{h['sigma_v_calibrated_gaussian_kms']:.1f}}}",
        rf"\newcommand{{\irSigThresh}}{{{h['sigma_v_calibrated_threshold_kms']:.1f}}}",
        rf"\newcommand{{\irTabTwoAbs}}{{{h['table2_comparison_R>2kpc']['mean_abs_dv_kms']:.1f}}}",
        rf"\newcommand{{\irTabTwoMed}}{{{h['table2_comparison_R>2kpc']['median_dv_kms']:.1f}}}",
        rf"\newcommand{{\irTabTwoN}}{{{h['table2_comparison_R>2kpc']['n_bins']}}}",
        rf"\newcommand{{\irEwAmp}}{{{ew['amp_kms']:.1f}}}",
        rf"\newcommand{{\irEwPeriod}}{{{ew['period_pc'] / 1000:.1f}}}",
        rf"\newcommand{{\irEwPhase}}{{{ew['phase_pc'] / 1000:.1f}}}",
        rf"\newcommand{{\irEwRms}}{{{ew['mid_disc_rms_kms']:.1f}}}",
    ]
    path = Path(out) / "papers" / "innerrc" / "generated" / "macros.tex"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Sofue & Kohno inner-RC replication (plan 86).")
    p.add_argument("--out", default=".")
    p.add_argument("--anchor", action="store_true", help="run the offline anchor leg")
    p.add_argument("--hi4pi", action="store_true", help="run the HI4PI real-data leg (network)")
    p.add_argument("--paper", action="store_true", help="macros from the committed evidence")
    p.add_argument("--lon-step", type=int, default=2)
    p.add_argument("--delete-after", action="store_true")
    args = p.parse_args(argv)
    if args.anchor:
        m = run_anchor(args.out)
        slim = {k: v for k, v in m.items() if k != "sensitivity_variants"}
        print(json.dumps(slim, indent=2))
        return 0
    if args.paper:
        print(paper_macros(args.out))
        return 0
    if args.hi4pi:
        m = run_hi4pi(args.out, lon_step=args.lon_step, delete_after=args.delete_after)
        slim = {k: v for k, v in m.items() if not str(k).startswith("rc_")}
        print(json.dumps(slim, indent=2))
        return 0
    p.error("choose a mode: --anchor, --hi4pi, or --paper")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
