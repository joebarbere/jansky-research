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
    "parse_paper_tables",
    "rc_from_terminal",
    "rho_dm_local_gev",
    "rotation_curve_weighted",
    "run_anchor",
    "sensitivity_scan",
    "synthetic_spectrum",
    "threshold_tvm",
]

# Galactic constants used by the paper (their eq. 6-9): keep identical for the anchor.
R0_PC = 8178.0
V0_KMS = 235.1
G_KPC = 4.30091e-6  # G in kpc (km/s)^2 / Msun
GEV_PER_MSUN_PC3 = 38.0  # 1 Msun/pc^3 = 38.0 GeV/cm^3 (m_p c^2 conversion)


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
) -> float:
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
        refined = [(popt[j + 1], popt[j]) for j in range(0, len(popt), 3) if popt[j] >= min_amp_k]
        if refined:
            centers = [(c, a, 0.0) for c, a in refined]
    except RuntimeError:
        pass  # keep greedy centers if refinement fails to converge
    best = max(centers, key=lambda t: sign * t[0])
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
    ell = np.deg2rad(np.abs(np.asarray(longitudes_deg, float)))
    r = r0_pc * np.sin(ell)
    v = (np.abs(np.asarray(vterm_kms, float)) - sigma_v_kms) + v0_kms * np.sin(ell)
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
        bounds=([50, 50, 50, 1000, 10, 3000], [800, 2000, 600, 20000, 500, 100000]),
        maxfev=20000,
    )
    perr = np.sqrt(np.diag(pcov))
    names = ["v_bulge", "a_bulge", "v_disc", "a_disc", "v_halo", "h_halo"]
    out: dict = {n: float(v) for n, v in zip(names, popt, strict=True)}
    out.update({f"d{n}": float(e) for n, e in zip(names, perr, strict=True)})
    out["halo"] = halo
    resid = v_kms[ok] - f(r_pc[ok], *popt)
    out["rms_kms"] = float(np.sqrt(np.mean(resid**2)))
    out["rho_dm_gev"] = rho_dm_local_gev(out["v_halo"], out["h_halo"], halo=halo)
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
    rhos = [v["rho_dm_gev"] for v in variants.values() if "rho_dm_gev" in v]
    return {
        "variants": variants,
        "rho_dm_min_gev": float(min(rhos)),
        "rho_dm_max_gev": float(max(rhos)),
        "n_converged": len(rhos),
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
    paper_rho = rho_dm_local_gev(paper["v_halo"], paper["h_halo"])  # == their 0.107 if consistent
    metrics = {
        "source": "Sofue & Kohno 2025 published RC tables (arXiv:2509.23581 source, vendored)",
        "n_inner_rows": int(tables["inner"]["R_pc"].size),
        "n_unified_rows": int(uni["R_pc"].size),
        "anchor_fit": fit,
        "paper_table1_rms_kms": paper_rms,
        "paper_table1_rho_dm_gev": paper_rho,
        "paper_rho_dm_gev": 0.107,
        "rho_dm_ratio_vs_paper": fit["rho_dm_gev"] / 0.107,
        "r_max_pc": float(uni["R_pc"].max()),
        "sensitivity": {k: v for k, v in scan.items() if k != "variants"},
        "sensitivity_variants": {
            k: (
                {kk: vv for kk, vv in v.items() if kk in ("rho_dm_gev", "rms_kms", "halo")}
                if "rho_dm_gev" in v
                else v
            )
            for k, v in scan["variants"].items()
        },
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "innerrc_anchor.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Sofue & Kohno inner-RC replication (plan 86).")
    p.add_argument("--out", default=".")
    p.add_argument("--anchor", action="store_true", help="run the offline anchor leg")
    args = p.parse_args(argv)
    if args.anchor:
        m = run_anchor(args.out)
        slim = {k: v for k, v in m.items() if k != "sensitivity_variants"}
        print(json.dumps(slim, indent=2))
        return 0
    p.error("choose a mode (--anchor; the HI4PI leg lands in a later increment)")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
