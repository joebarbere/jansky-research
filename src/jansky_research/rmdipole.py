"""The first rotation-measure dipole/anisotropy test --- SPICE-RACS DR2 (plan 38).

The cosmic-dipole-anomaly literature is entirely source-count/flux based (Boehme+ arXiv:2509.16732;
arXiv:2509.18689; RMP colloquium arXiv:2505.23526; Mittal & Lewis arXiv:2605.27520); no published
test uses Faraday rotation measures as the tracer (GATE-0 full-text pass, 2026-07-05 --- closest
prior art is the uniform-cosmological-field RM dipole of Kronberg 1976 / Vallee 1990 / Kolatt 1998
on samples of ~10^2 RMs). SPICE-RACS DR2 (arXiv:2605.16917; ~3.3x10^5 goodRM components,
~2.5x10^5 at S/N>=8, 87.5% of sky) is the first catalogue large enough. The kinematic
expectation for RM is a Doppler rescaling RM_obs = RM(1+beta cos theta)^2 --- a signed dipole of
amplitude 2*beta ~ 2.5e-3, two orders below this test's sensitivity and largely absorbed by the
local nn subtraction --- so this is framed strictly as an ISOTROPY test with honest nulls: fit a
monopole+dipole model to the per-source extragalactic-RM residual power (rm^2 - sigma^2, noise-
debiased as in `rmstructure`) and to |residual|, with the significance calibrated by footprint-
preserving scrambles (the Dec <= +49 deg footprint is the dominant systematic --- the scramble is
load-bearing, not decorative). The scramble preserves Dec-band means, so it tests the
RA-PROJECTED dipole: conservative for a general dipole, blind to a polar-axis-aligned one
(the CMB apex at Dec -6.9 deg is nearly fully in-band).

Residuals: DR2 ships nearest-neighbour local-RM columns (`nn_rm_med` from the survey's own GRM
machinery), so the per-source extragalactic residual is rm - nn_rm_med at |b| >= b_min --- the
same local-subtraction convention as the DR2 team's CGM work (arXiv:2605.16924) --- with a
latitude-band median subtraction as the cross-check path. NOTE: local subtraction absorbs any
SIGNED mean-RM dipole by construction; the statistics tested here are dipoles in the residual
*amplitude/power*, which survive it. Plan deviation from healpy: the fit is per-source weighted
least squares on the monopole+dipole design matrix (no pixelisation dependency); coarse Dec-ring
maps are built only for the figure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "dipole_fit",
    "extragalactic_residuals",
    "fit_dipole",
    "footprint_scramble_null",
    "synthetic_dipole_catalogue",
    "compare_directions",
    "run",
]

DR2_LOCAL = Path("data/spice-racs.dr2.fits")
DR2_DAP = "https://data.csiro.au/collection/csiro:64891"
# Planck 2018 CMB dipole apex, equatorial (RA, Dec) deg (l=264.021, b=48.253)
CMB_DIPOLE_RA_DEC = (167.942, -6.944)


def _unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(np.asarray(ra_deg, float))
    dec = np.radians(np.asarray(dec_deg, float))
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=1)


def _radec_of(vec: np.ndarray) -> tuple[float, float]:
    x, y, z = vec / np.linalg.norm(vec)
    return float(np.degrees(np.arctan2(y, x)) % 360.0), float(np.degrees(np.arcsin(z)))


def dipole_fit(
    ra_deg: np.ndarray, dec_deg: np.ndarray, y: np.ndarray, *, n_boot: int = 0, seed: int = 0
) -> dict:
    r"""Weighted-equal least-squares fit of :math:`y_i = m + \vec p\cdot\hat n_i` per source.

    Returns the monopole ``m``, the dipole vector, the FRACTIONAL amplitude ``|p|/m``, and the
    apex direction (RA, Dec). With ``n_boot`` > 0, source-level bootstrap gives the amplitude SE
    and the angular scatter of the apex (deg). Equal weights by design: the y-variance of the
    power statistic is dominated by the intrinsic RM-power distribution, not the catalogued
    measurement errors, and the significance calibration is delegated to the scramble null.
    """
    n = _unit_vectors(ra_deg, dec_deg)
    a = np.column_stack([np.ones(len(y)), n])
    coef, *_ = np.linalg.lstsq(a, np.asarray(y, float), rcond=None)
    m, p = float(coef[0]), coef[1:]
    amp = float(np.linalg.norm(p) / abs(m)) if m != 0 else float("inf")
    ra_apex, dec_apex = _radec_of(p)
    out = {
        "monopole": m,
        "dipole_vec": p,
        "amp": amp,
        "apex_ra": ra_apex,
        "apex_dec": dec_apex,
        "n_sources": int(len(y)),
    }
    if n_boot:
        rng = np.random.default_rng(seed)
        amps, vecs = np.empty(n_boot), np.empty((n_boot, 3))
        for k in range(n_boot):
            pick = rng.integers(0, len(y), len(y))
            c, *_ = np.linalg.lstsq(a[pick], np.asarray(y, float)[pick], rcond=None)
            amps[k] = np.linalg.norm(c[1:]) / abs(c[0])
            vecs[k] = c[1:] / np.linalg.norm(c[1:])
        out["amp_se"] = float(amps.std())
        mean_vec = vecs.mean(axis=0)
        mean_vec /= np.linalg.norm(mean_vec)
        out["apex_scatter_deg"] = float(
            np.degrees(np.arccos(np.clip(vecs @ mean_vec, -1, 1)).mean())
        )
    return out


def extragalactic_residuals(
    s: dict, *, b_min: float = 45.0, method: str = "nn", nn_min: int = 5
) -> dict:
    """Per-source extragalactic-RM residuals at |b| >= b_min.

    ``method='nn'`` subtracts the catalogue's own nearest-neighbour local RM (`nn_rm_med`,
    requiring >= ``nn_min`` neighbours) --- the DR2-native GRM estimate. ``method='latitude'``
    subtracts the median RM in 5-deg |b| bands (the `rmstructure`-style Galactic-floor path,
    kept as the systematics cross-check). Any SIGNED mean-RM dipole is absorbed by either
    subtraction; downstream statistics must be amplitude/power dipoles.
    """
    ab = np.abs(np.asarray(s["gal_b"], float))
    m = ab >= b_min
    rm = np.asarray(s["rm"], float)
    if method == "nn":
        nn = np.asarray(s["nn_rm_med"], float)
        cnt = np.asarray(s["nn_rm_count"], float)
        m &= np.isfinite(nn) & (cnt >= nn_min)
        resid = rm - nn
    elif method == "latitude":
        resid = np.full(rm.size, np.nan)
        for lo in np.arange(0.0, 90.0, 5.0):
            band = (ab >= lo) & (ab < lo + 5.0)
            if band.sum() >= 10:
                resid[band] = rm[band] - np.median(rm[band])
        m &= np.isfinite(resid)
    else:
        raise ValueError(f"unknown method {method!r}")
    return {
        "ra": np.asarray(s["ra"], float)[m],
        "dec": np.asarray(s["dec"], float)[m],
        "gal_b": np.asarray(s["gal_b"], float)[m],
        "resid": resid[m],
        "rm_err": np.asarray(s["rm_err"], float)[m],
        "method": method,
        "b_min": b_min,
    }


def fit_dipole(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    resid: np.ndarray,
    rm_err: np.ndarray,
    *,
    stat: str = "power",
    clip_quantile: float | None = None,
    n_boot: int = 100,
    seed: int = 0,
) -> dict:
    """Dipole fit in residual-RM POWER (rm^2 - sigma^2, noise-debiased), |residual|, or NOISE.

    The power statistic is the primary one: its monopole is the extragalactic RM variance and
    its noise debias is exact (E[r_obs^2] = r^2 + sigma^2, mirroring the `rmstructure` SF
    convention). |residual| is the robustness companion (outlier-resistant, but its noise bias
    is not analytically removable --- interpret only relative to its own scramble null).
    ``stat='noise'`` fits the dipole of the catalogued noise power sigma^2 itself --- the
    systematics probe: a "power dipole" sharing its apex with the noise map is mis-modelled
    noise, not sky. ``clip_quantile`` drops sources with |residual| above that quantile before
    the fit (tail-drivenness diagnostic; the clip fraction is recorded).
    """
    r = np.asarray(resid, float)
    err = np.asarray(rm_err, float)
    ra = np.asarray(ra_deg, float)
    dec = np.asarray(dec_deg, float)
    n_clipped = 0
    if clip_quantile is not None:
        keep = np.abs(r) <= np.quantile(np.abs(r), clip_quantile)
        n_clipped = int((~keep).sum())
        r, err, ra, dec = r[keep], err[keep], ra[keep], dec[keep]
    if stat == "power":
        y = r**2 - err**2
    elif stat == "abs":
        y = np.abs(r)
    elif stat == "noise":
        y = err**2
    else:
        raise ValueError(f"unknown stat {stat!r}")
    out = dipole_fit(ra, dec, y, n_boot=n_boot, seed=seed)
    out["stat"] = stat
    out["y"] = y
    out["ra_used"], out["dec_used"] = ra, dec
    out["n_clipped"] = n_clipped
    return out


def footprint_scramble_null(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    y: np.ndarray,
    *,
    amp_obs: float,
    n_scramble: int = 200,
    dec_band_deg: float = 5.0,
    seed: int = 0,
) -> dict:
    """Footprint-preserving scramble null for the dipole amplitude.

    Permutes y AMONG the real source positions within Dec bands of width ``dec_band_deg``:
    the footprint (Dec <= +49 deg) and any Dec-dependent systematic (leakage vs beam distance,
    scan pattern) are preserved exactly, while all RA structure --- including a real dipole's
    RA projection --- is destroyed. The p-value is the fraction of scrambles whose fitted
    amplitude reaches ``amp_obs``. This is the load-bearing significance in the slice.
    """
    rng = np.random.default_rng(seed)
    dec = np.asarray(dec_deg, float)
    bands = np.floor((dec + 90.0) / dec_band_deg).astype(int)
    order = np.argsort(bands, kind="stable")
    sorted_bands = bands[order]
    amps = np.empty(n_scramble)
    y_arr = np.asarray(y, float)
    for k in range(n_scramble):
        perm = order.copy()
        # permute indices within each contiguous same-band run
        starts = np.flatnonzero(np.r_[True, sorted_bands[1:] != sorted_bands[:-1]])
        for s0, s1 in zip(starts, np.r_[starts[1:], len(perm)], strict=True):
            seg = perm[s0:s1]
            perm[s0:s1] = rng.permutation(seg)
        y_scr = np.empty_like(y_arr)
        y_scr[order] = y_arr[perm]
        amps[k] = dipole_fit(ra_deg, dec, y_scr)["amp"]
    p = float((amps >= amp_obs).sum() + 1) / (n_scramble + 1)
    return {"amps": amps, "amp_obs": float(amp_obs), "p_value": p, "n_scramble": n_scramble}


def synthetic_dipole_catalogue(
    n_sources: int = 20000,
    *,
    amp: float = 0.3,
    apex_ra: float = 170.0,
    apex_dec: float = -7.0,
    sigma0: float = 13.0,
    noise: float = 2.0,
    dec_max: float = 49.0,
    ra_deg: np.ndarray | None = None,
    dec_deg: np.ndarray | None = None,
    seed: int = 0,
) -> dict:
    r"""Mock residual catalogue with a KNOWN power dipole on a DR2-like (or the real) footprint.

    Residuals are drawn with position-dependent variance
    :math:`\sigma^2(\hat n) = \sigma_0^2\,(1 + A\,\hat n\cdot\hat d)` --- so the noise-debiased
    power statistic has fractional dipole amplitude exactly ``amp`` toward the apex --- plus
    Gaussian measurement noise recorded in ``rm_err`` (the debiasing target). Pass the real DR2
    ``ra_deg``/``dec_deg`` to inject on the true footprint; otherwise positions are uniform on
    the sphere south of ``dec_max`` (the DR2 declination limit).
    """
    rng = np.random.default_rng(seed)
    if ra_deg is None or dec_deg is None:
        ra = rng.uniform(0.0, 360.0, n_sources)
        smax = np.sin(np.radians(dec_max))
        dec = np.degrees(np.arcsin(rng.uniform(-1.0, smax, n_sources)))
    else:
        ra, dec = np.asarray(ra_deg, float), np.asarray(dec_deg, float)
        n_sources = ra.size
    d_hat = _unit_vectors(np.array([apex_ra]), np.array([apex_dec]))[0]
    mu = _unit_vectors(ra, dec) @ d_hat
    var = sigma0**2 * (1.0 + amp * mu)
    resid = rng.normal(0.0, np.sqrt(var)) + rng.normal(0.0, noise, n_sources)
    return {
        "ra": ra,
        "dec": dec,
        "resid": resid,
        "rm_err": np.full(n_sources, noise),
        "true_amp": amp,
        "true_apex": (apex_ra, apex_dec),
    }


def compare_directions(apex_ra: float, apex_dec: float) -> dict:
    """Angular separation of a fitted apex from the CMB dipole apex (deg)."""
    v = _unit_vectors(np.array([apex_ra]), np.array([apex_dec]))[0]
    cmb = _unit_vectors(*(np.array([x]) for x in CMB_DIPOLE_RA_DEC))[0]
    return {"sep_cmb_deg": float(np.degrees(np.arccos(np.clip(v @ cmb, -1.0, 1.0))))}


def load_dr2_for_dipole(
    path: str | Path = DR2_LOCAL, *, snr_min: float = 8.0
) -> dict:  # pragma: no cover - needs the 9.3M-row DAP file
    """Load DR2 goodRM rows with the nn_rm_* local-GRM columns (CSIRO DAP csiro:64891)."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        good = np.asarray(d["goodRM_flag"], bool) & (np.asarray(d["snr_polint"], float) >= snr_min)
        cols = ("ra", "dec", "l", "b", "rm", "rm_err", "nn_rm_med", "nn_rm_count")
        out = {c: np.asarray(d[c], float)[good] for c in cols}
    out["gal_l"], out["gal_b"] = out.pop("l"), out.pop("b")
    m = np.isfinite(out["rm"]) & (out["rm_err"] > 0)
    return {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in out.items()}


def _dipole_leg(
    res: dict, *, stat: str, n_scramble: int, clip_quantile: float | None = None, seed: int = 0
) -> dict:
    fit = fit_dipole(
        res["ra"],
        res["dec"],
        res["resid"],
        res["rm_err"],
        stat=stat,
        clip_quantile=clip_quantile,
        seed=seed,
    )
    null = footprint_scramble_null(
        fit["ra_used"],
        fit["dec_used"],
        fit["y"],
        amp_obs=fit["amp"],
        n_scramble=n_scramble,
        seed=seed,
    )
    sep = compare_directions(fit["apex_ra"], fit["apex_dec"])
    return {
        "stat": stat,
        "clip_quantile": clip_quantile,
        "n_clipped": fit["n_clipped"],
        "n_sources": fit["n_sources"],
        "amp": round(fit["amp"], 4),
        "amp_se": round(fit.get("amp_se", float("nan")), 4),
        "apex_ra": round(fit["apex_ra"], 1),
        "apex_dec": round(fit["apex_dec"], 1),
        "apex_scatter_deg": round(fit.get("apex_scatter_deg", float("nan")), 1),
        "sep_cmb_deg": round(sep["sep_cmb_deg"], 1),
        "p_scramble": round(null["p_value"], 4),
        "null_amps": null["amps"],
    }


def run(out: str = ".", *, offline: bool = True, n_scramble: int = 200) -> dict:
    """Offline: injected-dipole recover-a-known on a DR2-like footprint; real: DR2 both stats."""
    import json

    if offline:
        syn = synthetic_dipole_catalogue(seed=0)
        res = {k: syn[k] for k in ("ra", "dec", "resid", "rm_err")}
        legs = [_dipole_leg(res, stat="power", n_scramble=n_scramble)]
        # a no-dipole control on the same footprint: the scramble p must NOT be small
        ctrl = synthetic_dipole_catalogue(amp=0.0, seed=1)
        res0 = {k: ctrl[k] for k in ("ra", "dec", "resid", "rm_err")}
        legs.append(_dipole_leg(res0, stat="power", n_scramble=n_scramble, seed=1))
        source = "synthetic dipole catalogue (DR2-like footprint)"
        extra = {
            "true_amp": syn["true_amp"],
            "true_apex_ra": syn["true_apex"][0],
            "true_apex_dec": syn["true_apex"][1],
            "sep_true_deg": round(
                float(
                    np.degrees(
                        np.arccos(
                            np.clip(
                                _unit_vectors(
                                    np.array([legs[0]["apex_ra"]]),
                                    np.array([legs[0]["apex_dec"]]),
                                )[0]
                                @ _unit_vectors(
                                    np.array([syn["true_apex"][0]]),
                                    np.array([syn["true_apex"][1]]),
                                )[0],
                                -1,
                                1,
                            )
                        )
                    )
                ),
                1,
            ),
        }
    else:  # pragma: no cover - needs the local 9.3M-row DAP file
        s = load_dr2_for_dipole()
        legs = []
        # primary first: b>=45, catalogue-native nn subtraction (legs[0] feeds macros/figure)
        for b_min in (45.0, 30.0):
            for method in ("nn", "latitude"):
                res = extragalactic_residuals(s, b_min=b_min, method=method)
                for stat in ("power", "abs"):
                    leg = _dipole_leg(res, stat=stat, n_scramble=n_scramble)
                    leg.update({"b_min": b_min, "method": method})
                    legs.append(leg)
        # tail-drivenness + systematics diagnostics on the primary sample (b>=45, nn):
        # clipped power (is the power dipole carried by the |residual| tail?) and the
        # noise-power dipole (does the catalogued sigma^2 map share the apex?)
        prim = extragalactic_residuals(s, b_min=45.0, method="nn")
        for stat, clip in (("power", 0.99), ("noise", None)):
            leg = _dipole_leg(prim, stat=stat, n_scramble=n_scramble, clip_quantile=clip)
            leg.update({"b_min": 45.0, "method": "nn"})
            legs.append(leg)
        # recover-a-known ON THE REAL FOOTPRINT: inject a known power dipole at the actual
        # source positions of the primary sample and confirm amplitude+direction come back
        inj = synthetic_dipole_catalogue(amp=0.3, ra_deg=prim["ra"], dec_deg=prim["dec"], seed=0)
        leg = _dipole_leg(
            {k: inj[k] for k in ("ra", "dec", "resid", "rm_err")},
            stat="power",
            n_scramble=n_scramble,
        )
        leg.update({"b_min": 45.0, "method": "injection", "true_amp": inj["true_amp"]})
        legs.append(leg)
        source = f"SPICE-RACS DR2 goodRM ({DR2_DAP})"
        extra = {}

    metrics = {
        "source": source,
        "is_real": not offline,
        "legs": [{k: v for k, v in leg.items() if k != "null_amps"} for leg in legs],
        **extra,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "rmdipole_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(legs, op / "papers" / "rmdipole" / "figures")
    _write_macros(metrics, op / "papers" / "rmdipole" / "generated" / "macros.tex")
    _write_legs_table(metrics, op / "papers" / "rmdipole" / "generated" / "legs_table.tex")
    return metrics


def _write_legs_table(m: dict, path: str | Path) -> None:
    """Emit the per-leg results as LaTeX tabular rows (every paper number is pipeline-made)."""
    rows = [
        "% Auto-generated by jansky_research.rmdipole._write_legs_table -- do not edit.",
        "% Columns: sample & subtraction & statistic & N & amp+-se & apex & sep_CMB & p_scramble",
    ]
    for leg in m["legs"]:
        method = str(leg.get("method", "synthetic"))
        b_min = leg.get("b_min")
        sample = rf"$|b|\ge{b_min:.0f}\arcdeg$" if b_min is not None else "mock footprint"
        stat = str(leg["stat"])
        if leg.get("clip_quantile"):
            stat += rf" (clip {100 * leg['clip_quantile']:.0f}\%)"
        if method == "injection":
            stat += rf" [inj.\ {leg.get('true_amp')}]"
        amp = f"${leg['amp']:.3f} \\pm {leg['amp_se']:.3f}$"
        apex = rf"$({leg['apex_ra']:.0f}\arcdeg, {leg['apex_dec']:+.0f}\arcdeg)$"
        sep = rf"${leg['sep_cmb_deg']:.0f}\arcdeg$"
        rows.append(
            f"{sample} & {method} & {stat} & {leg['n_sources']} & {amp} & {apex} & {sep} & "
            f"${leg['p_scramble']:.4g}$ \\\\"
        )
    # the closing \hline lives HERE, not in main.tex: an \hline straight after \input inside
    # tabular is "Misplaced \noalign" (the file boundary breaks the post-\\ lookahead)
    rows.append(r"\hline")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(rows) + "\n")


def _figure(legs: list[dict], out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.4, 3.9))
    lead = legs[0]
    ax.hist(lead["null_amps"], bins=30, color="0.7", label="footprint-scramble null")
    ax.axvline(lead["amp"], color="C3", lw=2, label=f"observed (p={lead['p_scramble']})")
    ax.set(
        xlabel="fractional dipole amplitude |p|/m",
        ylabel="scrambles",
        title=f"RM-power dipole vs null ({lead['stat']})",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "rmdipole.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(leg: dict | None, key: str) -> str:
        if leg is None or leg.get(key) is None:
            return "--"
        v = leg.get(key)
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "rmdReal" if m.get("is_real") else "rmdSyn"
    lead = m["legs"][0] if m.get("legs") else None
    # the injection control: the dedicated injection leg (real runs) or the lead leg itself,
    # which IS the injected-dipole recovery in the offline/synthetic run
    inj = next((x for x in m.get("legs", []) if x.get("method") == "injection"), None)
    if inj is None and not m.get("is_real"):
        inj = dict(lead or {}, true_amp=m.get("true_amp"))
    legs = m.get("legs", [])
    clip = next((x for x in legs if x.get("clip_quantile")), None)
    noise = next((x for x in legs if x.get("stat") == "noise"), None)
    wabs = next((x for x in legs if x.get("stat") == "abs" and x.get("b_min") == 30.0), None)
    lines = [
        "% Auto-generated by jansky_research.rmdipole._write_macros -- do not edit.",
        "% Synthetic (rmdSyn*) and real (rmdReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under rmdReal* (an offline rebuild resets rmdReal* to placeholders by design).",
        rf"\newcommand{{\rmdSource}}{{{m['source']}}}",
    ]
    for ns in ("rmdSyn", "rmdReal"):
        leg = lead if ns == pref else None
        ileg = inj if ns == pref else None
        lines += [
            rf"\newcommand{{\{ns}N}}{{{g(leg, 'n_sources')}}}",
            rf"\newcommand{{\{ns}Amp}}{{{g(leg, 'amp')}}}",
            rf"\newcommand{{\{ns}AmpSe}}{{{g(leg, 'amp_se')}}}",
            rf"\newcommand{{\{ns}ApexRa}}{{{g(leg, 'apex_ra')}}}",
            rf"\newcommand{{\{ns}ApexDec}}{{{g(leg, 'apex_dec')}}}",
            rf"\newcommand{{\{ns}SepCmb}}{{{g(leg, 'sep_cmb_deg')}}}",
            rf"\newcommand{{\{ns}PScramble}}{{{g(leg, 'p_scramble')}}}",
            rf"\newcommand{{\{ns}InjAmp}}{{{g(ileg, 'amp')}}}",
            rf"\newcommand{{\{ns}InjAmpSe}}{{{g(ileg, 'amp_se')}}}",
            rf"\newcommand{{\{ns}InjTrue}}{{{g(ileg, 'true_amp')}}}",
            rf"\newcommand{{\{ns}ClipAmp}}{{{g(clip if ns == pref else None, 'amp')}}}",
            rf"\newcommand{{\{ns}ClipP}}{{{g(clip if ns == pref else None, 'p_scramble')}}}",
            rf"\newcommand{{\{ns}ClipN}}{{{g(clip if ns == pref else None, 'n_clipped')}}}",
            rf"\newcommand{{\{ns}NoiseAmp}}{{{g(noise if ns == pref else None, 'amp')}}}",
            rf"\newcommand{{\{ns}NoiseP}}{{{g(noise if ns == pref else None, 'p_scramble')}}}",
            rf"\newcommand{{\{ns}WideAbsP}}{{{g(wabs if ns == pref else None, 'p_scramble')}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="RM dipole/anisotropy test on SPICE-RACS DR2.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--n-scramble", type=int, default=200)
    args = p.parse_args(argv)
    m = run(args.out, offline=args.offline, n_scramble=args.n_scramble)
    print(json.dumps(m, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
