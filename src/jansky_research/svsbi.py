"""SBI population inference for the RACS Stokes-V coherent-emitter class (plan 40).

Neural simulation-based inference (SBI / neural posterior estimation) has been applied to pulsar
populations (arXiv:2312.14848, 2412.04070), magnetars (arXiv:2503.11875), FRB selection functions
(arXiv:2606.26334), and per-source QU-fitting (VROOM-SBI, arXiv:2605.27538) --- but no SBI
POPULATION inference exists for circularly-polarized stellar/coherent emitters (GATE-0 full-text
pass, 2026-07-08). The classical prior art is qualitative: Callingham+2021, Pritchard+2021, and
Driessen+2024 report detection statistics, but there is no calibrated posterior on the population
parameters. This slice infers three: the coherent radio luminosity-function slope and break, and
the **beaming fraction** (what fraction of coherent emitters beam toward us) --- conditioned on the
merged `stokesv_discovery` census (60 nearby M dwarfs x 2 RACS-mid epochs; 39 with usable
two-epoch Stokes-V measurements, of which 2 are confident detections).

Two repo assets are the two hard inputs: the validated forced-V measurement
(`stokesv.measure_circular_pol`) and the per-field leakage floor (`stokesv.leakage_floor`, ~5.7%),
which is effectively the selection function. The forward model folds a drawn population through
each REAL target's distance (Gaia), local V rms, and leakage floor to predict which stars are
detected; SBI infers the population from the observed detection pattern + the detections'
fluxes. With N=39 and 2 detections the posterior is WIDE and carries a beaming--luminosity
degeneracy (few luminous beamers vs many faint ones): the deliverable is the FIRST calibrated
posterior with SBC-validated coverage, not a tight number --- stated from the outset.

The physics (`draw_population`, `forward_model`, `summary_stats`) is pure NumPy and tested in
core CI; only `train_npe`/`build_posterior` touch the `sbi` package (ROCm-verified: sbi 0.26.1
trains NPE on gfx1102 unchanged). Parameter vector theta = (lf_slope, log10 L_break, f_beam),
with L in erg/s/Hz and f_beam in [0, 1].
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "THETA_NAMES",
    "prior_bounds",
    "draw_population",
    "forward_model",
    "summary_stats",
    "simulate",
    "parent_from_census",
    "sbc_ranks",
    "run",
]

THETA_NAMES = ("lf_slope", "log_Lbreak", "f_beam")
# prior ranges: slope of the (differential) luminosity function dN/dL ~ L^-slope over the
# coherent-emitter luminosity range; log10 break luminosity (erg/s/Hz); beaming fraction
PRIOR_LOW = np.array([1.0, 12.0, 0.02])
PRIOR_HIGH = np.array([3.0, 15.0, 0.60])
CENSUS_CSV = Path("results/stokesv_discovery_realtargets.csv")
PC_CM = 3.0856775814913673e18  # cm per parsec
DET_NSIGMA = 5.0  # a Stokes-V detection: |V| above DET_NSIGMA * rms AND above the leakage floor


def prior_bounds() -> tuple[np.ndarray, np.ndarray]:
    """The (low, high) box-prior bounds on theta = (lf_slope, log_Lbreak, f_beam)."""
    return PRIOR_LOW.copy(), PRIOR_HIGH.copy()


def _sample_luminosity(
    n: int, slope: float, log_lbreak: float, rng: np.random.Generator, *, log_lmin: float = 12.0
) -> np.ndarray:
    r"""Draw coherent radio luminosities (erg/s/Hz) from a truncated power law with a break.

    :math:`dN/dL \propto L^{-\mathrm{slope}}` for :math:`L_\min \le L \le L_\mathrm{break}`; the
    break is an exponential cutoff above :math:`L_\mathrm{break}` (a Schechter-like taper), so the
    slope controls the faint end and ``log_lbreak`` the bright cutoff. Inverse-CDF over a log
    grid --- the :math:`dL` measure carries the ``np.gradient(grid)`` Jacobian so the REALISED
    differential slope is exactly ``slope`` (without it the log-grid CDF would realise
    :math:`L^{-(\mathrm{slope}+1)}`).
    """
    grid = np.logspace(log_lmin, log_lbreak + 1.0, 512)
    lbreak = 10.0**log_lbreak
    pdf = grid ** (-slope) * np.exp(-grid / lbreak) * np.gradient(grid)  # dN = (dN/dL) dL
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def draw_population(theta: np.ndarray, distances_pc: np.ndarray, rng: np.random.Generator) -> dict:
    """Draw a coherent-emitter realisation for the real parent sample.

    Each of the ``len(distances_pc)`` stars is assigned a coherent luminosity from the
    luminosity function and beams toward us with probability ``f_beam``; a non-beaming star has
    zero observed flux. The observed Stokes-V flux density (mJy) is
    :math:`S = L / (4\\pi d^2)` converted to mJy (1 erg/s/Hz/cm^2 = 1e26 mJy). Returns the
    per-star observed flux and the beaming mask.
    """
    slope, log_lbreak, f_beam = float(theta[0]), float(theta[1]), float(theta[2])
    d_cm = np.asarray(distances_pc, float) * PC_CM
    n = d_cm.size
    lum = _sample_luminosity(n, slope, log_lbreak, rng)
    beams = rng.random(n) < f_beam
    flux_mjy = np.where(beams, lum / (4.0 * np.pi * d_cm**2) * 1.0e26, 0.0)
    return {"flux_mjy": flux_mjy, "beams": beams}


def forward_model(
    theta: np.ndarray,
    parent: dict,
    rng: np.random.Generator,
) -> dict:
    """Fold a drawn population through the real per-target selection: which stars are detected.

    Adds Gaussian V noise at each target's ``v_rms`` (two epochs; the best epoch is what a
    variability census keeps), and calls a source detected if the noisy |V| clears both
    ``DET_NSIGMA`` * rms and the target's leakage floor (``leakage_frac`` * I-flux proxy). The
    leakage floor is the selection function reused from `stokesv`. Returns per-target detection
    flags and observed |V|.
    """
    pop = draw_population(theta, parent["distance_pc"], rng)
    v_true = pop["flux_mjy"]
    v_rms = np.asarray(parent["v_rms"], float)
    n_epochs = int(parent.get("n_epochs", 2))
    # observe over n_epochs, keep the max-|V| epoch (the two-epoch variability census logic)
    obs = v_true[None, :] + rng.normal(0.0, v_rms[None, :], (n_epochs, v_rms.size))
    best = obs[np.argmax(np.abs(obs), axis=0), np.arange(v_rms.size)]
    floor = np.asarray(parent["leakage_floor_mjy"], float)
    detected = (np.abs(best) >= DET_NSIGMA * v_rms) & (np.abs(best) >= floor)
    return {"detected": detected, "v_obs": best, "v_true": v_true}


def summary_stats(fm: dict) -> np.ndarray:
    """Reduce a forward-model realisation to a fixed-length summary vector for SBI.

    (n_detections; log10 of the brightest detected |V| or a floor if none; log10 of the summed
    detected |V|; the detected fraction) --- enough to constrain the beaming fraction (via the
    count) and the luminosity function (via the brightest-flux distribution).
    """
    det = fm["detected"]
    n_det = int(det.sum())
    if n_det:
        vmax = float(np.max(np.abs(fm["v_obs"][det])))
        vsum = float(np.sum(np.abs(fm["v_obs"][det])))
    else:
        vmax = vsum = 1e-3
    n = det.size
    return np.array(
        [n_det, np.log10(max(vmax, 1e-3)), np.log10(max(vsum, 1e-3)), n_det / max(n, 1)],
        dtype=np.float32,
    )


def simulate(theta: np.ndarray, parent: dict, *, seed: int = 0) -> np.ndarray:
    """One theta -> summary vector (draw a population, observe it, summarise). Deterministic in seed."""
    rng = np.random.default_rng(seed)
    return summary_stats(forward_model(np.asarray(theta, float), parent, rng))


def _fetch_gaia_distances(gaia_ids: list[str]) -> dict:  # pragma: no cover - network
    """Gaia DR3 distances (pc) for a list of source_ids, via VizieR (1000/parallax_mas)."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["Source", "Plx"])
    v.ROW_LIMIT = -1
    out: dict[str, float] = {}
    ids = [g for g in gaia_ids if g and g not in ("", "nan")]
    if not ids:
        return out
    tab = v.query_constraints(catalog="I/355/gaiadr3", Source="=,".join(ids))
    if tab:
        t = tab[0]
        for src, plx in zip(t["Source"], t["Plx"], strict=True):
            if np.isfinite(plx) and plx > 0:
                out[str(int(src))] = 1000.0 / float(plx)
    return out


def parent_from_census(
    csv_path: str | Path = CENSUS_CSV,
    *,
    default_distance_pc: float = 10.0,
    fetch_distances: bool = True,
) -> dict:
    """Build the SBI parent sample from the merged `stokesv_discovery` census CSV.

    One row per PHYSICAL target: per-target V rms (median over epochs, mJy), a leakage-floor flux
    (leakage fraction * the I-flux proxy), a distance (Gaia DR3 parallax by ``gaia_id`` when
    ``fetch_distances``, else ``default_distance_pc``), and the observed best-epoch |V|.

    Unresolved binaries are DEDUPLICATED: two catalogue targets whose forced photometry returns
    byte-identical V and I at both epochs are the same radio source (the RACS beam does not
    resolve a ~2" pair, e.g. GJ 65 = CNS5 424/425) and are collapsed to one entry --- otherwise
    the detection count that drives the beaming-fraction inference is inflated.
    """
    import csv

    rows = list(csv.DictReader(open(csv_path)))
    by: dict[str, list] = {}
    for r in rows:
        by.setdefault(r["name"], []).append(r)
    per_target: list[dict[str, Any]] = []
    for name, rs in by.items():
        evs = [float(r["e_v"]) for r in rs if r["e_v"] not in ("nan", "")]
        vvs = [float(r["v_mjy"]) for r in rs if r["v_mjy"] not in ("nan", "")]
        ivs = [float(r["i_mjy"]) for r in rs if r["i_mjy"] not in ("nan", "")]
        if not evs or not vvs:
            continue
        gaia = next((r["gaia_id"] for r in rs if r.get("gaia_id")), "")
        per_target.append(
            {
                "name": name,
                "gaia_id": gaia,
                "v_rms": float(np.median(evs)),
                "floor": 0.057 * max(np.nanmax(ivs) if ivs else 0.0, 0.0),
                "v_best": max(abs(v) for v in vvs),
                "flux_key": tuple(round(v, 4) for v in sorted(vvs)),  # unresolved-pair fingerprint
            }
        )
    # dedupe unresolved binaries: keep the first target of each identical-flux group
    seen: set = set()
    uniq = []
    for t in per_target:
        # only dedupe bright (detected) dupes -- faint targets sharing a rounded flux are noise
        if float(t["v_best"]) > 1.0 and t["flux_key"] in seen:
            continue
        seen.add(t["flux_key"])
        uniq.append(t)

    dist_map = _fetch_gaia_distances([str(t["gaia_id"]) for t in uniq]) if fetch_distances else {}
    dist = [dist_map.get(str(t["gaia_id"]).split(".")[0], default_distance_pc) for t in uniq]
    return {
        "name": np.array([t["name"] for t in uniq]),
        "v_rms": np.array([t["v_rms"] for t in uniq]),
        "leakage_floor_mjy": np.array([t["floor"] for t in uniq]),
        "distance_pc": np.array(dist, float),
        "v_best_obs": np.array([t["v_best"] for t in uniq]),
        "n_epochs": 2,
    }


def observed_summary(parent: dict) -> np.ndarray:
    """The real data's summary vector: detections in the actual census (not a simulation)."""
    v = np.asarray(parent["v_best_obs"], float)
    rms = np.asarray(parent["v_rms"], float)
    floor = np.asarray(parent["leakage_floor_mjy"], float)
    detected = (v >= DET_NSIGMA * rms) & (v >= floor)
    return summary_stats({"detected": detected, "v_obs": v})


# --------------------------------------------------------------------------------------------
# SBI: NPE training + simulation-based calibration (the `sbi`-package leg)
# --------------------------------------------------------------------------------------------


def _require_sbi() -> tuple[Any, Any]:
    try:
        import torch
        from sbi.inference import NPE
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "jansky_research.svsbi NPE training needs the `sbi` extra: "
            "`uv sync --extra sbi` (or a ROCm venv with sbi installed)"
        ) from exc
    return torch, NPE


def train_npe(
    parent: dict,
    *,
    n_sims: int = 5000,
    device: str = "cpu",
    seed: int = 0,
    max_epochs: int = 200,
):  # pragma: no cover - needs the sbi extra + is a GPU training job
    """Train an NPE posterior estimator on ``n_sims`` (theta, summary) pairs. Returns the posterior."""
    torch, NPE = _require_sbi()
    from sbi.utils import BoxUniform

    rng = np.random.default_rng(seed)
    low = torch.as_tensor(PRIOR_LOW, dtype=torch.float32, device=device)
    high = torch.as_tensor(PRIOR_HIGH, dtype=torch.float32, device=device)
    prior = BoxUniform(low=low, high=high, device=device)
    thetas = rng.uniform(PRIOR_LOW, PRIOR_HIGH, size=(n_sims, 3))
    x = np.stack([simulate(t, parent, seed=int(rng.integers(1 << 31))) for t in thetas])
    inf = NPE(prior=prior, device=device)
    inf.append_simulations(
        torch.as_tensor(thetas, dtype=torch.float32, device=device),
        torch.as_tensor(x, dtype=torch.float32, device=device),
    ).train(max_num_epochs=max_epochs)
    return inf.build_posterior()


def sbc_ranks(
    parent: dict,
    posterior,  # noqa: ANN001 - an sbi posterior
    *,
    n_trials: int = 200,
    n_post: int = 200,
    device: str = "cpu",
    seed: int = 0,
):  # pragma: no cover - needs a trained posterior
    """Simulation-based calibration: rank of each true theta within its posterior samples.

    For uniform coverage the ranks are uniform on [0, n_post]; a KS test against uniform is the
    calibration check. Returns the per-parameter rank arrays.
    """
    import torch

    rng = np.random.default_rng(seed)
    ranks = np.zeros((n_trials, 3), int)
    for i in range(n_trials):
        theta_true = rng.uniform(PRIOR_LOW, PRIOR_HIGH)
        x = simulate(theta_true, parent, seed=int(rng.integers(1 << 31)))
        samples = posterior.sample(
            (n_post,),
            x=torch.as_tensor(x, dtype=torch.float32, device=device),
            show_progress_bars=False,
        )
        s = np.asarray(samples.cpu())
        ranks[i] = (s < theta_true[None, :]).sum(axis=0)
    return ranks


def _ks_uniform(ranks_1d: np.ndarray, n_post: int) -> float:
    """KS distance of a rank distribution from uniform [0, n_post] (0 = perfectly calibrated)."""
    r = np.sort(ranks_1d / n_post)
    n = r.size
    cdf_emp = np.arange(1, n + 1) / n
    return float(np.max(np.abs(cdf_emp - r)))


def run(out: str = ".", *, offline: bool = True, device: str = "cpu", n_sims: int = 3000) -> dict:
    """Offline: physics recover-a-known + prior-predictive checks; real: NPE + SBC on the census."""
    import json

    # a synthetic parent sample (offline) or the real stokesv_discovery census (real leg)
    if offline:
        parent = _synthetic_parent(seed=0)
        source = "synthetic parent sample + prior-predictive forward model"
    else:  # pragma: no cover - needs the census CSV
        parent = parent_from_census()
        source = f"stokesv_discovery census ({parent['name'].size} M dwarfs, RACS-mid V)"

    # physics recover-a-known (always, pure NumPy): the forward model must respond monotonically
    # to beaming fraction and luminosity. This validates the MACHINERY, so it runs on a
    # detection-rich synthetic parent regardless of leg (the real 39-star census has too few
    # detections to average the monotonic trend out of Poisson noise)
    checks = _physics_checks(_synthetic_parent(seed=99))

    metrics: dict = {
        "source": source,
        "is_real": not offline,
        "n_targets": int(parent["v_rms"].size),
        "median_v_rms_mjy": round(float(np.median(parent["v_rms"])), 3),
        **checks,
    }

    if not offline:  # pragma: no cover - needs the sbi extra (GPU/CPU training)
        metrics.update(_real_leg(parent, device=device, n_sims=n_sims))
        n_det = int(observed_summary(parent)[0])
        metrics["n_v_detections_observed"] = n_det

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "svsbi_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(parent, metrics, op / "papers" / "svsbi" / "figures")
    _write_macros(metrics, op / "papers" / "svsbi" / "generated" / "macros.tex")
    return metrics


def _synthetic_parent(*, n: int = 400, seed: int = 0) -> dict:
    """A stokesv_discovery-like synthetic parent sample (distances, rms, leakage floors).

    Larger than the real 39-target census so the pure-NumPy physics checks average over O(10)
    detections rather than the census's O(1); the real leg uses the actual census.
    """
    rng = np.random.default_rng(seed)
    return {
        "name": np.array([f"SYN-{i}" for i in range(n)]),
        "v_rms": rng.uniform(0.10, 0.30, n),  # mJy, RACS-mid-like
        "leakage_floor_mjy": rng.uniform(0.02, 0.20, n),
        "distance_pc": np.exp(rng.uniform(np.log(2.0), np.log(15.0), n)),  # nearby M dwarfs
        "n_epochs": 2,
    }


def _physics_checks(parent: dict) -> dict:
    """Pure-NumPy recover-a-known: the forward model + summary behave as physics requires."""
    rng = np.random.default_rng(1)

    def mean_ndet(theta, k=120):
        return float(
            np.mean([summary_stats(forward_model(theta, parent, rng))[0] for _ in range(k)])
        )

    # evaluated in the detectable regime (log L_break ~ 14+, where a nearby-M-dwarf population
    # produces detections at the RACS 5-sigma floor; fainter breaks yield ~0 detections)
    lo_beam = mean_ndet(np.array([2.0, 14.0, 0.05]))
    hi_beam = mean_ndet(np.array([2.0, 14.0, 0.50]))
    lo_lum = mean_ndet(np.array([2.0, 12.5, 0.50]))
    hi_lum = mean_ndet(np.array([2.0, 14.5, 0.50]))
    return {
        "ndet_low_beaming": round(lo_beam, 2),
        "ndet_high_beaming": round(hi_beam, 2),
        "ndet_faint_lf": round(lo_lum, 2),
        "ndet_bright_lf": round(hi_lum, 2),
        "beaming_monotonic": bool(hi_beam > lo_beam),
        "luminosity_monotonic": bool(hi_lum > lo_lum),
    }


def _real_leg(parent: dict, *, device: str, n_sims: int) -> dict:  # pragma: no cover - sbi extra
    """NPE posterior on the real census + SBC coverage."""
    import torch

    posterior = train_npe(parent, n_sims=n_sims, device=device)
    x_o = torch.as_tensor(observed_summary(parent), dtype=torch.float32, device=device)
    post_samples = np.asarray(posterior.sample((4000,), x=x_o, show_progress_bars=False).cpu())
    n_sbc = 150
    ranks = sbc_ranks(parent, posterior, n_trials=n_sbc, n_post=n_sbc, device=device)
    ks = {THETA_NAMES[j]: round(_ks_uniform(ranks[:, j], n_sbc), 3) for j in range(3)}
    # the formal one-sample KS 95% critical value at n trials (Kolmogorov): 1.358/sqrt(n)
    ks_crit = 1.358 / np.sqrt(n_sbc)
    q = {
        THETA_NAMES[j]: [round(float(np.percentile(post_samples[:, j], p)), 3) for p in (5, 50, 95)]
        for j in range(3)
    }
    # posterior/prior width ratio per parameter: ~1 means the data barely update the prior
    prior_w = PRIOR_HIGH - PRIOR_LOW  # 90% of a uniform ~= 0.9 of the full range
    width_ratio = {
        THETA_NAMES[j]: round((q[THETA_NAMES[j]][2] - q[THETA_NAMES[j]][0]) / (0.9 * prior_w[j]), 2)
        for j in range(3)
    }
    return {
        "posterior_median": {k: v[1] for k, v in q.items()},
        "posterior_ci90": {k: [v[0], v[2]] for k, v in q.items()},
        "posterior_prior_width_ratio": width_ratio,
        "sbc_ks": ks,
        "sbc_ks_max": round(float(max(ks.values())), 3),
        "sbc_ks_crit": round(float(ks_crit), 3),
        "sbc_calibrated": bool(max(ks.values()) < ks_crit),  # formal KS 95% critical value
        "n_sims": n_sims,
    }


def _figure(parent: dict, m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    # prior-predictive: n_detections vs beaming fraction at fixed LF
    rng = np.random.default_rng(3)
    fbeams = np.linspace(PRIOR_LOW[2], PRIOR_HIGH[2], 12)
    nd = [
        np.mean(
            [
                summary_stats(forward_model(np.array([2.0, 13.5, fb]), parent, rng))[0]
                for _ in range(30)
            ]
        )
        for fb in fbeams
    ]
    ax1.plot(fbeams, nd, "o-", color="C0")
    ax1.set(xlabel="beaming fraction", ylabel="mean # detections", title="Prior-predictive")
    post = m.get("posterior_median")
    if post is not None:
        ci = m["posterior_ci90"]
        ys = list(range(3))
        for j, k in enumerate(THETA_NAMES):
            ax2.plot(ci[k], [j, j], "-", color="C3", lw=2)
            ax2.plot(post[k], j, "o", color="C3")
        ax2.set(yticks=ys, yticklabels=list(THETA_NAMES), title="Posterior (90% CI)")
    else:
        ax2.text(0.5, 0.5, "real NPE leg pending", ha="center", va="center")
        ax2.set_axis_off()
    fig.tight_layout()
    fig.savefig(out / "svsbi.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(d: str, key: str | None = None) -> str:
        sub = m.get(d)
        v = sub.get(key) if isinstance(sub, dict) and key is not None else sub
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    ci = m.get("posterior_ci90") or {}

    def ci_bound(param: str, idx: int) -> str:
        v = ci.get(param)
        return str(v[idx]) if isinstance(v, list) and len(v) == 2 else "--"

    pref = "svbReal" if m.get("is_real") else "svbSyn"
    lines = [
        "% Auto-generated by jansky_research.svsbi._write_macros -- do not edit.",
        "% Synthetic (svbSyn*) and real (svbReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under svbReal* (an offline rebuild resets svbReal* to placeholders by design).",
        rf"\newcommand{{\svbSource}}{{{m['source']}}}",
        rf"\newcommand{{\svbNTargets}}{{{m['n_targets']}}}",
    ]
    for ns in ("svbSyn", "svbReal"):
        live = ns == pref
        for macro, d, key in (
            ("NDet", "n_v_detections_observed", None),
            ("Fbeam", "posterior_median", "f_beam"),
            ("Slope", "posterior_median", "lf_slope"),
            ("LogLbreak", "posterior_median", "log_Lbreak"),
            ("Calibrated", "sbc_calibrated", None),
            ("KsMax", "sbc_ks_max", None),
            ("KsCrit", "sbc_ks_crit", None),
            ("FbeamWidthRatio", "posterior_prior_width_ratio", "f_beam"),
            ("NSims", "n_sims", None),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(d, key) if live else '--'}}}")
        for macro, param, idx in (
            ("FbeamLo", "f_beam", 0),
            ("FbeamHi", "f_beam", 1),
            ("LogLbreakLo", "log_Lbreak", 0),
            ("LogLbreakHi", "log_Lbreak", 1),
            ("SlopeLo", "lf_slope", 0),
            ("SlopeHi", "lf_slope", 1),
        ):
            lines.append(
                rf"\newcommand{{\{ns}{macro}}}{{{ci_bound(param, idx) if live else '--'}}}"
            )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="SBI population inference for RACS Stokes-V emitters.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--device", default="cpu")
    p.add_argument("--n-sims", type=int, default=3000)
    args = p.parse_args(argv)
    print(
        json.dumps(
            run(args.out, offline=args.offline, device=args.device, n_sims=args.n_sims), indent=2
        )
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
