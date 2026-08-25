"""Tests for jansky_research.vlbi — multi-decade VLBI flux variability. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import vlbi


def test_lightcurve_metrics_steady_vs_variable():
    # one steady source (constant) and one with a clear single-epoch flare
    flux = np.array(
        [
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 5.0, 1.0, 1.0],
        ]
    )
    err = np.full_like(flux, 0.05)
    eta, v, pval, nep, mean, md = vlbi.lightcurve_metrics(flux, err)
    assert nep.tolist() == [5, 5]
    assert eta[1] > eta[0]  # the flaring source is far more significant
    assert v[1] > v[0]  # and higher amplitude
    assert pval[0] > 0.01 and pval[1] < 1e-3  # steady consistent with constant; variable not


def test_lightcurve_metrics_too_few_epochs():
    flux = np.array([[1.0, np.nan, np.nan, 1.0]])  # only 2 finite < MIN_EPOCHS
    err = np.array([[0.1, np.nan, np.nan, 0.1]])
    eta, v, pval, nep, mean, md = vlbi.lightcurve_metrics(flux, err)
    assert nep[0] == 2
    assert np.isnan(eta[0]) and np.isnan(v[0])


def test_sx_index_flat_and_steep():
    # source 0: equal S and X flux -> alpha = 0; source 1: X brighter -> inverted (alpha > 0)
    flux_s = np.array([[1.0, 1.0], [1.0, 1.0]])
    flux_x = np.array([[1.0, 1.0], [2.0, 2.0]])
    err = np.full_like(flux_s, 0.05)
    alpha, aerr = vlbi.sx_index(flux_s, flux_x, err, err)
    assert abs(alpha[0]) < 1e-9
    assert alpha[1] > 0  # rising to X band
    assert np.all(np.isfinite(aerr))


def test_select_variable_excludes_short_curves():
    # a realistic steady population (the sigma-clip needs a real distribution to set the threshold)
    rng = np.random.default_rng(0)
    n = 200
    eta = np.abs(rng.normal(1.0, 0.3, n))
    v = np.abs(rng.normal(0.05, 0.01, n))
    nep = np.full(n, 10)
    # two clear, well-sampled outliers
    eta[0], v[0] = 80.0, 0.6
    eta[1], v[1] = 120.0, 0.9
    # one outlier just as extreme but with too few epochs -> must be excluded
    eta[2], v[2], nep[2] = 200.0, 1.0, 2
    mask, eta_thr, v_thr = vlbi.select_variable(eta, v, nep)
    assert mask[0] and mask[1]  # the well-sampled outliers are flagged
    assert not mask[2]  # the short curve is excluded despite being the most extreme
    assert np.isfinite(eta_thr) and np.isfinite(v_thr)


def test_variability_floor_from_controls():
    # controls (steady) have low V; their median sets the floor; non-controls above it are flagged
    v = np.array([0.20, 0.18, 0.22, 0.55, 0.40, 0.10])
    nep = np.array([10, 10, 10, 10, 10, 10])
    is_control = np.array([True, True, True, False, False, False])
    floor, above = vlbi.variability_floor(v, nep, is_control)
    assert abs(floor - 0.20) < 1e-9  # median of the three controls
    assert above.tolist() == [False, False, False, True, True, False]  # 0.55/0.40 above, 0.10 below
    # controls themselves are never flagged, and short curves are excluded
    nep2 = np.array([10, 10, 10, 2, 10, 10])
    _, above2 = vlbi.variability_floor(v, nep2, is_control)
    assert not above2[3]  # too few epochs


def test_variability_floor_no_controls_returns_nan():
    v = np.array([0.2, 0.5, 0.3])
    floor, above = vlbi.variability_floor(
        v, np.array([10, 10, 10]), np.array([False, False, False])
    )
    assert np.isnan(floor) and not above.any()


def test_synthetic_population_recovers_injected_variables():
    pop = vlbi.synthetic_lightcurves(n_sources=500, seed=1)
    assert pop["flux_x"].shape[0] == 500
    eta, v, pval, nep, mean, md = vlbi.lightcurve_metrics(pop["flux_x"], pop["err_x"])
    mask, _, _ = vlbi.select_variable(eta, v, nep)
    truth = pop["is_variable"]
    # most injected variables (that have enough epochs) are recovered, with high purity
    completeness = (mask & truth).sum() / truth.sum()
    purity = (mask & truth).sum() / max(mask.sum(), 1)
    assert completeness > 0.6
    assert purity > 0.6


def test_floor_diagnostics_known_values():
    # controls 0.18/0.19/0.20/0.25 -> floor = 0.195; drop-one medians are 0.19, 0.20, 0.19, 0.19
    v = np.array([0.18, 0.19, 0.20, 0.25, 0.30, 0.196, 0.10])
    nep = np.full(7, 10)
    ctrl = np.array([True, True, True, True, False, False, False])
    d = vlbi.floor_diagnostics(v, nep, ctrl)
    assert abs(d["floor"] - 0.195) < 1e-9
    assert d["n_above"] == 2  # 0.30 and 0.196
    # a median threshold sits above half its own controls by construction
    assert d["n_controls_above_floor"] == 2
    assert len(d["jackknife"]) == 4
    assert d["jk_floor_lo"] == 0.19 and d["jk_floor_hi"] == 0.2
    # dropping a LOW control raises the floor and can remove a marginal source
    assert d["jk_n_above_lo"] == 1 and d["jk_n_above_hi"] == 2
    # the most conservative floor the controls supply: only 0.30 clears 0.25
    assert d["n_above_at_control_max"] == 1


def test_epoch_confound_detects_injected_trend():
    # V rises with log n by construction; controls follow the same trend at a lower level
    rng = np.random.default_rng(0)
    nep = np.array([4, 5, 8, 10, 15, 20, 30, 50, 80, 120] * 2)
    ctrl = np.zeros(20, dtype=bool)
    ctrl[:4] = True
    v = 0.05 + 0.15 * np.log10(nep) + rng.normal(0, 0.01, 20)
    v[ctrl] -= 0.04
    ec = vlbi.epoch_confound(v, nep, ctrl, n_perm=500)
    assert ec["spearman_rho"] > 0.8
    assert ec["spearman_p_perm"] < 0.01
    assert ec["ols_slope"] > 0.1
    # an epoch-matched floor exists and counts non-controls against their own-n threshold
    assert "n_above_epoch_matched" in ec
    assert 0 <= ec["n_above_epoch_matched"] <= 16
    # with no confound the correlation must not fire
    flat = 0.2 + rng.normal(0, 0.01, 20)
    ec0 = vlbi.epoch_confound(flat, nep, ctrl, n_perm=500)
    assert ec0["spearman_p_perm"] > 0.05


def test_v_sampling_scatter_shrinks_with_epochs():
    rng = np.random.default_rng(1)
    short = 1.0 + rng.normal(0, 0.2, (1, 6))
    long = 1.0 + rng.normal(0, 0.2, (1, 96))
    f = np.full((2, 96), np.nan)
    f[0, :6] = short
    f[1] = long
    e = np.where(np.isfinite(f), 0.05, np.nan)
    se = vlbi.v_sampling_scatter(f, e, n_boot=500)
    assert np.isfinite(se).all()
    assert se[0] > 2.0 * se[1]  # 16x the epochs -> much smaller sampling scatter
    # too few epochs -> nan
    f2 = f.copy()
    f2[0, 3:] = np.nan
    se2 = vlbi.v_sampling_scatter(f2, np.where(np.isfinite(f2), 0.05, np.nan), n_boot=100)
    assert np.isnan(se2[0])


def test_synthetic_controls_are_steady_and_flagged():
    pop = vlbi.synthetic_lightcurves(n_sources=300, n_controls=4, seed=2)
    assert pop["is_control"].sum() == 4
    assert not (pop["is_control"] & pop["is_variable"]).any()
    # the default fixture (n_controls=0) is unchanged
    pop0 = vlbi.synthetic_lightcurves(n_sources=300, seed=2)
    assert pop0["is_control"].sum() == 0
    assert np.array_equal(pop0["is_variable"], pop["is_variable"])


def test_floor_fixture_metrics_and_amplitude_ladder():
    strong = vlbi.floor_fixture_metrics(var_amp=2.0)
    assert strong["completeness"] > 0.8
    # the selection function, measured: a median floor admits ~half of any steady population, so
    # on a steady-dominated blind population the floor selector's purity is LOW. It is a
    # consistency-with-steadiness threshold, not a blind-survey classifier.
    assert strong["purity"] < 0.5
    assert strong["n_selected"] > 2 * strong["n_injected"]
    # completeness must degrade as the injected amplitude approaches the floor
    weak = vlbi.floor_fixture_metrics(var_amp=0.25)
    assert weak["completeness"] < strong["completeness"]


def test_run_offline(tmp_path):
    m = vlbi.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_sources"] == 400
    assert m["n_candidates"] >= 1
    assert m["completeness"] > 0.5 and m["purity"] > 0.5
    assert m["median_alpha_sx"] is not None
    # the control-floor path now runs offline too (4 synthetic controls)
    assert m["n_controls"] == 4
    assert m["v_floor"] > 0
    assert m["floor_diagnostics"]["jk_floor_lo"] <= m["floor_diagnostics"]["jk_floor_hi"]
    assert m["epoch_confound"]["spearman_rho"] is not None
    assert m["median_v_noncontrol"] is not None
    assert m["syn_floor_validation"]["completeness"] > 0.8
    amps = [s["var_amp"] for s in m["syn_floor_amp_sweep"]]
    assert amps == [0.25, 0.5, 1.0, 2.0]
    assert (tmp_path / "results" / "vlbi_metrics.json").exists()
    assert (tmp_path / "papers" / "vlbi" / "figures" / "etav.pdf").exists()
    macros = (tmp_path / "papers" / "vlbi" / "generated" / "macros.tex").read_text()
    assert r"\viNcand" in macros and r"\viCompleteness" in macros
    assert r"\viSynFloorCompleteness" in macros and r"\viFloorJkLo" in macros
    assert r"\viMedVnonctrl" in macros and r"\viNctrlAbove" in macros


def test_synthetic_run_does_not_clobber_real_figure(tmp_path):
    import json

    # plant a real-marked results JSON and a fake real figure
    (tmp_path / "results").mkdir(parents=True)
    (tmp_path / "results" / "vlbi_metrics.json").write_text(
        json.dumps({"source": "Astrogeo VLBI (18 sources)", "n_sources": 18})
    )
    figdir = tmp_path / "papers" / "vlbi" / "figures"
    figdir.mkdir(parents=True)
    (figdir / "etav.pdf").write_bytes(b"REAL")
    vlbi.run(out=str(tmp_path), offline=True)
    assert (figdir / "etav.pdf").read_bytes() == b"REAL"  # synthetic run left the real figure alone
    # and the real metrics were not gutted by the synthetic run
    kept = json.loads((tmp_path / "results" / "vlbi_metrics.json").read_text())
    assert kept["n_sources"] == 18
