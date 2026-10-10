"""Tests for the hiblend two-beam blending test (plan 97)."""

from __future__ import annotations

import numpy as np

from jansky_research import hiblend as h


def test_beam_response_is_half_at_half_fwhm():
    assert h.beam_response(np.array([0.0]), 2.9)[0] == 1.0
    assert np.isclose(h.beam_response(np.array([1.45]), 2.9)[0], 0.5)
    # the differential response the test relies on peaks near 2 arcmin
    s = np.linspace(0.1, 6, 600)
    d = h.beam_response(s, h.ALFA_FWHM_ARCMIN) - h.beam_response(s, h.FAST_FWHM_ARCMIN)
    assert 1.8 < s[np.argmax(d)] < 2.2


def test_crossmatch_is_mutual_and_velocity_gated():
    d = 1.0 / 60.0
    ra1, dec1, v1 = np.array([180.0, 180.0]), np.array([30.0, 31.0]), np.array([5000.0, 5000.0])
    # source 0 matches; source 1's only positional counterpart is 500 km/s away
    ra2 = np.array([180.0, 180.0, 180.0])
    dec2 = np.array([30.0 + 0.3 * d, 30.0 + 1.2 * d, 31.0])
    v2 = np.array([5020.0, 5010.0, 5500.0])
    i1, i2 = h.crossmatch(ra1, dec1, v1, ra2, dec2, v2)
    assert i1.tolist() == [0] and i2.tolist() == [0]  # nearest of two candidates, not the farther


def test_find_neighbours_excludes_self_and_flags_overlap():
    d = 1.0 / 60.0
    t = (np.array([180.0]), np.array([30.0]), np.array([5000.0]), np.array([200.0]))
    n_ra = np.array([180.0, 180.0, 180.0])
    n_dec = np.array([30.0, 30.0 + 2 * d, 30.0 + 3 * d])
    n_v = np.array([5000.0, 5100.0, 5900.0])
    nb = h.find_neighbours(
        *t, n_ra, n_dec, n_v, np.full(3, 200.0), np.ones(3), self_index=np.array([0])
    )
    assert nb.target.tolist() == [0, 0]
    assert np.allclose(np.sort(nb.sep_arcmin), [2.0, 3.0], atol=1e-6)
    assert nb.overlap.tolist() == [True, False]  # 100 km/s overlaps, 900 does not


def test_predicted_log_ratio_single_neighbour():
    nb = h.Neighbours(
        target=np.array([0]), sep_arcmin=np.array([2.0]), dv_kms=np.array([0.0]),
        flux=np.array([1.0]), overlap=np.array([True]),
    )  # fmt: skip
    r = h.predicted_log_ratio(np.array([1.0, 1.0]), nb)
    ba, bf = h.beam_response(np.array([2.0]), 3.5)[0], h.beam_response(np.array([2.0]), 2.9)[0]
    assert np.isclose(r[0], np.log10((1 + ba) / (1 + bf))) and r[1] == 0.0
    nb.overlap[:] = False
    assert h.predicted_log_ratio(np.array([1.0]), nb)[0] == 0.0
    assert h.predicted_log_ratio(np.array([1.0]), nb, force_overlap=True)[0] > 0


def test_sample_masks_partition_the_three_samples():
    nb = h.Neighbours(
        target=np.array([1, 2, 3, 3]),
        sep_arcmin=np.array([2.0, 3.0, 0.5, 2.0]),
        dv_kms=np.array([0.0, 900.0, 0.0, 0.0]),
        flux=np.ones(4),
        overlap=np.array([True, False, True, True]),
    )
    m = h.sample_masks(4, nb)
    assert m["isolated"].tolist() == [True, False, False, False]
    assert m["primary"].tolist() == [False, True, False, False]  # target 3 has a <1' pair
    assert m["null"].tolist() == [False, False, True, False]
    assert h.null_pairs(nb).tolist() == [False, True, False, False]


def test_target_flux_estimate_noise_is_uncorrelated_with_the_response():
    rng = np.random.default_rng(1)
    n = 200_000
    sa, sf = 0.12, 0.06  # unequal survey errors (dex)
    ea, ef = rng.normal(0, sa, n), rng.normal(0, sf, n)
    fa, ff = 10**ea, 10**ef  # true flux 1, no neighbours
    empty = h.Neighbours(*(np.empty(0) for _ in range(4)), overlap=np.zeros(0, bool))  # type: ignore[arg-type]
    empty.target = empty.target.astype(int)
    est = h.target_flux_estimate(fa, ff, np.full(n, sa), np.full(n, sf), empty)
    resp = np.log10(fa / ff)
    assert abs(np.corrcoef(np.log10(est), resp)[0, 1]) < 0.01
    # using one survey's flux instead would be strongly correlated with the response
    assert np.corrcoef(np.log10(ff), resp)[0, 1] < -0.3


def test_fit_beta_recovers_a_planted_slope_with_cluster_bootstrap():
    rng = np.random.default_rng(2)
    r = rng.uniform(0, 0.1, 2000)
    y = 0.01 + 0.8 * r + rng.normal(0, 0.02, 2000)
    out = h.fit_beta(y, r, np.repeat(np.arange(1000), 2), rng, n_boot=200)
    assert out["n_clusters"] == 1000 and abs(out["beta"] - 0.8) < 3 * out["beta_se"]
    assert h.fit_beta(y[:2], r[:2], np.arange(2), rng) == {"n": 2}


def test_target_clusters_links_close_targets_only():
    d = 1.0 / 60.0
    labels = h.target_clusters(
        np.array([180.0, 180.0, 181.0]), np.array([30.0, 30.0 + 3 * d, 30.0])
    )
    assert labels[0] == labels[1] != labels[2]
    assert h.target_clusters(np.array([1.0]), np.array([1.0])).tolist() == [0]


def test_analyse_recovers_planted_blending_and_a_null_without_it():
    """C1 on the fixture: blending planted at the model's strength comes back as beta ~ 1; with
    none planted beta ~ 0; the spectral null control is ~ 0 either way."""
    rng = np.random.default_rng(3)
    on = h.analyse(h.synthetic_field(3000, strength=1.0, seed=11), rng, n_boot=150)
    off = h.analyse(h.synthetic_field(3000, strength=0.0, seed=11), rng, n_boot=150)
    assert on["n_primary"] > 500 and on["n_isolated"] > 1000 and on["n_null"] > 100
    b, b0 = on["primary"], off["primary"]
    assert abs(b["beta"] - 1.0) < 3 * b["beta_se"] and b["beta_sigma"] > 3
    assert abs(b0["beta"]) < 3 * b0["beta_se"]
    assert abs(on["null"]["beta"]) < 3 * on["null"]["beta_se"]


def test_power_reports_detection_fraction():
    p = h.power(2000, strength=1.0, n_real=3, seed=5, n_boot=60)
    assert set(p) == {"n_real", "strength", "detect_frac", "beta_mean", "beta_std"}
    assert p["detect_frac"] == 1.0  # a planted beta = 1 is easy at this sample size
    assert h.power(2000, strength=0.0, n_real=3, seed=5, n_boot=60)["detect_frac"] <= 1 / 3


def _two_survey_sky(n=4000, *, strength=1.0, seed=7):
    """One galaxy population 'observed' by FAST and Arecibo beams, blending at ``strength``."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(150, 200, n)
    dec = rng.uniform(5, 25, n)
    # companions for 40% of galaxies at 1-5', mostly in the same velocity window
    k = int(0.4 * n)
    host = rng.choice(n, k, replace=False)
    sep, ang = rng.uniform(1.2, 5.0, k), rng.uniform(0, 2 * np.pi, k)
    ra = np.concatenate([ra, ra[host] + sep / 60 * np.cos(ang) / np.cos(np.radians(dec[host]))])
    dec = np.concatenate([dec, dec[host] + sep / 60 * np.sin(ang)])
    m = ra.size
    v = rng.uniform(3000, 12000, m)
    v[n:] = v[host] + rng.normal(0, 50, k)
    w50 = rng.uniform(80, 300, m)
    s = 10 ** rng.uniform(-0.3, 1.0, m)
    gal_nb = h.find_neighbours(ra, dec, v, w50, ra, dec, v, w50, s, self_index=np.arange(m))
    sa, sf = h.blended_fluxes(s, gal_nb, strength=strength)
    sig = 0.06
    fashi = {"ra": ra, "dec": dec, "cz": v, "w50": w50,
             "flux": sf * 10 ** rng.normal(0, sig, m), "flux_err": sf * sig * np.log(10)}  # fmt: skip
    jit = rng.normal(0, 0.2 / 60, (2, m))
    alfalfa = {"ra": ra + jit[0], "dec": dec + jit[1], "v": v + rng.normal(0, 10, m), "w50": w50,
               "flux": sa * 10 ** rng.normal(0, sig, m), "flux_err": sa * sig * np.log(10),
               "code": np.ones(m, int)}  # fmt: skip
    return fashi, alfalfa


def test_build_field_matches_and_counts_neighbours_once():
    fashi, alfalfa = _two_survey_sky(1500)
    f = h.build_field(fashi, alfalfa)
    assert len(f["flux_f"]) > 0.9 * len(fashi["ra"])  # nearly every galaxy matched
    assert f["n_alfalfa_only_neighbours"] < 0.1 * len(alfalfa["ra"])  # matched ones not duplicated
    assert np.all(f["sig_a_dex"] > 0) and f["nb"].target.max() < len(f["flux_f"])
    shifted = h.build_field(fashi, alfalfa, shift_dec_arcmin=h.SHIFT_ARCMIN)
    assert len(shifted["flux_f"]) < 0.01 * len(f["flux_f"])  # C4: chance matches are rare


def test_injection_field_plants_into_isolated_targets_only():
    fashi, alfalfa = _two_survey_sky(3000)
    f = h.build_field(fashi, alfalfa)
    rng = np.random.default_rng(4)
    inj = h.injection_field(f, rng, strength=1.0)
    masks = h.sample_masks(len(inj["flux_f"]), inj["nb"])
    assert masks["primary"].sum() > 100 and masks["isolated"].sum() > 100
    r1 = h.analyse(inj, rng, n_boot=150)["primary"]
    r0 = h.analyse(h.injection_field(f, rng, strength=0.0), rng, n_boot=150)["primary"]
    assert abs(r1["beta"] - 1.0) < 0.35 and abs(r0["beta"]) < 3 * r0["beta_se"]


def test_run_gated_follows_the_frozen_order():
    fashi, alfalfa = _two_survey_sky(3000, strength=1.0)
    f = h.build_field(fashi, alfalfa)
    s = h.build_field(fashi, alfalfa, shift_dec_arcmin=h.SHIFT_ARCMIN)
    out = h.run_gated(f, s, np.random.default_rng(5), n_power=4, n_inj=3, n_boot=150)
    assert list(out["gates"]) == ["C0_power", "C1_planted", "C4_match", "C2_spectral_null"]
    assert all(g["pass"] for g in out["gates"].values())
    assert out["outcome"].startswith("blending supported")
    # an underpowered sample stops at C0 and never computes beta
    small = h._subset(f, np.arange(len(f["flux_f"])) < 60, f["nb"])
    small["nb"] = h.find_neighbours(small["ra"], small["dec"], small["v"], small["w50"],
                                    small["ra"], small["dec"], small["v"], small["w50"],
                                    small["flux_f"], self_index=np.arange(60))  # fmt: skip
    stop = h.run_gated(small, s, np.random.default_rng(6), n_power=3, n_inj=2, n_boot=50)
    assert stop["stopped_at"].startswith("C0") and "primary" not in stop


def test_fit_ols_recovers_coefficients_and_honours_weights():
    rng = np.random.default_rng(8)
    x = rng.normal(size=(3000, 2))
    y = 0.1 + 0.5 * x[:, 0] - 0.3 * x[:, 1] + rng.normal(0, 0.05, 3000)
    out = h.fit_ols(y, x, np.arange(3000), rng, n_boot=100)
    assert np.allclose(out["coef"], [0.1, 0.5, -0.3], atol=0.01) and all(s > 0 for s in out["se"])
    # zero weights drop rows: fitting only the first half must match an unweighted half-fit
    w = (np.arange(3000) < 1500).astype(float)
    half = h.fit_ols(y[:1500], x[:1500], np.arange(1500), rng, n_boot=20)
    assert np.allclose(
        h.fit_ols(y, x, np.arange(3000), rng, weights=w, n_boot=20)["coef"], half["coef"]
    )
    assert h.fit_ols(y[:3], x[:3], np.arange(3), rng) == {"n": 3}


def test_diagnostics_run_and_report_every_block():
    fashi, alfalfa = _two_survey_sky(3000, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    out = h.diagnostics(
        f, np.random.default_rng(9), n_cat_ra=f["n_cat_ra"], n_cat_dec=f["n_cat_dec"], n_boot=50
    )
    assert len(out["D0_heldout_isolated_by_snr"]) == 5
    assert len(out["D1_null_by_snr_tercile"]) == 3 and "coef" in out["D1_null_snr_matched"]
    assert (
        out["D2_null_terms"]["terms"][3] == "log S_n,max" and len(out["D2_null_terms"]["coef"]) == 4
    )
    assert len(out["D3_null_density"]["coef"]) == 3 and len(out["D3_isolated_density"]["coef"]) == 2
    # no blending and no systematic planted: the held-out calibration residual is flat in S/N
    assert all(
        abs(b["sigma"]) < 3 for b in out["D0_heldout_isolated_by_snr"] if b["sigma"] is not None
    )


def test_heldout_check_fails_plan97_form_and_passes_v2_on_a_d0_shaped_hinge():
    """Plan 98's required test: a top-quintile S/N hinge (the shape D0 found) defeats the
    quadratic calibration out of sample; the decile linear spline absorbs it. With no hinge,
    both forms pass, so C3' is not failing everything."""
    rng = np.random.default_rng(1)
    hinge = h.synthetic_field(30000, pair_frac=0.1, scatter_dex=0.04, seed=3, snr_hinge_dex=0.6)
    v1 = h.heldout_calibration_check({**hinge, "calib": "v1"}, rng, n_boot=100)
    v2 = h.heldout_calibration_check({**hinge, "calib": "v2"}, rng, n_boot=100)
    assert not v1["pass"] and v2["pass"]
    assert [d["train_parity"] for d in v2["directions"]] == [0, 1]
    assert all(len(d["bins"]) == 10 for d in v2["directions"])
    flat = h.synthetic_field(30000, pair_frac=0.1, scatter_dex=0.04, seed=3)
    assert h.heldout_calibration_check({**flat, "calib": "v1"}, rng, n_boot=100)["pass"]


def test_calibration_design_v2_has_one_hinge_per_knot():
    x = np.linspace(1, 3, 50)
    d = h.calibration_design_v2(x, x, x, x, x, h.snr_knots(x))
    assert d.shape == (50, 2 + 9 + 4)
    assert np.all(d[:, 2] >= 0) and d[0, 2] == 0 and d[-1, 2] > 0


def test_run_gated_v2_runs_c3prime_first_and_checks_c2_per_tercile(monkeypatch):
    """The v2 wrapper's own logic: C3' gates everything, then run_gated (tested above) runs on
    the v2 calibration, then C2 must also pass in every null S/N tercile. run_gated is stubbed
    because C1 is not stable on a sky this small (injected blends push 16-24% of targets past
    the calibration's S/N range; survey/hiblend-findings.md step 5)."""
    fashi, alfalfa = _two_survey_sky(3000, strength=1.0)
    f = h.build_field(fashi, alfalfa)
    s = h.build_field(fashi, alfalfa, shift_dec_arcmin=h.SHIFT_ARCMIN)
    assert f["log_n15"].shape == f["flux_f"].shape
    seen: list = []

    def stub_gated(field, shifted, rng, **kw):
        seen.append((field["calib"], shifted["calib"]))
        res = h.analyse(field, rng, n_boot=50)
        null = {**res["null"], "pass": True}
        return {"gates": {"C0_power": {"pass": True}, "C1_planted": {"pass": True},
                          "C4_match": {"pass": True}, "C2_spectral_null": null},
                "C3_calibration_coef": res["calibration_coef"],
                "samples": {"n_null": res["n_null"]}, "primary": res["primary"]}  # fmt: skip

    passing = {"pass": True, "directions": [{"max_abs_median": 0.0, "p": 1.0}]}
    monkeypatch.setattr(h, "heldout_calibration_check", lambda field, rng: passing)
    monkeypatch.setattr(h, "run_gated", stub_gated)
    out = h.run_gated_v2(f, s, np.random.default_rng(5), n_boot=100)
    assert seen == [("v2", "v2")]
    assert list(out["gates"]) == [
        "C3prime_heldout", "C0_power", "C1_planted", "C4_match", "C2_spectral_null",
    ]  # fmt: skip
    c2 = out["gates"]["C2_spectral_null"]
    assert len(c2["by_snr_tercile"]) == 3 and c2["pass_overall"] is True
    assert c2["pass"] == all(h._within_2sigma(t) for t in c2["by_snr_tercile"])
    assert len(out["secondary_primary_by_snr_tercile"]) == 3 and "outcome" in out
    # a tercile that fails turns an overall C2 pass into a fail
    monkeypatch.setattr(h, "_within_2sigma", lambda b: False)
    bad = h.run_gated_v2(f, s, np.random.default_rng(5), n_boot=50)
    assert not bad["gates"]["C2_spectral_null"]["pass"] and bad["outcome"].startswith("ambiguous")
    # a control failing inside run_gated stops the run as ambiguous
    monkeypatch.setattr(h, "run_gated", lambda *a, **k: {"gates": {"C0_power": {"pass": False}},
                                                         "stopped_at": "C0: underpowered"})  # fmt: skip
    c0 = h.run_gated_v2(f, s, np.random.default_rng(5))
    assert c0["outcome"] == "ambiguous: C0: underpowered" and "primary" not in c0
    # a failed C3' stops the run before any control or beta
    monkeypatch.setattr(
        h, "heldout_calibration_check", lambda field, rng: {**passing, "pass": False}
    )
    stop = h.run_gated_v2(f, s, np.random.default_rng(5))
    assert stop["stopped_at"].startswith("C3'") and list(stop["gates"]) == ["C3prime_heldout"]
    assert "primary" not in stop


def test_write_paper_fills_every_macro_from_the_committed_results(tmp_path):
    import shutil
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "results"
    (tmp_path / "results").mkdir()
    for name in (
        "hiblend_metrics.json",
        "hiblend_diagnostics.json",
        "hiblend_v2_metrics.json",
        "hiblend_referee1.json",
        "hiblend_referee2.json",
        "hiblend_referee3.json",
    ):
        shutil.copy(root / name, tmp_path / "results" / name)
    macros, fig = h.write_paper(tmp_path)
    text = macros.read_text()
    assert fig.stat().st_size > 0 and r"\hbRealVtwoMaxB" in text
    assert "{--}" not in text  # every quoted number has a committed value
    assert abs(h.equal_flux_signal_dex() - 0.0445) < 0.001  # plan 97: "at most ~0.045 dex"
    # a perfect calibration with these bin errors passes the 0.010 dex rule about 1/3 of the time
    one = [{"bins": [{"se": 0.005}] * 10}]
    assert 0.6 < h.amplitude_pass_probability(one) ** (1 / 10) < 0.97


def test_injection_can_freeze_the_calibration_covariates():
    fashi, alfalfa = _two_survey_sky(1500, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    sub = h.injection_field(f, np.random.default_rng(3), strength=1.0, freeze_covariates=True)
    assert np.all(sub["cov_flux_f"] <= sub["flux_f"]) and np.any(sub["cov_flux_f"] < sub["flux_f"])
    assert "cov_flux_a" not in h.injection_field(f, np.random.default_rng(3), strength=1.0)


def test_referee1_checks_reproduce_their_splits_and_report_every_block():
    fashi, alfalfa = _two_survey_sky(3000, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    out = h.referee1_checks(f, n_inj=2, n_boot=50)
    assert len(out["R1_D0_plan97"]["bins"]) == 5 and len(out["R1_C3prime_v2"]) == 2
    b = out["R1_C3prime_v2"][0]["bins"][0]
    assert b["lo"] <= b["x_median"] <= b["hi"] and "mean" in b and "se_mean" in b
    assert set(out["R3_c1_covariates"]) == {
        "production_strength_1", "frozen_strength_1", "production_strength_0", "frozen_strength_0",
    }  # fmt: skip
    assert set(out["R4_strip_bootstrap"]) == {
        "cluster_6arcmin",
        "strip_2deg",
        "strip_4deg",
        "strip_8deg",
    }
    # the C3' medians R1 recomputes are the ones heldout_calibration_check reports
    c3 = h.heldout_calibration_check({**f, "calib": "v2"}, np.random.default_rng(0), n_boot=20)
    assert [x["median"] for x in out["R1_C3prime_v2"][0]["bins"]] == [
        x["median"] for x in c3["directions"][0]["bins"]
    ]


def test_injection_alfa_scale_multiplies_only_the_alfalfa_addition():
    fashi, alfalfa = _two_survey_sky(1500, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    n = len(f["flux_f"])
    zero = h.injection_field(f, np.random.default_rng(3), strength=0.0)  # same draws, no flux
    base = h.injection_field(f, np.random.default_rng(3), strength=1.0)
    two = h.injection_field(f, np.random.default_rng(3), strength=1.0, alfa_scale=np.full(n, 2.0))
    assert np.allclose(two["flux_f"], base["flux_f"])
    assert np.allclose(two["flux_a"] - zero["flux_a"], 2 * (base["flux_a"] - zero["flux_a"]))
    assert np.any(base["flux_a"] > zero["flux_a"])
    out = h.referee2_checks(f, n_inj=2)["R5_injection_on_survey_scale"]
    assert {"frozen_strength_1", "production_strength_0", "alfa_scale_isolated"} <= set(out)


def test_injection_alfa_map_replaces_the_alfalfa_increment():
    fashi, alfalfa = _two_survey_sky(1500, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    zero = h.injection_field(f, np.random.default_rng(3), strength=0.0)
    base = h.injection_field(f, np.random.default_rng(3), strength=1.0)
    tripled = h.injection_field(
        f, np.random.default_rng(3), strength=1.0, alfa_map=lambda idx, add: 3 * add
    )
    assert np.allclose(tripled["flux_a"] - zero["flux_a"], 3 * (base["flux_a"] - zero["flux_a"]))
    assert np.allclose(tripled["flux_f"], base["flux_f"])
    out = h.referee3_checks(f, n_inj=2)["R6_injection_through_curve"]
    assert {"frozen_strength_1", "effective_factor_10pct_isolated"} <= set(out)


def test_referee4_checks_report_absorption_per_u():
    fashi, alfalfa = _two_survey_sky(1500, strength=0.0)
    f = h.build_field(fashi, alfalfa)
    out = h.referee4_checks(f, n_inj=2, us=(1.0,))["R7_constant_u"]["u_1"]
    assert out["absorption"] == round(out["production"]["mean"] - out["frozen"]["mean"], 4)
