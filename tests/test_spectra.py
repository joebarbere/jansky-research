"""Tests for jansky_research.spectra — spectral index, classify, cross-match, USS. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import spectra


def test_spectral_index_recovers_injected_alpha():
    nu_lo, nu_hi = 147.5, 1400.0
    alpha_true = -1.4
    s_lo = np.array([100.0])
    s_hi = s_lo * (nu_hi / nu_lo) ** alpha_true
    alpha, e = spectra.spectral_index(s_lo, nu_lo, s_hi, nu_hi, 0.1 * s_lo, 0.05 * s_hi)
    assert np.isclose(alpha[0], alpha_true, atol=1e-6)
    assert e[0] > 0 and np.isfinite(e[0])


def test_classify_thresholds():
    assert spectra.classify(-1.6) == "uss"
    assert spectra.classify(-0.8) == "steep"
    assert spectra.classify(-0.2) == "flat"
    assert spectra.classify(0.4) == "inverted"
    assert spectra.classify(spectra.USS_THRESHOLD - 0.01) == "uss"


def test_crossmatch_respects_radius():
    ra_lo = np.array([180.0, 181.0])
    dec_lo = np.array([30.0, 30.0])
    # one counterpart 3" away, one 60" away (outside a 15" radius)
    ra_hi = np.array([180.0 + 3 / 3600.0, 181.0 + 60 / 3600.0])
    dec_hi = np.array([30.0, 30.0])
    i_lo, i_hi, sep = spectra.crossmatch(ra_lo, dec_lo, ra_hi, dec_hi, radius_arcsec=15.0)
    assert i_lo.tolist() == [0]  # only the close pair survives
    assert sep[0] < 15.0


def test_find_uss_recovers_injected_population():
    low, high = spectra.synthetic_field(n=300, f_uss=0.05, f_inverted=0.05, seed=1)
    res = spectra.find_uss(low, high)
    assert res["alpha"].size > 250  # nearly all cross-match (small jitter)
    n_uss = int(res["is_uss"].sum())
    assert 8 <= n_uss <= 25  # ~15 injected USS recovered
    # the USS sources really are steep
    assert np.median(res["alpha"][res["is_uss"]]) < spectra.USS_THRESHOLD
    assert set(np.unique(res["cls"])) <= {"uss", "steep", "flat", "inverted"}


def test_synthetic_field_shapes():
    low, high = spectra.synthetic_field(n=50, seed=0)
    for d in (low, high):
        assert {"ra", "dec", "flux", "eflux"} <= set(d)
        assert d["ra"].size == 50


def test_analyze_and_run_offline(tmp_path):
    low, high = spectra.synthetic_field(n=300, seed=2)
    m = spectra.analyze(spectra.find_uss(low, high), source="synthetic")
    assert m["n_matched"] > 250 and m["n_uss"] >= 1
    assert m["alpha_min"] < m["alpha_median"]

    metrics = spectra.run(out=str(tmp_path), offline=True)
    assert metrics["source"] == "synthetic"
    assert (tmp_path / "results" / "uss_metrics.json").exists()
    assert (tmp_path / "results" / "uss_candidates.csv").exists()
    figs = {p.name for p in (tmp_path / "papers" / "spectra" / "figures").glob("*.pdf")}
    assert {"alpha_hist.pdf", "alpha_vs_flux.pdf"} <= figs
    macros = (tmp_path / "papers" / "spectra" / "generated" / "macros.tex").read_text()
    # the fixture's boundary measurement is computed and emitted in every mode
    assert metrics["syn_cut_purity"] is not None and metrics["syn_cut_completeness"] is not None
    assert metrics["syn_selection_bias_mc"]["predicted_offset"] < 0  # selection biases steep
    assert r"\usSynCutPurity" in macros and r"\usRealDeltaAlpha" in macros
    assert metrics["chance_matches_expected"] > 0
    assert metrics["matched_sensitivity"]["alpha_match"] < 0


def test_reference_crossmatch_aligns_and_flags_misses():
    res = {"ra": np.array([180.0, 181.0]), "dec": np.array([30.0, 30.0])}
    ref = {
        "ra": np.array([180.0 + 2 / 3600.0]),
        "dec": np.array([30.0]),
        "spindex": np.array([-1.5]),
        "e_spindex": np.array([0.05]),
        "scode": np.array(["S"]),
    }
    rc = spectra.reference_crossmatch(res, ref)
    assert np.isclose(rc["ref_alpha"][0], -1.5) and rc["ref_scode"][0] == "S"
    assert np.isnan(rc["ref_alpha"][1]) and rc["ref_scode"][1] == ""


def test_population_offset_excludes_limits_and_splits_dec():
    n = 400
    rng = np.random.default_rng(3)
    ref_alpha = rng.normal(-0.8, 0.15, n)
    alpha = ref_alpha + rng.normal(0, 0.05, n)  # no injected offset
    e_alpha = np.full(n, 0.05)
    scode = np.array(["S"] * n, dtype=object)
    scode[:50] = "L"
    ref_alpha_l = ref_alpha.copy()
    ref_alpha_l[:50] -= 5.0  # absurd limit rows that would wreck the mean if included
    rc = {"ref_alpha": ref_alpha_l, "ref_e": np.full(n, 0.03), "ref_scode": scode}
    dec = np.concatenate([np.full(200, 29.5), np.full(200, 30.5)])
    po = spectra.population_offset(alpha, e_alpha, rc, dec)
    assert po["n_pairs"] == 350  # the 50 limit rows are out
    assert abs(po["mean_offset"]) < 3 * po["se_offset"]  # consistent with zero
    assert po["dec_le_n"] + po["dec_gt_n"] == 350
    # scatter should be near the combined formal errors, not wildly above
    assert 0.5 < po["scatter"] / po["expected_scatter"] < 2.0


def test_selection_bias_mc_predicts_a_negative_offset():
    rng = np.random.default_rng(4)
    truth = rng.normal(-0.8, 0.2, 2000)
    e = np.full(2000, 0.15)
    sb = spectra.selection_bias_mc(truth, e, n_mc=300)
    # sources selected on noisy alpha < -1.3 are preferentially scattered-down: offset < 0
    assert sb["predicted_offset"] < -0.05
    assert sb["predicted_offset_sd"] > 0
    # with tiny errors the selection bias vanishes
    sb0 = spectra.selection_bias_mc(truth, np.full(2000, 1e-4), n_mc=100)
    assert sb0["n_realizations"] == 0 or abs(sb0["predicted_offset"]) < 0.01


def test_matched_sensitivity_arithmetic():
    ms = spectra.matched_sensitivity(alpha_probe=-1.65)
    # equally-sensitive index from the documented floors: ln(2.5/24.5)/ln(1400/147.5)
    assert np.isclose(ms["alpha_match"], np.log(2.5 / 24.5) / np.log(1400 / 147.5), atol=1e-3)
    # at the probe index the required 150 MHz flux is ~100 mJy (the truncation the paper cites)
    assert 80 < ms["s150_needed_at_probe_mjy"] < 130
    # steeper probe -> more required flux
    assert ms["s150_needed_at_probe_mjy"] > ms["s150_needed_at_threshold_mjy"]


def test_uss_confusion_scores_both_directions():
    # our cut flags 2: one real USS, one noise-scattered impostor; ref holds 3 USS, we recover 1
    res = {
        "ra": np.array([10.0, 10.01, 10.02, 10.03]),
        "dec": np.array([0.0, 0.0, 0.0, 0.0]),
        "alpha": np.array([-1.5, -1.4, -0.8, -0.7]),
        "is_uss": np.array([True, True, False, False]),
    }
    ref = {
        "ra": np.array([10.0, 10.01, 10.05, 10.06]),
        "dec": np.array([0.0, 0.0, 0.0, 0.0]),
        "spindex": np.array([-1.6, -0.9, -1.5, -1.7]),
        "e_spindex": np.array([0.05, 0.05, 0.05, 0.05]),
        "scode": np.array(["S", "S", "S", "S"]),
    }
    cf = spectra.uss_confusion(res, ref)
    assert cf["n_flagged"] == 2 and cf["n_flagged_ref_uss"] == 1
    assert cf["purity"] == 0.5
    assert cf["n_ref_uss"] == 3 and cf["n_ref_uss_recovered"] == 1
    assert np.isclose(cf["completeness"], 0.333, atol=0.001)


def test_chance_matches_arithmetic():
    # 100 targets x 1000 background over 10 deg^2 with a 15" disc
    exp = spectra.chance_matches(100, 1000, 10.0)
    assert np.isclose(exp, 100 * 1000 * np.pi * (15 / 3600) ** 2 / 10.0)


def test_synthetic_field_injects_at_the_threshold():
    low, _high = spectra.synthetic_field(n=1000, f_uss=0.1, seed=5)
    at = low["alpha_true"]
    uss_true = at[at < spectra.USS_THRESHOLD]
    # both the deep and the near-threshold populations exist
    assert (uss_true < -1.5).any() and (uss_true > -1.45).any()


def test_synthetic_run_does_not_clobber_real_artifacts(tmp_path):
    import json

    (tmp_path / "results").mkdir(parents=True)
    (tmp_path / "results" / "uss_metrics.json").write_text(
        json.dumps({"source": "tgss-x-nvss @ (180.00, +30.00) r=2 deg", "n_matched": 456})
    )
    (tmp_path / "results" / "uss_candidates.csv").write_text("REAL")
    figdir = tmp_path / "papers" / "spectra" / "figures"
    figdir.mkdir(parents=True)
    (figdir / "alpha_hist.pdf").write_bytes(b"REAL")
    spectra.run(out=str(tmp_path), offline=True)
    assert (tmp_path / "results" / "uss_candidates.csv").read_text() == "REAL"
    assert (figdir / "alpha_hist.pdf").read_bytes() == b"REAL"
    kept = json.loads((tmp_path / "results" / "uss_metrics.json").read_text())
    assert kept["n_matched"] == 456
