"""Tests for jansky_research.rmstructure -- RM structure functions by latitude. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import rmstructure as rms


def test_structure_function_debiases_pure_noise():
    # a pure-noise RM sky must debias to SF ~ 0 (within errors), not 2*sigma^2
    rng = np.random.default_rng(0)
    n = 400
    ra, dec = rng.uniform(0, 20, n), rng.uniform(-10, 10, n)
    noise = 3.0
    rm = rng.normal(0, noise, n)
    out = rms.structure_function(ra, dec, rm, np.full(n, noise), n_boot=40)
    good = np.isfinite(out["sf"])
    assert good.sum() >= 4
    # undebiased would be ~2*9=18; debiased must hug zero
    assert np.nanmedian(np.abs(out["sf"][good])) < 4.0


def test_structure_function_recovers_coherence_scale():
    s = rms.synthetic_rm_screen(seed=1)
    hi = np.abs(s["gal_b"]) > 10.0
    out = rms.structure_function(s["ra"][hi], s["dec"][hi], s["rm"][hi], s["rm_err"][hi], n_boot=30)
    brk = rms._sf_break(out["sep_deg"], out["sf"])
    # Gaussian-blob ACF: half-plateau crossing at ~1.7 sigma; allow a factor-2 window
    assert 1.0 * s["coherence_deg"] < brk < 3.0 * s["coherence_deg"]
    # SF rises: small-sep bins well below the plateau
    good = np.isfinite(out["sf"])
    assert out["sf"][good][0] < 0.5 * np.nanmedian(out["sf"][good][-3:])


def test_plane_enhancement_shows_in_sf_amplitude():
    s = rms.synthetic_rm_screen(seed=2)
    lo, hi = np.abs(s["gal_b"]) < 10.0, np.abs(s["gal_b"]) > 10.0
    sf_lo = rms.structure_function(
        s["ra"][lo], s["dec"][lo], s["rm"][lo], s["rm_err"][lo], n_boot=20
    )
    sf_hi = rms.structure_function(
        s["ra"][hi], s["dec"][hi], s["rm"][hi], s["rm_err"][hi], n_boot=20
    )
    assert np.nanmedian(sf_lo["sf"][-3:]) > 4.0 * np.nanmedian(sf_hi["sf"][-3:])


def test_pair_subsampling_records_fraction():
    s = rms.synthetic_rm_screen(n_sources=300, seed=3)
    out = rms.structure_function(s["ra"], s["dec"], s["rm"], s["rm_err"], max_pairs=5000, n_boot=10)
    assert 0.0 < out["pair_fraction"] < 1.0
    assert np.isfinite(out["sf"]).any()


def test_run_offline_writes_artifacts(tmp_path):
    m = rms.run(str(tmp_path), offline=True)
    assert m["source"] == "synthetic RM screen"
    # Assert on the ENSEMBLE, not the single default seed: 17 of the 30 committed seed ratios
    # fall outside the old (3, 7) window, so a per-seed assertion only passed because seed 0
    # sits high -- a test that locks the outlier in (the svbSynNTargets lesson shape).
    assert 2.0 < m["enhancement_ratio_ensemble"] < 5.0  # band median of an injected 5x boost
    assert 0.3 < m["enhancement_ratio_ensemble_sd"] < 3.0
    assert (
        abs(m["enhancement_ratio"] - m["enhancement_ratio_ensemble"])
        < 3.0 * m["enhancement_ratio_ensemble_sd"]
    )
    assert m["sf_plateau_low_b"] > m["sf_plateau_high_b"]
    assert m["sf_plateau_ratio"] > 4.0
    assert 1.0 < m["sf_break_high_b_deg"] < 6.0
    assert m["true_coherence_deg"] == 2.0
    assert isinstance(m["jackknife_se_by_block_deg"], dict)
    assert m["n_jackknife_blocks_effective"] <= m["n_jackknife_blocks"]
    saved = json.loads((tmp_path / "results" / "rmstructure_metrics.json").read_text())
    assert saved == m
    fig = tmp_path / "papers" / "rmstructure" / "figures" / "rmstructure_syn.pdf"
    assert fig.stat().st_size > 0
    # per-leg figure name: an offline run must never overwrite the real leg's figure
    assert not (tmp_path / "papers" / "rmstructure" / "figures" / "rmstructure_real.pdf").exists()
    macros = (tmp_path / "papers" / "rmstructure" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\rmsSynRatio}" in macros and r"\newcommand{\rmsRealRatio}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    rms._write_macros({"source": "x", "sf_break_low_b_deg": None}, p)
    assert r"\newcommand{\rmsSynBreakLo}{--}" in p.read_text()


def test_latitude_ladder_recovers_profile_shape():
    s = rms.synthetic_rm_screen(n_sources=2500, seed=4)
    lad = rms.latitude_ladder(s, b_edges=(0.0, 5.0, 10.0, 20.0), max_pairs=100_000, n_boot=10)
    fin = np.isfinite(lad["sigma_rm"])
    assert fin.sum() >= 2
    # injected plane boost -> monotone-decreasing sigma_RM with |b|
    vals = lad["sigma_rm"][fin]
    assert vals[0] > vals[-1]
    # floor subtraction: sigma_gal <= sigma_rm everywhere, and 0/NaN at the floor bin
    assert np.all(lad["sigma_gal"][fin] <= lad["sigma_rm"][fin] + 1e-9)
    assert lad["floor_sigma"] > 0


def test_latitude_ladder_thin_bins_are_nan():
    s = rms.synthetic_rm_screen(n_sources=300, seed=5)
    lad = rms.latitude_ladder(s, b_edges=(0.0, 1.0, 90.0), max_pairs=50_000, n_boot=5)
    assert np.isnan(lad["sigma_rm"][0])  # ~no sources in |b|<1 -> honest NaN


def test_spatial_block_jackknife_beats_iid_bootstrap_on_correlated_field():
    # On a blob-correlated screen the i.i.d. bootstrap resamples within patches; the block
    # jackknife treats a block as the exchangeable unit and must report a LARGER error.
    s = rms.synthetic_rm_screen(n_sources=2000, seed=7)
    from jansky_research.rmsky import _ratio_bootstrap_se, enhancement_ratio

    jk = rms.spatial_block_jackknife(
        s["gal_l"],
        s["gal_b"],
        lambda m: enhancement_ratio(s["rm"][m], s["gal_b"][m], pole_deg=15.0),
        block_deg=10.0,
    )
    boot = _ratio_bootstrap_se(s["rm"], s["gal_b"], pole_deg=15.0)
    assert np.isfinite(jk["se"]) and jk["se"] > 0
    assert jk["n_blocks"] >= 5
    assert np.isclose(jk["stat"], enhancement_ratio(s["rm"], s["gal_b"], pole_deg=15.0), atol=1e-12)
    assert jk["se"] > boot  # the whole point of the fix


def test_spatial_block_jackknife_too_few_blocks_is_nan():
    rng = np.random.default_rng(0)
    gl = rng.uniform(0.0, 5.0, 100)  # everything inside one 10-deg block
    gb = rng.uniform(0.0, 5.0, 100)
    jk = rms.spatial_block_jackknife(gl, gb, lambda m: float(m.sum()), block_deg=10.0)
    assert np.isnan(jk["se"])


def test_latitude_ladder_commits_uncertainties_and_pair_accounting():
    s = rms.synthetic_rm_screen(n_sources=2500, seed=4)
    lad = rms.latitude_ladder(s, b_edges=(0.0, 5.0, 10.0, 20.0), max_pairs=100_000, n_boot=10)
    fin = np.isfinite(lad["sigma_rm"])
    assert np.all(np.isfinite(lad["plateau_err"][fin]))
    assert np.all(np.isfinite(lad["sigma_rm_err"][fin]))
    assert np.all(lad["n_pairs"][fin] > 0)
    assert np.all(lad["pair_fraction"][fin] > 0)
    # the floor's own uncertainty and the alternative-floor sensitivity are quoted
    assert lad["floor_sigma_err"] > 0
    assert np.isfinite(lad["sigma_gal_plane"])
    assert lad["sigma_gal_plane_floor_lo"] <= lad["sigma_gal_plane_floor_hi"] or np.isnan(
        lad["sigma_gal_plane_floor_lo"]
    )


def test_run_offline_commits_seed_ensemble(tmp_path):
    m = rms.run(out=str(tmp_path), offline=True)
    import json

    ens = json.loads((tmp_path / "results" / "rmstructure_synthetic.json").read_text())
    assert len(ens["ratios"]) == rms.N_SYNTHETIC_REALIZATIONS
    assert ens["injected_plane_boost"] == 5.0
    assert "synthetic" in ens["source"]
    assert np.isclose(np.mean(ens["ratios"]), m["enhancement_ratio_ensemble"], atol=0.01)
    assert m["enhancement_ratio_ensemble_sem"] is not None
    macros = (tmp_path / "papers" / "rmstructure" / "generated" / "macros.tex").read_text()
    # per-namespace provenance, so a rebuild of one leg cannot invert the other's marker
    assert r"\newcommand{\rmsSynSource}" in macros
    assert r"\newcommand{\rmsRealSource}{--}" in macros
    assert r"\rmsRealRatioJkSe" in macros and r"\rmsSynRatioEnsSem" in macros
    assert r"\rmsRealNAfterGoodrm" in macros
    assert r"\rmsRealBreakHiFlagM" in macros and r"\rmsRealBreakFracSmall" in macros
    assert r"\rmsSynPlatRatio" in macros


def test_matched_flag_sf_identity_when_all_flagged():
    # with every source flagged good, the two curves share the exact same pairs
    s = rms.synthetic_rm_screen(n_sources=500, seed=8)
    out = rms.matched_flag_sf(
        s["ra"], s["dec"], s["rm"], s["rm_err"], np.ones(s["rm"].size, bool), n_boot=20
    )
    f, u = out["flagged"], out["unflagged"]
    assert np.allclose(f["sf"], u["sf"], equal_nan=True)
    assert np.array_equal(f["n_pairs"], u["n_pairs"])
    assert f["break_deg"] == u["break_deg"]
    assert out["n_flagged"] == out["n_total"]


def test_matched_flag_sf_detects_flagged_contamination():
    # contaminate 25% of sources with large independent noise NOT recorded in rm_err (the
    # leakage failure mode); flagging them out must lower the small-separation SF, and the
    # flagged curve must use only pairs of clean sources
    s = rms.synthetic_rm_screen(n_sources=800, seed=9)
    rng = np.random.default_rng(10)
    n = s["rm"].size
    bad = rng.random(n) < 0.25
    rm = s["rm"].copy()
    rm[bad] += rng.normal(0.0, 60.0, int(bad.sum()))
    out = rms.matched_flag_sf(s["ra"], s["dec"], rm, s["rm_err"], ~bad, n_boot=40)
    f, u = out["flagged"], out["unflagged"]
    assert np.all(f["n_pairs"] <= u["n_pairs"])
    g = np.isfinite(f["sf"]) & np.isfinite(u["sf"])
    # unrecorded white noise raises the unflagged SF at (at least) the smallest separations
    assert u["sf"][g][0] > f["sf"][g][0]
    # bootstrap accounting is committed, not asserted from a single draw
    assert out["n_boot_finite"] > 0
    assert 0.0 <= out["frac_boot_unflagged_break_smaller"] <= 1.0
    for leg in (f, u):
        assert np.isfinite(leg["break_deg"])
        assert leg["break_p16"] <= leg["break_p84"]


def test_structure_function_block_jackknife_plateau():
    s = rms.synthetic_rm_screen(n_sources=1200, seed=11)
    blocks = np.floor(s["gal_l"] / 10.0) * 1000 + np.floor((s["gal_b"] + 90.0) / 10.0)
    out = rms.structure_function(
        s["ra"], s["dec"], s["rm"], s["rm_err"], n_boot=10, block_ids=blocks
    )
    assert out["plateau_jk_blocks"] >= 5
    assert np.isfinite(out["plateau_jk_se"]) and out["plateau_jk_se"] > 0


def test_latitude_ladder_carries_jackknife_errors():
    s = rms.synthetic_rm_screen(n_sources=2500, seed=4)
    lad = rms.latitude_ladder(s, b_edges=(0.0, 10.0, 20.0), max_pairs=100_000, n_boot=10)
    fin = np.isfinite(lad["sigma_rm"])
    assert np.all(np.isfinite(lad["plateau_jk_err"][fin]))
    assert np.all(np.isfinite(lad["sigma_rm_jk_err"][fin]))
    assert np.all(lad["jk_blocks"][fin] >= 5)
    assert np.isfinite(lad["floor_sigma_jk_err"])


def test_spatial_block_jackknife_effective_blocks():
    # a statistic supported only at |b| > 15 is untouched by low-|b| blocks: the effective
    # count must exclude them (the "601 blocks, ~330 contribute zero" honesty item)
    rng = np.random.default_rng(3)
    n = 2000
    gl = rng.uniform(0.0, 40.0, n)
    gb = rng.uniform(-20.0, 20.0, n)
    vals = rng.normal(0.0, 1.0, n)
    pole = np.abs(gb) > 15.0

    def stat(m):
        keep = m & pole
        return float(np.median(vals[keep])) if keep.any() else float("nan")

    jk = rms.spatial_block_jackknife(gl, gb, stat, block_deg=10.0)
    assert jk["n_blocks_effective"] < jk["n_blocks"]
    assert jk["n_blocks_effective"] > 0
