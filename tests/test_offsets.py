"""Tests for jansky_research.offsets — radio-optical position offsets. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import offsets


def test_radio_optical_offset():
    ra_r, dec_r = np.array([10.0]), np.array([20.0])
    mas = 100.0
    # a pure-North optical offset of 100 mas
    north = offsets.radio_optical_offset(ra_r, dec_r, ra_r, dec_r + mas / offsets.DEG_TO_MAS)
    assert abs(north["ddec_mas"][0] - mas) < 1e-3 and abs(north["dra_mas"][0]) < 1e-6
    assert abs(north["offset_mas"][0] - mas) < 1e-3 and abs(north["pa_deg"][0]) < 1e-3
    # a pure-East offset (RA increases; account for cos dec)
    cosd = np.cos(np.radians(20.0))
    east = offsets.radio_optical_offset(ra_r, dec_r, ra_r + mas / offsets.DEG_TO_MAS / cosd, dec_r)
    assert abs(east["dra_mas"][0] - mas) < 1e-3 and abs(east["pa_deg"][0] - 90.0) < 1e-3


def test_normalised_offset():
    x = offsets.normalised_offset(
        np.array([3.0]), np.array([4.0]), np.array([1.0]), np.array([1.0])
    )
    assert np.isclose(x[0], 5.0)  # sqrt(3^2 + 4^2)


def test_offset_statistics():
    # a clean Gaussian-noise sample has an X>3 fraction near the Rayleigh expectation (~1.1%)
    rng = np.random.default_rng(0)
    x = np.hypot(rng.normal(0, 1, 100000), rng.normal(0, 1, 100000))
    s = offsets.offset_statistics(x)
    assert np.isclose(s["rayleigh_expectation"], np.exp(-4.5), rtol=1e-6)
    assert abs(s["frac_x_gt_cut"] - s["rayleigh_expectation"]) < 0.003  # noise ~ Rayleigh
    assert 0.5 < s["excess_ratio"] < 1.6


def test_synthetic_field_recovers_excess():
    radio, optical, truth = offsets.synthetic_field(
        n_sources=5000, structured_fraction=0.15, seed=1
    )
    assert truth.sum() >= 1
    off = offsets.radio_optical_offset(radio["ra"], radio["dec"], optical["ra"], optical["dec"])
    sig_a = np.hypot(radio["e_a"], optical["e_a"])
    sig_d = np.hypot(radio["e_d"], optical["e_d"])
    x = offsets.normalised_offset(off["dra_mas"], off["ddec_mas"], sig_a, sig_d)
    s = offsets.offset_statistics(x, off["offset_mas"])
    # injected structure makes the X>3 fraction far exceed the Rayleigh expectation
    assert s["excess_ratio"] > 3.0
    # and the tail is dominated by the injected structured sources
    assert truth[x > 3.0].mean() > 0.7


def test_jet_axis_angles_known():
    # offset along the jet -> 0; perpendicular -> 90; anti-jet -> downstream 180, axis 0
    down, axis = offsets.jet_axis_angles(
        np.array([30.0, 120.0, 210.0]), np.array([30.0, 30.0, 30.0])
    )
    assert np.allclose(down, [0.0, 90.0, 180.0])
    assert np.allclose(axis, [0.0, 90.0, 0.0])


def test_alignment_stats_recovers_injected_alignment():
    off_pa, jet_pa, x = offsets.synthetic_alignment(seed=2)
    a = offsets.alignment_stats(off_pa, jet_pa, x)
    # the full sample is non-uniform (aligned) and the significant subset is tighter
    assert a["median_axis_deg"] < 45.0  # random would be ~45
    assert a["frac_axis_lt30"] > 0.33  # random 1/3
    assert a["ks_p"] < 1e-3
    assert a["median_axis_signif_deg"] < a["median_axis_deg"]  # strengthens with significance
    # a purely random sample is consistent with uniform
    rng = np.random.default_rng(0)
    rand = offsets.alignment_stats(
        rng.uniform(0, 360, 400), rng.uniform(0, 360, 400), rng.uniform(0, 5, 400)
    )
    assert rand["ks_p"] > 0.01
    assert abs(rand["median_axis_deg"] - 45.0) < 8.0


def test_match_jets_positional():
    radio_ra = np.array([10.0, 200.0])
    radio_dec = np.array([5.0, -20.0])
    jets = {
        "ra": np.array([10.00005, 350.0]),  # first matches within ~0.2"
        "dec": np.array([5.0, 40.0]),
        "jet_pa": np.array([77.0, 12.0]),
        "delpa": np.array([10.0, 20.0]),
    }
    mask, jet_pa, _delpa = offsets.match_jets(radio_ra, radio_dec, jets, max_arcsec=1.0)
    assert mask.tolist() == [True, False]
    assert jet_pa[0] == 77.0


def test_mahalanobis_restores_the_rayleigh_null_under_correlation():
    """Correlated RA/Dec errors make the naive X's noise tail heavier than Rayleigh — the
    round-8 referee's finding that moved the published 24.0x. The Mahalanobis statistic
    whitens with the full covariance, so its pure-noise tail is Rayleigh by construction."""
    radio, optical, _ = offsets.synthetic_field(
        n_sources=60000, structured_fraction=0.0, rho=0.5, seed=3
    )
    off = offsets.radio_optical_offset(radio["ra"], radio["dec"], optical["ra"], optical["dec"])
    sig_a = np.hypot(radio["e_a"], optical["e_a"])
    sig_d = np.hypot(radio["e_d"], optical["e_d"])
    x_naive = offsets.normalised_offset(off["dra_mas"], off["ddec_mas"], sig_a, sig_d)
    x_mahal = offsets.mahalanobis_offset(off["dra_mas"], off["ddec_mas"], radio, optical)
    rayleigh = float(np.exp(-4.5))
    frac_naive = float((x_naive > 3.0).mean())
    frac_mahal = float((x_mahal > 3.0).mean())
    assert frac_naive > 1.3 * rayleigh  # the wrong null: visibly heavy on pure noise
    assert abs(frac_mahal - rayleigh) < 0.35 * rayleigh  # the right null: Rayleigh restored


def test_noise_core_scale_recovers_an_injected_inflation():
    rng = np.random.default_rng(4)
    x = np.hypot(rng.normal(0, 1.3, 40000), rng.normal(0, 1.3, 40000))  # errors 30% too small
    fit = offsets.noise_core_scale(x, x_max=1.5)
    assert abs(fit["scale"] - 1.3) < 0.08
    x1 = np.hypot(rng.normal(0, 1.0, 40000), rng.normal(0, 1.0, 40000))
    fit1 = offsets.noise_core_scale(x1, x_max=1.5)
    assert abs(fit1["scale"] - 1.0) < 0.05


def test_inflation_sweep_monotone():
    rng = np.random.default_rng(5)
    x = np.hypot(rng.normal(0, 1, 20000), rng.normal(0, 1, 20000))
    x[:2000] += 8.0  # a heavy structural tail
    sweep = offsets.inflation_sweep(x)
    ex = [sweep[k]["excess"] for k in ("1.2", "1.5", "2.0", "3.0")]
    assert all(a >= b for a, b in zip(ex, ex[1:], strict=False))
    assert ex[-1] > 1.0  # a genuinely heavy tail survives 3x inflation


def test_directional_split_uses_the_sign_symmetric_null():
    off_pa, jet_pa, _x = offsets.synthetic_alignment(
        n=2000, aligned_fraction=1.0, downstream_fraction=0.7, seed=6
    )
    d = offsets.directional_split(off_pa, jet_pa)
    assert d["n_axis_aligned"] > 1500
    assert d["n_downstream"] > d["n_upstream"]  # the injected downstream majority
    assert d["binom_p_two_sided"] < 1e-6
    # a sign-symmetric axis-aligned population must NOT fire (this is what "random 25%"
    # got wrong: it tested alignment again, not direction)
    off_sym, jet_sym, _ = offsets.synthetic_alignment(
        n=2000, aligned_fraction=1.0, downstream_fraction=0.5, seed=7
    )
    d_sym = offsets.directional_split(off_sym, jet_sym)
    assert d_sym["binom_p_two_sided"] > 0.01


def test_run_offline(tmp_path):
    m = offsets.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n"] > 1000
    assert m["excess_ratio"] > 1.0  # reproduces the offset excess tail
    assert m["frac_struct_in_tail"] > 0.6
    # the jet-alignment analysis is present and recovers the injected alignment
    assert m["n_jet"] > 100
    assert m["median_axis_deg"] < 45.0
    assert m["frac_axis_signif"] > m["frac_axis_lt30"]  # tighter in the significant subset
    # the round-8 additions: measured error model, correct directional null, sweeps
    assert m["x_statistic"].startswith("mahalanobis")
    assert 0.9 < m["noise_core_scale"] < 1.2  # fixture errors are exact by construction
    assert "1.5" in m["inflation_sweep"] and "3.0" in m["inflation_sweep"]
    assert m["n_downstream"] + m["n_upstream"] == m["n_axis_aligned"]
    assert m["downstream_binom_p"] < 0.05  # the fixture's 80% downstream majority
    assert len(m["frac_x_gt3_by_err_quartile"]) == 4
    assert m["frac_x_gt3_dec_jackknife_se_pct"] >= 0.0
    assert (tmp_path / "results" / "offsets_metrics.json").exists()
    assert (tmp_path / "papers" / "offsets" / "figures" / "xnorm.pdf").exists()
    macros = (tmp_path / "papers" / "offsets" / "generated" / "macros.tex").read_text()
    assert r"\offFracTail" in macros and r"\offExcess" in macros
    assert r"\offJetMedAxis" in macros and r"\offJetKsP" in macros
    assert r"\offCoreScale" in macros and r"\offNdown" in macros and r"\offInflTwo" in macros
