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
