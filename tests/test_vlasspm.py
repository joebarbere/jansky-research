"""Tests for the blind VLASS proper-motion search (plan 64)."""

import json

import numpy as np
import pytest

from jansky_research import vlasspm as v


def _cat(ra, dec, t, flux=None, err=None, ident=None):
    ra, dec, t = (np.atleast_1d(np.asarray(x, float)) for x in (ra, dec, t))
    n = ra.size
    return v.EpochCatalog(
        ra,
        dec,
        np.broadcast_to(t, (n,)).copy(),
        np.full(n, 3.0) if flux is None else np.asarray(flux, float),
        np.full(n, 0.2) if err is None else np.asarray(err, float),
        np.arange(n) if ident is None else np.asarray(ident),
    )


def test_tangent_offsets_wrap_and_scale():
    dx, dy = v.tangent_offsets_arcsec(359.9999, 0.0, 0.0001, 0.0)
    assert dx == pytest.approx(0.72, rel=1e-3)  # across RA=0, not -359.9998 deg
    assert dy == 0.0
    dx, _ = v.tangent_offsets_arcsec(10.0, 60.0, 10.0 + 1 / 3600, 60.0)
    assert dx == pytest.approx(0.5, rel=1e-3)  # cos(60 deg)


def test_orphan_masks():
    a = _cat([10.0, 20.0], [0.0, 0.0], 2018.0)
    b = _cat([10.0 + 1 / 3600, 30.0], [0.0, 0.0], 2021.0)
    oa, ob = v.orphan_masks(a, b)
    assert oa.tolist() == [False, True]  # 1" apart = static; 20 deg has no partner
    assert ob.tolist() == [False, True]


def test_link_pairs_rate_window_and_significance():
    ra0, dec0 = 150.0, 10.0
    # a mover at 2"/yr in +Dec over 3 yr, and a pair implying 20"/yr (out of the window)
    a = _cat([ra0, ra0 + 1.0], [dec0, dec0], 2018.0)
    b = _cat([ra0, ra0 + 1.0], [dec0 + 6 / 3600, dec0 + 60 / 3600], 2021.0)
    p = v.link_pairs(a, b)
    assert len(p) == 1 and p.i[0] == 0 and p.j[0] == 0
    assert p.mu_dec[0] == pytest.approx(2.0, rel=1e-6)
    assert p.mu[0] == pytest.approx(2.0, rel=1e-6)
    # the same geometry with huge positional errors is not significant
    noisy = v.link_pairs(
        _cat(a.ra, a.dec, 2018.0, err=[5, 5]), _cat(b.ra, b.dec, 2021.0, err=[5, 5])
    )
    assert len(noisy) == 0


def test_link_pairs_empty_and_nonpositive_dt():
    a = _cat([1.0], [1.0], 2021.0)
    assert len(v.link_pairs(a, _cat([], [], 2021.0))) == 0
    assert len(v.link_pairs(a, _cat([1.0], [1.0 + 6 / 3600], 2018.0))) == 0  # b before a


def test_collinearity_accepts_line_rejects_kink():
    ra0, dec0 = 40.0, -5.0
    a = _cat([ra0], [dec0], 2018.0)
    b = _cat([ra0], [dec0 + 6 / 3600], 2021.0)  # 2"/yr north
    on_line = _cat([ra0], [dec0 + 11 / 3600], 2023.5)  # continues: +5" in 2.5 yr
    kinked = _cat([ra0 + 5 / 3600], [dec0 + 6 / 3600], 2023.5)  # turned east
    p = v.link_pairs(a, b)
    assert len(v.collinearity_test(a, b, on_line, p)) == 1
    assert len(v.collinearity_test(a, b, kinked, p)) == 0
    assert len(v.collinearity_test(a, b, _cat([], [], 2023.5), p)) == 0


def test_flux_consistent():
    fa, fb, fc = np.array([1.0, 1.0, 1.0]), np.array([2.0, 5.0, 0.0]), np.array([1.5, 1.0, 1.0])
    assert v.flux_consistent(fa, fb, fc, max_ratio=3.0).tolist() == [True, False, False]
    # off by default: flare stars vary by ~10x between epochs (UV Ceti)
    assert v.flux_consistent(fa, fb, fc).tolist() == [True, True, True]


def test_search_recovers_planted_movers_without_false_candidates():
    e1, e2, e3, mu = v.synthetic_epochs(seed=1)
    res = v.search(e1, e2, e3)
    a, b, c = res.orphans
    cand = res.candidates
    same = (
        (a.ident[cand.i] >= 0)
        & (a.ident[cand.i] == b.ident[cand.j])
        & (b.ident[cand.j] == c.ident[cand.k])
    )
    assert same.all()  # no false candidates in the fixture
    # Above BOTH static-radius floors: 2.5"/(E1->E2 ~3 yr) and the tighter 2.5"/(E2->E3 ~2 yr),
    # and isolated in every epoch (a chance neighbour within 30" costs completeness by design).
    iso = [{int(x) for x in e.ident[v.isolated_mask(e)] if x >= 0} for e in (e1, e2, e3)]
    fast = [k for k in np.flatnonzero(mu > 1.5) if all(k in s_ for s_ in iso)]
    assert len(fast) >= 10
    assert set(fast) <= set(a.ident[cand.i].tolist())
    # recovered rates match the planted ones
    for n in range(len(cand)):
        assert cand.mu[n] == pytest.approx(mu[a.ident[cand.i[n]]], rel=0.25)


def test_scramble_null_destroys_real_movers():
    e1, e2, e3, _ = v.synthetic_epochs(seed=2)
    res = v.search(e1, e2, e3)
    assert len(res.candidates) > 5
    null = v.scramble_null(res.orphans, n_reps=5, seed=3)
    assert null["n_reps"] == 5
    assert null["candidates_mean"] < 1.0
    assert len(null["candidates_per_rep"]) == 5


def test_completeness_rises_with_rate_and_has_a_floor():
    e1, e2, e3, _ = v.synthetic_epochs(n_movers=0, seed=4)
    comp = v.completeness(e1, e2, e3, n=600, seed=5)
    per = np.asarray(comp["per_bin"])
    assert per[0] < 0.2  # below ~0.5"/yr a mover stays inside the static-match radius
    assert per[-1] > 0.8
    assert 0.0 < comp["overall"] < 1.0
    assert len(comp["bin_edges"]) == len(per) + 1


def test_inject_movers_appends_and_tags():
    e1, e2, e3, _ = v.synthetic_epochs(n_movers=0, n_static=500, n_variable=0, seed=6)
    i1, i2, i3, mu = v.inject_movers(e1, e2, e3, n=50, seed=7)
    assert len(i1) == len(e1) + 50 and len(i3) == len(e3) + 50
    assert (i1.ident[-50:] >= 10_000_000).all()
    assert mu.min() >= v.MU_MIN_ARCSEC_YR and mu.max() <= v.MU_MAX_ARCSEC_YR


def test_surface_density_limit_divides_by_effective_area():
    assert v.surface_density_limit(0, 100.0, 1.0) == pytest.approx(0.02996)
    assert v.surface_density_limit(0, 100.0, 0.25) == pytest.approx(4 * 0.02996)
    assert v.surface_density_limit(3, 100.0, 0.5) == pytest.approx(0.06)
    assert v.surface_density_limit(0, 100.0, 0.0) == float("inf")


def test_epochcatalog_concat_subset_len():
    a = _cat([1.0, 2.0], [0.0, 0.0], 2018.0)
    b = _cat([3.0], [0.0], 2018.0)
    c = v.EpochCatalog.concat([a, b])
    assert len(c) == 3 and c.ra.tolist() == [1.0, 2.0, 3.0]
    assert len(c.subset(np.array([True, False, True]))) == 2
    assert (
        v.EpochCatalog(np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)).ident == -1
    ).all()


def test_run_offline_writes_synthetic_results(tmp_path):
    m = v.run(str(tmp_path), n_null=3, n_inject=300)
    assert m["is_real"] is False
    assert m["syn_n_recovered"] >= 10 and m["syn_n_false"] == 0
    saved = json.loads((tmp_path / "results" / "vlasspm_metrics.json").read_text())
    assert saved["syn_n_movers"] == m["syn_n_movers"]


def test_e3_residuals_are_calibrated():
    """True movers' E3 residuals, in sigma units, must follow a Rayleigh(1) distribution.

    Guards the error propagation: E2's error enters the prediction twice (rate and anchor), and
    treating those as independent under-disperses the tolerance (median ~1.4 instead of 1.18).
    """
    res = []
    for seed in range(4):
        e1, e2, e3, _ = v.synthetic_epochs(seed=seed, n_movers=200, n_static=5000, n_variable=0)
        r = v.search(e1, e2, e3)
        a, b, c = r.orphans
        t = r.triplets
        same = (a.ident[t.i] >= 0) & (a.ident[t.i] == b.ident[t.j]) & (b.ident[t.j] == c.ident[t.k])
        res += t.resid_sigma[same].tolist()
    assert len(res) > 300
    assert np.median(res) == pytest.approx(np.sqrt(2 * np.log(2)), abs=0.08)


def test_epoch_triples():
    assert v.epoch_triples(3) == [(0, 1, 2)]
    assert v.epoch_triples(4) == [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


def test_synthetic_epochs_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        v.synthetic_epochs(pos_err=(0.5, 0.3), epochs_t=(2018.0, 2021.0, 2024.0))


def test_mover_absent_from_first_epoch_found_through_later_triple():
    """The UV Ceti case: undetected in E1, present in E2-E4, recovered via the E2-E3-E4 triple."""
    *cats, mu = v.synthetic_epochs(
        n_movers=10,
        n_static=3000,
        n_variable=0,
        seed=8,
        pos_err=(0.5, 0.3, 0.2, 0.2),
        epochs_t=(2018.5, 2021.5, 2024.0, 2026.0),
    )
    k = int(np.argmax(mu))  # the fastest mover, so the recovery assertion always applies
    assert mu[k] > 1.5
    cats[0] = cats[0].subset(cats[0].ident != k)  # drop it from E1 entirely
    res = v.search_multi(cats)
    found = {t: set(r.orphans[0].ident[r.candidates.i].tolist()) for t, r in res.items()}
    assert all(k not in found[t] for t in found if t[0] == 0)  # no triple using E1 has it
    assert k in found[(1, 2, 3)]


def test_completeness_reports_per_triple():
    *cats, _ = v.synthetic_epochs(
        n_movers=0,
        n_static=3000,
        n_variable=0,
        seed=9,
        pos_err=(0.5, 0.3, 0.2, 0.2),
        epochs_t=(2018.5, 2021.5, 2024.0, 2026.0),
    )
    comp = v.completeness(*cats, n=300, seed=10)
    assert set(comp["per_triple"]) == {"E1-E2-E3", "E1-E2-E4", "E1-E3-E4", "E2-E3-E4"}
    assert comp["overall"] >= max(comp["per_triple"].values())  # union beats any one triple


def test_calibrate_floors_recovers_injected_astrometric_errors():
    *cats, _ = v.synthetic_epochs(
        n_movers=0,
        n_variable=0,
        n_static=20000,
        seed=11,
        pos_err=(0.6, 0.35, 0.15),
        epochs_t=(2018.5, 2021.5, 2024.0),
    )
    cal = v.calibrate_floors(cats, flux_min_mjy=0.0)
    assert cal["floor_arcsec"] == pytest.approx([0.6, 0.35, 0.15], abs=0.04)
    assert set(cal["pair_sigma_arcsec"]) == {"E1-E2", "E1-E3", "E2-E3"}


def test_calibrate_floors_underdetermined():
    a = _cat([1.0], [0.0], 2018.0)
    assert v.calibrate_floors([a, a])["floor_arcsec"] is None


def test_isolated_mask():
    c = _cat([10.0, 10.0 + 10 / 3600, 20.0], [0.0, 0.0, 0.0], 2020.0)
    assert v.isolated_mask(c).tolist() == [False, False, True]  # 10" pair is not isolated
    assert v.isolated_mask(c, radius_arcsec=5.0).tolist() == [True, True, True]
    assert v.isolated_mask(_cat([], [], 2020.0)).size == 0


def test_search_isolation_removes_split_extended_sources():
    """A double source whose two components shift between epochs mimics a mover; isolation stops it."""
    ra0, dec0 = 120.0, 20.0
    a = _cat([ra0, ra0 + 20 / 3600], [dec0, dec0], 2018.0)
    b = _cat([ra0, ra0 + 20 / 3600], [dec0 + 6 / 3600, dec0 + 6 / 3600], 2021.0)
    c = _cat([ra0, ra0 + 20 / 3600], [dec0 + 11 / 3600, dec0 + 11 / 3600], 2023.5)
    assert len(v.search(a, b, c, isolation_arcsec=None).candidates) == 2
    assert len(v.search(a, b, c).candidates) == 0
