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


def test_calibrate_floors_banded_sees_a_worse_southern_epoch():
    """E1 twice as noisy south of Dec -20: the banded fit must see it; the all-sky fit cannot."""
    rng = np.random.default_rng(12)
    n = 30000
    ra = rng.uniform(0, 60, n)
    dec = rng.uniform(-38, 28, n)
    south = dec < -20
    cats = []
    for e, (fn, fs) in enumerate([(0.3, 0.9), (0.2, 0.2), (0.15, 0.15)]):
        err = np.where(south, fs, fn)
        cats.append(
            v.EpochCatalog(
                ra + rng.normal(0, 1, n) * err / 3600 / np.cos(np.radians(dec)),
                dec + rng.normal(0, 1, n) * err / 3600,
                np.full(n, 2018.0 + 3 * e),
                np.full(n, 20.0),
                err,
            )
        )
    band = v.calibrate_floors_banded(cats, flux_min_mjy=0.0)
    south_band = next(b for b in band["bands"] if b["lo"] == -40.0)
    north_band = next(b for b in band["bands"] if b["lo"] == 0.0)
    assert south_band["floor_arcsec"][0] == pytest.approx(0.9, abs=0.08)
    assert north_band["floor_arcsec"][0] == pytest.approx(0.3, abs=0.05)
    f = v.banded_floor(np.array([-30.0, 10.0, 80.0]), 0, band, fallback=0.5)
    assert f[0] == pytest.approx(south_band["floor_arcsec"][0])
    assert f[1] == pytest.approx(north_band["floor_arcsec"][0])
    assert f[2] == 0.5  # empty band (Dec >= 30 here) falls back


# ------------------------------------------------------------ shape-aware error model (run 4)


def test_condon_point_source_limit():
    """Point source in a round beam: rho^2 = 2 SNR^2, so sigma = (FWHM / sqrt(8 ln2)) / SNR."""
    fwhm, snr = 2.5, 20.0
    sM, sm = v.condon_position_errors(fwhm, fwhm, fwhm, fwhm, snr)
    expect = fwhm / np.sqrt(8 * np.log(2)) / snr
    assert float(sM) == pytest.approx(expect, rel=1e-12)
    assert float(sm) == pytest.approx(expect, rel=1e-12)
    # and it scales as 1/SNR
    s2, _ = v.condon_position_errors(fwhm, fwhm, fwhm, fwhm, 2 * snr)
    assert float(s2) == pytest.approx(expect / 2, rel=1e-12)


def test_condon_resolved_source_has_larger_error():
    beam = (3.0, 2.0)
    point = v.condon_position_errors(3.0, 2.0, *beam, 10.0)
    resolved = v.condon_position_errors(6.0, 2.5, *beam, 10.0)  # same peak S/N, extended
    assert resolved[0] > point[0] and resolved[1] > point[1]
    # an elongated fit is less certain along its long axis
    assert resolved[0] > resolved[1]


def test_condon_reproduces_the_vlass_catalogue_errors():
    """Row 0 of QL3.1 (VLASS3QL J083036.46-241605.9): the catalogue's E_RA/E_DEC (deg) are this
    formula rotated by PA, and E_RA is on-sky (no cos dec) --- the run 1-3 loader got that wrong."""
    cov = v.condon_cov(
        3.6592156526759996,
        2.50531610598,
        152.8896549865025,
        2.43607401847836,
        1.77422451972948,
        10.97892573861 / 0.16255225636999998,
    )
    # 0.25%: PyBDSF used the rms map at the source, the catalogue quotes the island rms
    assert np.sqrt(cov[0, 0]) == pytest.approx(5.05267868e-06 * 3600, rel=5e-3)
    assert np.sqrt(cov[0, 1]) == pytest.approx(7.09375924e-06 * 3600, rel=5e-3)
    assert np.sqrt(cov[0, 0]) != pytest.approx(
        5.05267868e-06 * 3600 * np.cos(np.radians(24.27)), rel=5e-2
    )


def test_ellipse_cov_rotation_by_pa():
    # PA 0: major axis north (y); PA 90: east (x); PA 45: correlated
    c0 = v.ellipse_cov(2.0, 1.0, 0.0)[0]
    c90 = v.ellipse_cov(2.0, 1.0, 90.0)[0]
    c45 = v.ellipse_cov(2.0, 1.0, 45.0)[0]
    assert c0 == pytest.approx([1.0, 4.0, 0.0], abs=1e-12)
    assert c90 == pytest.approx([4.0, 1.0, 0.0], abs=1e-12)
    assert c45 == pytest.approx([2.5, 2.5, 1.5], abs=1e-12)
    a, b, pa = v.cov_axes(v.ellipse_cov([2.0, 3.0], [1.0, 0.5], [30.0, -60.0]))
    assert a == pytest.approx([2.0, 3.0]) and b == pytest.approx([1.0, 0.5])
    assert pa == pytest.approx([30.0, -60.0])


def test_abs_cov_and_chi2():
    m = np.array([[3.0, -1.0, 0.0]])  # indefinite: grew in x, shrank in y
    assert v.abs_cov(m)[0] == pytest.approx([3.0, 1.0, 0.0])
    r = v.ellipse_cov(2.0, 1.0, 30.0) - v.ellipse_cov(1.0, 2.0, 30.0)  # eigen +3/-3
    ax = v.cov_axes(v.abs_cov(r))
    assert ax[0][0] == pytest.approx(np.sqrt(3.0)) and ax[1][0] == pytest.approx(np.sqrt(3.0))
    iso = np.array([[0.04, 0.04, 0.0]])
    assert v.chi2_2d(0.3, 0.4, iso)[0] == pytest.approx(0.25 / 0.04)
    # an offset along an elongated error ellipse is less significant than across it
    el = v.ellipse_cov(1.0, 0.1, 0.0)
    assert v.chi2_2d(0.0, 0.5, el)[0] < v.chi2_2d(0.5, 0.0, el)[0]


def test_error_model_sys_cov():
    shape = v.ellipse_cov(4.0, 2.0, 0.0)
    beam_a, beam_b = v.ellipse_cov(4.8, 1.7, 0.0), v.ellipse_cov(3.5, 2.2, 0.0)
    ref = 0.5 * (beam_a + beam_b)
    assert v.NO_MODEL.sys_cov(shape, beam_a, ref) == pytest.approx(np.zeros((1, 3)))
    s = v.ErrorModel(0.1, 0.0).sys_cov(shape, beam_a, ref)
    assert s[0] == pytest.approx([0.04, 0.16, 0.0])  # 0.1 x the FWHM along each source axis
    b = v.ErrorModel(0.0, 0.1).sys_cov(np.zeros((1, 3)), beam_a, ref)
    assert (b[0, :2] > 0).all()  # a beam change moves a centroid on both axes


def test_link_significance_uses_the_error_ellipse():
    """A 1.2" shift along a static source's long, uncertain axis is not significant; across it,
    it is. The run 1-3 circularised error could not tell these apart."""
    ra0, dec0 = 100.0, 10.0
    cov = v.ellipse_cov(0.5, 0.1, 0.0)  # long axis north
    a = v.EpochCatalog([ra0], [dec0], [2018.0], [3.0], [0.36], cov=cov)
    b_n = v.EpochCatalog([ra0], [dec0 + 1.2 / 3600], [2021.0], [3.0], [0.36], cov=cov)
    east = ra0 + 1.2 / 3600 / np.cos(np.radians(dec0))
    b_e = v.EpochCatalog([east], [dec0], [2021.0], [3.0], [0.36], cov=cov)
    assert len(v.link_pairs(a, b_n)) == 0
    assert len(v.link_pairs(a, b_e)) == 1


def test_fit_error_model_recovers_injected_systematics():
    a, b = v.synthetic_statics(seed=3)
    terms = v.static_pair_terms(a, b, v.match_statics(a, b))
    model, surface = v.fit_error_model(terms)
    assert model.k_struct == pytest.approx(0.1, abs=0.02)
    assert model.q_beam == pytest.approx(0.03, abs=0.03)
    assert np.asarray(surface["loss"]).shape == (31, 21)
    rows = v.calibration_table(terms, model, by="size", edges=v.SIZE_BIN_EDGES)
    for r in rows:
        if r["n"] > 1000:
            assert r["f_gt3"] == pytest.approx(v.RAYLEIGH_TAIL_3SIGMA, abs=0.006)
            assert r["median_z"] == pytest.approx(v.RAYLEIGH_MEDIAN, abs=0.08)


def test_measurement_only_errors_are_caught_being_optimistic_for_resolved_sources():
    """The run 1-3 failure mode: faint, resolved statics scatter more than the catalogue says.
    The calibration table must flag it (tail far above Rayleigh), and the fitted model fix it."""
    a, b = v.synthetic_statics(seed=4)
    terms = v.static_pair_terms(a, b, v.match_statics(a, b))
    big = v.calibration_table(terms, v.NO_MODEL, by="size", edges=(3.0, np.inf))[0]
    assert big["f_gt3"] > 10 * v.RAYLEIGH_TAIL_3SIGMA  # flagged
    model, _ = v.fit_error_model(terms)
    fixed = v.calibration_table(terms, model, by="size", edges=(3.0, np.inf))[0]
    assert fixed["f_gt3"] < 2 * v.RAYLEIGH_TAIL_3SIGMA  # handled
    # ... and it is what makes static sources look like movers: link significance >= 3
    resolved = np.flatnonzero(v.cov_axes(a.shape)[0] > 3.0)[:300]
    z_old = [v.triplet_significance(a, i, b, i, b, i)["z_12"] for i in resolved]
    z_new = [v.triplet_significance(a, i, b, i, b, i, model)["z_12"] for i in resolved]
    assert np.mean(np.array(z_old) >= 3) > 0.15
    assert np.mean(np.array(z_new) >= 3) < 0.04


def test_triplet_significance_on_a_true_mover():
    ra0, dec0 = 40.0, -5.0
    a = _cat([ra0], [dec0], 2018.0)
    b = _cat([ra0], [dec0 + 6 / 3600], 2021.0)
    c = _cat([ra0], [dec0 + 11 / 3600], 2023.5)
    s = v.triplet_significance(a, 0, b, 0, c, 0)
    assert s["z_12"] == pytest.approx(6 / np.hypot(0.2, 0.2))
    assert s["z_13"] == pytest.approx(11 / np.hypot(0.2, 0.2))
    assert s["resid"] == pytest.approx(0.0, abs=1e-6)


def _shaped_field(seed=0, n=4000, n_epochs=3, crowd=False):
    """A shaped static field (for injections): beams, rms and floors on every component."""
    cats = v.synthetic_statics(n=n, n_epochs=n_epochs, seed=seed, area_side_deg=4.0)
    if crowd:  # add a 2nd component 10-20" from every source: non-isolated pairs everywhere
        out = []
        for c in cats:
            off = np.random.default_rng(seed).uniform(10, 20, len(c)) / 3600
            extra = c.subset(np.ones(len(c), bool))
            extra.dec = extra.dec + off
            extra.ident = extra.ident + 10**6
            out.append(v.EpochCatalog.concat([c, extra]))
        cats = out
    return cats


def test_injected_scatter_follows_the_quoted_errors_plus_model():
    cats = _shaped_field(seed=5)
    assert all(c.has_shape for c in cats)
    model = v.ErrorModel(0.1, 0.03)
    out, mu, _, _ = v._inject(
        *cats, n=3000, seed=6, size_arcsec=3.0, model=model, mu_range=(1.0, 1.0001)
    )
    # E1 -> E2 offset minus the true motion, normalised by the pair covariance (with the model)
    a, b = out[0], out[1]
    ia = np.flatnonzero(a.ident >= 10_000_000)
    ib = np.flatnonzero(b.ident >= 10_000_000)
    assert (a.ident[ia] == b.ident[ib]).all()
    dx, dy = v.tangent_offsets_arcsec(a.ra[ia], a.dec[ia], b.ra[ib], b.dec[ib])
    dt = b.t_yr[ib] - a.t_yr[ia]
    sep = np.hypot(dx, dy)
    assert sep.mean() == pytest.approx(float(np.mean(mu * dt)), rel=0.05)
    # the quoted per-detection covariance is the Condon ellipse for a 3" source + floor
    maj = v.cov_axes(a.cov[ia])[0]
    point_out, _, _, _ = v._inject(*cats, n=3000, seed=6, size_arcsec=0.0)
    maj_pt = v.cov_axes(point_out[0].cov[point_out[0].ident >= 10_000_000])[0]
    assert np.median(maj) > np.median(maj_pt)  # a resolved injection is quoted a larger error


def test_injected_static_offsets_are_rayleigh_under_the_model():
    """Zero-motion injections: their epoch-to-epoch offsets, normalised by exactly the covariance
    the search uses (measurement + model with the pair-mean beam), must be Rayleigh(1)."""
    cats = _shaped_field(seed=7, n_epochs=2)
    model = v.ErrorModel(0.1, 0.03)
    out, _, _, _ = v._inject(
        *cats, n=4000, seed=8, size_arcsec=2.0, model=model, mu_range=(1e-9, 1.1e-9)
    )
    a, b = out
    ia = np.flatnonzero(a.ident >= 10_000_000)
    ib = np.flatnonzero(b.ident >= 10_000_000)
    dx, dy = v.tangent_offsets_arcsec(a.ra[ia], a.dec[ia], b.ra[ib], b.dec[ib])
    z = np.sqrt(v.chi2_2d(dx, dy, v._pair_cov(a, ia, b, ib, model)))
    # with a two-epoch field the injection's beam reference IS the pair mean, as in linkage
    assert np.median(z) == pytest.approx(v.RAYLEIGH_MEDIAN, abs=0.06)
    assert np.mean(z > 3) == pytest.approx(v.RAYLEIGH_TAIL_3SIGMA, abs=0.006)
    zo = np.sqrt(v.chi2_2d(dx, dy, a.cov[ia] + b.cov[ib]))  # measurement-only: too optimistic
    assert np.mean(zo > 3) > 3 * v.RAYLEIGH_TAIL_3SIGMA


def test_realistic_placement_pays_the_isolation_cost():
    """Realistic injections sit at random sky positions, so their nearest-neighbour distances
    follow the field's (not the 60-120" the run 1-3 injections were given), and the ones that
    land within 30" of a real component are lost to the isolation cut, as a real mover would be."""
    cats = _shaped_field(seed=9, n=6000, crowd=True)
    first = 10_000_000
    out, mu, realistic, nn = v._inject(
        *cats, n=3000, seed=10, realistic_frac=0.5, mu_range=(2.0, 3.0)
    )
    assert realistic.mean() == pytest.approx(0.5, abs=0.05)
    rng = np.random.default_rng(11)  # random points in the (shrunk) field for comparison
    e1 = cats[0]
    pts_ra = rng.uniform(e1.ra.min() + 0.3, e1.ra.max() - 0.3, 20000)
    pts_dec = rng.uniform(e1.dec.min() + 0.3, e1.dec.max() - 0.3, 20000)
    d, _ = v.cKDTree(v._xyz(e1.ra, e1.dec)).query(v._xyz(pts_ra, pts_dec), k=1)
    nn_random = np.degrees(2 * np.arcsin(d / 2)) * 3600
    assert np.median(nn[realistic]) == pytest.approx(np.median(nn_random), rel=0.1)
    res = v.search_multi(out)
    found = set()
    for r in res.values():
        found |= {int(x) for x in v._recovered_idents(r, first)}
    found_arr = np.array(sorted(found))
    assert found_arr.size > 100
    # a mover with a neighbour inside 30" is never found (1" slack: detections scatter)
    assert (nn[found_arr] >= v.ISOLATION_ARCSEC - 1.0).all()
    lost_near = realistic & (nn < v.ISOLATION_ARCSEC)
    assert lost_near.sum() > 20  # the realistic class DOES land near sources
    res2 = v.completeness(*cats, n=600, seed=10, realistic_frac=0.5)
    assert res2["realistic"]["n"] + res2["isolated"]["n"] == 600
    assert res2["realistic_frac"] == 0.5 and res2["size_arcsec"] == 0.0


def test_scramble_and_search_accept_a_model():
    e1, e2, e3, _ = v.synthetic_epochs(seed=2)
    model = v.ErrorModel(0.1, 0.0)  # shapes are zero in this fixture: identical to no model
    r0 = v.search(e1, e2, e3)
    r1 = v.search(e1, e2, e3, model=model)
    assert len(r0.candidates) == len(r1.candidates)
    null = v.scramble_null(r1.orphans, n_reps=2, seed=3, model=model)
    assert null["n_reps"] == 2


def test_epochcatalog_carries_shape_fields_through_subset_and_concat():
    cats = v.synthetic_statics(n=50, seed=1)
    c = v.EpochCatalog.concat([cats[0], cats[1]])
    assert c.cov.shape == (100, 3) and c.beam.shape == (100, 3) and c.rms.shape == (100,)
    s = c.subset(np.arange(100) < 10)
    assert s.shape.shape == (10, 3) and s.floor.shape == (10,)
    plain = _cat([1.0], [0.0], 2018.0)
    assert not plain.has_shape
    assert plain.cov[0] == pytest.approx([0.04, 0.04, 0.0])
