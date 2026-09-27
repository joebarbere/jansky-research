"""Tests for jansky_research.fashienv -- environment-split FASHI HI mass function. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import fashienv as fe


def test_comoving_xyz_radius_matches_distance():
    z = np.array([0.01, 0.03, 0.06])
    xyz = fe.comoving_xyz(np.array([10.0, 200.0, 300.0]), np.array([0.0, 30.0, -10.0]), z)
    r = np.linalg.norm(xyz, axis=1)
    # |xyz| must equal the comoving distance for that z
    assert np.allclose(r, fe._comoving_distance_mpc(z), rtol=1e-6)


def test_void_membership_inside_and_outside():
    spheres = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
    radii = np.array([10.0, 5.0])
    gal = np.array([[3.0, 0.0, 0.0], [50.0, 0.0, 0.0], [101.0, 2.0, 0.0]])
    m = fe.void_membership(gal, spheres, radii)
    assert list(m) == [True, False, True]  # in sphere 0, in the gap, in sphere 1


def test_assign_groups_and_radius():
    # one group at (RA,Dec)=(150,2), cz=6000, R200=1 Mpc; a member and a field galaxy
    grp_ra, grp_dec, grp_cz, grp_r200 = (
        np.array([150.0]),
        np.array([2.0]),
        np.array([6000.0]),
        np.array([1.0]),
    )
    # member: 0.3 Mpc projected at cz 6000 -> ~0.2 deg; field: 5 deg away
    gal_ra = np.array([150.2, 155.0])
    gal_dec = np.array([2.0, 2.0])
    gal_cz = np.array([6050.0, 6000.0])
    idx = fe.assign_groups(gal_ra, gal_dec, gal_cz, grp_ra, grp_dec, grp_cz, grp_r200)
    assert idx[0] == 0 and idx[1] == -1
    rr = fe.clustercentric_radius(gal_ra, gal_dec, gal_cz, idx, grp_ra, grp_dec, grp_cz, grp_r200)
    assert 0.0 < rr[0] < 1.0 and np.isnan(rr[1])


def test_assign_groups_velocity_cut():
    grp = (np.array([150.0]), np.array([2.0]), np.array([6000.0]), np.array([1.0]))
    # spatially coincident but 3000 km/s away -> not a member
    idx = fe.assign_groups(np.array([150.0]), np.array([2.0]), np.array([9000.0]), *grp)
    assert idx[0] == -1


def test_vmax_scales_with_flux():
    # a brighter source is detectable to larger distance -> larger Vmax
    v = fe.vmax_1vmax(
        np.array([9.0, 9.0]),
        np.array([100.0, 100.0]),
        np.array([1.0, 4.0]),
        flux_limit=0.3,
    )
    assert v[1] > v[0]


def test_schechter_shape():
    lm = np.linspace(7, 11, 50)
    phi = fe.schechter(lm, -2.5, 9.9, -1.3)
    assert np.all(phi > 0)
    # the knee: phi turns over above logM*
    assert phi[np.argmin(np.abs(lm - 10.5))] < phi[np.argmin(np.abs(lm - 9.9))]


def test_himf_recovers_injected_schechter():
    # a flux-limited mock drawn from a SINGLE Schechter must fit back to it
    cat = fe.synthetic_environment_catalogue(n=60000, void_frac=0.0, seed=1)
    vmax = fe.vmax_1vmax(cat["log_mhi"], cat["dist_mpc"], cat["flux"])
    h = fe.himf(cat["log_mhi"], vmax, area_sr=cat["area_sr"])
    fit = fe.fit_schechter(h)
    assert abs(fit["log_m_star"] - cat["truth"]["wall"][0]) < 0.2
    assert abs(fit["alpha"] - cat["truth"]["wall"][1]) < 0.2


def test_environment_split_recovers_both_himfs():
    cat = fe.synthetic_environment_catalogue(n=120000, seed=2)
    for env, mask in (("void", cat["is_void"]), ("wall", ~cat["is_void"])):
        vmax = fe.vmax_1vmax(cat["log_mhi"][mask], cat["dist_mpc"][mask], cat["flux"][mask])
        fit = fe.fit_schechter(fe.himf(cat["log_mhi"][mask], vmax, area_sr=cat["area_sr"]))
        t_lms, t_a = cat["truth"][env]
        assert abs(fit["log_m_star"] - t_lms) < 0.25, env
        assert abs(fit["alpha"] - t_a) < 0.25, env
    # the injected signal: void knee is below wall knee
    assert cat["truth"]["void"][0] < cat["truth"]["wall"][0]


def test_fit_schechter_insufficient_bins():
    h = {
        "logm": np.array([9.0, 9.5]),
        "phi": np.array([1e-3, 1e-4]),
        "phi_err": np.array([1e-4, 1e-5]),
        "counts": np.array([5, 5]),
        "dlog": 0.5,
    }
    fit = fe.fit_schechter(h)
    assert np.isnan(fit["log_m_star"])


def test_run_offline_writes_artifacts(tmp_path):
    m = fe.run(str(tmp_path), offline=True)
    assert m["source"].startswith("synthetic")
    # recover-a-known: fitted void/wall knees near the injected truth, void below wall
    assert abs(m["himf_void"]["log_m_star"] - m["true_void_logmstar"]) < 0.3
    assert abs(m["himf_wall"]["log_m_star"] - m["true_wall_logmstar"]) < 0.3
    assert m["himf_void"]["log_m_star"] < m["himf_wall"]["log_m_star"]
    # the knee offset recovers the injected sign and rough magnitude
    assert m["void_knee_offset"] < 0 and abs(m["void_knee_offset"] - m["true_knee_offset"]) < 0.2
    saved = json.loads((tmp_path / "results" / "fashienv_metrics.json").read_text())
    assert saved["n_sources"] == m["n_sources"]
    assert (tmp_path / "papers" / "fashienv" / "figures" / "fashienv.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "fashienv" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\feSynVoidLogMStar}" in macros
    assert r"\newcommand{\feRealVoidLogMStar}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    fe._write_macros(
        {"source": "x", "is_real": True, "n_sources": 5, "himf_void": {"log_m_star": None}}, p
    )
    txt = p.read_text()
    assert r"\newcommand{\feRealVoidLogMStar}{--}" in txt
    assert r"\newcommand{\feSynVoidLogMStar}{--}" in txt


def test_void_jackknife_measures_sample_variance_and_is_committed():
    """The error the fit errors cannot see -- and a test that could have gone either way.

    The quoted knee-offset error is Poisson counting noise within one realisation of
    large-scale structure. For a void statistic the usual worry is that void-to-void variance
    dominates it, which counting noise cannot see (the `rmstructure` failure mode). Deleting
    each occupied void in turn and refitting measures it directly. Here it comes out SMALLER
    than the fit error, so the offset survives -- but the number has to be committed either
    way, because the paper's argument now depends on it.
    """
    from pathlib import Path

    path = Path("results/fashienv_metrics.json")
    if not path.exists():  # pragma: no cover - absent in a bare checkout
        pytest.skip("committed fashienv results not present")
    m = json.loads(path.read_text())
    jk = m.get("void_jackknife")
    assert jk and jk.get("B_jackknife_err") is not None, "the paper quotes this; commit it"
    assert jk.get("B_minus_optA_same_jackknife_err") is not None  # the paired weighting error
    assert 0 < jk["n_ok"] == jk["n_voids_occupied"]
    # the jackknife mean must sit on the full-sample offset, or the resampling is wrong
    assert jk["B_mean_offset"] == pytest.approx(m["void_knee_offset"], abs=0.01)
    assert jk["B_min_offset"] <= m["void_knee_offset"] <= jk["B_max_offset"]


def test_comparison_bin_size_is_committed():
    """The comparison ("wall") bin's size is committed, and it is the in-footprint pool.

    n_wall and n_field were absent from the evidence file, so a reader could not learn how big
    the comparison sample was -- and it is 58% of all of DR1, with a knee 0.010 dex from the
    all-sky fit. That is the single most important number for reading this result.
    """
    from pathlib import Path

    path = Path("results/fashienv_metrics.json")
    if not path.exists():  # pragma: no cover - absent in a bare checkout
        pytest.skip("committed fashienv results not present")
    m = json.loads(path.read_text())
    assert "n_wall" in m and "n_field" in m
    assert m["n_wall"] + m["n_in_void"] == m["n_classifiable_void"]
    # Sixth referee round: the comparison sample is the in-SDSS-footprint pool, not the padded
    # box (which put 24.5% unclassifiable galaxies into wall AND field). This used to assert the
    # wall knee matched the global fit to 0.05 dex -- true only BECAUSE the box pool padded the
    # wall with out-of-footprint galaxies, so it locked the defect in. The robustness leg must
    # use the same pool as the headline split.
    assert m["n_classifiable_void"] == m["robustness"]["n_pool"]
    assert m["wall_pct_of_sample"] < 50.0


def test_vmax_from_catalogue_is_completeness_weighted_per_sr():
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    v = fe.vmax_from_catalogue(np.array([1e6, 2e6, np.nan, 5e5]), np.array([0.5, 1.0, 0.9, 0.0]))
    assert v[0] == pytest.approx(0.5e6 / omega)
    assert v[1] == pytest.approx(2e6 / omega)
    assert v[2] == 0.0 and v[3] == 0.0  # non-finite / zero completeness drop out of himf
    # multiplied back by the survey area it is C * Vmax: the DR2 1/(C * Vmax) weight
    assert v[0] * omega == pytest.approx(0.5e6)


def test_load_fashi_dr2_maps_columns_and_units(tmp_path):
    p = tmp_path / "t2.csv"
    p.write_text(
        "ra,dec,v_opt,z_opt,W_50,S_sum,distance,mass,completeness,Vmax\n"
        "10.0,5.0,3000.0,0.01,120.0,450.0,43.0,9.1,0.8,2.5e6\n"
        "11.0,6.0,bad,0.02,100.0,,86.0,9.5,,\n"
    )
    d = fe.load_fashi_dr2(p)
    assert d["ra"].tolist() == [10.0, 11.0]
    assert d["flux"][0] == pytest.approx(0.45)  # mJy km/s -> Jy km/s
    assert d["cz"][0] == 3000.0 and np.isnan(d["cz"][1])
    assert d["completeness"][0] == 0.8 and np.isnan(d["completeness"][1])
    assert d["vmax_mpc3"][0] == 2.5e6 and np.isnan(d["vmax_mpc3"][1])
    assert set(d) >= {"ra", "dec", "cz", "z", "w50", "flux", "dist_mpc", "log_mhi"}


def test_himf_and_fit_uses_supplied_vmax():
    cat = fe.synthetic_environment_catalogue()
    lm, dd, ff = cat["log_mhi"], cat["dist_mpc"], cat["flux"]
    v1 = fe.vmax_1vmax(lm, dd, ff)
    h_default, _ = fe._himf_and_fit(lm, dd, ff, cat["area_sr"])
    h_same, _ = fe._himf_and_fit(lm, dd, ff, cat["area_sr"], vmax=v1)
    h_half, _ = fe._himf_and_fit(lm, dd, ff, cat["area_sr"], vmax=2.0 * v1)
    np.testing.assert_allclose(h_same["phi"], h_default["phi"])
    np.testing.assert_allclose(h_half["phi"], 0.5 * h_default["phi"])  # twice the volume
    mask = cat["is_void"]
    h_m, _ = fe._himf_and_fit(lm, dd, ff, cat["area_sr"], mask=mask, vmax=v1)
    lm_m = lm[mask & (v1 > 0)]
    in_bins = (lm_m >= 6.5) & (lm_m < 11.0)  # himf's default bin range
    assert h_m["counts"].sum() == int(in_bins.sum())  # exactly the masked sources, no others


def test_dr2_macros_are_emitted_from_nested_metrics(tmp_path):
    m = {
        "source": "x", "is_real": True, "release": "DR2", "n_sources": 10, "n_dr2_catalogue": 156411,
        "c_min": 0.5, "himf_global": {"log_m_star": 9.9, "log_m_star_err": 0.02, "alpha": -1.3},
        "optA_single_flux_cut": {"void_knee_offset": -0.25, "himf_global": {"alpha": -1.78}},
        "dr1_optA_single_flux_cut": {"n_sources": 41741},
    }  # fmt: skip
    p = tmp_path / "macros.tex"
    fe._write_macros(m, p)
    t = p.read_text()
    assert r"\newcommand{\feRealNDRTwo}{156411}" in t
    assert r"\newcommand{\feRealOptAVoidKneeOffset}{-0.250}" in t  # floats to 3 decimals
    assert r"\newcommand{\feRealOptAGlobalAlpha}{-1.780}" in t
    assert r"\newcommand{\feRealDROneN}{41741}" in t
    assert r"\newcommand{\feRealEdsVoidKneeOffset}{--}" in t  # absent key -> placeholder
    # a DR1 or synthetic run never emits the DR2 block
    fe._write_macros({**m, "release": "DR1"}, tmp_path / "m1.tex")
    assert "feRealNDRTwo" not in (tmp_path / "m1.tex").read_text()


def test_void_membership_holes_matches_the_brute_force_test():
    rng = np.random.default_rng(3)
    gal = rng.uniform(-50, 50, (4000, 3))
    holes = rng.uniform(-40, 40, (30, 3))
    radii = rng.uniform(3, 12, 30)
    brute = fe.void_membership(gal, holes, radii)
    np.testing.assert_array_equal(fe.void_membership_holes(gal, holes, radii), brute)
    assert fe.void_membership_holes(gal[:0], holes, radii).size == 0


def test_cmb_to_helio_cz_follows_the_dipole():
    # Toward the CMB apex the Sun approaches: heliocentric cz is 369.82 km/s LOWER; anti-apex higher.
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    apex = SkyCoord(
        l=fe.CMB_DIPOLE_LB[0] * u.deg, b=fe.CMB_DIPOLE_LB[1] * u.deg, frame="galactic"
    ).icrs
    anti = SkyCoord(
        l=fe.CMB_DIPOLE_LB[0] * u.deg + 180 * u.deg,
        b=-fe.CMB_DIPOLE_LB[1] * u.deg,
        frame="galactic",
    ).icrs
    out = fe.cmb_to_helio_cz(
        np.array([apex.ra.deg, anti.ra.deg]),
        np.array([apex.dec.deg, anti.dec.deg]),
        np.array([10000.0, 10000.0]),
    )
    assert out[0] == pytest.approx(10000.0 - 369.82, abs=0.05)
    assert out[1] == pytest.approx(10000.0 + 369.82, abs=0.05)


def test_random_void_positions_is_rigid_and_stays_in_footprint():
    rng = np.random.default_rng(5)
    # two voids of three holes each, at distances ~100 and ~200
    base = np.array(
        [[100.0, 0, 0], [102, 3, 0], [99, -2, 4], [0, 200.0, 0], [3, 203, 1], [-2, 198, 5]]
    )
    vid = np.array([1, 1, 1, 2, 2, 2])
    fra = rng.uniform(150, 200, 5000)  # footprint: a patch of sky
    fdec = rng.uniform(10, 40, 5000)
    moved, diag = fe.random_void_positions(base, vid, fra, fdec, rng)
    assert diag["n_unplaced"] == 0
    for v in (1, 2):
        a, b = base[vid == v], moved[vid == v]
        np.testing.assert_allclose(np.linalg.norm(b, axis=1), np.linalg.norm(a, axis=1), rtol=1e-9)
        # internal geometry preserved: all pairwise hole separations unchanged
        da = np.linalg.norm(a[:, None] - a[None], axis=2)
        db = np.linalg.norm(b[:, None] - b[None], axis=2)
        np.testing.assert_allclose(db, da, atol=1e-9)
        c = b.mean(axis=0)
        ra = np.degrees(np.arctan2(c[1], c[0])) % 360
        dec = np.degrees(np.arcsin(c[2] / np.linalg.norm(c)))
        assert 146 <= ra <= 204 and 6 <= dec <= 44  # inside the patch (+ one cell of slack)


def test_void_members_and_paired_jackknife():
    cat = fe.synthetic_environment_catalogue()
    n = cat["log_mhi"].size
    rng = np.random.default_rng(9)
    xyz = rng.uniform(-100, 100, (n, 3))
    holes = rng.uniform(-80, 80, (40, 3))
    radii = np.full(40, 18.0)
    vid = np.repeat(np.arange(20), 2)
    vm = fe.void_members(xyz, holes, radii, vid)
    in_void = fe.void_membership_holes(xyz, holes, radii)
    assert (vm["n_voids"] > 0).sum() == in_void.sum()
    base = np.ones(n, bool)
    v1 = fe.vmax_1vmax(cat["log_mhi"], cat["dist_mpc"], cat["flux"])
    out = fe.void_jackknife(
        cat, cat["area_sr"], in_void, np.ones(n, bool), vm, {"a": (base, v1), "b": (base, 2 * v1)}
    )
    assert out["n_voids_occupied"] > 5 and out["n_ok"] > 5
    # scaling the volume by a constant shifts phi, not the knee: the paired difference is ~0
    assert out["b_minus_a_jackknife_err"] == pytest.approx(0.0, abs=1e-6)
    assert out["a_jackknife_err"] > 0


def test_fit_schechter_reports_fit_quality():
    cat = fe.synthetic_environment_catalogue()
    h, fit = fe._himf_and_fit(cat["log_mhi"], cat["dist_mpc"], cat["flux"], cat["area_sr"])
    assert fit["n_bins"] >= 4 and np.isfinite(fit["red_chi2"]) and fit["red_chi2"] > 0
    assert -1.0 <= fit["corr_mstar_alpha"] <= 1.0


def _assign_groups_bruteforce(gra, gdec, gcz, pra, pdec, pcz, r200, dv_max=1000.0, h0=70.0):
    """The original O(N_gal x N_grp) loop, kept as the oracle for the KD-tree version."""
    ra, dec = np.radians(gra), np.radians(gdec)
    qra, qdec = np.radians(pra), np.radians(pdec)
    out = np.full(len(gcz), -1, int)
    for i in range(len(gcz)):
        near = np.abs(gcz[i] - pcz) <= dv_max
        if not near.any():
            continue
        h = (
            np.sin((qdec - dec[i]) / 2) ** 2
            + np.cos(dec[i]) * np.cos(qdec) * np.sin((qra - ra[i]) / 2) ** 2
        )
        sep = 2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))) * (pcz / h0)
        ratio = np.where(near, sep / np.maximum(r200, 1e-6), np.inf)
        j = int(np.argmin(ratio))
        if ratio[j] <= 1.0:
            out[i] = j
    return out


@pytest.mark.parametrize("h0", [70.0, 67.8])
def test_assign_groups_kdtree_matches_the_bruteforce_oracle(h0):
    rng = np.random.default_rng(21)
    n, m = 3000, 400
    gra, gdec = rng.uniform(150, 170, n), rng.uniform(10, 25, n)
    gcz = rng.uniform(3000, 15000, n)
    pra, pdec = rng.uniform(150, 170, m), rng.uniform(10, 25, m)
    pcz, r200 = rng.uniform(3000, 15000, m), rng.uniform(0.3, 1.5, m)
    fast = fe.assign_groups(gra, gdec, gcz, pra, pdec, pcz, r200, h0=h0)
    slow = _assign_groups_bruteforce(gra, gdec, gcz, pra, pdec, pcz, r200, h0=h0)
    np.testing.assert_array_equal(fast, slow)
    assert (fast >= 0).sum() > 50  # the test exercises real assignments, not just -1s


def test_random_void_positions_constrained_keeps_holes_in_footprint_and_apart():
    rng = np.random.default_rng(31)
    # 12 voids of 3 holes each, radius 5, at distances 80-120
    base, vid = [], []
    for k in range(12):
        c = rng.normal(0, 1, 3)
        c = c / np.linalg.norm(c) * rng.uniform(80, 120)
        for _ in range(3):
            base.append(c + rng.normal(0, 3, 3))
            vid.append(k)
    base, vid = np.asarray(base), np.asarray(vid)
    rad = np.full(len(base), 5.0)
    fra, fdec = rng.uniform(140, 220, 20000), rng.uniform(0, 50, 20000)
    moved, diag = fe.random_void_positions(
        base, vid, fra, fdec, rng, hole_radius=rad, constrained=True, no_overlap=True
    )
    assert diag["n_unplaced"] == 0 and diag["spill_frac"] == 0.0
    ra = np.degrees(np.arctan2(moved[:, 1], moved[:, 0])) % 360
    dec = np.degrees(np.arcsin(moved[:, 2] / np.linalg.norm(moved, axis=1)))
    assert np.all(
        (ra >= 138) & (ra <= 222) & (dec >= -2) & (dec <= 52)
    )  # every hole, not just centre
    for a in range(12):
        for b in range(a + 1, 12):
            d = np.linalg.norm(moved[vid == a][:, None] - moved[vid == b][None], axis=2)
            assert d.min() >= 10.0 - 1e-9  # radius 5 + radius 5: no overlap between voids


def test_random_group_positions_land_in_footprint():
    rng = np.random.default_rng(41)
    fra, fdec = rng.uniform(150, 200, 5000), rng.uniform(10, 40, 5000)
    ra, dec = fe.random_group_positions(np.zeros(300), np.zeros(300), fra, fdec, rng)
    assert ra.size == 300 and np.all((ra >= 148) & (ra <= 202) & (dec >= 8) & (dec <= 42))


def test_knee_offset_common_bins_uses_only_shared_bins():
    cat = fe.synthetic_environment_catalogue()
    v, w = cat["is_void"], ~cat["is_void"]
    out = fe.knee_offset_common_bins(cat, cat["area_sr"], v, w, None)
    assert out["n_common_bins"] >= 4
    assert np.isfinite(out["common_knee_offset"])
    assert out["common_knee_offset"] < 0  # the mock's void knee is injected below the wall's


def test_void_null_runs_both_variants_and_reports_fairness_diagnostics():
    cat = fe.synthetic_environment_catalogue()
    n = cat["log_mhi"].size
    rng = np.random.default_rng(51)
    # galaxies on a sky patch, at their catalogue distances (Mpc/h-ish frame is irrelevant here)
    ra, dec = rng.uniform(150, 210, n), rng.uniform(5, 45, n)
    d = np.asarray(cat["dist_mpc"], float) * 0.7
    r_, d_ = np.radians(ra), np.radians(dec)
    xyz = (
        np.column_stack([np.cos(d_) * np.cos(r_), np.cos(d_) * np.sin(r_), np.sin(d_)]) * d[:, None]
    )
    holes, vid = [], []
    for k in range(15):
        c = xyz[rng.integers(n)]
        for _ in range(2):
            holes.append(c + rng.normal(0, 2, 3))
            vid.append(k)
    voids = {
        "sphere_xyz": np.asarray(holes),
        "sphere_radius": np.full(30, 8.0),
        "void_id": np.asarray(vid),
    }
    in_void = fe.void_membership_holes(xyz, voids["sphere_xyz"], voids["sphere_radius"])
    env = {"xyz": xyz, "in_void": in_void, "classifiable": np.ones(n, bool)}
    v1 = fe.vmax_1vmax(cat["log_mhi"], cat["dist_mpc"], cat["flux"])
    wts = {"B": (np.ones(n, bool), v1), "optA_same": (np.ones(n, bool), 2 * v1)}
    for constrained in (False, True):
        out = fe.void_null(
            cat, cat["area_sr"], env, voids, (ra, dec), wts,
            n=6, seed=3, constrained=constrained, n_jackknife=1,
        )  # fmt: skip
        assert out["n_reps"] == 6 and len(out["rows"]) == 6
        assert out["B"]["n_ok"] >= 3
        row = out["rows"][0]
        for key in ("n_in_void", "n_bins_void", "real_overlap_frac", "spill_frac", "B_offset"):
            assert key in row
        assert len(out["B_jackknife_err_placements"]) <= 1
        assert "B_regression" not in out  # needs > 20 placements
        # scaling Vmax by a constant cannot move a knee: where both fits are well posed
        # (sensible offsets), the two weightings agree placement by placement
        sane = [
            r for r in out["rows"]
            if abs(r["B_offset"]) < 1.0 and abs(r["optA_same_offset"]) < 1.0
        ]  # fmt: skip
        assert sane, "no well-posed placements -- the test would be vacuous"
        for r in sane:
            assert r["B_offset"] == pytest.approx(r["optA_same_offset"], abs=1e-3)


def test_void_null_regression_block_with_enough_placements():
    cat = fe.synthetic_environment_catalogue()
    n = cat["log_mhi"].size
    rng = np.random.default_rng(52)
    ra, dec = rng.uniform(150, 210, n), rng.uniform(5, 45, n)
    d = np.asarray(cat["dist_mpc"], float) * 0.7
    r_, d_ = np.radians(ra), np.radians(dec)
    xyz = (
        np.column_stack([np.cos(d_) * np.cos(r_), np.cos(d_) * np.sin(r_), np.sin(d_)]) * d[:, None]
    )
    holes = xyz[rng.integers(n, size=20)]
    voids = {"sphere_xyz": holes, "sphere_radius": np.full(20, 25.0), "void_id": np.arange(20)}
    env = {
        "xyz": xyz,
        "in_void": fe.void_membership_holes(xyz, holes, voids["sphere_radius"]),
        "classifiable": np.ones(n, bool),
    }
    v1 = fe.vmax_1vmax(cat["log_mhi"], cat["dist_mpc"], cat["flux"])
    out = fe.void_null(
        cat, cat["area_sr"], env, voids, (ra, dec), {"B": (np.ones(n, bool), v1)},
        n=25, seed=4, constrained=False,
    )  # fmt: skip
    reg = out["B_regression"]
    assert 0.0 <= reg["r2_overlap_only"] <= reg["r2_occupancy_z_overlap"] <= 1.0
    assert reg["r2_occupancy_z"] <= reg["r2_occupancy_z_overlap"] + 1e-12
    lo, hi = reg["overlap_range_minmax"]
    assert 0.0 <= lo <= hi <= 1.0


def test_dmax_from_vmax_roundtrip():
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    d = np.array([10.0, 100.0, 350.0])
    assert np.allclose(fe.dmax_from_vmax(omega / 3.0 * d**3), d)


def test_environment_volume_fraction_inner_sphere():
    rng = np.random.default_rng(3)
    rd = 300.0 * np.cbrt(rng.uniform(0, 1, 400_000))
    g = fe.environment_volume_fraction(rd, rd < 100.0, np.array([50.0, 100.0, 200.0, 300.0]))
    assert np.allclose(g, [1.0, 1.0, (100 / 200) ** 3, (100 / 300) ** 3], atol=0.01)
    assert fe.environment_volume_fraction(rd, rd < 100.0, np.array([0.0]))[0] == 0.0


def test_survey_randoms_in_footprint_and_uniform_in_volume():
    rng = np.random.default_rng(4)
    fra = rng.uniform(150, 210, 5000)
    fdec = np.degrees(np.arcsin(rng.uniform(0, np.sin(np.radians(50)), 5000)))
    r = fe.survey_randoms(fra, fdec, 20_000, rng, d_max_mpc=300.0)
    cells = fe._sky_cells(fra, fdec, 1.0)
    keys = zip(*fe._cells_arrays(r["ra"], r["dec"], 1.0), strict=True)
    assert all(c in cells for c in keys)
    assert abs(np.median(r["d_mpc"]) - 300.0 * 0.5 ** (1 / 3)) < 3.0
    # z inverts the distance relation; the Mpc/h frame is d * h
    assert np.allclose(fe._comoving_distance_mpc(r["z"], fe.H0), r["d_mpc"], rtol=1e-6)
    assert np.allclose(np.linalg.norm(r["xyz_h"], axis=1), r["d_mpc"] * fe.H0 / 100.0, rtol=1e-6)


def _distance_split_mock(n=60_000, seed=5):
    """One Schechter HIMF everywhere; the 'environment' is simply the far half of the volume."""
    rng = np.random.default_rng(seed)
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    d = 300.0 * np.cbrt(rng.uniform(0, 1, n))
    lm = rng.uniform(7.0, 10.9, 20 * n)
    w = fe.schechter(lm, 0.0, 9.9, -1.3)
    lm = rng.choice(lm, size=n, p=w / w.sum())
    s_lim = 0.3
    flux = 10**lm / (2.356e5 * d**2)
    dmax = np.minimum(np.sqrt(10**lm / (2.356e5 * 0.3)), 300.0)  # depth-blind Vmax
    keep = flux > s_lim
    cat = {"log_mhi": lm[keep], "dist_mpc": d[keep], "flux": flux[keep]}
    vcat = omega / 3.0 * dmax[keep] ** 3
    return cat, vcat, omega, rng


def test_env_vmax_offset_removes_distance_selection_bias():
    cat, vcat, omega, rng = _distance_split_mock()
    n = cat["log_mhi"].size
    far = cat["dist_mpc"] > 150.0
    ones, comp = np.ones(n, bool), np.ones(n)
    rd = 300.0 * np.cbrt(rng.uniform(0, 1, 400_000))
    env = fe.env_vmax_offset(cat, omega, ones, ones, far, vcat, comp, rd, rd >= 0, rd > 150.0)
    vm = fe.vmax_from_catalogue(vcat, comp)
    _a, fi = fe._himf_and_fit(cat["log_mhi"], None, None, omega, mask=far, vmax=vm)
    _b, fo = fe._himf_and_fit(cat["log_mhi"], None, None, omega, mask=~far, vmax=vm)
    survey = fi["log_m_star"] - fo["log_m_star"]
    # Same HIMF in both 'environments': the survey-wide Vmax manufactures an offset, the
    # environment-restricted one does not.
    assert abs(survey) > 0.1
    assert abs(env["offset"]) < 0.05
    assert set(env) >= {"offset", "offset_err", "offset_sigma", "fit_in", "fit_out"}


def test_env_vmax_leg_runs_and_pairs_placements():
    rng = np.random.default_rng(6)
    n = 8000
    ra = rng.uniform(150, 210, n)
    dec = np.degrees(np.arcsin(rng.uniform(0, np.sin(np.radians(50)), n)))
    d = 250.0 * np.cbrt(rng.uniform(0, 1, n))
    a = 0.5 * (1.0 + fe._Q0)
    x = d * fe.H0 / fe.C_KM_S
    z = (1.0 - np.sqrt(1.0 - 4.0 * a * x)) / (2.0 * a)
    grid = rng.uniform(7.5, 10.9, 20 * n)
    w = fe.schechter(grid, 0.0, 9.9, -1.3)
    lm = rng.choice(grid, size=n, p=w / w.sum())
    cat = {"ra": ra, "dec": dec, "z": z, "cz": fe.C_KM_S * z, "log_mhi": lm, "dist_mpc": d,
           "flux": 10**lm / (2.356e5 * d**2)}  # fmt: skip
    xyz = fe.comoving_xyz(ra, dec, z, h0=100.0)
    holes = xyz[rng.choice(n, 30, replace=False)]
    voids = {"sphere_xyz": holes, "sphere_radius": np.full(30, 12.0), "void_id": np.arange(30) // 3}
    gi = rng.choice(n, 150, replace=False)
    grp = {"grp_ra": ra[gi], "grp_dec": dec[gi], "grp_cz": cat["cz"][gi],
           "grp_r200": np.full(150, 6.0), "gal_ra": ra, "gal_dec": dec}  # fmt: skip
    env = {
        "xyz": xyz,
        "in_void": fe.void_membership_holes(xyz, holes, voids["sphere_radius"]),
        "classifiable": np.ones(n, bool),
        "in_group": fe.assign_groups(ra, dec, cat["cz"], grp["grp_ra"], grp["grp_dec"],
                                     grp["grp_cz"], grp["grp_r200"], h0=fe.TEMPEL_H0) >= 0,
    }  # fmt: skip
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    vcat = omega / 3.0 * np.full(n, 250.0) ** 3
    out = fe.env_vmax_leg(
        cat, omega, env, voids, grp, np.ones(n, bool), vcat, np.ones(n), (ra, dec), (ra, dec),
        n_rand=50_000, n_void=3, n_group=3,
    )  # fmt: skip
    assert len(out["g_void_cumulative"]) == len(fe.ENV_VMAX_D_GRID)
    # dist_mpc here IS the h70 comoving distance, so the frame ratio is exactly 1 + z
    assert out["frame_correct"] and out["frame_ratio_median"] > 1.0
    assert {"void_frame_uncorrected", "group_frame_uncorrected"} <= set(out)
    for k in ("void_null_constrained", "group_null"):
        assert "excess_sigma_quadrature" in out[k] and "null_env_mean_se" in out[k]
    assert 0.0 < out["rand_classifiable_frac"] <= 1.0
    for k in ("void_null_constrained", "group_null"):
        assert out[k]["n_reps"] == 3 and len(out[k]["rows"]) == 3
        assert set(out[k]["rows"][0]) == {"survey_vmax", "env_vmax", "n_in", "median_z"}


def test_env_vmax_offset_unit_frame_ratio_is_identity():
    cat, vcat, omega, rng = _distance_split_mock(n=20_000)
    n = cat["log_mhi"].size
    far = cat["dist_mpc"] > 150.0
    ones, comp = np.ones(n, bool), np.ones(n)
    rd = 300.0 * np.cbrt(rng.uniform(0, 1, 200_000))
    args = (cat, omega, ones, ones, far, vcat, comp, rd, rd >= 0, rd > 150.0)
    a = fe.env_vmax_offset(*args)
    b = fe.env_vmax_offset(*args, frame_ratio=np.ones(n))
    c = fe.env_vmax_offset(*args, frame_ratio=np.full(n, 0.5))
    assert a["offset"] == b["offset"]
    assert c["offset"] != a["offset"]  # the lookup distance matters


def test_label_shuffle_null_keeps_redshift_occupancy_and_flags_mass_difference():
    cat, vcat, omega, rng = _distance_split_mock(n=300_000, seed=8)
    n = cat["log_mhi"].size
    z = cat["dist_mpc"] * fe.H0 / fe.C_KM_S
    vm = fe.vmax_from_catalogue(vcat, np.ones(n))
    pool = np.ones(n, bool)
    rand_in = rng.uniform(size=n) < 0.3
    offs = fe.label_shuffle_null(cat["log_mhi"], z, pool, rand_in, vm, vm, omega, rng, n=15)
    assert np.all(np.isfinite(offs)) and abs(np.mean(offs)) < 0.05
    # members biased to high mass at fixed z: the measured offset sits far outside the null
    rank = np.argsort(np.argsort(cat["log_mhi"] - 0.5 * np.log10(cat["dist_mpc"] ** 2)))
    heavy = rank > 0.8 * n
    _a, fi = fe._himf_and_fit(cat["log_mhi"], None, None, omega, mask=heavy, vmax=vm)
    _b, fo = fe._himf_and_fit(cat["log_mhi"], None, None, omega, mask=~heavy, vmax=vm)
    null = fe.label_shuffle_null(cat["log_mhi"], z, pool, heavy, vm, vm, omega, rng, n=30)
    measured = fi["log_m_star"] - fo["log_m_star"]
    assert measured > np.max(null)
    assert measured - np.mean(null) > 3 * np.std(null)


def test_confusion_counts_beam_and_velocity_window():
    ra, dec = np.array([180.0]), np.array([30.0])
    cz, w50 = np.array([5000.0]), np.array([200.0])
    d = 1.0 / 60.0  # 1 arcmin in Dec
    opt_ra = np.array([180.0, 180.0, 180.0, 180.0])
    opt_dec = np.array([30.0, 30.0 + d, 30.0 + 2 * d, 30.0 + d])
    opt_cz = np.array([5000.0, 5150.0, 5000.0, 5300.0])
    # counterpart + a neighbour at 1' inside W50/2+100 = 200 km/s; the 2' one is inside one
    # FWHM (2.9') too; the 300 km/s one is outside the line window
    n = fe.confusion_counts(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz)
    assert n.tolist() == [3]
    assert fe.confusion_counts(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz, beam_arcmin=1.5)[0] == 2
    assert fe.confusion_counts(ra + 10, dec, cz, w50, opt_ra, opt_dec, opt_cz)[0] == 0
    # a fixed +/-350 km/s window also admits the 300 km/s neighbour
    kw = {"half_window_kms": 350.0}
    assert fe.confusion_counts(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz, **kw)[0] == 4


def test_robustness_leg_runs_on_mock():
    rng = np.random.default_rng(9)
    n = 8000
    ra = rng.uniform(150, 210, n)
    dec = np.degrees(np.arcsin(rng.uniform(0, np.sin(np.radians(50)), n)))
    d = 250.0 * np.cbrt(rng.uniform(0, 1, n))
    a = 0.5 * (1.0 + fe._Q0)
    z = (1.0 - np.sqrt(1.0 - 4.0 * a * d * fe.H0 / fe.C_KM_S)) / (2.0 * a)
    grid = rng.uniform(7.5, 10.9, 20 * n)
    w = fe.schechter(grid, 0.0, 9.9, -1.3)
    lm = rng.choice(grid, size=n, p=w / w.sum())
    cz = fe.C_KM_S * z
    cat = {"ra": ra, "dec": dec, "z": z, "cz": cz, "log_mhi": lm, "dist_mpc": d,
           "flux": 10**lm / (2.356e5 * d**2), "w50": 10 ** (0.25 * (lm - 9.0) + 2.3),
           "w20": 1.2 * 10 ** (0.25 * (lm - 9.0) + 2.3), "rms": rng.uniform(0.5, 2.5, n),
           "rms_beam": rng.uniform(0.3, 1.5, n), "ell_maj": rng.uniform(2.9, 5.0, n)}  # fmt: skip
    xyz = fe.comoving_xyz(ra, dec, z, h0=100.0)
    holes = xyz[rng.choice(n, 30, replace=False)]
    voids = {"sphere_xyz": holes, "sphere_radius": np.full(30, 12.0), "void_id": np.arange(30) // 3}
    gi = rng.choice(n, 150, replace=False)
    # optical catalogue = the galaxies themselves plus a close companion for every 10th one
    comp_i = np.arange(0, n, 10)
    grp = {"grp_ra": ra[gi], "grp_dec": dec[gi], "grp_cz": cz[gi], "grp_r200": np.full(150, 6.0),
           "gal_ra": np.concatenate([ra, ra[comp_i]]),
           "gal_dec": np.concatenate([dec, dec[comp_i] + 0.5 / 60.0]),
           "gal_cz": np.concatenate([cz, cz[comp_i] + 50.0]),
           "gal_group_id": np.concatenate([np.arange(n), comp_i]),
           "gal_ngal": np.full(n + comp_i.size, 2),
           "gal_rmag_abs": rng.uniform(-22, -17, n + comp_i.size)}  # fmt: skip
    env = {
        "xyz": xyz,
        "in_void": fe.void_membership_holes(xyz, holes, voids["sphere_radius"]),
        "classifiable": np.ones(n, bool),
        "in_group": fe.assign_groups(ra, dec, cz, grp["grp_ra"], grp["grp_dec"], grp["grp_cz"],
                                     grp["grp_r200"], h0=fe.TEMPEL_H0) >= 0,
    }  # fmt: skip
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    vcat = omega / 3.0 * np.full(n, 250.0) ** 3
    meas = {k: {"survey": 0.0, "survey_err": 0.03, "env": 0.0, "env_err": 0.03}
            for k in ("void", "group")}  # fmt: skip
    out = fe.robustness_leg(
        cat, omega, env, voids, grp, np.ones(n, bool), vcat, np.ones(n), meas,
        n_rand=50_000, n_shuffle=4,
    )  # fmt: skip
    for name in ("void", "group"):
        for wk in ("survey", "env"):
            s = out["shuffle"][name][wk]
            assert s["n_ok"] == 4 and np.isfinite(s["mean"])
        assert out["shuffle_z_only"][name]["env"]["n_ok"] == 2  # comparison runs use n // 2
        assert out["shuffle_rms_beam_terciles"][name]["env"]["n_ok"] == 2
        for fk in ("w50", "fixed300"):
            bl = out["blending"][name][fk]
            # roughly every 10th source has a planted companion
            assert 0.05 < bl["flagged_frac_out"] < 0.2
            assert np.isfinite(bl["survey_offset_unflagged"])
            assert bl["survey_offset_unflagged_err"] > 0
            assert "n_ok" in bl["env_shuffle_unflagged"]
        b = out["blending"][name]
        # every planted companion shares its counterpart's group id
        assert b["flagged_same_group_frac"] is None or b["flagged_same_group_frac"] > 0.9
        lw = b["linewidth"]
        assert set(lw["classes"]) == {"inner", "ring", "isolated", "no_counterpart"}
        assert "out_inner_matched_fine" in lw and "in_inner_dvlt100_w20w50" in lw
        assert "in_inner_logm_at_fixed_rmag" in b["hi_at_fixed_optical"]
        assert set(b["common_alpha"]) == {"alpha", "survey", "env"}
        assert set(out["footprint_strict"][name]) == {"survey", "env"}
    assert "median_flagged" in out["blending"]["group"]["r_over_r200"]
    assert "n_cells" in out["blending"]["group"]["segregation_unblendable"]


def test_measured_offsets_reads_both_weightings():
    m = {"void_knee_offset": -0.15, "void_knee_offset_err": 0.04, "group_knee_offset": 0.17,
         "group_knee_offset_err": 0.04,
         "env_vmax": {"void": {"offset": -0.09, "offset_err": 0.035},
                      "group": {"offset": 0.09, "offset_err": 0.032}}}  # fmt: skip
    got = fe.measured_offsets(m)
    assert got["void"] == {"survey": -0.15, "survey_err": 0.04, "env": -0.09, "env_err": 0.035}
    assert got["group"]["env"] == 0.09
    assert fe.measured_offsets({**m, "env_vmax": {}})["group"]["env"] is None


def test_robustness_macros_values_and_placeholders():
    rob = {
        "n_shuffle": 1000,
        "shuffle": {"void": {"env": {"mean": 0.027, "std": 0.014, "excess": -0.119,
                                     "excess_sigma_quadrature": 3.16, "n_reaching_measured": 0}}},
        "blending": {"n_pool": 55893,
                     "group": {"env_offset_all": 0.123, "flagged_same_group_frac": 0.8,
                               "w50": {"flagged_frac_in": 0.2553, "flagged_frac_out": 0.0155,
                                       "env_offset_unflagged": 0.038},
                               "linewidth": {"in_inner_matched": {"median_diff": 0.021, "se": 0.005,
                                                                  "sigma": 4.2,
                                                                  "n_flagged_matched": 900}}}},
    }  # fmt: skip
    text = "\n".join(fe._robustness_macros(rob))
    assert r"\feRealShufVoidEnvExcess}{-0.119}" in text
    assert r"\feRealShufVoidEnvSigma}{3.2}" in text
    assert r"\feRealBlendGroupWfiftyFlagIn}{25.5}" in text
    assert r"\feRealBlendGroupWfiftyFlagOut}{1.6}" in text
    assert r"\feRealBlendGroupSameGroupPct}{80.0}" in text
    assert r"\feRealBlendGroupLwInInnerDiff}{0.021}" in text
    assert r"\feRealBlendGroupLwInInnerN}{900}" in text
    assert r"\feRealShufGroupSurveyMean}{--}" in text  # absent -> placeholder, never a crash
    assert r"\feRealBlendNPool}{55893}" in text


def test_label_shuffle_depth_strata_remove_a_pure_depth_confound():
    """Members sit in a shallow region with NO intrinsic mass difference, and the weights use the
    survey-wide flux limit (a Vmax that does not know about local depth -- the case that
    matters). A z-only shuffle is fooled (it scrambles depth into the null); a (z, depth)
    shuffle is not."""
    rng = np.random.default_rng(11)
    n = 300_000
    omega = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    d = 300.0 * np.cbrt(rng.uniform(0, 1, n))
    grid = rng.uniform(7.0, 10.9, 20 * n)
    w = fe.schechter(grid, 0.0, 9.9, -1.3)
    lm = rng.choice(grid, size=n, p=w / w.sum())
    shallow = rng.uniform(size=n) < 0.3  # the "environment" is simply a shallower sky region
    s_lim = np.where(shallow, 0.9, 0.3)  # 3x worse flux limit there
    keep = 10**lm / (2.356e5 * d**2) > s_lim
    lm, d, shallow, s_lim = lm[keep], d[keep], shallow[keep], s_lim[keep]
    dmax = np.minimum(np.sqrt(10**lm / (2.356e5 * 0.3)), 300.0)  # depth-blind Vmax
    vm = fe.vmax_from_catalogue(omega / 3.0 * dmax**3, np.ones(lm.size))
    z = d * fe.H0 / fe.C_KM_S
    pool = np.ones(lm.size, bool)
    _a, fi = fe._himf_and_fit(lm, None, None, omega, mask=shallow, vmax=vm)
    _b, fo = fe._himf_and_fit(lm, None, None, omega, mask=~shallow, vmax=vm)
    measured = fi["log_m_star"] - fo["log_m_star"]
    z_only = fe.label_shuffle_null(lm, z, pool, shallow, vm, vm, omega, rng, n=20)
    depth = fe.label_shuffle_null(
        lm, z, pool, shallow, vm, vm, omega, rng, n=20, strata=shallow.astype(int)
    )
    # z-only: the measured offset looks significant; (z, depth): it is inside the null
    assert abs(measured - z_only.mean()) > 3 * z_only.std()
    assert abs(measured - depth.mean()) < 3 * max(depth.std(), 1e-3)


def test_linewidth_residual_flags_broadened_profiles_at_fixed_mass():
    rng = np.random.default_rng(12)
    n = 5000
    lm = rng.uniform(8, 10.5, n)
    z = rng.uniform(0.01, 0.05, n)
    lw = 2.3 + 0.25 * (lm - 9.0) + rng.normal(0, 0.05, n)
    blended = rng.uniform(size=n) < 0.2
    lw_b = lw + np.where(blended, 0.1, 0.0)  # blending broadens by 0.1 dex
    res = fe.linewidth_residual(lw_b, lm, z, ~blended)
    assert abs(np.median(res[~blended])) < 0.01
    assert abs(np.median(res[blended]) - 0.1) < 0.01
    assert np.isnan(fe.linewidth_residual(np.array([np.nan, 2.0, 2.1, 2.2]), np.array([9.0, 9.0, 9.5, 10.0]),
                                          np.array([0.02] * 4), np.ones(4, bool))[0])  # fmt: skip


def test_confused_same_group_needs_two_members_of_one_group():
    ra, dec = np.array([180.0, 180.0]), np.array([30.0, 40.0])
    cz, w50 = np.array([5000.0, 5000.0]), np.array([200.0, 200.0])
    d = 1.0 / 60.0
    opt_ra = np.array([180.0, 180.0, 180.0, 180.0])
    opt_dec = np.array([30.0, 30.0 + d, 40.0, 40.0 + d])
    opt_cz = np.array([5000.0, 5100.0, 5000.0, 5100.0])
    gid = np.array([7, 7, 8, 9])  # source 0: same group; source 1: two different groups
    ngal = np.array([3, 3, 1, 1])
    got = fe.confused_same_group(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz, gid, ngal)
    assert got.tolist() == [True, False]


def test_matched_median_diff_removes_a_mass_mismatch():
    rng = np.random.default_rng(13)
    n = 20_000
    lm = rng.uniform(8.0, 10.5, n)
    z = rng.uniform(0.01, 0.05, n)
    x = 0.25 * lm + rng.normal(0, 0.05, n)  # x rises with mass; no class effect at all
    flagged = rng.uniform(size=n) < 1 / (1 + np.exp(-(lm - 9.8) * 4))  # flagged are massive
    naive = np.median(x[flagged]) - np.median(x[~flagged])
    m = fe.matched_median_diff(x, flagged, ~flagged, lm, z)
    assert naive > 0.1  # the raw comparison is dominated by the mass mismatch
    assert abs(m["median_diff"]) < 3 * m["se"] and m["n_cells"] > 5
    assert fe.matched_median_diff(x, np.zeros(n, bool), ~flagged, lm, z)["n_cells"] == 0


def test_confusion_counts_annulus_and_per_source_radius():
    ra, dec = np.array([180.0, 180.0]), np.array([30.0, 30.0])
    cz, w50 = np.array([5000.0, 5000.0]), np.array([200.0, 200.0])
    d = 1.0 / 60.0
    opt_ra = np.full(3, 180.0)
    opt_dec = np.array([30.0, 30.0 + 2 * d, 30.0 + 4 * d])  # counterpart, 2', 4'
    opt_cz = np.full(3, 5000.0)
    kw = {"half_window_kms": 300.0}
    # per-source radii: 3' sees two, 5' sees three
    got = fe.confusion_counts(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz,
                              beam_arcmin=np.array([3.0, 5.0]), **kw)  # fmt: skip
    assert got.tolist() == [2, 3]
    # an annulus 2.5-5' excludes the counterpart and the 2' neighbour
    ring = fe.confusion_counts(ra, dec, cz, w50, opt_ra, opt_dec, opt_cz, beam_arcmin=5.0,
                               inner_arcmin=2.5, **kw)  # fmt: skip
    assert ring.tolist() == [1, 1]


def test_shuffle_on_an_outcome_tracking_stratifier_absorbs_a_real_offset():
    """Seventh round: DR2's rms scales with source size, hence with mass. A true mass offset
    between members and the rest survives a shuffle stratified on a mass-independent depth, and
    is absorbed by one stratified on a mass-tracking 'depth'."""
    cat, vcat, omega, rng = _distance_split_mock(n=300_000, seed=14)
    lm = cat["log_mhi"]
    n = lm.size
    z = cat["dist_mpc"] * fe.H0 / fe.C_KM_S
    vm = fe.vmax_from_catalogue(vcat, np.ones(n))
    rank = np.argsort(np.argsort(lm - np.log10(cat["dist_mpc"])))
    members = rank > 0.8 * n  # members truly more massive at fixed distance
    _a, fi = fe._himf_and_fit(lm, None, None, omega, mask=members, vmax=vm)
    _b, fo = fe._himf_and_fit(lm, None, None, omega, mask=~members, vmax=vm)
    measured = fi["log_m_star"] - fo["log_m_star"]
    pool = np.ones(n, bool)
    honest = rng.integers(0, 3, n)  # depth independent of mass
    tracking = np.digitize(lm - np.log10(cat["dist_mpc"]),
                           np.quantile(lm - np.log10(cat["dist_mpc"]), [0.2, 0.4, 0.6, 0.8]))  # fmt: skip
    null_h = fe.label_shuffle_null(lm, z, pool, members, vm, vm, omega, rng, n=20, strata=honest)
    null_t = fe.label_shuffle_null(lm, z, pool, members, vm, vm, omega, rng, n=20, strata=tracking)
    assert measured - null_h.mean() > 2 * null_h.std()
    assert abs(measured - null_t.mean()) < abs(measured - null_h.mean()) / 2


def test_strict_max_shift():
    m = {"void_knee_offset": -0.171, "group_knee_offset": 0.2,
         "env_vmax": {"void": {"offset": -0.088}, "group": {"offset": 0.123}},
         "robustness": {"footprint_strict": {
             "void": {"survey": {"offset": -0.172}, "env": {"offset": -0.088}},
             "group": {"survey": {"offset": 0.2}, "env": {"offset": 0.124}}}}}  # fmt: skip
    assert fe._strict_max_shift(m) == "0.001"
    assert fe._strict_max_shift({}) == "--"


def test_optical_classes_are_disjoint_and_match_the_counterpart():
    d = 1.0 / 60.0
    # four sources at Dec 30, 35, 40, 45 (well separated)
    ra = np.full(4, 180.0)
    dec = np.array([30.0, 35.0, 40.0, 45.0])
    cz = np.full(4, 5000.0)
    opt_ra = np.full(6, 180.0)
    opt_dec = np.array(
        [
            30.0,
            30.0 + 2 * d,  # src0: counterpart + inner neighbour (2')
            35.0,
            35.0 + 6 * d,  # src1: counterpart + ring neighbour (6')
            40.0,  # src2: counterpart only -> isolated
            45.0 + 3 * d,
        ]
    )  # src3: no counterpart within 1.45'
    opt_cz = np.array([5000.0, 5050.0, 5000.0, 5000.0, 5000.0, 5000.0])
    got = fe.optical_classes(ra, dec, cz, opt_ra, opt_dec, opt_cz, np.full(4, 4.0))
    oc = fe.OPT_CLASS
    assert got["cls"].tolist() == [oc["inner"], oc["ring"], oc["isolated"], oc["no_counterpart"]]
    assert got["counterpart"].tolist() == [0, 2, 4, -1]
    assert got["min_dv_inner"][0] == 50.0 and np.isnan(got["min_dv_inner"][1])


def test_fit_schechter_fixed_alpha_recovers_the_knee():
    lm = np.arange(7.0, 11.0, 0.25) + 0.125
    phi = fe.schechter(lm, -2.5, 9.9, -1.3)
    h = {"logm": lm, "phi": phi, "phi_err": 0.05 * phi, "counts": np.full(lm.size, 50)}
    f = fe.fit_schechter(h, alpha_fixed=-1.3)
    assert abs(f["log_m_star"] - 9.9) < 1e-3 and f["alpha"] == -1.3
    g = fe.fit_schechter(h, alpha_fixed=-1.0)  # a wrong slope moves the knee
    assert abs(g["log_m_star"] - 9.9) > 0.01
    assert np.isnan(
        fe.fit_schechter({**h, "counts": np.zeros(lm.size)}, alpha_fixed=-1.3)["log_m_star"]
    )


def test_placement_occupancy_macros_fit_the_trend():
    n_in = np.linspace(1500, 2300, 50)
    rows = [{"n_in": float(x), "env_vmax": 0.0002 * (x - 1500) - 0.1} for x in n_in]
    m = {"env_vmax": {"group_null": {"rows": rows}, "void_null_constrained": {"rows": []}}}
    text = "\n".join(fe._placement_occupancy_macros(m))
    assert r"\feRealEnvPlaceGroupSlopeHundred}{0.020}" in text
    assert r"\feRealEnvPlaceGroupMeanAtLo}{-0.100}" in text
    assert r"\feRealEnvPlaceGroupMeanAtHi}{0.060}" in text
    assert r"\feRealEnvPlaceVoidOccLo}{--}" in text
