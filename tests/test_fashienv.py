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
    assert jk and jk.get("jackknife_err") is not None, "the paper quotes this; commit it"
    assert 0 < jk["n_occupied"] <= jk["n_voids"]
    assert jk["n_ok"] == jk["n_occupied"]
    # the jackknife mean must sit on the full-sample offset, or the resampling is wrong
    assert jk["mean_offset"] == pytest.approx(m["void_knee_offset"], abs=0.01)
    assert jk["min_offset"] <= m["void_knee_offset"] <= jk["max_offset"]


def test_comparison_bin_size_is_committed():
    """The "wall" bin is a bounding box, not a footprint, and holds most of the catalogue.

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
    # if the wall bin ever stops matching the global fit, the "void versus everything"
    # framing in the paper needs revisiting
    assert abs(m["himf_wall"]["log_m_star"] - m["himf_global"]["log_m_star"]) < 0.05


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
    assert r"\newcommand{\feRealOptAVoidKneeOffset}{-0.25}" in t
    assert r"\newcommand{\feRealOptAGlobalAlpha}{-1.78}" in t
    assert r"\newcommand{\feRealDROneN}{41741}" in t
    assert r"\newcommand{\feRealEdsVoidKneeOffset}{--}" in t  # absent key -> placeholder
    # a DR1 or synthetic run never emits the DR2 block
    fe._write_macros({**m, "release": "DR1"}, tmp_path / "m1.tex")
    assert "feRealNDRTwo" not in (tmp_path / "m1.tex").read_text()
