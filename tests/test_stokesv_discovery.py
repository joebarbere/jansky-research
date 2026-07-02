"""Tests for jansky_research.stokesv_discovery -- epoch-pair V selection + variability. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import stokesv_discovery as svd


def test_epoch_position_propagates_gaia_pm():
    # 1000 mas/yr in pmra* over ~10 yr at dec=0 -> 10 arcsec of RA coordinate motion
    obs_mjd = 51544.5 + 365.25 * 26.0  # 2026.0
    ra, dec = svd.epoch_position(
        np.array([100.0]), np.array([0.0]), np.array([1000.0]), np.array([0.0]), 2016.0, obs_mjd
    )
    assert abs((ra[0] - 100.0) * 3600.0 - 10.0) < 0.05
    assert dec[0] == 0.0
    # at dec=60 the same pmra* gives a LARGER RA coordinate shift (1/cos d)
    ra60, _ = svd.epoch_position(
        np.array([100.0]), np.array([60.0]), np.array([1000.0]), np.array([0.0]), 2016.0, obs_mjd
    )
    assert (ra60[0] - 100.0) > (ra[0] - 100.0)


def test_two_epoch_variability_flags():
    v1 = np.array([10.0, 0.1, 5.0, 0.0])
    v2 = np.array([10.1, 8.0, -5.0, 0.1])
    e = np.full(4, 0.25)
    sel1 = np.array([True, False, True, False])
    sel2 = np.array([True, True, True, False])
    var = svd.two_epoch_variability(v1, e, v2, e, sel1, sel2)
    assert not var["variable"][0]  # steady emitter: tiny dv
    assert var["appeared"][1] and var["variable"][1]  # flare turning on
    assert var["variable"][2]  # handedness flip = strong signed change
    assert not (var["appeared"][3] or var["variable"][3])  # never selected


def test_synthetic_pair_selection_and_variability_recover_truth():
    s = svd.synthetic_epoch_pair(n_stars=800, seed=3)
    sel1, floor1 = svd.select_epoch_candidates(s["i_flux"], s["v1"], s["e_i"], s["e_v"])
    sel2, floor2 = svd.select_epoch_candidates(s["i_flux"], s["v2"], s["e_i"], s["e_v"])
    assert 0.0 < floor1 < 0.2 and 0.0 < floor2 < 0.2
    emitters = s["is_steady"] | s["is_flare"]
    union = sel1 | sel2
    assert union[emitters].mean() > 0.9  # completeness of the epoch union
    assert emitters[union].mean() > 0.9  # purity vs leakage contaminants
    # flares are selected in exactly one epoch -> flagged by appear/disappear
    var = svd.two_epoch_variability(s["v1"], s["e_v"], s["v2"], s["e_v"], sel1, sel2)
    flagged = var["appeared"] | var["disappeared"]
    assert flagged[s["is_flare"]].mean() > 0.9
    assert not flagged[s["is_steady"]].any()


def test_single_epoch_misses_flares_the_pair_catches():
    # the slice's raison d'etre: one epoch undercounts the emitter population
    s = svd.synthetic_epoch_pair(n_stars=800, seed=5)
    sel1, _ = svd.select_epoch_candidates(s["i_flux"], s["v1"], s["e_i"], s["e_v"])
    sel2, _ = svd.select_epoch_candidates(s["i_flux"], s["v2"], s["e_i"], s["e_v"])
    emitters = s["is_steady"] | s["is_flare"]
    assert sel1[emitters].mean() < (sel1 | sel2)[emitters].mean()


def test_run_offline_writes_artifacts(tmp_path):
    m = svd.run(str(tmp_path), offline=True)
    assert m["n_targets"] == 600
    assert m["completeness"] > 0.9 and m["purity"] > 0.9
    assert m["var_completeness"] > 0.9
    assert 0.0 < m["single_epoch_miss_frac"] < 1.0
    saved = json.loads((tmp_path / "results" / "stokesv_discovery_metrics.json").read_text())
    assert saved == m
    fig = tmp_path / "papers" / "stokesv_discovery" / "figures" / "stokesv_discovery.pdf"
    assert fig.exists() and fig.stat().st_size > 0
    macros = (tmp_path / "papers" / "stokesv_discovery" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\svdComp}" in macros and r"\newcommand{\svdMissFrac}" in macros


def test_write_macros_placeholder(tmp_path):
    path = tmp_path / "macros.tex"
    svd._write_macros({"source": "x", "completeness": None}, path)
    assert r"\newcommand{\svdComp}{--}" in path.read_text()
