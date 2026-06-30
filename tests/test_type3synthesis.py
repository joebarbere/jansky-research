"""Tests for jansky_research.type3synthesis -- the corona->0.4 AU type III synthesis. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import type3synthesis as syn


def test_crosscheck_track_offline_correlates():
    track = syn.crosscheck_track(offline=True)
    rg = np.asarray(track["r_geom"], float)
    rp = np.asarray(track["r_plasma"], float)
    assert rg.size >= 8
    # the geometric and plasma-frequency distances track each other (the centrepiece claim)
    assert np.corrcoef(rg, rp)[0, 1] > 0.8


def test_run_offline_spans_corona_to_interplanetary(tmp_path):
    m = syn.run(out=str(tmp_path), offline=True)
    assert m["n_instruments"] == 4
    # the ladder spans from the corona (few R_sun) to genuinely interplanetary (>0.1 AU)
    assert m["corona_r_lo"] < 5.0
    assert m["overall_r_hi_au"] > 0.1
    # the high-frequency (corona) and low-frequency (interplanetary) ends are present
    assert m["f_hi_mhz"] > 20.0 and m["f_lo_mhz"] < 0.5
    # the geometric cross-check correlation is reported
    assert m["geom_corr"] is not None and m["geom_corr"] > 0.8
    # outputs land where the paper expects them
    assert (tmp_path / "results" / "type3synthesis_metrics.json").exists()
    assert (tmp_path / "papers" / "type3synthesis" / "figures" / "type3synthesis.pdf").exists()
    macros = (tmp_path / "papers" / "type3synthesis" / "generated" / "macros.tex").read_text()
    assert r"\synGeomCorr" in macros and r"\synOverallRhiAU" in macros


def test_model_curves_monotone():
    c = syn._model_curves()
    # both density models give radius decreasing with frequency (higher freq -> deeper -> smaller r)
    assert c["r_corona"][0] > c["r_corona"][-1]  # f_corona ascending -> r descending
    assert c["r_helio"][0] > c["r_helio"][-1]
