"""Tests for jansky_research.type3synthesis -- the corona->0.4 AU type III synthesis. No network."""

from __future__ import annotations

import numpy as np
import pytest

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
    # Namespaced: an offline run fills the SYNTHETIC side only. Before this, its ladder
    # (corona speed 0.3002) overwrote the real one (0.1347) under shared names, and this
    # slice emitted no provenance macro so neither merge guard could even see the mode.
    assert r"\synSynGeomCorr" in macros and r"\synSynOverallRhiAU" in macros
    assert r"\newcommand{\synRealGeomCorr}{--}" in macros
    assert r"\synSource" in macros, "provenance marker missing; both merge guards are blind"


def test_model_curves_monotone():
    c = syn._model_curves()
    # both density models give radius decreasing with frequency (higher freq -> deeper -> smaller r)
    assert c["r_corona"][0] > c["r_corona"][-1]  # f_corona ascending -> r descending
    assert c["r_helio"][0] > c["r_helio"][-1]


def test_committed_ladder_matches_the_committed_siblings():
    """The round-8 blocker: the synthesis shipped pre-revision sibling vintages for five
    per-leg values, and no guard could see a stale-but-real file. This one can: every per-leg
    value in the committed synthesis JSON must equal the corresponding key in the sibling's
    own committed JSON. Runs only on the committed real artifacts."""
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    syn_path = root / "results" / "type3synthesis_metrics.json"
    if not syn_path.exists():
        pytest.skip("no committed synthesis metrics")
    syn = json.loads(syn_path.read_text())
    if str(syn.get("source", "")).lower().startswith("synthetic"):
        pytest.skip("committed synthesis metrics are synthetic")
    pairs = [
        ("corona_speed_c", "solarbursts", "speed_c"),
        ("corona_r_lo", "solarbursts", "r_lo_rsun"),
        ("corona_r_hi", "solarbursts", "r_hi_rsun"),
        ("helio_speed_c", "windwaves", "speed_c"),
        ("helio_r_hi", "windwaves", "r_hi_rsun"),
        ("ip_speed_c", "swaves", "speed_c"),
        ("ip_r_hi_rsun", "swaves", "r_hi_rsun"),
        ("geom_r_hi_rsun", "triangulate", "r_hi_rsun"),
        ("geom_ratio", "triangulate", "ratio_geom_plasma"),
    ]
    for syn_key, sibling, sib_key in pairs:
        sib_path = root / "results" / f"{sibling}_metrics.json"
        if not sib_path.exists():
            continue
        sib = json.loads(sib_path.read_text())
        if str(sib.get("source", "")).lower().startswith("synthetic"):
            continue
        assert syn.get(syn_key) == sib.get(sib_key), (
            f"{syn_key} = {syn.get(syn_key)} is stale: {sibling}.{sib_key} = {sib.get(sib_key)}"
        )


def test_offline_crosscheck_loglog_slope_can_fail():
    """The corr > 0.8 assertion clears for any monotone curve; the log-log slope of r_geom on
    r_plasma is the statistic that can fail (it IS ~0.65 on the real data, the additive-offset
    signature). On a low-noise zero-bias fixture it must be near 1 — and direction noise alone
    must flatten it below 1, the mechanism behind the real value."""
    from jansky_research import triangulate

    quiet = triangulate.synthetic_event(noise_deg=2.0, seed=3)
    tq = triangulate.triangulate_track(quiet["spec_a"], quiet["spec_b"])
    slope_q = float(np.polyfit(np.log10(tq["r_plasma"]), np.log10(tq["r_geom"]), 1)[0])
    assert 0.9 < slope_q < 1.1
    noisy = syn.crosscheck_track(offline=True)  # default 9-degree scatter
    slope_n = float(np.polyfit(np.log10(noisy["r_plasma"]), np.log10(noisy["r_geom"]), 1)[0])
    assert slope_n < slope_q  # noise flattens the slope (the outward additive bias)
