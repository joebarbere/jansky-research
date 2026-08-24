"""Tests for jansky_research.junodam -- DAM occurrence census. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import junodam as jdm


def test_io_mean_longitude_rate_and_range():
    jd = jdm.J2000_JD + np.array([0.0, 1.0])
    lon = jdm.io_mean_longitude(jd)
    assert 0 <= lon[0] < 360 and abs(lon[0] - 106.07719) < 1e-6
    # one day advances by the Lieske rate mod 360
    assert abs(((lon[1] - lon[0]) % 360.0) - (203.488955790 % 360.0)) < 1e-6


def test_io_phase_wraps():
    jd = np.array([jdm.J2000_JD])
    assert 0 <= jdm.io_phase(jd, np.array([350.0]))[0] < 360


def test_detect_active_threshold():
    af = np.array([0.0, 0.05, 0.1, 0.5])
    assert jdm.detect_active(af).tolist() == [False, False, True, True]


def test_sensitivity_corrected_active_flattens_pure_1r2_trend():
    # a range-INDEPENDENT intrinsic emitter seen through 1/r^2: raw detection rises near-in
    # (snr clears the floor more often at small range), but distance-correction removes it
    rng = np.random.default_rng(0)
    n = 4000
    dist = np.linspace(0.02, 0.6, n)  # AU, perijove-ish to far
    rng.shuffle(dist)
    ref = np.median(dist)
    intrinsic_snr = 10.0 ** rng.normal(-0.3, 0.15, n)  # range-independent, straddles 1
    snr_obs = intrinsic_snr * (ref / dist) ** 2  # observed: brighter (higher snr) when closer
    raw_active = snr_obs >= 1.0
    corr_active = jdm.sensitivity_corrected_active(snr_obs, dist)
    near = dist <= np.quantile(dist, 0.25)
    far = dist > np.quantile(dist, 0.75)
    raw_ratio = raw_active[near].mean() / max(raw_active[far].mean(), 1e-9)
    corr_ratio = corr_active[near].mean() / max(corr_active[far].mean(), 1e-9)
    assert raw_ratio > 3.0  # 1/r^2 sensitivity manufactures a strong proximity trend
    assert corr_ratio < 1.5  # the null collapses it to ~flat (intrinsic is range-independent)


def test_occurrence_map_masks_low_exposure():
    cml = np.array([10.0] * 5 + [200.0])
    pha = np.array([10.0] * 5 + [200.0])
    act = np.array([True, True, False, False, False, True])
    m = jdm.occurrence_map(cml, pha, act, n_bins=18, min_exposure=3)
    assert np.isfinite(m["occ"][0, 0]) and abs(m["occ"][0, 0] - 0.4) < 1e-9
    assert np.isnan(m["occ"][10, 10])  # single visit -> masked


def test_synthetic_orbit_recovers_io_contrast():
    s = jdm.synthetic_orbit(seed=1)
    m = jdm.occurrence_map(s["cml"], s["io_phase"], s["active"])
    con = jdm.io_region_contrast(m)
    expected = s["p_in"] / s["p_out"]
    assert con["contrast"] > 4.0  # strong recovered Io-region enhancement
    assert con["contrast"] < 1.5 * expected  # and not inflated
    assert con["cells_used"] > 200  # a month covers most of the plane


def test_run_offline_writes_artifacts(tmp_path):
    m = jdm.run(str(tmp_path), offline=True)
    assert m["io_contrast"] and m["io_contrast"] > 4.0
    saved = json.loads((tmp_path / "results" / "junodam_metrics.json").read_text())
    assert saved == m
    assert (tmp_path / "papers" / "junodam" / "figures" / "junodam.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "junodam" / "generated" / "macros.tex").read_text()
    # An offline run fills the SYNTHETIC namespace only. These macros were once shared
    # between modes, so an offline rebuild wrote the synthetic recovery (~7) into the macro
    # the paper uses for the real measured contrast (1.12).
    assert r"\newcommand{\jdSynContrast}" in macros
    assert r"\newcommand{\jdRealContrast}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    jdm._write_macros({"source": "x", "io_contrast": None}, p)
    assert r"\newcommand{\jdRealContrast}{--}" in p.read_text()


def test_offline_rebuild_cannot_clobber_the_real_contrast(tmp_path):
    """The mode-dependent macros are namespaced, so `make figures` cannot overwrite the
    measured Io-region contrast with the synthetic recovery of an injected one."""
    p = tmp_path / "m.tex"
    jdm._write_macros({"source": "Juno/Waves v02 (real)", "io_contrast": 1.12}, p)
    assert r"\newcommand{\jdRealContrast}{1.12}" in p.read_text()
    jdm._write_macros(
        {"source": "synthetic orbit", "io_contrast": 6.95, "expected_contrast": 8.75}, p
    )
    text = p.read_text()
    assert r"\newcommand{\jdRealContrast}{1.12}" in text, "offline run clobbered the real value"
    assert r"\newcommand{\jdSynContrast}{6.95}" in text
    assert r"\newcommand{\jdSynExpContrast}{8.75}" in text


def test_censored_census_never_promotes_and_only_removes():
    rng = np.random.default_rng(0)
    snr = 10.0 ** rng.normal(0.0, 0.5, 5000)
    dist = np.linspace(0.001, 0.02, 5000)
    raw = snr >= 1.0
    cens = jdm.sensitivity_censored_active(snr, dist)
    assert not (cens & ~raw).any()  # nothing sub-threshold is promoted
    corr = jdm.sensitivity_corrected_active(snr, dist)
    assert (corr & ~raw).any()  # ...which the upward-rescale variant does do (the defect)


def test_monthly_contrast_test_matches_hand_arithmetic():
    # the referee's reconstruction from the seven committed contrasts
    c = [1.56, 2.22, 0.87, 0.35, 0.93, 0.70, 0.84]
    r = jdm.monthly_contrast_test(c)
    assert r["n"] == 7
    assert r["p_sign_two_sided"] == pytest.approx(2 * 29 / 128, abs=0.01)  # 2-above/5-below
    lo, hi = r["ci95"]
    assert lo < 1.0 < hi  # does not reject unity...
    assert hi > 1.5  # ...and does not reject a 1.5x enhancement either -- the point


def test_box_shift_scan_peaks_at_zero_on_injected_boxes():
    s = jdm.synthetic_orbit(p_in=0.5, p_out=0.03, seed=2)
    scan = jdm.box_shift_scan(s["cml"], s["io_phase"], s["active"])
    best = max((r for r in scan if r["contrast"] is not None), key=lambda r: r["contrast"])
    assert best["shift_deg"] == 0.0  # the frame convention is only right if zero shift wins


def test_episode_stats_run_lengths():
    a = np.array([0, 1, 1, 0, 1, 0, 0, 1, 1, 1], bool)
    e = jdm.episode_stats(a)
    assert e["n_episodes"] == 3 and e["max_episode_bins"] == 3


def test_per_region_contrast_keys():
    s = jdm.synthetic_orbit(seed=3)
    m = jdm.occurrence_map(s["cml"], s["io_phase"], s["active"])
    reg = jdm.per_region_contrast(m)
    assert set(reg) == {"Io-A", "Io-B", "Io-C", "Io-D"}
    assert all(v is None or v > 1.0 for v in reg.values())  # every injected box is enhanced


def test_run_offline_commits_the_recovery_curve(tmp_path):
    m = jdm.run(out=str(tmp_path), offline=True)
    inj = [c["injected"] for c in m["recovery_curve"]]
    assert inj == [1.25, 1.5, 2.0, 8.75]
    # every recovered point is contracted toward unity (boundary-cell dilution), so the
    # measured real 1.12 must be read against this curve, not against the injected values
    for c in m["recovery_curve"]:
        assert c["recovered"] is not None and c["recovered"] <= c["injected"] + 0.4
