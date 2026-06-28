"""Tests for jansky_research.driftsearch — injection-recovery benchmark. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import driftsearch


def test_completeness_snr_interpolates():
    snrs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    p = np.array([0.0, 0.2, 0.5, 0.8, 1.0])
    assert abs(driftsearch.completeness_snr(snrs, p, 0.5) - 2.0) < 1e-9
    # never reaches the level -> nan
    assert np.isnan(driftsearch.completeness_snr(snrs, np.full(5, 0.3), 0.5))


def test_false_positive_rate_low_at_threshold():
    # Noise-only best S/N peaks near ~5; threshold 10 should almost never fire.
    fpr = driftsearch.false_positive_rate(n_trials=150, threshold=10.0, seed=0)
    assert fpr < 0.05


def test_injection_recovery_separates_weak_from_strong():
    res = driftsearch.injection_recovery(
        np.array([0.25, 3.0]),
        drift_rates=np.array([0.3]),
        n_trials=15,
        fpr_trials=100,
        seed=0,
    )
    pm = res.p_detect.mean(axis=1)
    assert pm[0] < 0.5  # a very weak tone is rarely recovered
    assert pm[1] > 0.9  # a strong tone is almost always recovered
    assert res.false_positive_rate < 0.05
    assert res.inj_snrs[0] <= res.completeness_snr_50 <= res.inj_snrs[-1]


def test_run_offline_writes_artifacts(tmp_path):
    m = driftsearch.run(out=str(tmp_path), n_trials=10)
    assert 0.0 <= m["false_positive_rate"] < 0.1
    assert np.isfinite(m["completeness_snr_50"])  # the grid spans the transition
    assert (tmp_path / "results" / "drift_metrics.json").exists()
    assert (tmp_path / "papers" / "driftsearch" / "figures" / "drift_recovery.pdf").exists()
    assert (tmp_path / "papers" / "driftsearch" / "generated" / "macros.tex").exists()
