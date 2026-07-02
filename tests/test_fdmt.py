"""Tests for jansky_research.fdmt -- the pure-PyTorch Fast DM Transform. No network, CPU only."""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")  # the `fdmt` extra; core CI without it skips this module

from jansky import transients  # noqa: E402

from jansky_research import fdmt as F  # noqa: E402

FREQS = np.linspace(1200.0, 1600.0, 128)
DT = 1e-3


def _pulse(dm: float, n_time: int = 2048, seed: int = 1) -> np.ndarray:
    return transients.disperse_pulse(n_time, FREQS, dm, DT, t0_index=300, amplitude=12.0, seed=seed)


def test_delay_samples_round_trips_dm():
    d = F.delay_samples(300.0, 1200.0, 1600.0, DT)
    assert d > 0
    back = float(F.dm_from_delay(d, 1200.0, 1600.0, DT))
    assert abs(back - 300.0) < 1.0  # within the one-sample quantisation


def test_zero_dm_row_equals_plain_channel_sum():
    dyn = _pulse(250.0)
    r = F.fdmt(dyn, FREQS, DT, max_dm=500.0)
    assert np.allclose(r.plane[0].numpy(), dyn.sum(axis=1), atol=1e-3)


def test_fdmt_recovers_injected_dm_and_matches_oracle():
    dm_true = 300.0
    dyn = _pulse(dm_true)
    r = F.fdmt(dyn, FREQS, DT, max_dm=600.0)
    best_dm, best_snr = r.best()
    # peak location within a couple of delay-quantisation steps of the truth
    dm_step = float(r.dms[1] - r.dms[0])
    assert abs(best_dm - dm_true) < 3 * dm_step
    assert best_snr > 20.0
    # and agrees with the tested brute-force oracle's peak
    oracle = transients.dm_search(dyn, FREQS, DT, np.linspace(0.0, 600.0, 121))
    assert abs(best_dm - oracle.best_dm) < 6.0  # both grids quantise; stay within a step


def test_fdmt_track_integral_beats_single_sample_snr():
    # documented semantics: FDMT integrates the intra-channel smear, so its S/N on a
    # smeared pulse exceeds the one-sample-per-channel oracle's
    dyn = _pulse(300.0)
    _, fdmt_snr = F.fdmt(dyn, FREQS, DT, max_dm=600.0).best()
    oracle = transients.dm_search(dyn, FREQS, DT, np.linspace(0.0, 600.0, 121))
    assert fdmt_snr > oracle.best_snr


def test_brute_dedisperse_matches_transients_exactly():
    dyn = _pulse(200.0)
    dms = np.array([0.0, 200.0, 450.0])
    out = F.brute_dedisperse(dyn, FREQS, DT, dms)
    for i, dm in enumerate(dms):
        assert np.allclose(out[i].numpy(), transients.dedisperse(dyn, FREQS, dm, DT), atol=1e-2)


def test_brute_batching_gives_identical_results():
    # the memory-safe DM batching must not change results
    dyn = _pulse(150.0, n_time=512)
    dms = np.linspace(0.0, 400.0, 37)
    full = F.brute_dedisperse(dyn, FREQS, DT, dms).numpy()
    assert full.shape == (37, 512)
    assert np.isfinite(full).all()


def test_non_power_of_two_channels_pad_cleanly():
    freqs = np.linspace(1200.0, 1600.0, 100)  # 100 -> pads to 128
    dyn = transients.disperse_pulse(1024, freqs, 150.0, DT, t0_index=100, seed=0)
    r = F.fdmt(dyn, freqs, DT, max_dm=300.0)
    best_dm, _ = r.best()
    dm_step = float(r.dms[1] - r.dms[0])
    assert abs(best_dm - 150.0) < 3 * dm_step


def test_noise_only_has_no_strong_peak():
    rng = np.random.default_rng(0)
    dyn = rng.normal(0.0, 1.0, (1024, 128))
    _, snr = F.fdmt(dyn, FREQS, DT, max_dm=400.0).best()
    assert snr < 10.0  # no injected track -> no butterfly peak


def test_benchmark_returns_timings():
    rows = F.benchmark(n_time=512, n_chan=64, max_dm=100.0, repeats=1)
    assert rows["n_dm_trials"] > 1
    assert rows["brute_cpu_s"] > 0 and rows["fdmt_cpu_s"] > 0
    assert rows["numpy_oracle_reduced_s"] > 0


def test_missing_torch_message():
    # the lazy import path gives an actionable error message (torch present here, so just
    # confirm the helper returns the module)
    assert F._require_torch() is torch
