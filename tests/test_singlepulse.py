"""Tests for jansky_research.singlepulse -- FDMT single-pulse recover-a-known. Offline."""

from __future__ import annotations

import json
import struct

import numpy as np
import pytest

pytest.importorskip("torch")

from jansky_research import singlepulse as sp  # noqa: E402


def _write_fil(path, dyn8, fch1=1600.0, foff=-1.5625, tsamp=5.12e-4):
    """Write a minimal 8-bit SIGPROC file the reader must round-trip."""

    def s(key):
        return struct.pack("<i", len(key)) + key.encode()

    hdr = s("HEADER_START")
    hdr += s("source_name") + s("J0534+2200")
    for k, v in [("nchans", dyn8.shape[1]), ("nbits", 8), ("nifs", 1)]:
        hdr += s(k) + struct.pack("<i", v)
    for k, v in [("fch1", fch1), ("foff", foff), ("tsamp", tsamp), ("tstart", 58543.0)]:
        hdr += s(k) + struct.pack("<d", v)
    hdr += s("HEADER_END")
    path.write_bytes(hdr + dyn8.astype(np.uint8).tobytes())


def test_read_sigproc_round_trips(tmp_path):
    rng = np.random.default_rng(0)
    dyn = rng.integers(0, 255, (128, 64)).astype(np.uint8)
    f = tmp_path / "t.fil"
    _write_fil(f, dyn)
    out, freqs, hdr = sp.read_sigproc(f)
    assert out.shape == (128, 64)
    assert np.array_equal(out, dyn.astype(np.float32))
    assert hdr["source_name"] == "J0534+2200"
    assert freqs[0] == 1600.0 and freqs[1] < freqs[0]  # foff negative


def test_read_sigproc_rejects_non_filterbank(tmp_path):
    f = tmp_path / "bad.fil"
    f.write_bytes(struct.pack("<i", 3) + b"NOP" + b"x" * 50)
    with pytest.raises(ValueError, match="not a SIGPROC"):
        sp.read_sigproc(f)


def test_search_recovers_crab_like_dm():
    dyn, freqs, dt = sp.synthetic_observation(seed=2)
    r = sp.search(dyn, freqs, dt, max_dm=120.0)
    assert abs(r["best_dm"] - sp.CRAB_DM) < 2.0
    assert r["best_snr"] > 20.0
    assert r["sp_snr"] > 5.0  # boxcar finds individual pulses


def test_run_offline_recovers_dm_and_period(tmp_path):
    m = sp.run(str(tmp_path), offline=True)
    assert abs(m["recovered_dm"] - m["true_dm"]) < 2.0
    assert abs(m["recovered_period_ms"] - m["true_period_ms"]) < 0.5
    saved = json.loads((tmp_path / "results" / "singlepulse_metrics.json").read_text())
    assert saved == m
    fig = tmp_path / "papers" / "torchfdmt" / "figures" / "singlepulse.pdf"
    assert fig.exists() and fig.stat().st_size > 0
    macros = (tmp_path / "papers" / "torchfdmt" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\spRecDm}" in macros
    assert r"\newcommand{\spRealDm}{--}" in macros  # offline: real macros are placeholders


def test_write_macros_placeholders(tmp_path):
    p = tmp_path / "m.tex"
    sp._write_macros({"source": "x"}, p)
    t = p.read_text()
    assert r"\newcommand{\spTrueDm}{--}" in t and r"\newcommand{\spRealDmErr}{--}" in t


def test_speedup_macros_are_derived_not_typed(tmp_path):
    """The paper once carried a hand-typed "~24x" beside the two macros making it 29x."""
    mac = tmp_path / "m.tex"
    sp._write_macros(
        {
            "source": "x",
            "bench_brute_cpu_s": 44.12,
            "bench_brute_gpu_s": 1.5,
            "bench_fdmt_cpu_s": 0.45,
            "bench_fdmt_gpu_s": 0.44,
        },
        mac,
    )
    text = mac.read_text()
    assert r"\newcommand{\spBruteSpeedup}{29}" in text
    assert r"\newcommand{\spFdmtSpeedup}{1}" in text


def test_speedup_macro_is_a_placeholder_without_a_gpu_leg(tmp_path):
    """A CPU-only run must not emit a ratio; preserve_live_macros then keeps any real one."""
    mac = tmp_path / "m.tex"
    sp._write_macros({"source": "x", "bench_brute_cpu_s": 44.12, "bench_brute_gpu_s": None}, mac)
    assert r"\newcommand{\spBruteSpeedup}{--}" in mac.read_text()


def test_committed_speedup_matches_the_committed_timings():
    """The macro on disk must equal the ratio of the two numbers in the evidence file."""
    import json
    import re
    from pathlib import Path

    metrics = Path("results/singlepulse_metrics.json")
    macros = Path("papers/torchfdmt/generated/macros.tex")
    if not (metrics.exists() and macros.exists()):  # pragma: no cover
        pytest.skip("real results not present")
    m = json.loads(metrics.read_text())
    got = re.search(r"\\newcommand\{\\spBruteSpeedup\}\{([^}]*)\}", macros.read_text())
    assert got is not None
    assert got.group(1) == f"{m['bench_brute_cpu_s'] / m['bench_brute_gpu_s']:.0f}"


def test_dm_step_is_one_delay_sample():
    """The FDMT quantises DM at one sample of band-crossing delay; 0.3% must be readable
    against that grid rather than as a bare percentage."""
    freqs = np.linspace(702.0, 4030.0, 1024)
    step = sp._dm_step(freqs, 5.12e-4)
    # 4.148808e3 * (702^-2 - 4030^-2) s per unit DM = 8.16 ms; one 0.512 ms sample is ~0.063.
    assert step == pytest.approx(0.0627, abs=5e-4)
    # and the committed real offset is a few trials, not a fitted precision
    assert abs(56.59 - 56.77) / step == pytest.approx(2.9, abs=0.3)


def test_dm_curve_is_the_normalised_statistic():
    # The figure must plot the statistic best() maximises: the curve's maximum sits at the
    # recovered DM with the reported height, not at the highest delay row (the raw-sum trend).
    dyn, freqs, dt = sp.synthetic_observation(n_time=2048)
    s = sp.search(dyn, freqs, dt, max_dm=120.0)
    i = int(np.argmax(s["dm_curve"]))
    assert abs(float(s["dms"][i]) - s["best_dm"]) < 1e-9
    assert abs(float(s["dm_curve"][i]) - s["best_snr"]) < 1e-6
    # the raw curve is monotone-trending upward with row count and peaks elsewhere; it is
    # committed only as the diagnostic
    assert "dm_curve_raw" in s and s["dm_curve_raw"].shape == s["dm_curve"].shape
    assert s["n_time"] == 2048


def test_shift_null_bounds_the_noise_maximum():
    rng = np.random.default_rng(1)
    dyn = rng.normal(0.0, 1.0, (1024, 64)).astype(np.float32)
    freqs = np.linspace(1200.0, 1600.0, 64)
    null = sp.shift_null(dyn, freqs, 1e-3, max_dm=60.0, n_reps=15, seed=0)
    assert null["n_reps"] == 15 and null["best_snrs"].shape == (15,)
    assert 3.0 < null["p50"] < 9.0  # a noise plane's maximum is itself several sigma
    assert null["max"] >= null["p99"] >= null["p50"] >= null["mean"] - 2.0
    # an injected strong pulse beats the null decisively
    s = sp.search(*sp.synthetic_observation(n_time=1024), max_dm=120.0)
    assert s["best_snr"] > null["max"]
