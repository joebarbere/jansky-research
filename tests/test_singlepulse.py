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
