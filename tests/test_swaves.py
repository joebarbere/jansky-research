"""Tests for jansky_research.swaves -- STEREO/WAVES interplanetary type III. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import swaves


def test_parse_swaves_ascii():
    text = "\n".join(
        [
            "  125.0   175.0   225.0",  # frequency axis (kHz)
            "  1.0   1.0   1.0",  # background row (skipped)
            "0   5.0   1.0   1.0",  # minute 0: index + 3 intensities
            "1   1.0   5.0   1.0",  # minute 1
            "2   1.0   1.0   5.0",  # minute 2
        ]
    )
    d = swaves.parse_swaves_ascii(text)
    assert d["data"].shape == (3, 3)  # (n_freq, n_time)
    assert np.allclose(d["freqs"], [0.225, 0.175, 0.125])  # MHz, descending
    assert np.allclose(d["times"], [0.0, 60.0, 120.0])  # minutes -> seconds


def test_synthetic_burst_over_hfr_band():
    b = swaves.synthetic_ip_burst(seed=0)
    assert b["data"].shape[0] == 319  # HFR channels
    assert b["freqs"][0] > b["freqs"][-1]  # descending
    assert b["freqs"].min() < 0.2 and b["freqs"].max() > 15.0  # spans ~0.125-16 MHz


def test_run_offline_recovers_speed_and_reaches_interplanetary(tmp_path):
    m = swaves.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_ridge"] > 30
    assert 0.05 < m["speed_c"] < 0.4
    assert 0.85 < m["recovery_ratio"] < 1.15
    # the HFR band reaches genuinely interplanetary distances (tens of R_sun)
    assert m["r_hi_rsun"] > 20.0 and m["r_hi_au"] > 0.1
    assert (tmp_path / "results" / "swaves_metrics.json").exists()
    assert (tmp_path / "papers" / "swaves" / "figures" / "swaves.pdf").exists()
    macros = (tmp_path / "papers" / "swaves" / "generated" / "macros.tex").read_text()
    assert r"\swSpeedC" in macros and r"\swRhiAU" in macros
