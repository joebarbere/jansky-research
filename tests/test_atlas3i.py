"""Tests for jansky_research.atlas3i — BL 3I/ATLAS reproduction. No network."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from jansky_research import atlas3i

FIXTURE = Path(__file__).parent / "data" / "gb_atlas_index_sample.html"


def test_parse_index_real_listing():
    scans = atlas3i.parse_index(FIXTURE.read_text())
    # 6 L-band scans + 1 S-band scan for blc25, x3 products each
    assert len(scans) == 21
    fine = [s for s in scans if s.product == "0000"]
    assert len(fine) == 7
    assert all(s.node == "blc25" for s in fine)
    first = fine[0]
    assert (first.sec, first.on, first.scan_id) == (17226, True, "0013")
    assert first.filename == "blc25_guppi_61027_17226_DIAG_3I_ATLAS_0013.rawspec.0000.h5"
    # each filename appears in 2 places per index line (href + label prefix) -> deduplicated
    assert len({s.filename for s in scans}) == len(scans)


def test_cadence_selects_abacad_block():
    scans = atlas3i.parse_index(FIXTURE.read_text())
    blk = atlas3i.cadence(scans, "blc25", band="L")
    assert [s.sec for s in blk] == [17226, 17543, 17861, 18179, 18497, 18814]
    assert [s.on for s in blk] == [True, False, True, False, True, False]
    # the S-band scan at 21817 is not part of the L cadence, so blc25 S is incomplete here
    with pytest.raises(ValueError, match="expected 6"):
        atlas3i.cadence(scans, "blc25", band="S")
    with pytest.raises(ValueError, match="expected 6"):
        atlas3i.cadence(scans, "blc99", band="L")


def test_cadence_rejects_broken_onoff_pattern():
    scans = atlas3i.parse_index(FIXTURE.read_text())
    broken = [
        atlas3i.Scan(s.node, s.sec, not s.on, s.scan_id, s.product, s.filename)
        for s in scans
        if s.product == "0000" and s.sec < 20000
    ]
    with pytest.raises(ValueError, match="not an ABACAD"):
        atlas3i.cadence(broken, "blc25", band="L")


def test_normalize_flattens_bandpass_and_excises_dc_spike():
    rng = np.random.default_rng(0)
    n_time, n_freq = 16, 1024
    bandpass = 1.0 + 0.5 * np.sin(np.linspace(0, 3, n_freq))
    wf = rng.normal(10.0, 1.0, (n_time, n_freq)) * bandpass
    wf[:, n_freq // 2] += 1e3
    out, bad = atlas3i.normalize(wf)
    assert bad[n_freq // 2]
    assert bad.sum() < 5  # only the spike (and at most stray outliers) excised
    assert abs(np.median(out) - 1.0) < 0.05  # bandpass flattened to ~unity
    assert out[:, n_freq // 2].max() <= 1.0  # spike gone


def test_dedoppler_concentrates_a_drifting_tone():
    tsamp, foff = 18.6, -2.79
    n_time, n_freq, drift = 16, 512, 0.4
    wf = np.zeros((n_time, n_freq))
    t = np.arange(n_time) * tsamp
    chans = np.rint(200 + drift * t / foff).astype(int)
    wf[np.arange(n_time), chans] = 1.0
    dd = atlas3i.dedoppler(wf, tsamp_s=tsamp, foff_hz=foff, drift_rates=np.array([0.0, drift]))
    # at the true drift all power stacks into the start channel; at zero drift it smears
    assert dd[1].max() == pytest.approx(1.0)
    assert int(dd[1].argmax()) == 200
    assert dd[0].max() < 0.5


def test_find_hits_recovers_tone_and_suppresses_neighbours():
    wfs, _, _ = atlas3i.synthetic_cadence(rfi_chan=None, dc_spike=False, seed=1)
    wfn, _ = atlas3i.normalize(wfs[0])
    hits = atlas3i.find_hits(
        wfn, tsamp_s=18.6, foff_hz=-2.79, drift_rates=np.linspace(-1, 1, 41), threshold=10.0
    )
    assert len(hits) == 1  # neighbour suppression: one tone -> one hit
    assert abs(hits[0].chan - 1024) < 8
    assert hits[0].drift_hz_s == pytest.approx(0.4, abs=0.11)


def test_onoff_filter_keeps_sky_signal_drops_rfi():
    tsamp, foff = 18.6, -2.79
    drifts = np.linspace(-1, 1, 41)
    wfs, secs, on = atlas3i.synthetic_cadence(tsamp_s=tsamp, foff_hz=foff, seed=0)
    hits = []
    for wf in wfs:
        wfn, _ = atlas3i.normalize(wf)
        hits.append(
            atlas3i.find_hits(wfn, tsamp_s=tsamp, foff_hz=foff, drift_rates=drifts, threshold=10.0)
        )
    survivors = atlas3i.onoff_filter(hits, secs, on, foff_hz=foff)
    assert any(abs(h.chan - 1024) < 32 for h in survivors)  # injected tone survives
    assert not any(abs(h.chan - 3000) < 32 for h in survivors)  # always-on RFI rejected


def test_onoff_filter_validates_inputs():
    with pytest.raises(ValueError, match="equal length"):
        atlas3i.onoff_filter([[]], [1, 2], [True], foff_hz=-2.79)
    with pytest.raises(ValueError, match="start with an ON"):
        atlas3i.onoff_filter([[], []], [1, 2], [False, True], foff_hz=-2.79)


def test_eirp_limit_reproduces_paper_headline():
    # Paper params: 16 sigma, GBT L-band SEFD ~10 Jy, 5-min scans, ~2.79 Hz channels, 1.8 au.
    eirp = atlas3i.eirp_limit_w(snr=16.0)
    assert 0.03 < eirp < 0.3  # ~100 mW headline figure
    # scaling sanity: EIRP ~ d^2 and ~ snr
    assert atlas3i.eirp_limit_w(snr=16.0, distance_au=3.6) == pytest.approx(4 * eirp)
    assert atlas3i.eirp_limit_w(snr=32.0) == pytest.approx(2 * eirp)


def test_run_offline_round_trip(tmp_path):
    m = atlas3i.run(out=str(tmp_path))
    assert m["recovered_injected"] is True
    assert m["rfi_rejected"] is True
    assert m["dc_spike_excised"] is True
    assert m["survivor_drift_hz_s"] == pytest.approx(0.4, abs=0.11)
    assert (tmp_path / "results" / "atlas3i_metrics.json").exists()
    assert (tmp_path / "papers" / "atlas3i" / "figures" / "atlas3i_dedoppler.pdf").exists()
    macros = (tmp_path / "papers" / "atlas3i" / "generated" / "macros.tex").read_text()
    assert r"\aiEirpMw" in macros
