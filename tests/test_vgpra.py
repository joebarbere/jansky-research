"""Tests for jansky_research.vgpra -- Voyager 2 PRA Uranus/Neptune rotation periods. Offline."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from jansky_research import vgpra as vg

SAMPLE = Path(__file__).parent / "data" / "vg2_nep_pra_6sec_sample.tab"


def test_read_pra_series_parses_vendored_real_sample():
    ds = vg.read_pra_series(SAMPLE)
    # 8 major frames x 8 sweeps = 64 spectra, 70 channels each
    assert ds["spectra"].shape == (64, 70)
    assert ds["freqs_khz"].shape == (70,)
    # channel grid: 1326.0 kHz down to 1.2 kHz
    assert abs(ds["freqs_khz"][0] - 1326.0) < 1e-6 and abs(ds["freqs_khz"][-1] - 1.2) < 1e-6
    # times monotone non-decreasing, starting at 0, spanning ~8*48 s
    t = ds["times_hr"]
    assert t[0] == 0.0 and np.all(np.diff(t) >= 0)
    assert 0.08 < t[-1] < 0.12  # ~ 7*48s + 7*6s in hours
    # millibell values are physical (the vendored rows sit ~1000-6000 mB)
    assert np.all(ds["spectra"] > 0) and ds["spectra"].max() < 10000


def test_band_flux_is_positive_linear_power():
    ds = vg.read_pra_series(SAMPLE)
    flux = vg.band_flux(ds["spectra"], ds["freqs_khz"])
    assert flux.shape == (64,)
    assert np.all(flux > 0)  # linear power sum
    # out-of-coverage band -> zeros
    assert np.all(vg.band_flux(ds["spectra"], ds["freqs_khz"], band_khz=(5000.0, 6000.0)) == 0)


def test_detect_bursts_finds_injected_episodes():
    s = vg.synthetic_flyby(period_hr=16.0, n_rot=12, seed=3)
    bursts = vg.detect_bursts(s["times_hr"], s["flux"])
    # roughly one episode per rotation (allow slack for merges/misses)
    assert 6 <= bursts.size <= 18
    # bursts cluster at one rotation phase -> low circular spread
    phase = (bursts / 16.0) % 1.0
    r = np.hypot(np.cos(2 * np.pi * phase).mean(), np.sin(2 * np.pi * phase).mean())
    assert r > 0.5  # concentrated, not uniform


def test_period_posterior_recovers_known_period():
    # clean periodic epochs at 17.0 h over 20 cycles
    rng = np.random.default_rng(0)
    epochs = np.sort((np.arange(20) + rng.uniform(0, 0.1, 20)) * 17.0)
    post = vg.period_posterior(epochs, p_lo=16.6, p_hi=17.4, n_boot=100, seed=0)
    assert abs(post["best_period_hr"] - 17.0) < 0.05
    assert post["z2"] > 10 and post["n_bursts"] == 20
    assert post["boot_sigma_hr"] >= 0.0
    # too few epochs -> NaN posterior, no crash
    assert np.isnan(
        vg.period_posterior(np.array([1.0, 2.0]), p_lo=16.0, p_hi=18.0)["best_period_hr"]
    )


def test_bin_series_block_averages():
    t = np.array([0.0, 0.02, 0.05, 0.11, 0.13])
    y = np.array([1.0, 3.0, 2.0, 10.0, 20.0])
    tb, yb = vg.bin_series(t, y, 0.1)
    assert tb.size == 2  # bins [0,0.1) and [0.1,0.2)
    assert abs(yb[0] - 2.0) < 1e-9 and abs(yb[1] - 15.0) < 1e-9


def test_flux_period_posterior_recovers_known_period():
    # PRIMARY method: Lomb-Scargle on the continuous flux series recovers the injected period to
    # within a few minutes (the bootstrap band on a clean synthetic is far tighter than real data)
    s = vg.synthetic_flyby(period_hr=16.5, n_rot=25, seed=2)
    post = vg.flux_period_posterior(
        s["times_hr"], s["flux"], p_lo=15.8, p_hi=17.2, n_boot=80, seed=0
    )
    assert abs(post["best_period_hr"] - 16.5) < 0.05  # ~ 3 min
    assert post["n_binned"] > 50 and post["ls_power"] > 0.01
    assert np.isfinite(post["boot_sigma_hr"]) and post["boot_hi_hr"] > post["boot_lo_hr"]
    # too-short series -> NaN, no crash
    short = vg.flux_period_posterior(np.arange(3.0), np.ones(3), p_lo=16.0, p_hi=18.0)
    assert np.isnan(short["best_period_hr"])


def test_band_stability_flags_achromatic_vs_bandlimited():
    # a rotation signal present across ALL channels -> stable peak across sub-bands (small spread)
    s = vg.synthetic_flyby(period_hr=16.4, n_rot=25, seed=1)
    n = s["times_hr"].size
    spectra = np.tile(s["flux"][:, None], (1, 70)) * 100.0  # same modulation in every channel
    freqs = vg.FREQ_KHZ
    stab = vg.band_stability(s["times_hr"], spectra, freqs, p_lo=15.5, p_hi=17.5)
    assert stab["band_spread_hr"] < 0.2  # achromatic -> consistent across bands
    assert len(stab["band_peaks_hr"]) == len(vg.STABILITY_BANDS_KHZ)
    # a signal in ONE narrow channel range only -> other bands see noise -> large spread
    rng = np.random.default_rng(0)
    spectra2 = rng.normal(2000, 50, (n, 70))
    lowband = vg.FREQ_KHZ < 200  # inject only in the lowest-frequency channels
    spectra2[:, lowband] += (s["flux"][:, None] - s["flux"].mean()) * 30
    stab2 = vg.band_stability(s["times_hr"], spectra2, freqs, p_lo=15.5, p_hi=17.5)
    assert stab2["band_spread_hr"] > stab["band_spread_hr"]


def test_synthetic_flyby_recover_a_known():
    s = vg.synthetic_flyby(period_hr=17.24, n_rot=17, seed=0)
    bursts = vg.detect_bursts(s["times_hr"], s["flux"])
    post = vg.period_posterior(bursts, p_lo=16.6, p_hi=18.0, n_boot=120, seed=1)
    assert vg._consistent(post, 17.24)  # secondary (Rayleigh) also inside its bootstrap band


def test_compare_periods_three_way():
    tw = {"best_period_hr": 17.25, "boot_sigma_hr": 0.1}
    rows = vg.compare_periods("URANUS", tw)
    srcs = {r["source"] for r in rows}
    assert {"radio_1986", "lamy_2025", "this_work"} <= srcs
    assert next(r for r in rows if r["source"] == "this_work")["period_hr"] == 17.25


def test_run_offline_recovers_and_writes(tmp_path):
    m = vg.run(str(tmp_path), offline=True)
    assert m["recovered_injected"]  # synthetic period recovered within bootstrap band
    assert abs(m["best_period_hr"] - 17.24) < 0.2
    saved = json.loads((tmp_path / "results" / "vgpra_metrics.json").read_text())
    assert saved["injected_period_hr"] == m["injected_period_hr"]
    assert (tmp_path / "papers" / "vgpra" / "figures" / "vgpra.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "vgpra" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\vgSynRecoveredOK}{yes}" in macros
    assert r"\newcommand{\vgRealUPeriod}{--}" in macros  # real namespace placeholder offline


def test_write_macros_real_namespace(tmp_path):
    real = {
        "source": "PDS-PPI VG2-PRA",
        "is_real": True,
        "planets": {
            "URANUS": {  # band-wanders + offset -> NOT recovered; consistent only within wide 2 sigma
                "best_period_hr": 18.44,
                "boot_sigma_hr": 1.94,
                "band_spread_hr": 1.76,
                "total_unc_hr": 1.94,
                "rayleigh_period_hr": 14.72,
                "span_hr": 312.0,
                "recovers_hist": False,
                "consistent_hist": True,
                "consistent_lamy": True,
            },
            "NEPTUNE": {  # railed to window edge -> NOT recovered, NOT consistent
                "best_period_hr": 20.0,
                "boot_sigma_hr": 1.94,
                "band_spread_hr": 0.70,
                "total_unc_hr": 1.94,
                "rayleigh_period_hr": 15.95,
                "span_hr": 504.0,
                "recovers_hist": False,
                "consistent_hist": False,
            },
        },
    }
    p = tmp_path / "m.tex"
    vg._write_macros(real, p)
    txt = p.read_text()
    assert r"\newcommand{\vgRealUPeriod}{18.44}" in txt
    assert (
        r"\newcommand{\vgRealURecovers}{no}" in txt
    )  # neither planet recovered by the blind method
    assert r"\newcommand{\vgRealNRecovers}{no}" in txt
    assert r"\newcommand{\vgRealNHistOK}{no}" in txt  # Neptune not even consistent (railed edge)
    # synthetic namespace still live in a real build (the method works on a clean injected signal)
    assert r"\newcommand{\vgSynRecoveredOK}{yes}" in txt
