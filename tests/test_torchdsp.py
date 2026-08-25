"""Tests for jansky_research.torchdsp -- the pure-PyTorch coherent-DSP suite. Offline.

Skipped without the ``fdmt`` extra (torch), like test_fdmt.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

pytest.importorskip("torch")

from jansky import rfi as rfi_cpu  # noqa: E402

from jansky_research import torchdsp as td  # noqa: E402


def test_chirp_round_trip_recollapses_impulse():
    syn = td.synthetic_dispersed_voltage(dm=150.0, seed=1)
    # dispersed: energy spread far beyond any 9 samples
    p_before = np.abs(syn["voltage"]) ** 2
    assert np.sort(p_before)[-9:].sum() / p_before.sum() < 0.5
    v = td.coherent_dedisperse(
        syn["voltage"], syn["dm"], syn["f0_mhz"], chan_bw_mhz=syn["chan_bw_mhz"]
    )
    p_after = np.asarray((v.abs() ** 2).cpu())
    peak = int(np.argmax(p_after))
    assert abs(peak - syn["impulse_index"]) <= 1
    w = 4
    assert p_after[peak - w : peak + w + 1].sum() / p_after.sum() > 0.85


def test_chirp_wrong_dm_does_not_recollapse():
    syn = td.synthetic_dispersed_voltage(dm=150.0, seed=2)
    v = td.coherent_dedisperse(syn["voltage"], 75.0, syn["f0_mhz"], chan_bw_mhz=syn["chan_bw_mhz"])
    p = np.asarray((v.abs() ** 2).cpu())
    peak = int(np.argmax(p))
    w = 4
    assert p[max(0, peak - w) : peak + w + 1].sum() / p.sum() < 0.5


def test_dedisperse_channelized_aligns_channels():
    # two channels, same impulse: after inter-channel correction the peaks coincide
    n = 2048
    freqs = np.array([500.0, 700.0])
    chans = []
    for f0 in freqs:
        syn = td.synthetic_dispersed_voltage(
            n_time=n, dm=30.0, f0_mhz=float(f0), noise=0.01, seed=3
        )
        chans.append(syn["voltage"])
    v = np.stack(chans)
    # emulate the inter-channel delay the ISM adds (relative to the top channel)
    from jansky.constants import DM_CONST

    dt_s = 1.0e-6 / td.CHIME_CHAN_BW_MHZ
    lag = int(round(DM_CONST * 30.0 * (freqs[0] ** -2 - freqs[1] ** -2) / dt_s))
    v[0] = np.roll(v[0], lag)
    ded = td.dedisperse_channelized(v, 30.0, freqs, chan_bw_mhz=td.CHIME_CHAN_BW_MHZ)
    p = np.asarray((ded.abs() ** 2).cpu())
    assert abs(int(np.argmax(p[0])) - int(np.argmax(p[1]))) <= 1


def test_spectral_kurtosis_matches_cpu_oracle():
    rng = np.random.default_rng(4)
    power = rng.chisquare(2, size=(256, 32))
    sk_t = np.asarray(td.spectral_kurtosis(power, axis=0).cpu())
    sk_cpu = rfi_cpu.spectral_kurtosis(power, axis=0)
    assert np.allclose(sk_t, sk_cpu, atol=1e-10)


def test_sumthreshold_sequential_equals_cpu_oracle():
    rng = np.random.default_rng(5)
    x = rng.standard_normal(512)
    x[100:130] += 2.5  # faint extended RFI: the SumThreshold specialty
    x[300] += 12.0
    m_t = td.sumthreshold(x, sequential=True)
    m_cpu = rfi_cpu.sumthreshold(x)
    assert np.array_equal(m_t, m_cpu)
    # with an initial mask threaded through
    seed = np.zeros(512, bool)
    seed[400:410] = True
    assert np.array_equal(
        td.sumthreshold(x, mask=seed, sequential=True), rfi_cpu.sumthreshold(x, mask=seed)
    )


def test_sumthreshold_parallel_agrees_and_catches_rfi():
    rng = np.random.default_rng(6)
    x = rng.standard_normal(1024)
    x[500:540] += 2.5
    m_seq = td.sumthreshold(x, sequential=True)
    m_par = td.sumthreshold(x)
    # the parallel variant must catch the injected block and agree closely overall
    assert m_par[500:540].mean() > 0.9
    inter = np.logical_and(m_seq, m_par).sum()
    union = np.logical_or(m_seq, m_par).sum()
    assert inter / union > 0.8


def test_sumthreshold2d_sequential_equals_cpu_oracle():
    rng = np.random.default_rng(7)
    dyn = rng.standard_normal((96, 24))
    dyn[:, 5] += 3.0
    dyn[40:44, :] += 2.5
    m_t = td.sumthreshold2d(dyn, sequential=True)
    m_cpu = rfi_cpu.sumthreshold2d(dyn)
    assert np.array_equal(m_t, m_cpu)


def test_ffa_recovers_injected_period_and_matches_fold_oracle():
    rng = np.random.default_rng(8)
    n = 1 << 15
    p_true = 217.4
    ts = rng.standard_normal(n).astype(np.float32)
    idx = np.arange(0, n, p_true).astype(int)
    ts[idx[idx < n]] += 5.0
    out = td.ffa_search(ts, pmin_samples=190, pmax_samples=250)
    assert abs(out["period_samples"] - p_true) < 1.0
    oracle = td.fold_snr(ts, p_true)
    assert out["snr"] > 0.5 * oracle  # FFA S/N within reach of the exact-period fold
    assert oracle > 10.0


def test_ffa_flat_noise_finds_nothing_loud():
    rng = np.random.default_rng(9)
    ts = rng.standard_normal(1 << 14).astype(np.float32)
    out = td.ffa_search(ts, pmin_samples=100, pmax_samples=140)
    assert out["snr"] < 8.0


def test_chirp_group_delay_matches_cold_plasma_law():
    """The check the self-composing round trip cannot perform: the chirp's numerical group
    delay must equal the cold-plasma law. A wrong constant or exponent fails loudly."""
    gd = td.chirp_group_delay_residual(dm=556.0, f0_mhz=400.0)
    assert gd["max_resid_samples"] < 1e-6
    # sanity that the check CAN fail: a chirp for a different DM is not the law for this one
    f = np.linspace(-td.CHIME_CHAN_BW_MHZ / 2, td.CHIME_CHAN_BW_MHZ / 2, 101)
    from jansky.constants import DM_CONST

    tau_law = -1.0e6 * DM_CONST * 556.0 * ((400.0 + f) ** -2 - 400.0**-2)
    tau_wrong = -1.0e6 * DM_CONST * 500.0 * ((400.0 + f) ** -2 - 400.0**-2)
    assert np.max(np.abs(tau_wrong - tau_law)) * td.CHIME_CHAN_BW_MHZ > 1e-3


def test_dedisperse_channelized_conjugate_flag_changes_the_kernel():
    syn = td.synthetic_dispersed_voltage(dm=150.0, seed=11)
    v = syn["voltage"][None, :]
    freqs = np.array([syn["f0_mhz"]])
    right = td.dedisperse_channelized(v, 150.0, freqs, chan_bw_mhz=syn["chan_bw_mhz"])
    wrong = td.dedisperse_channelized(
        v, 150.0, freqs, chan_bw_mhz=syn["chan_bw_mhz"], conjugate=True
    )
    p_right = np.asarray((right[0].abs() ** 2).cpu())
    p_wrong = np.asarray((wrong[0].abs() ** 2).cpu())
    w = 4
    pk = int(np.argmax(p_right))
    frac_right = p_right[pk - w : pk + w + 1].sum() / p_right.sum()
    pk_w = int(np.argmax(p_wrong))
    frac_wrong = p_wrong[max(0, pk_w - w) : pk_w + w + 1].sum() / p_wrong.sum()
    assert frac_right > 0.85  # the correct sign re-collapses the impulse
    assert frac_wrong < 0.5  # the conjugate doubles the smear instead


def test_mask_agreement_reports_direction():
    a = np.zeros((4, 4), bool)
    b = np.zeros((4, 4), bool)
    a[0, :2] = True  # parallel over-flags two
    b[1, :3] = True  # oracle-only three
    out = td._mask_agreement(a, b)
    assert out["n_par_only"] == 2 and out["n_cpu_only"] == 3
    assert out["jaccard"] == 0.0


def test_run_offline_writes_artifacts(tmp_path):
    m = td.run(str(tmp_path), offline=True)
    assert m["dedisp"]["peak_offset_samples"] <= 1
    assert m["dedisp"]["reconcentrated_energy_frac"] > 0.85
    # the validated path must be the shipped one: dedisperse_channelized round-trips too
    assert m["dedisp"]["peak_offset_samples_channelized"] <= 1
    assert m["dedisp"]["reconcentrated_energy_frac_channelized"] > 0.85
    assert m["dedisp"]["group_delay_max_resid_samples"] < 1e-6
    assert m["sk_max_diff"] < 1e-9
    assert m["sumthreshold_sequential_equals_oracle"]
    assert m["sumthreshold_parallel_jaccard"] > 0.8
    # the divergence direction is committed, and the sweep bounds it off the default config
    assert m["sumthreshold_par_only_flags"] >= 0 and m["sumthreshold_cpu_only_flags"] >= 0
    sweep = m["sumthreshold_agreement_sweep"]
    assert "thr3.0_iter2" in sweep and all(0.0 <= v["jaccard"] <= 1.0 for v in sweep.values())
    assert m["rfi_line_caught"] and m["rfi_burst_caught"]
    assert m["ffa_period_err_samples"] < 1.0
    # the injected period is not on the drift grid; the honest resolution is committed
    assert 0 < m["ffa_drift_resolution_samples"] < 0.1
    assert m["ffa_snr_ratio_predicted"] == pytest.approx(m["ffa_snr_ratio_measured"], abs=0.15)
    saved = json.loads((tmp_path / "results" / "torchdsp_metrics.json").read_text())
    assert saved["ffa_period_found"] == m["ffa_period_found"]
    assert (tmp_path / "papers" / "torchdsp" / "figures" / "torchdsp.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "torchdsp" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\tdSynFfaErr}" in macros and r"\newcommand{\tdRealFfaErr}{--}" in macros
    # derived speedups exist as macros (placeholders here: no benchmark in an offline run)
    assert r"\newcommand{\tdRealSpeedupFfa}{--}" in macros


def test_write_macros_derives_speedups_from_the_two_columns(tmp_path):
    p = tmp_path / "m.tex"
    td._write_macros(
        {
            "source": "x",
            "device": "cpu",
            "is_real": True,
            "benchmark": {"chirp_s": 2.0, "sumthreshold_s": 8.0, "ffa_s": 0.65},
            "benchmark_cpu": {"chirp_s": 2.0, "sumthreshold_s": 3.0, "ffa_s": 8.42},
            "benchmark_device": "cuda",
        },
        p,
    )
    txt = p.read_text()
    # 8.42 / 0.65 = 12.95... -> 13.0, derived, never hand-typed
    assert r"\newcommand{\tdRealSpeedupFfa}{13.0}" in txt
    assert r"\newcommand{\tdRealBenchFfaGpu}{0.65}" in txt
    assert r"\newcommand{\tdRealBenchFfaCpu}{8.42}" in txt


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    td._write_macros({"source": "x", "device": "cpu", "is_real": True, "sk_max_diff": None}, p)
    txt = p.read_text()
    assert r"\newcommand{\tdRealSkMaxDiff}{--}" in txt
    assert r"\newcommand{\tdSynSkMaxDiff}{--}" in txt
