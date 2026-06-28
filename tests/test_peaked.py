"""Tests for jansky_research.peaked — peaked-spectrum (GPS/CSS) selection. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import peaked


def test_classify_sed():
    assert peaked.classify_sed(0.5, -0.6) == "peaked"  # rising then falling -> turnover in band
    assert peaked.classify_sed(-0.8, -0.8) == "steep"  # falling throughout
    assert peaked.classify_sed(0.4, 0.3) == "inverted"  # still rising at 3 GHz
    assert peaked.classify_sed(-0.2, -0.2) == "flat"
    assert peaked.classify_sed(np.nan, -0.6) == "nan"


def test_two_point_indices_textbook():
    # a flat power law S ~ nu^-0.7 gives alpha_low = alpha_high = -0.7
    nu = peaked.NU_GHZ
    s_t = (nu["tgss"] / nu["nvss"]) ** -0.7
    s_n = 1.0
    s_v = (nu["vlass"] / nu["nvss"]) ** -0.7
    al, ah = peaked.two_point_indices(np.array([s_t]), np.array([s_n]), np.array([s_v]))
    assert np.isclose(al[0], -0.7, atol=1e-6) and np.isclose(ah[0], -0.7, atol=1e-6)


def test_peak_frequency_recovers_turnover():
    # a parabola in log-log peaking at 1 GHz
    nu = np.array([0.15, 1.0, 5.0])
    flux = 10.0 ** (1.0 - 0.8 * (np.log10(nu / 1.0)) ** 2)
    nu_peak, is_peak = peaked.peak_frequency(flux, nu)
    assert is_peak and abs(nu_peak - 1.0) < 0.05
    # a monotonic power law is not a peak
    _, is_peak2 = peaked.peak_frequency(nu**-0.7, nu)
    assert not is_peak2


def test_find_peaked_upper_limit_recovers_injected():
    from jansky_research.spectra import crossmatch

    tgss, nvss, vlass, truth = peaked.synthetic_field(n_sources=2000, peaked_fraction=0.05, seed=1)
    # the injected peaked sources are faint at 150 MHz -> below the TGSS limit
    assert tgss["ra"].size < nvss["ra"].size  # TGSS is shallower than NVSS
    res = peaked.find_peaked(tgss, nvss, vlass)
    pk = res["is_peaked"]
    assert pk.sum() >= 1
    assert res["alpha_low_is_limit"][pk].mean() > 0.9  # recovered chiefly via the TGSS upper limit
    # recovery: injected-peaked positions that land on a flagged candidate
    rp = np.flatnonzero(truth)
    i, _, _ = crossmatch(nvss["ra"][rp], nvss["dec"][rp], res["ra"][pk], res["dec"][pk], 5.0)
    assert i.size / truth.sum() > 0.6  # most bright injected peaked sources recovered
    # purity: flagged candidates are mostly truly peaked
    j, _, _ = crossmatch(res["ra"][pk], res["dec"][pk], nvss["ra"][rp], nvss["dec"][rp], 5.0)
    assert j.size / pk.sum() > 0.7


def test_run_offline(tmp_path):
    m = peaked.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_nvss_vlass"] > 1000  # most of the field is NVSS+VLASS detected
    assert m["n_peaked"] >= 1
    assert m["n_peaked_recovered"] >= 1
    assert m["n_peaked_recovered"] <= m["n_injected_peaked"]
    assert (tmp_path / "results" / "peaked_metrics.json").exists()
    assert (tmp_path / "papers" / "peaked" / "figures" / "curvature.pdf").exists()
    macros = (tmp_path / "papers" / "peaked" / "generated" / "macros.tex").read_text()
    assert r"\pkNpeaked" in macros
