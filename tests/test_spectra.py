"""Tests for jansky_research.spectra — spectral index, classify, cross-match, USS. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import spectra


def test_spectral_index_recovers_injected_alpha():
    nu_lo, nu_hi = 147.5, 1400.0
    alpha_true = -1.4
    s_lo = np.array([100.0])
    s_hi = s_lo * (nu_hi / nu_lo) ** alpha_true
    alpha, e = spectra.spectral_index(s_lo, nu_lo, s_hi, nu_hi, 0.1 * s_lo, 0.05 * s_hi)
    assert np.isclose(alpha[0], alpha_true, atol=1e-6)
    assert e[0] > 0 and np.isfinite(e[0])


def test_classify_thresholds():
    assert spectra.classify(-1.6) == "uss"
    assert spectra.classify(-0.8) == "steep"
    assert spectra.classify(-0.2) == "flat"
    assert spectra.classify(0.4) == "inverted"
    assert spectra.classify(spectra.USS_THRESHOLD - 0.01) == "uss"


def test_crossmatch_respects_radius():
    ra_lo = np.array([180.0, 181.0])
    dec_lo = np.array([30.0, 30.0])
    # one counterpart 3" away, one 60" away (outside a 15" radius)
    ra_hi = np.array([180.0 + 3 / 3600.0, 181.0 + 60 / 3600.0])
    dec_hi = np.array([30.0, 30.0])
    i_lo, i_hi, sep = spectra.crossmatch(ra_lo, dec_lo, ra_hi, dec_hi, radius_arcsec=15.0)
    assert i_lo.tolist() == [0]  # only the close pair survives
    assert sep[0] < 15.0


def test_find_uss_recovers_injected_population():
    low, high = spectra.synthetic_field(n=300, f_uss=0.05, f_inverted=0.05, seed=1)
    res = spectra.find_uss(low, high)
    assert res["alpha"].size > 250  # nearly all cross-match (small jitter)
    n_uss = int(res["is_uss"].sum())
    assert 8 <= n_uss <= 25  # ~15 injected USS recovered
    # the USS sources really are steep
    assert np.median(res["alpha"][res["is_uss"]]) < spectra.USS_THRESHOLD
    assert set(np.unique(res["cls"])) <= {"uss", "steep", "flat", "inverted"}


def test_synthetic_field_shapes():
    low, high = spectra.synthetic_field(n=50, seed=0)
    for d in (low, high):
        assert {"ra", "dec", "flux", "eflux"} <= set(d)
        assert d["ra"].size == 50


def test_analyze_and_run_offline(tmp_path):
    low, high = spectra.synthetic_field(n=300, seed=2)
    m = spectra.analyze(spectra.find_uss(low, high), source="synthetic")
    assert m["n_matched"] > 250 and m["n_uss"] >= 1
    assert m["alpha_min"] < m["alpha_median"]

    metrics = spectra.run(out=str(tmp_path), offline=True)
    assert metrics["source"] == "synthetic"
    assert (tmp_path / "results" / "uss_metrics.json").exists()
    assert (tmp_path / "results" / "uss_candidates.csv").exists()
    figs = {p.name for p in (tmp_path / "papers" / "spectra" / "figures").glob("*.pdf")}
    assert {"alpha_hist.pdf", "alpha_vs_flux.pdf"} <= figs
    assert (tmp_path / "papers" / "spectra" / "generated" / "macros.tex").exists()
