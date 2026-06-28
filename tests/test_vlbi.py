"""Tests for jansky_research.vlbi — multi-decade VLBI flux variability. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import vlbi


def test_lightcurve_metrics_steady_vs_variable():
    # one steady source (constant) and one with a clear single-epoch flare
    flux = np.array(
        [
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 5.0, 1.0, 1.0],
        ]
    )
    err = np.full_like(flux, 0.05)
    eta, v, pval, nep, mean = vlbi.lightcurve_metrics(flux, err)
    assert nep.tolist() == [5, 5]
    assert eta[1] > eta[0]  # the flaring source is far more significant
    assert v[1] > v[0]  # and higher amplitude
    assert pval[0] > 0.01 and pval[1] < 1e-3  # steady consistent with constant; variable not


def test_lightcurve_metrics_too_few_epochs():
    flux = np.array([[1.0, np.nan, np.nan, 1.0]])  # only 2 finite < MIN_EPOCHS
    err = np.array([[0.1, np.nan, np.nan, 0.1]])
    eta, v, pval, nep, mean = vlbi.lightcurve_metrics(flux, err)
    assert nep[0] == 2
    assert np.isnan(eta[0]) and np.isnan(v[0])


def test_sx_index_flat_and_steep():
    # source 0: equal S and X flux -> alpha = 0; source 1: X brighter -> inverted (alpha > 0)
    flux_s = np.array([[1.0, 1.0], [1.0, 1.0]])
    flux_x = np.array([[1.0, 1.0], [2.0, 2.0]])
    err = np.full_like(flux_s, 0.05)
    alpha, aerr = vlbi.sx_index(flux_s, flux_x, err, err)
    assert abs(alpha[0]) < 1e-9
    assert alpha[1] > 0  # rising to X band
    assert np.all(np.isfinite(aerr))


def test_select_variable_excludes_short_curves():
    # a realistic steady population (the sigma-clip needs a real distribution to set the threshold)
    rng = np.random.default_rng(0)
    n = 200
    eta = np.abs(rng.normal(1.0, 0.3, n))
    v = np.abs(rng.normal(0.05, 0.01, n))
    nep = np.full(n, 10)
    # two clear, well-sampled outliers
    eta[0], v[0] = 80.0, 0.6
    eta[1], v[1] = 120.0, 0.9
    # one outlier just as extreme but with too few epochs -> must be excluded
    eta[2], v[2], nep[2] = 200.0, 1.0, 2
    mask, eta_thr, v_thr = vlbi.select_variable(eta, v, nep)
    assert mask[0] and mask[1]  # the well-sampled outliers are flagged
    assert not mask[2]  # the short curve is excluded despite being the most extreme
    assert np.isfinite(eta_thr) and np.isfinite(v_thr)


def test_synthetic_population_recovers_injected_variables():
    pop = vlbi.synthetic_lightcurves(n_sources=500, seed=1)
    assert pop["flux_x"].shape[0] == 500
    eta, v, pval, nep, mean = vlbi.lightcurve_metrics(pop["flux_x"], pop["err_x"])
    mask, _, _ = vlbi.select_variable(eta, v, nep)
    truth = pop["is_variable"]
    # most injected variables (that have enough epochs) are recovered, with high purity
    completeness = (mask & truth).sum() / truth.sum()
    purity = (mask & truth).sum() / max(mask.sum(), 1)
    assert completeness > 0.6
    assert purity > 0.6


def test_run_offline(tmp_path):
    m = vlbi.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_sources"] == 400
    assert m["n_candidates"] >= 1
    assert m["completeness"] > 0.5 and m["purity"] > 0.5
    assert m["median_alpha_sx"] is not None
    assert (tmp_path / "results" / "vlbi_metrics.json").exists()
    assert (tmp_path / "papers" / "vlbi" / "figures" / "etav.pdf").exists()
    macros = (tmp_path / "papers" / "vlbi" / "generated" / "macros.tex").read_text()
    assert r"\viNcand" in macros and r"\viCompleteness" in macros
