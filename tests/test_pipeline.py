"""Tests for jansky_research.pipeline — offline end-to-end + CSV loader. No network."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import frbstats, pipeline


def test_analyze_returns_metrics():
    m = pipeline.analyze(frbstats.synthetic_catalog(seed=0), source="synthetic")
    assert m["source"] == "synthetic"
    assert "weibull" in m and "energy" in m and "ks" in m
    assert m["weibull"]["k"] > 0


def test_run_offline_writes_all_artifacts(tmp_path):
    metrics = pipeline.run(tmp_path, offline=True)
    assert metrics["source"] == "synthetic"
    # metrics.json
    saved = json.loads((tmp_path / "results" / "metrics.json").read_text())
    assert saved["weibull"]["k"] == metrics["weibull"]["k"]
    # three figures
    figs = sorted((tmp_path / "papers" / "frbstats" / "figures").glob("*.pdf"))
    assert {f.name for f in figs} == {"wait_time_cdf.pdf", "fluence_ccdf.pdf", "dm_by_class.pdf"}
    # macros, with a headline command present
    macros = (tmp_path / "papers" / "frbstats" / "generated" / "macros.tex").read_text()
    assert r"\weibullK" in macros and r"\energyGamma" in macros


def test_build_catalog_offline():
    cat, source = pipeline.build_catalog(offline=True)
    assert source == "synthetic"
    assert cat["repeater"].any()


def test_load_catalog_csv(tmp_path):
    csv_path = tmp_path / "cat.csv"
    csv_path.write_text(
        "tns_name,repeater_name,mjd_400,fluence,dm_fitb,width_fitb\n"
        "FRB1,-,58000.0,1.2,500.0,1.0\n"
        "FRB2,FRB20180916B,58001.0,3.4,520.0,2.0\n"
        "FRB3,FRB20180916B,58001.5,2.1,515.0,1.5\n"
    )
    cat = pipeline.load_catalog_csv(csv_path)
    assert cat["repeater"].tolist() == [False, True, True]
    assert np.allclose(cat["mjd"], [58000.0, 58001.0, 58001.5])
    assert np.allclose(cat["dm"], [500.0, 520.0, 515.0])


def test_cli_offline(tmp_path, capsys):
    assert pipeline._main(["--offline", "--out", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "weibull" in out
