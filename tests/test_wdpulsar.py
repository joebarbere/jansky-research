"""Tests for jansky_research.wdpulsar -- the WD-pulsar candidate radio survey. Offline."""

from __future__ import annotations

import csv
import json

import numpy as np
import pytest

from jansky_research import wdpulsar as wd


def test_candidate_table_loads_and_validates():
    cat = wd.load_candidate_table()
    assert cat["ra_deg"].size == 56
    assert (cat["type"] == "pulsar").sum() == 1
    # the paper's two distinct counts: 26 with no prior Simbad characterisation (abstract),
    # 17 whose classification was determined in this work (Table 2's asterisk)
    assert int(cat["previously_uncharacterized"].sum()) == 26
    assert int(cat["class_this_work"].sum()) == 17
    j = np.where(cat["short_name"] == "J1912-4410")[0][0]
    assert abs(cat["ra_deg"][j] - 288.057) < 0.01
    assert abs(cat["dec_deg"][j] + 44.179) < 0.01


def test_candidate_table_validation_catches_corruption(tmp_path):
    cat_rows = list(csv.DictReader(open(wd.CANDIDATES_LOCAL)))
    bad = tmp_path / "bad.csv"
    with open(bad, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cat_rows[0]))
        w.writeheader()
        w.writerows(cat_rows[:-1])  # drop a row
    with pytest.raises(ValueError, match="55 rows"):
        wd.load_candidate_table(bad)


def test_targets_include_control():
    cat = wd.load_candidate_table()
    targets = wd.candidate_targets(cat)
    assert len(targets) == 57
    assert targets[-1]["name"] == "AR_Sco" and targets[-1]["type"] == "control"
    assert len(wd.candidate_targets(cat, include_control=False)) == 56


def test_injection_roundtrip_recovers_iv_and_blank_is_null():
    rt = wd.injection_roundtrip()
    assert abs(rt["i_out"] - rt["i_in"]) < 5 * 0.25  # within a few noise sigma
    assert abs(rt["v_out"] - rt["v_in"]) < 5 * 0.25
    assert rt["class"] in ("circular", "highly_circular")
    assert rt["blank_i_sig"] < 5.0  # a blank field must not fake a detection


def _sweep_csv(tmp_path, rows):
    p = tmp_path / "sweep.csv"
    fields = ["name", "obs_id", "i_mjy", "e_i", "v_mjy", "e_v"]
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return p


def test_summarize_sweep_detections_limits_and_leakage_vetting(tmp_path):
    targets = [
        {"name": "DET", "ra_deg": 0.0, "dec_deg": -30.0, "type": "polar"},
        {"name": "LIM", "ra_deg": 1.0, "dec_deg": -30.0, "type": "YSO"},
        {"name": "LEAKY", "ra_deg": 2.0, "dec_deg": -30.0, "type": "CV"},
        {"name": "AR_Sco", "ra_deg": 245.447, "dec_deg": -22.886, "type": "control"},
    ]
    rows = [
        # DET: bright I + genuine deep V in epoch 2
        {"name": "DET", "obs_id": "a", "i_mjy": 6.0, "e_i": 0.3, "v_mjy": -0.2, "e_v": 0.3},
        {"name": "DET", "obs_id": "b", "i_mjy": 6.5, "e_i": 0.3, "v_mjy": -2.4, "e_v": 0.3},
        # LIM: nothing at any epoch
        {"name": "LIM", "obs_id": "a", "i_mjy": 0.2, "e_i": 0.25, "v_mjy": 0.1, "e_v": 0.25},
        # LEAKY: enormous I, V formally significant but BELOW 0.6% leakage of I
        {"name": "LEAKY", "obs_id": "a", "i_mjy": 900.0, "e_i": 0.5, "v_mjy": 2.5, "e_v": 0.3},
        # control detected
        {"name": "AR_Sco", "obs_id": "a", "i_mjy": 8.6, "e_i": 0.3, "v_mjy": -2.2, "e_v": 0.3},
    ]
    s = wd.summarize_sweep(_sweep_csv(tmp_path, rows), targets)
    by = {e["name"]: e for e in s["per_target"]}
    assert by["DET"]["i_det"] and by["DET"]["v_det"] and by["DET"]["n_epochs"] == 2
    assert by["DET"]["class"] in ("circular", "highly_circular")
    assert not by["LIM"]["i_det"] and by["LIM"]["v_limit_mjy"] == pytest.approx(0.75)
    assert by["LEAKY"]["i_det"] and not by["LEAKY"]["v_det"]  # leakage-vetted away
    assert s["control_i_det"] and s["control_i_mjy"] == pytest.approx(8.6)
    assert s["n_i_detections"] == 3 and s["n_v_detections"] == 2


def test_summarize_sweep_uncovered_target(tmp_path):
    targets = [{"name": "NOCOV", "ra_deg": 0.0, "dec_deg": 80.0, "type": "YSO"}]
    s = wd.summarize_sweep(_sweep_csv(tmp_path, []), targets)
    assert s["per_target"][0]["n_epochs"] == 0
    assert s["n_measured"] == 0


def test_run_offline_writes_artifacts(tmp_path):
    m = wd.run(str(tmp_path), offline=True)
    assert m["n_candidates"] == 56 and m["n_uncharacterized"] == 26
    assert m["injection"]["class"] in ("circular", "highly_circular")
    saved = json.loads((tmp_path / "results" / "wdpulsar_metrics.json").read_text())
    assert saved["n_candidates"] == 56
    assert (tmp_path / "papers" / "wdpulsar" / "figures" / "wdpulsar.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "wdpulsar" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\wdNCand}{56}" in macros
    assert r"\newcommand{\wdRealNCandIDet}{--}" in macros
    table = (tmp_path / "papers" / "wdpulsar" / "generated" / "limits_table.tex").read_text()
    assert "Auto-generated" in table


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    wd._write_macros(
        {
            "source": "x",
            "is_real": True,
            "n_candidates": 56,
            "n_uncharacterized": 26,
            "n_class_this_work": 17,
            "n_candidate_i_det": None,
        },
        p,
    )
    txt = p.read_text()
    assert (
        r"\newcommand{\wdRealNCandIDet}{--}" in txt and r"\newcommand{\wdSynNCandIDet}{--}" in txt
    )


def test_limits_table_rows(tmp_path):
    m = {
        "per_target": [
            {  # a detection: I shown, V as a limit
                "name": "A",
                "type": "polar",
                "n_epochs": 2,
                "i_mjy": 5.0,
                "v_mjy": None,
                "i_limit_mjy": 0.7,
                "v_limit_mjy": 0.8,
                "class": "weak",
            },
            {  # a non-detection: I shown as a 3sigma limit
                "name": "C",
                "type": "CV",
                "n_epochs": 3,
                "i_mjy": None,
                "v_mjy": None,
                "i_limit_mjy": 0.6,
                "v_limit_mjy": 0.5,
                "class": "nan",
            },
        ]
    }
    p = tmp_path / "t.tex"
    wd._write_limits_table(m, p)
    txt = p.read_text()
    assert "A & polar & 2 & 5.00 & $<$0.80 & weak" in txt
    assert "C & CV & 3 & $<$0.60 & $<$0.50 & nan" in txt
    # the 6-column form: 5 ampersands per data row
    assert all(r.count("&") == 5 for r in txt.splitlines() if r.endswith(r"\\"))
