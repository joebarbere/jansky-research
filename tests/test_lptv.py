"""Tests for jansky_research.lptv -- LPT v3 catalogue + Stokes-V forced photometry. Offline."""

from __future__ import annotations

import csv
import json

import pytest

from jansky_research import lptv as lv


def test_v3_catalogue_has_sixteen_with_2026_rows():
    pos = lv.lpt_positions()
    assert len(pos) == 16
    names = {p["name"] for p in pos}
    for n in ("ASKAP J142431.2-612611", "ASKAP J165130.3-450520", "ASKAP J170036.6-445758"):
        assert n in names


def test_new_row_coordinates_decode_from_name():
    # the lpt provenance discipline: RA/Dec must match the source-name sexagesimal
    pos = {p["name"]: p for p in lv.lpt_positions()}
    j1424 = pos["ASKAP J142431.2-612611"]
    ra_name = 15 * (14 + 24 / 60 + 31.2 / 3600)
    dec_name = -(61 + 26 / 60 + 11 / 3600)
    assert abs(j1424["ra_deg"] - ra_name) * 3600 < 2.0
    assert abs(j1424["dec_deg"] - dec_name) * 3600 < 2.0


def test_catalogue_stats_binary_boundary_not_significant():
    c = lv.catalogue_stats()
    assert c["n_lpt"] == 16
    # the plan's headline question: does the ~78-min WD-binary period boundary move at N=16?
    assert not c["binary_boundary_significant"]  # p >= 0.05, still not significant
    assert 0.0 <= c["period_split_p"] <= 1.0


def test_injection_roundtrip_recovers_v_and_blank_is_null():
    rt = lv.injection_roundtrip()
    assert abs(rt["v_out"] - rt["v_in"]) < 5 * 0.2  # within a few noise sigma
    assert rt["class"] in ("circular", "highly_circular")
    assert rt["handedness"] == "LCP"  # injected V<0
    assert abs(rt["blank_v_sig"]) < 5.0  # a blank field is not a fake detection


def test_handedness_changes_detects_flip():
    rows = [
        {"v_det": True, "v_mjy": "2.5"},
        {"v_det": True, "v_mjy": "-3.0"},  # sign flip -> RCP to LCP
    ]
    assert lv.handedness_changes(rows) == "flip"
    same = [{"v_det": True, "v_mjy": "2.5"}, {"v_det": True, "v_mjy": "3.1"}]
    assert lv.handedness_changes(same) is None
    # non-detections don't count toward a flip
    mixed = [{"v_det": True, "v_mjy": "2.5"}, {"v_det": False, "v_mjy": "-3.0"}]
    assert lv.handedness_changes(mixed) is None


def _sweep_csv(tmp_path, rows):
    p = tmp_path / "sweep.csv"
    fields = ["name", "epoch", "i_mjy", "e_i", "v_mjy", "e_v", "offset_arcsec"]
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({"offset_arcsec": "1.0", **r})  # small offset unless the row overrides
    return p


def test_confusion_veto_flags_offset_outlier(tmp_path):
    # a source detected only in one epoch that is BOTH far off-centre AND >>brighter than its
    # other epochs is a nearby confusing source, not the target -> suspect, not a detection
    # (the real ASKAP J183950 case: 240 mJy at 5.3" vs a ~0.5 mJy source)
    targets = [{"name": "CONF", "ra_deg": 0.0, "dec_deg": -50.0, "period_s": 6000.0}]
    rows = [
        {
            "name": "CONF",
            "epoch": "low",
            "i_mjy": "0.5",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        {
            "name": "CONF",
            "epoch": "low",
            "i_mjy": "0.6",
            "e_i": "0.2",
            "v_mjy": "-0.1",
            "e_v": "0.2",
        },
        {
            "name": "CONF",
            "epoch": "low",
            "i_mjy": "240.0",
            "e_i": "0.3",
            "v_mjy": "96.0",
            "e_v": "0.3",
            "offset_arcsec": "5.3",
        },
    ]
    s = lv.summarize_v_sweep(_sweep_csv(tmp_path, rows), targets)
    c = s["per_target"][0]
    assert c["v_det"] and c["suspect_confusion"] and not c["believable"]
    assert s["n_v_detections"] == 0 and s["n_suspect_confusion"] == 1
    assert c["class"] == "nan"  # a suspect is not classified


def test_confusion_veto_keeps_oncentre_bright_burst(tmp_path):
    # a bright burst that is well-centred (small offset) is a believable detection even as a large
    # I outlier -- a real flaring source (the accreting-CV LPT at 0.7", 15% circular)
    targets = [{"name": "BURST", "ra_deg": 0.0, "dec_deg": -50.0, "period_s": 5000.0}]
    rows = [
        {
            "name": "BURST",
            "epoch": "mid",
            "i_mjy": "0.4",
            "e_i": "0.2",
            "v_mjy": "0.0",
            "e_v": "0.2",
        },
        {
            "name": "BURST",
            "epoch": "mid",
            "i_mjy": "21.6",
            "e_i": "0.2",
            "v_mjy": "-3.2",
            "e_v": "0.15",
            "offset_arcsec": "0.7",
        },
    ]
    s = lv.summarize_v_sweep(_sweep_csv(tmp_path, rows), targets)
    b = s["per_target"][0]
    assert b["believable"] and not b["suspect_confusion"] and b["secure"]
    assert s["n_v_detections"] == 1 and s["n_v_secure"] == 1


def test_offset_detection_is_candidate_not_secure(tmp_path):
    # the delicate real case (J1651): a believable detection at 3.2" offset, 8.6x I outlier --
    # inside neither confusion threshold, so kept, but off-centre -> candidate, not secure
    targets = [{"name": "CAND", "ra_deg": 0.0, "dec_deg": -50.0, "period_s": 23000.0}]
    rows = [
        {
            "name": "CAND",
            "epoch": "mid",
            "i_mjy": "0.5",
            "e_i": "0.2",
            "v_mjy": "0.0",
            "e_v": "0.2",
        },
        {
            "name": "CAND",
            "epoch": "mid",
            "i_mjy": "0.5",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        {
            "name": "CAND",
            "epoch": "mid",
            "i_mjy": "4.37",
            "e_i": "0.24",
            "v_mjy": "2.56",
            "e_v": "0.20",
            "offset_arcsec": "3.2",
        },  # 8.6x median, 3.2" -> kept but candidate
    ]
    s = lv.summarize_v_sweep(_sweep_csv(tmp_path, rows), targets)
    c = s["per_target"][0]
    assert c["believable"] and not c["suspect_confusion"] and not c["secure"]
    assert s["n_v_detections"] == 1 and s["n_v_secure"] == 0 and s["n_v_candidate"] == 1


def test_summarize_v_sweep_detections_limits_and_leakage(tmp_path):
    targets = [
        {"name": "DET", "ra_deg": 0.0, "dec_deg": -50.0, "period_s": 1200.0},
        {"name": "LIM", "ra_deg": 1.0, "dec_deg": -50.0, "period_s": 3600.0},
        {"name": "LEAKY", "ra_deg": 2.0, "dec_deg": -50.0, "period_s": 600.0},
    ]
    rows = [
        # DET: genuine deep V well above the leakage floor of its I
        {
            "name": "DET",
            "epoch": "low1",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "-2.0",
            "e_v": "0.2",
        },
        {
            "name": "DET",
            "epoch": "low2",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        # LIM: nothing
        {
            "name": "LIM",
            "epoch": "low1",
            "i_mjy": "0.1",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        # LEAKY: huge I, V formally 5sigma but below 0.6% leakage of I
        {
            "name": "LEAKY",
            "epoch": "low1",
            "i_mjy": "500.0",
            "e_i": "0.3",
            "v_mjy": "1.5",
            "e_v": "0.2",
        },
    ]
    s = lv.summarize_v_sweep(_sweep_csv(tmp_path, rows), targets)
    by = {e["name"]: e for e in s["per_target"]}
    assert by["DET"]["v_det"] and by["DET"]["class"] in ("circular", "highly_circular")
    assert not by["LIM"]["v_det"] and by["LIM"]["v_limit_mjy"] == pytest.approx(0.6)
    assert not by["LEAKY"]["v_det"]  # leakage-vetted away
    assert s["n_v_detections"] == 1 and s["n_measured"] == 3


def test_summarize_v_sweep_handedness_flip_counted(tmp_path):
    targets = [{"name": "FLIP", "ra_deg": 0.0, "dec_deg": -50.0, "period_s": 900.0}]
    rows = [
        {
            "name": "FLIP",
            "epoch": "low1",
            "i_mjy": "4.0",
            "e_i": "0.2",
            "v_mjy": "2.5",
            "e_v": "0.2",
        },
        {
            "name": "FLIP",
            "epoch": "mid",
            "i_mjy": "4.0",
            "e_i": "0.2",
            "v_mjy": "-2.6",
            "e_v": "0.2",
        },
    ]
    s = lv.summarize_v_sweep(_sweep_csv(tmp_path, rows), targets)
    assert s["n_handedness_flips"] == 1


def test_write_v_table_all_row_types(tmp_path):
    p = tmp_path / "vt.tex"
    m = {
        "per_target": [
            {"name": "SECURE_x", "period_min": 80.7, "n_epochs": 5, "believable": True,
             "secure": True, "v_mjy": -3.25, "class": "circular", "offset_arcsec": 0.7,
             "v_limit_mjy": 0.4, "handedness_change": None},
            {"name": "CAND", "period_min": 388.6, "n_epochs": 7, "believable": True,
             "secure": False, "v_mjy": 2.56, "class": "highly_circular", "offset_arcsec": 3.2,
             "v_limit_mjy": 0.5, "handedness_change": None},
            {"name": "CONF", "period_min": 387.0, "n_epochs": 10, "believable": False,
             "secure": False, "suspect_confusion": True, "v_mjy": None, "v_limit_mjy": 0.49,
             "handedness_change": None},
            {"name": "LIM", "period_min": 22.0, "n_epochs": 8, "believable": False,
             "suspect_confusion": False, "v_mjy": None, "v_limit_mjy": 0.45,
             "handedness_change": None},
            {"name": "UNCOV", "period_min": 125.5, "n_epochs": 0},
        ]
    }
    lv._write_v_table(m, p)
    txt = p.read_text()
    assert "SECURE" in txt and r"\_" in txt  # underscore escaped in the name
    assert "circular" in txt and "_" not in txt.split("SECURE")[1].split("&")[4]  # class no underscore
    assert "cand., 3.2" in txt  # candidate flagged with offset
    assert "confused" in txt  # suspect reported as a flagged limit
    assert "uncovered" in txt  # zero-epoch target


def test_run_offline_writes_artifacts(tmp_path):
    m = lv.run(str(tmp_path), offline=True)
    assert m["n_lpt"] == 16 and not m["binary_boundary_significant"]
    assert m["injection"]["class"] in ("circular", "highly_circular")
    saved = json.loads((tmp_path / "results" / "lptv_metrics.json").read_text())
    assert saved["n_lpt"] == 16
    assert (tmp_path / "papers" / "lptv" / "figures" / "lptv.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "lptv" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\lvNLpt}{16}" in macros
    assert r"\newcommand{\lvRealNVDet}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    lv._write_macros(
        {
            "source": "x",
            "is_real": True,
            "n_lpt": 16,
            "n_wd_binary": 7,
            "median_period_min": 73.4,
            "period_split_p": 0.52,
            "injection": {"class": "circular"},
            "n_v_detections": None,
        },
        p,
    )
    txt = p.read_text()
    assert r"\newcommand{\lvRealNVDet}{--}" in txt and r"\newcommand{\lvSynNVDet}{--}" in txt
