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
            {
                "name": "SECURE_x",
                "period_min": 80.7,
                "n_epochs": 5,
                "believable": True,
                "secure": True,
                "v_mjy": -3.25,
                "class": "circular",
                "offset_arcsec": 0.7,
                "v_limit_mjy": 0.4,
                "handedness_change": None,
            },
            {
                "name": "CAND",
                "period_min": 388.6,
                "n_epochs": 7,
                "believable": True,
                "secure": False,
                "v_mjy": 2.56,
                "class": "highly_circular",
                "offset_arcsec": 3.2,
                "v_limit_mjy": 0.5,
                "handedness_change": None,
            },
            {
                "name": "CONF",
                "period_min": 387.0,
                "n_epochs": 10,
                "believable": False,
                "secure": False,
                "suspect_confusion": True,
                "v_mjy": None,
                "v_limit_mjy": 0.49,
                "handedness_change": None,
            },
            {
                "name": "LIM",
                "period_min": 22.0,
                "n_epochs": 8,
                "believable": False,
                "suspect_confusion": False,
                "v_mjy": None,
                "v_limit_mjy": 0.45,
                "handedness_change": None,
            },
            {"name": "UNCOV", "period_min": 125.5, "n_epochs": 0},
        ]
    }
    lv._write_v_table(m, p)
    txt = p.read_text()
    assert "SECURE" in txt and r"\_" in txt  # underscore escaped in the name
    assert (
        "circular" in txt and "_" not in txt.split("SECURE")[1].split("&")[4]
    )  # class no underscore
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


# --- VAST extension -------------------------------------------------------------------


def _vast_csv(tmp_path):
    import csv as _csv

    from jansky_research.lptv import J1839_PERIOD_S, J1839_T0_MJD

    p = tmp_path / "vast.csv"
    fields = [
        "name",
        "epoch",
        "obs_id",
        "band",
        "epoch_mjd",
        "duration_s",
        "i_mjy",
        "e_i",
        "v_mjy",
        "e_v",
        "offset_arcsec",
        "note",
    ]
    period_d = J1839_PERIOD_S / 86400.0
    # an epoch whose MIDPOINT lands exactly 100 cycles after T0 (main-pulse phase 0)
    dur = 726.0
    t_min = J1839_T0_MJD + 100 * period_d - dur / 2.0 / 86400.0
    rows = [
        # J1839 detection at main-pulse phase: 60 mJy, 40% circular
        dict(
            name="ASKAP J183950.5-075635",
            epoch="low",
            obs_id="A-1",
            band="low",
            epoch_mjd=f"{t_min:.5f}",
            duration_s="726.0",
            i_mjy="60.0",
            e_i="0.5",
            v_mjy="24.0",
            e_v="0.4",
            offset_arcsec="0.5",
            note="",
        ),
        # same source, quiet epoch
        dict(
            name="ASKAP J183950.5-075635",
            epoch="low",
            obs_id="A-2",
            band="low",
            epoch_mjd="60500.10000",
            duration_s="726.0",
            i_mjy="0.1",
            e_i="0.4",
            v_mjy="-0.2",
            e_v="0.35",
            offset_arcsec="1.0",
            note="",
        ),
        # bright-I, quiet-V epoch (GPM-like)
        dict(
            name="GPM J1839-10",
            epoch="low",
            obs_id="B-1",
            band="low",
            epoch_mjd="60600.20000",
            duration_s="726.0",
            i_mjy="6.0",
            e_i="0.4",
            v_mjy="0.1",
            e_v="0.35",
            offset_arcsec="0.8",
            note="",
        ),
        # unreleased row
        dict(
            name="GPM J1839-10",
            epoch="low",
            obs_id="B-2",
            band="low",
            epoch_mjd="60000.00000",
            duration_s="726.0",
            i_mjy="nan",
            e_i="nan",
            v_mjy="nan",
            e_v="nan",
            offset_arcsec="nan",
            note="unreleased: obscore obs_release_date is NULL",
        ),
        # off-mosaic nan row
        dict(
            name="GPM J1839-10",
            epoch="low",
            obs_id="B-3",
            band="low",
            epoch_mjd="60010.00000",
            duration_s="726.0",
            i_mjy="nan",
            e_i="nan",
            v_mjy="nan",
            e_v="nan",
            offset_arcsec="nan",
            note="",
        ),
    ]
    with p.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return p


def test_fold_phase_exact_cycles():
    from jansky_research.lptv import J1839_PERIOD_S, J1839_T0_MJD, fold_phase

    ph, pe = fold_phase(J1839_T0_MJD + 100 * J1839_PERIOD_S / 86400.0)
    assert ph == pytest.approx(0.0, abs=1e-6) or ph == pytest.approx(1.0, abs=1e-6)
    # period-error contribution: 100 cycles * 0.332 s / P
    assert pe == pytest.approx(100 * 0.332 / J1839_PERIOD_S, rel=1e-6)


def test_summarize_vast_reduces_and_folds(tmp_path):
    from jansky_research.lptv import summarize_vast

    m = summarize_vast(_vast_csv(tmp_path))
    assert m["n_rows"] == 5
    assert m["n_measured"] == 3
    assert m["n_unreleased"] == 1
    assert m["n_offmosaic_nan"] == 1
    assert m["n_failed"] == 0
    assert m["n_sources_covered"] == 2
    # exactly one leakage-vetted V detection, at main-pulse phase, with its phase attached
    assert len(m["v_detections"]) == 1
    det = m["v_detections"][0]
    assert det["obs_id"] == "A-1"
    assert det["v_over_i"] == pytest.approx(0.4, abs=0.01)
    assert min(det["phase"], 1 - det["phase"]) < 0.001
    # 100 cycles x 0.332 s / P = 0.00143, rounded to 4 dp in the record
    assert det["phase_err_period"] == pytest.approx(0.0014, abs=2e-4)
    # total error adds exposure smearing (726/P/sqrt(12)=0.009) and the 0.013 anchor
    assert det["phase_err_total"] == pytest.approx(0.0159, abs=5e-4)
    # the bright-I epoch is reported separately, not as a V detection
    assert [e["obs_id"] for e in m["i_bright_epochs"]] == ["B-1"]
    # per-target limits exist for the quiet source
    assert m["per_target"]["ASKAP J183950.5-075635"]["n_epochs"] == 2


def test_summarize_vast_decision_boundaries(tmp_path):
    """Rows that sit AT the thresholds: 4.9-sigma V excluded, 5.1 included; a V just under
    the leakage floor excluded; negative-I rows still leakage-guarded via |I|."""
    import csv as _csv

    from jansky_research.lptv import summarize_vast

    fields = [
        "name",
        "epoch",
        "obs_id",
        "band",
        "epoch_mjd",
        "duration_s",
        "i_mjy",
        "e_i",
        "v_mjy",
        "e_v",
        "offset_arcsec",
        "note",
    ]
    rows = [
        # 4.9 sigma: below det threshold
        dict(name="S1", obs_id="a", i_mjy="1.0", e_i="0.3", v_mjy="0.98", e_v="0.2"),
        # 5.1 sigma and far above leakage: detected
        dict(name="S1", obs_id="b", i_mjy="1.0", e_i="0.3", v_mjy="1.02", e_v="0.2"),
        # high sigma but |V| below 0.6% of I: leakage-vetoed
        dict(name="S2", obs_id="c", i_mjy="2000.0", e_i="0.3", v_mjy="10.0", e_v="0.2"),
        # negative I with significant V: leakage guard uses |I|, small |V| still passes
        dict(name="S3", obs_id="d", i_mjy="-1.0", e_i="0.3", v_mjy="1.5", e_v="0.2"),
    ]
    for r in rows:
        r.setdefault("epoch", "low")
        r.setdefault("band", "low")
        r.setdefault("epoch_mjd", "60500.00000")
        r.setdefault("duration_s", "726.0")
        r.setdefault("offset_arcsec", "0.5")
        r.setdefault("note", "")
    p = tmp_path / "b.csv"
    with p.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    m = summarize_vast(p)
    assert [d["obs_id"] for d in m["v_detections"]] == ["b", "d"]
    # the 5.1-sigma detection carries v_over_i; the negative-I one records None, not a crash
    assert m["v_detections"][1]["v_over_i"] is None


def test_fold_phase_noninteger_cycles():
    """A midpoint 100.25 cycles after T0 must fold to 0.25 — catches sign/unit errors in
    the duration/2 midpoint arithmetic feeding the paper's phases."""
    from jansky_research.lptv import J1839_PERIOD_S, J1839_T0_MJD, fold_phase

    ph, _ = fold_phase(J1839_T0_MJD + 100.25 * J1839_PERIOD_S / 86400.0)
    assert ph == pytest.approx(0.25, abs=1e-6)
