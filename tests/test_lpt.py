"""Tests for jansky_research.lpt -- the LPT population catalogue. Offline (vendored CSV)."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import lpt


def test_load_sample_shapes_and_flags():
    s = lpt.load_sample()
    assert s["period_s"].size == 16  # v3: +3 2026 discoveries (J1424, J1651, J1700)
    assert (s["period_s"] > 60).all()  # every period > 1 minute (unit sanity)
    assert s["pdot_is_measurement"].sum() == 2  # CHIME J0630+25 + CHIME/ILT J1634+44
    assert s["is_wd_binary"].sum() == 7
    # the one clear spin-up is negative
    neg = s["pdot"][s["pdot_is_measurement"]] < 0
    assert neg.sum() == 1


def test_population_table_death_line_headline():
    s = lpt.load_sample()
    pop = lpt.population_table(s)
    assert pop["n_lpt"] == 16  # v3
    # every Pdot-constrained object sits below the pulsar death line -- the class puzzle
    # (the 3 new rows carry no Pdot constraint, so the 9/9 headline is unchanged)
    assert pop["n_below_death_line"] == pop["n_pdot_constrained"] == 9
    assert pop["period_min_min"] == 7.0 and pop["period_max_hr"] > 6


def test_period_split_stat_honest_at_small_n():
    s = lpt.load_sample()
    out = lpt.period_split_stat(s["period_s"], s["is_wd_binary"])
    assert out["delta_log_median"] > 0  # WD binaries do sit at longer periods...
    assert out["p_perm"] > 0.05  # ...but NOT significantly at N=16 (the honest result)


def test_split_stat_round_trips_injected_split():
    p, wd = lpt.synthetic_lpt_population(seed=0)
    out = lpt.period_split_stat(p, wd)
    assert out["p_perm"] < 0.05  # a REAL split registers -> the real non-detection is informative


def test_split_stat_degenerate():
    out = lpt.period_split_stat(np.array([100.0, 200.0]), np.array([True, False]))
    assert np.isnan(out["p_perm"])


def test_run_writes_artifacts(tmp_path):
    m = lpt.run(str(tmp_path), offline=True)
    assert m["n_lpt"] == 16 and m["n_pdot_measurements"] == 2
    saved = json.loads((tmp_path / "results" / "lpt_metrics.json").read_text())
    assert saved == m
    assert (tmp_path / "papers" / "lpt" / "figures" / "lpt_ppdot.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "lpt" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\lptNdeath}{9}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    lpt._write_macros({"n_lpt": None}, p)
    assert r"\newcommand{\lptN}{--}" in p.read_text()


def test_split_is_exact_at_this_n():
    s = lpt.load_sample()
    out = lpt.period_split_stat(s["period_s"], s["is_wd_binary"])
    assert out["method"].startswith("exact")
    # 16 choose 7 partitions, enumerated -- the published 0.5219 was Monte Carlo seed noise
    assert "11440" in out["method"]


def test_label_sensitivity_unknowns_matter():
    s = lpt.load_sample()
    rows = {r["labelling"]: r for r in lpt.label_sensitivity(s)}
    prim = rows["unknowns_excluded"]
    pub = rows["as_published_unknowns_with_rest"]
    assert prim["n_wd"] == 7 and prim["n_rest"] == 6  # three unknowns are OUT of the primary
    assert pub["n_rest"] == 9
    # counting two of the three longest periods with the null class shrinks the contrast
    assert prim["delta_log_median"] > 2 * pub["delta_log_median"]
    # every labelling is still non-significant -- the conclusion survives, the effect size moves
    assert all(r["p"] > 0.05 for r in rows.values())


def test_death_line_margins_are_structural():
    s = lpt.load_sample()
    dm = lpt.death_line_margins(s)
    assert len(dm["objects"]) == 9
    # no achievable measurement could have landed above the line: min margin is tens
    assert dm["min_margin"] > 10
    # the claim-carrying subset: every no-companion constrained object is below
    assert dm["n_below_nonbinary"] == dm["n_nonbinary_constrained"] == 6
    # and the count is 9/9 across the whole death-valley constant range
    assert all(row["n_below"] == 9 for row in dm["valley_sweep"])
    assert min(row["min_margin"] for row in dm["valley_sweep"]) > 1


def test_split_power_is_size_at_observed_offset():
    # at the observed contrast the exact test has essentially no power at N=16 --
    # the honest statement replacing "the test has power"
    weak = lpt.split_power(0.18, 7, 9, n_sims=150)
    assert weak["power"] < 0.2
    # while a large disjoint-scale offset is comfortably detectable
    strong = lpt.split_power(1.5, 7, 9, n_sims=150)
    assert strong["power"] > 0.6


def test_csv_matches_pinned_literature():
    """Cell-by-cell pin of the audited values: any future CSV edit must change this test too.

    Values audited against the discovery papers during the 2026-08-21 lptduty ephemeris audit
    (data/lpt_sample.csv flags column records the audit); sources are the discovery_arxiv IDs.
    """
    pins = {
        # name: (period_s, pdot, pdot_type, discovery_arxiv)
        "GCRT J1745-3009": (4620.72, None, "none", "astro-ph/0503052"),
        "GLEAM-X J162759.5-523504.3": (1091.1690, 1.2e-9, "upper_limit", "2503.08033"),
        "GPM J1839-10": (1318.1957, 3.6e-13, "upper_limit", "2503.08036"),
        "ASKAP J1935+2148": (3225.313, 1.2e-10, "upper_limit", "2407.12266"),
        "CHIME J0630+25": (421.35584, 5.2e-12, "measurement_disputed", "2407.07480"),
        "ASKAP J1832-0911": (2656.247, 9.8e-10, "upper_limit", "2411.16606"),
        "ILT J1101+5521": (7531.78, 1.71e-11, "upper_limit", "2408.11536"),
        "GLEAM-X J0704-37": (10496.5575, 3.9e-11, "upper_limit", "2408.15757"),
        "ASKAP J183950.5-075635": (23221.74, 1.6e-7, "upper_limit", "2501.09133"),
        "ASKAP J144834-685644": (5631.07, 2.2e-8, "upper_limit", "2507.13453"),
        "CHIME J1634+44": (841.24, -9.03e-12, "measurement", "2507.05139"),
        "ASKAP J175534.9-252749.1": (4186.3285, -1.0e-11, "consistent_zero", "2507.14448"),
        "ASKAP J174508.9-505149": (4842.0, None, "none", "2606.04232"),
        "ASKAP J142431.2-612611": (2147.27, None, "none", "2603.07857"),
        "ASKAP J165130.3-450520": (23317.9, None, "none", "2606.20067"),
        "ASKAP J170036.6-445758": (16895.9, None, "none", "2606.20067"),
    }
    import csv as _csv

    rows = {r["name"]: r for r in _csv.DictReader(open(lpt.SAMPLE_CSV))}
    assert set(rows) == set(pins)
    for name, (p, pd, ptype, arx) in pins.items():
        r = rows[name]
        assert float(r["period_s"]) == p, name
        if pd is None:
            assert r["pdot_s_s"] == "", name
        else:
            assert float(r["pdot_s_s"]) == pd, name
        assert (r["pdot_type"] or "none") == ptype, name
        assert r["discovery_arxiv"] == arx, name


def test_lptxray_coordinate_caches_match_the_csv():
    """The lptxray JSONs embed per-object coordinates; they must not silently diverge."""
    import json as _json
    from pathlib import Path

    s = lpt.load_sample()
    by_name = {
        str(n): (float(r), float(d)) for n, r, d in zip(s["name"], s["ra"], s["dec"], strict=True)
    }
    root = Path(lpt.SAMPLE_CSV).parents[1]
    for fname in ("data/lptxray/coverage.json", "data/lptxray/cones.json"):
        f = root / fname
        if not f.exists():
            continue
        data = _json.loads(f.read_text())
        items = data.items() if isinstance(data, dict) else []
        for name, entry in items:
            if name not in by_name or not isinstance(entry, dict):
                continue
            meta = entry.get("meta", entry)
            ra, dec = meta.get("ra"), meta.get("dec")
            if ra is None or dec is None:
                continue
            assert abs(ra - by_name[name][0]) < 1e-4, (fname, name)
            assert abs(dec - by_name[name][1]) < 1e-4, (fname, name)


def test_catalogue_table_is_generated(tmp_path):
    m = lpt.run(str(tmp_path), offline=True)
    table = (tmp_path / "papers" / "lpt" / "generated" / "table.tex").read_text()
    assert table.count("\\\\") >= 16  # one row per object
    assert "unknown" in table  # unknown companion status is shown, not folded into "no"
    for key in ("hyman2005", "wang2026vaster", "rose2026", "caleb2024"):
        assert key in table
    assert m["split_method"].startswith("exact")


def test_lptv_ephemeris_constant_matches_the_csv():
    """lptv hard-codes J1839's period as its published-ephemeris anchor; it must equal the CSV."""
    from jansky_research.lptv import J1839_PERIOD_S

    s = lpt.load_sample()
    csv_p = float(s["period_s"][list(s["name"]).index("ASKAP J183950.5-075635")])
    assert J1839_PERIOD_S == csv_p
