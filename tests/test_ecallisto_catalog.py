"""Tests for jansky_research.ecallisto_catalog -- the e-Callisto day-scan worker. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import ecallisto_catalog as ec


def test_synthetic_day_shape():
    specs = ec.synthetic_day(n_stations=6, n_bursts=2, seed=0)
    assert len(specs) == 6
    for _name, spec in specs:
        assert {"data", "freqs", "times"} <= set(spec)
        assert spec["data"].ndim == 2


def test_scan_spectrum_flags_a_burst_and_not_noise():
    from jansky_research import solarbursts

    burst = ec.scan_spectrum(solarbursts.synthetic_burst(seed=1))
    assert burst["is_burst"] is True
    assert burst["drift_mhz_s"] < 0  # type III drifts high -> low
    rng = np.random.default_rng(0)
    t = solarbursts.synthetic_burst(seed=1)
    noise = {"data": rng.normal(0, 1, t["data"].shape), "freqs": t["freqs"], "times": t["times"]}
    assert ec.scan_spectrum(noise)["is_burst"] is False


def test_scan_day_recovers_injected_bursts():
    specs = ec.synthetic_day(n_stations=8, n_bursts=3, seed=0)
    rows = ec.scan_day_specs(specs)
    assert len(rows) == 8
    assert sum(r["is_burst"] for r in rows) == 3  # exactly the injected bursts
    # the flagged stations are the first three (the burst stations)
    flagged = {r["station"] for r in rows if r["is_burst"]}
    assert flagged == {"STATION00", "STATION01", "STATION02"}


def test_synthetic_coincident_day_and_coincidence():
    # a real burst at 4 stations (same UT) + 3 single-station RFI (distinct times) + 3 quiet
    specs = ec.synthetic_coincident_day(n_coincident=4, n_rfi=3, n_quiet=3, seed=0)
    assert len(specs) == 10
    rows = ec.scan_day_specs(specs)
    assert sum(r["is_burst"] for r in rows) == 7  # 4 real + 3 RFI detected as candidates
    events = ec.coincident_events(rows, dt_tol_s=60.0, min_stations=2)
    # coincidence confirms exactly the real burst and rejects the single-station RFI
    assert len(events) == 1
    assert events[0]["n_stations"] == 4
    assert events[0]["median_drift_mhz_s"] < 0


def test_coincident_events_on_rows():
    # two stations at ~the same time -> 1 event; an isolated single-station candidate -> rejected
    rows = [
        {"station": "A", "is_burst": True, "t_peak_s": 300.0, "drift_mhz_s": -6.0},
        {"station": "B", "is_burst": True, "t_peak_s": 320.0, "drift_mhz_s": -6.5},
        {"station": "C", "is_burst": True, "t_peak_s": 700.0, "drift_mhz_s": -7.0},  # lone RFI
        {"station": "D", "is_burst": False, "t_peak_s": 305.0, "drift_mhz_s": None},  # not a burst
    ]
    events = ec.coincident_events(rows, dt_tol_s=60.0, min_stations=2)
    assert len(events) == 1 and events[0]["n_stations"] == 2
    assert sorted(events[0]["stations"]) == ["A", "B"]


def test_run_offline_writes_catalogue(tmp_path):
    m = ec.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic-day"
    assert m["n_scanned"] == 10 and m["n_bursts"] == 7
    # the coincidence QC confirms one real burst and rejects the single-station RFI
    assert m["n_events"] == 1 and m["max_event_stations"] == 4 and m["n_rfi_rejected"] == 3
    assert (tmp_path / "results" / "ecallisto_catalog.csv").exists()
    assert (tmp_path / "results" / "ecallisto_metrics.json").exists()
    assert (tmp_path / "papers" / "ecallisto_pipeline" / "figures" / "ecallisto.pdf").exists()
    macros = (tmp_path / "papers" / "ecallisto_pipeline" / "generated" / "macros.tex").read_text()
    assert r"\ecNbursts" in macros and r"\ecNevents" in macros
