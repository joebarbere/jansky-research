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


def test_run_offline_writes_catalogue(tmp_path):
    m = ec.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic-day"
    assert m["n_scanned"] == 8 and m["n_bursts"] == 3
    assert (tmp_path / "results" / "ecallisto_catalog.csv").exists()
    assert (tmp_path / "results" / "ecallisto_metrics.json").exists()
    assert (tmp_path / "papers" / "ecallisto_pipeline" / "figures" / "ecallisto.pdf").exists()
    macros = (tmp_path / "papers" / "ecallisto_pipeline" / "generated" / "macros.tex").read_text()
    assert r"\ecNbursts" in macros
