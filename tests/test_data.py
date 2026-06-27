"""Tests for jansky_research.data — registry, cache, offline fallback. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import data


def test_registry_well_formed():
    assert data.list_datasets() == sorted(data.DATASETS)
    for name, spec in data.DATASETS.items():
        assert spec.name == name
        assert spec.url.startswith("http")
        assert spec.filename
        assert spec.category in {"small", "large"}


def test_data_dir_env_override(tmp_path, monkeypatch):
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path / "cache"))
    d = data.data_dir()
    assert d == (tmp_path / "cache")
    assert d.is_dir()


def test_fetch_unknown_raises():
    with pytest.raises(KeyError):
        data.fetch("does-not-exist")


def test_fetch_returns_cached_without_network(tmp_path, monkeypatch):
    # A pre-existing cached file is returned without any download attempt.
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path))
    name = data.list_datasets()[0]
    cached = tmp_path / data.DATASETS[name].filename
    cached.write_text("cached")

    def _boom(*a, **k):  # would fail if a download were attempted
        raise AssertionError("network should not be touched for a cached file")

    monkeypatch.setattr(data, "_download", _boom)
    assert data.fetch(name) == cached


def test_fetch_wraps_download_failure(tmp_path, monkeypatch):
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(
        data, "_download", lambda *a, **k: (_ for _ in ()).throw(OSError("offline"))
    )
    with pytest.raises(RuntimeError, match="synthetic"):
        data.fetch(data.list_datasets()[0])


def test_fetch_monkeypatched_download(tmp_path, monkeypatch):
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path))
    name = data.list_datasets()[0]

    def fake_download(url, target):
        target.write_text("downloaded")

    monkeypatch.setattr(data, "_download", fake_download)
    path = data.fetch(name, force=True)
    assert path.read_text() == "downloaded"


def test_synthetic_dynamic_spectrum_shape_and_determinism():
    a = data.synthetic_dynamic_spectrum(n_time=64, n_chan=32, seed=1)
    b = data.synthetic_dynamic_spectrum(n_time=64, n_chan=32, seed=1)
    assert a.shape == (64, 32)
    assert np.array_equal(a, b)  # deterministic for a fixed seed
    # the injected pulse row is brighter than the median noise row
    assert a[32].mean() > a[0].mean()


def test_cli_list_runs(capsys, tmp_path, monkeypatch):
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path))
    assert data._main(["--list"]) == 0
    out = capsys.readouterr().out
    assert "Cache directory" in out
    assert data.list_datasets()[0] in out
