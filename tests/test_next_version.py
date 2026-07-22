"""Tests for scripts/next_version.py — the changelog parser + bump recipe. No network, no git."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "next_version", Path(__file__).resolve().parent.parent / "scripts" / "next_version.py"
)
assert _SPEC and _SPEC.loader
nv = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(nv)


CHANGELOG = """\
# Changelog

## [Unreleased]

### Added
- a new `foo` slice
- a `--bar` flag

### Fixed
- a typo in the paper

## [1.2.3] - 2026-01-01

### Added
- the original release
"""


def test_parse_unreleased_collects_only_unreleased():
    sections = nv.parse_unreleased(CHANGELOG)
    assert sections["Added"] == ["a new `foo` slice", "a `--bar` flag"]
    assert sections["Fixed"] == ["a typo in the paper"]
    # The [1.2.3] block must not leak into Unreleased.
    assert "the original release" not in sections["Added"]


def test_parse_unreleased_ignores_subbullets():
    text = "## [Unreleased]\n### Added\n- top level\n  - nested detail\n"
    assert nv.parse_unreleased(text) == {"Added": ["top level"]}


def test_parse_unreleased_empty_when_no_section():
    assert nv.parse_unreleased("# Changelog\n\nnothing here\n") == {}


def test_recommend_initial_release_when_no_tag():
    version, reason = nv.recommend(None, {"Added": ["anything"]})
    assert version == "1.0.0"
    assert "initial release" in reason


def test_recommend_minor_for_added():
    version, reason = nv.recommend("v1.0.0", {"Added": ["a new slice"]})
    assert version == "1.1.0"
    assert reason.startswith("MINOR")


def test_recommend_patch_for_fixed_only():
    version, reason = nv.recommend("v1.4.2", {"Fixed": ["a bug"]})
    assert version == "1.4.3"
    assert reason.startswith("PATCH")


def test_recommend_major_for_removed():
    version, reason = nv.recommend("v1.4.2", {"Removed": ["dropped the old CLI"]})
    assert version == "2.0.0"
    assert reason.startswith("MAJOR")


def test_recommend_major_for_breaking_changed():
    sections = {"Changed": ["**BREAKING** renamed run() argument"], "Added": ["something"]}
    version, reason = nv.recommend("v1.4.2", sections)
    assert version == "2.0.0"
    assert reason.startswith("MAJOR")


def test_recommend_added_beats_fixed():
    version, _ = nv.recommend("v1.0.0", {"Added": ["x"], "Fixed": ["y"]})
    assert version == "1.1.0"


def test_recommend_raises_on_empty():
    with pytest.raises(ValueError):
        nv.recommend("v1.0.0", {})
