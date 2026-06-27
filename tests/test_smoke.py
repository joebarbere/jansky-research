"""Smoke tests: the package imports and the jansky dependency resolves."""

from __future__ import annotations

import jansky_research


def test_package_version():
    assert isinstance(jansky_research.__version__, str)
    assert jansky_research.__version__


def test_jansky_dependency_importable():
    # The whole point of the sibling repo: jansky's tested helpers are available.
    import jansky
    from jansky import transients  # a representative helper module

    assert hasattr(transients, "dispersion_delay")
    assert jansky  # imported successfully
