"""jansky-research — amateur radio-astronomy research tooling.

A sibling of the `jansky` teaching course. Where `jansky` *teaches* radio
astronomy, this package *does* a small slice of original amateur research:
new tooling that fills a surveyed open-source gap, an automated analysis of a
public dataset, and a reproducible write-up.

The package reuses `jansky`'s tested helpers (``jansky.transients``,
``jansky.rfi``, ``jansky.timing``, ``jansky.seti``, ``jansky.sourcecounts``,
``jansky.formats`` ...) rather than reimplementing them. The science focus and
the gap-filling tool module are chosen by the survey (see ``plans/``) and land
here after GATE 1.

Layers (mirroring jansky's tested-helper style):

* ``data``     — research-dataset registry + offline synthetic fallback.
* ``pipeline`` — fetch -> analyze -> metrics; the single entry point shared by
  the Makefile, the notebooks, and the Airflow DAG. (Added in P3.)
* ``report``   — figure/table emitters that write the paper's inputs. (Added in P3.)
"""

from __future__ import annotations

__version__ = "0.0.1"

__all__ = ["__version__"]
