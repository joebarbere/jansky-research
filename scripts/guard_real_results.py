"""Refuse to package papers when pipeline outputs are synthetic.

The offline figure DAG (``make figures``) regenerates every slice's metrics/figures/macros from
synthetic fixtures — correct for CI smoke builds, but silently substituting those outputs into
compiled papers mislabels synthetic validation numbers as real-data results (caught 2026-07-31:
the released ``hi`` paper carried the synthetic 231 km/s flat level captioned as LAB data).

This guard scans ``results/*.json`` for a ``source`` field containing ``synthetic`` and fails,
listing the offenders, unless the file is allowlisted in ``scripts/guard_real_allowlist.txt``
(one filename per line, ``#`` comments; an entry is only added after verifying the corresponding
paper's own text states the synthetic/offline provenance). ``make papers-zip`` runs this before
compiling anything.

Exit codes: 0 = all real (or allowlisted); 1 = synthetic-sourced results present.
"""

from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path("results")
ALLOWLIST = Path("scripts/guard_real_allowlist.txt")


def _allowlist() -> set[str]:
    if not ALLOWLIST.exists():
        return set()
    lines = ALLOWLIST.read_text().splitlines()
    return {ln.strip() for ln in lines if ln.strip() and not ln.lstrip().startswith("#")}


def _source_of(path: Path) -> str | None:
    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return None
    if not isinstance(data, dict):
        return None
    src = data.get("source")
    return src if isinstance(src, str) else None


def main() -> int:
    allow = _allowlist()
    offenders = []
    checked = 0
    for path in sorted(RESULTS.glob("*.json")):
        src = _source_of(path)
        if src is None:
            continue
        checked += 1
        if "synthetic" in src.lower() and path.name not in allow:
            offenders.append((path.name, src))
    if offenders:
        print("GUARD FAILED: synthetic-sourced results present — papers built from these would")
        print("mislabel synthetic validation output as real data. Re-run the real legs")
        print("(`make reproduce` / scripts/*_real.py) or allowlist with justification:")
        for name, src in offenders:
            print(f"  results/{name}  (source: {src!r})")
        return 1
    print(
        f"guard-real: OK — {checked} sourced results checked, none synthetic (allowlist: {len(allow)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
