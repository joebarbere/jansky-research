#!/usr/bin/env python3
"""Prose lint for the traditional-style skill (plan 89).

Two modes over a paper directory (papers/<slice>):

  default        score main.tex's AI-marker metrics against the pre-LLM corpus
                 percentiles in results/stylecorpus_fingerprints.json; findings
                 in the triage table style (HIGH/MED/LOW). Exit 1 on any HIGH.
  --diff-guard   compare the worktree main.tex against HEAD and FAIL unless the
                 edit was prose-only: the multisets of generated-macro
                 invocations, \\cite keys, numeric literals, and \\software{}
                 content must be unchanged. Exit 1 on any difference.

Usage:
  uv run python .claude/skills/traditional-style/prose_lint.py papers/<slice>
  uv run python .claude/skills/traditional-style/prose_lint.py papers/<slice> --diff-guard
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import stylecorpus as sc  # noqa: E402


def _macro_names(paper: Path) -> set[str]:
    macros = paper / "generated" / "macros.tex"
    return sc.macro_names_from_tex(macros.read_text()) if macros.exists() else set()


def run_lint(paper: Path) -> int:
    corpus = json.loads((REPO / "results" / "stylecorpus_fingerprints.json").read_text())
    fp = sc.fingerprint_latex((paper / "main.tex").read_text())
    findings = sc.lint_paper(fp, corpus["all"])
    if not findings:
        print(f"{paper.name}: no style findings (corpus p90 clean)")
        return 0
    order = {"HIGH": 0, "MED": 1, "LOW": 2}
    for sev, metric, msg in sorted(findings, key=lambda f: order[f[0]]):
        print(f"{sev:5} {metric:28} {msg}")
    n_high = sum(1 for sev, _, _ in findings if sev == "HIGH")
    print(f"{paper.name}: {len(findings)} findings ({n_high} HIGH)")
    return 1 if n_high else 0


def run_diff_guard(paper: Path) -> int:
    rel = (paper / "main.tex").resolve().relative_to(REPO)
    head = subprocess.run(
        ["git", "-C", str(REPO), "show", f"HEAD:{rel.as_posix()}"],
        capture_output=True, text=True, check=True,
    ).stdout
    names = _macro_names(paper)
    problems = sc.compare_signatures(
        sc.guard_signature(head, names),
        sc.guard_signature((paper / "main.tex").read_text(), names),
    )
    if not problems:
        print(f"{paper.name}: diff-guard clean (prose-only edit)")
        return 0
    print(f"{paper.name}: diff-guard FAILED — the edit changed protected content:")
    for p in problems:
        print(f"  {p}")
    return 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paper", type=Path, help="papers/<slice> directory")
    ap.add_argument("--diff-guard", action="store_true")
    args = ap.parse_args(argv)
    if not (args.paper / "main.tex").exists():
        ap.error(f"{args.paper}/main.tex not found")
    return run_diff_guard(args.paper) if args.diff_guard else run_lint(args.paper)


if __name__ == "__main__":
    raise SystemExit(main())
