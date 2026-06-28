#!/usr/bin/env python3
"""Check a research repo's readiness to publish to Zenodo, JOSS, RNAAS, and arXiv — in order.

A static pre-flight that inspects the repo's *artifacts* (no network, no accounts) and prints an
ordered readiness report mirroring the dependency chain (Zenodo DOI unblocks JOSS; RNAAS and arXiv are
independent). It flags what's present, what's missing, and the next manual step — so nothing in the
account-bound `TODO.md` flow is missed.

Marks: [x] ready · [ ] missing/action-needed · [~] present but verify. Exit 0 always (it's a report).

Part of the `research-publish` skill; see SKILL.md.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ARXIV_PAPERS = ("frbstats", "vlass", "peaked", "southern")  # the fresh-angle papers (see TODO.md)


def _exists(p: Path) -> bool:
    return p.exists() and (p.is_dir() or p.stat().st_size > 0)


def check(repo: Path) -> None:
    def mark(ok, msg, note=""):
        m = "[x]" if ok is True else ("[~]" if ok == "verify" else "[ ]")
        print(f"  {m} {msg}" + (f"  — {note}" if note else ""))

    print("\n=== 0. One-time setup ===")
    cff = repo / "CITATION.cff"
    cff_txt = cff.read_text() if cff.exists() else ""
    mark(_exists(repo / "LICENSE"), "OSI license (LICENSE)")
    mark(_exists(cff), "CITATION.cff present")
    mark(
        "0009-0008-3289-4447" in cff_txt or "orcid" in cff_txt.lower(),
        "ORCID in CITATION.cff",
        "keep author name consistent across services",
    )

    print("\n=== 1. Zenodo (do first — JOSS needs its DOI) ===")
    zj = repo / ".zenodo.json"
    zok = False
    if zj.exists():
        try:
            z = json.loads(zj.read_text())
            zok = bool(z.get("title") and z.get("creators"))
            mark(
                zok,
                ".zenodo.json valid (title + creators)",
                f"creators: {', '.join(c.get('name', '?') for c in z.get('creators', []))[:60]}",
            )
        except json.JSONDecodeError:
            mark(False, ".zenodo.json present but INVALID JSON")
    else:
        mark(False, ".zenodo.json present")
    readme = (repo / "README.md").read_text() if (repo / "README.md").exists() else ""
    has_badge = bool(re.search(r"zenodo\.org/badge|doi\.org/10\.5281", readme))
    mark(has_badge, "Zenodo DOI badge in README", "added after the first release is archived")
    print("      next: enable the repo on zenodo.org, cut a GitHub release, copy the concept DOI.")

    print("\n=== 2. JOSS (needs the step-1 DOI) ===")
    pm = repo / "joss" / "paper.md"
    pm_txt = pm.read_text() if pm.exists() else ""
    mark(_exists(pm), "joss/paper.md present")
    mark(_exists(repo / "joss" / "paper.bib"), "joss/paper.bib present")
    mark(bool(re.search(r"orcid", pm_txt, re.I)), "paper.md has an ORCID'd author")
    mark(
        "verify" if pm.exists() else False,
        "JOSS 'substantial scholarly effort' bar",
        "lead with the validated analyses + Airflow/Podman layer",
    )
    print("      next: archive on Zenodo (step 1), then submit at joss.theoj.org/papers/new.")

    print("\n=== 3. RNAAS (independent; quickest) ===")
    rn = repo / "papers" / "frbstats" / "rnaas.tex"
    mark(_exists(rn), "papers/frbstats/rnaas.tex present")
    if rn.exists():
        words = len(re.sub(r"[^A-Za-z0-9 ]", " ", rn.read_text()).split())
        mark(
            "verify", f"length sanity (~{words} raw words incl. LaTeX)", "RNAAS limit is 1000 words"
        )
    print("      next: build the note, confirm the current AAS fee, submit via Editorial Manager.")

    print("\n=== 4. arXiv (fresh-angle papers; endorsement has lead time) ===")
    for p in ARXIV_PAPERS:
        d = repo / "papers" / p
        main_ok = _exists(d / "main.tex")
        pkg_ok = _exists(d / "arxiv-submission" / "arxiv-source.tar.gz")
        meta_ok = _exists(d / "arxiv-submission" / "metadata.yaml")
        status = (
            "ready (run `make arxiv`)" if main_ok and not pkg_ok else ("packaged" if pkg_ok else "")
        )
        mark(main_ok, f"papers/{p}/main.tex", status)
        if pkg_ok:
            mark(
                "verify",
                f"  papers/{p}/arxiv-submission/ package",
                "fill metadata.yaml TODOs (primary_category, comments)"
                if meta_ok
                else "no metadata.yaml",
            )
    print(
        "      next: register + line up an endorser early; `make arxiv`; fill metadata.yaml; upload."
    )
    print("\nSee TODO.md for the full account-bound walkthrough.\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--repo", default=".", help="repo root to inspect (default: cwd)")
    args = ap.parse_args(argv)
    repo = Path(args.repo).resolve()
    print(f"Publishing readiness for {repo.name} (artifacts only; accounts/DOIs are manual):")
    check(repo)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
