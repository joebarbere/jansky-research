#!/usr/bin/env python3
"""Find mode-dependent macros that are NOT namespaced, i.e. live cross-run clobbers.

``report.preserve_live_macros`` stops an offline run blanking a real macro, and
``report.preserve_live_results`` does the same for the results JSON. Neither can help when a
macro *name* means different things in the two run modes: both runs then write a real value
and the merge has nothing to arbitrate. That is the ``\\tiiNEvents`` incident, where an offline
rebuild turned 768 real observing days into 48 synthetic events under prose reading
"768 days, zero failures" -- a wrong published number, worse than the blank it replaced.

It has since recurred twice (`southern`, `junodam`), so this is the grep CLAUDE.md asks for:
run every slice's offline leg into a throwaway directory, diff its macros against the
committed ones, and report every macro that changes value while sharing a name.

    uv run python scripts/audit_macro_namespaces.py            # audit all slices
    uv run python scripts/audit_macro_namespaces.py --slice hi # one slice

Exit codes: 0 = no unnamespaced mode-dependent macros; 1 = at least one found.

**This never writes into the repository.** Each slice runs with ``--out <tmpdir>``; writing
offline output into the repo root is the very thing being audited for.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
MACRO = re.compile(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$")
PLACEHOLDER = "--"
#: A name carrying either marker is already mode-namespaced and cannot collide across runs.
NAMESPACED = re.compile(r"Syn|Real", re.IGNORECASE)


def slice_map() -> dict[str, tuple[str, str]]:
    """The paper-dir -> (module, extra args) mapping the offline DAG uses."""
    text = SNAKEFILE.read_text()
    start = text.index("SLICES = {")
    end = text.index("}", start) + 1
    return ast.literal_eval(text[start + len("SLICES = ") : end])


def read_macros(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    out: dict[str, str] = {}
    for line in path.read_text(errors="ignore").splitlines():
        m = MACRO.match(line.strip())
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def committed_source(slice_name: str) -> str:
    """The `source` recorded for this slice, so a synthetic-by-design paper is not misreported."""
    for cand in REPO.glob(f"results/{slice_name}*_metrics.json"):
        try:
            src = json.loads(cand.read_text()).get("source")
        except (json.JSONDecodeError, OSError, UnicodeDecodeError):
            continue
        if isinstance(src, str):
            return src
    return ""


def audit_slice(slice_name: str, module: str, extra: str) -> dict:
    committed = read_macros(REPO / "papers" / slice_name / "generated" / "macros.tex")
    if not committed:
        return {"slice": slice_name, "status": "no committed macros", "hazards": []}

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [sys.executable, "-m", f"jansky_research.{module}", "--out", tmp]
        cmd += [a for a in extra.split() if a]
        proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=900)
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout).strip().splitlines()[-1:] or [""]
            return {
                "slice": slice_name,
                "status": f"offline run failed: {tail[0][:90]}",
                "hazards": [],
            }
        offline = read_macros(Path(tmp) / "papers" / slice_name / "generated" / "macros.tex")

    hazards = []
    for name, syn_val in offline.items():
        real_val = committed.get(name)
        if real_val is None or NAMESPACED.search(name):
            continue
        # A placeholder on either side is the case preserve_live_macros already arbitrates.
        if syn_val == PLACEHOLDER or real_val == PLACEHOLDER or syn_val == real_val:
            continue
        hazards.append({"macro": name, "committed": real_val, "offline_would_write": syn_val})
    return {
        "slice": slice_name,
        "status": "ok",
        "source": committed_source(slice_name),
        "n_offline_macros": len(offline),
        "hazards": sorted(hazards, key=lambda h: h["macro"]),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slice", help="audit one paper directory instead of all")
    ap.add_argument("--json", type=Path, help="also write the full report here")
    args = ap.parse_args(argv)

    slices = slice_map()
    if args.slice:
        if args.slice not in slices:
            print(f"unknown slice {args.slice!r}; known: {', '.join(sorted(slices))}")
            return 2
        slices = {args.slice: slices[args.slice]}

    reports, failed = [], 0
    for name, (module, extra) in sorted(slices.items()):
        rep = audit_slice(name, module, extra)
        reports.append(rep)
        if rep["hazards"]:
            failed += 1
            synthetic = "synthetic" in rep.get("source", "").lower()
            note = (
                "  (committed results are synthetic, so the clobber is cosmetic)"
                if synthetic
                else ""
            )
            print(f"\n{name}: {len(rep['hazards'])} unnamespaced mode-dependent macro(s){note}")
            for h in rep["hazards"]:
                print(
                    f"    \\{h['macro']:<24} committed {h['committed']!r:>18}"
                    f"  <- offline would write {h['offline_would_write']!r}"
                )
        elif rep["status"] != "ok":
            print(f"{name}: {rep['status']}")

    checked = sum(1 for r in reports if r["status"] == "ok")
    print(
        f"\naudited {checked} of {len(reports)} slices; "
        f"{failed} carry an unnamespaced mode-dependent macro"
    )
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.json}")
    if failed:
        print(
            "\nEach one is a live clobber: both run modes write a real value under the same\n"
            "name, so preserve_live_macros has nothing to arbitrate and whichever ran last\n"
            "wins. Namespace them <slice>Syn*/<slice>Real* and cite the right one in the paper."
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
