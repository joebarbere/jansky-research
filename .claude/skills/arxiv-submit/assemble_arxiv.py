#!/usr/bin/env python3
"""Assemble + validate an arXiv submission package (no upload — arXiv has no submit API).

Given a paper directory, finds the main LaTeX file, collects its dependencies (\\input files,
figures, the compiled .bbl, generated macros), builds an arxiv-source.tar.gz, extracts the
title/authors/abstract, writes a metadata.yaml with every arXiv submission property, validates
against arXiv's rules, and writes a CHECKLIST.md. Stdlib only.

    python assemble_arxiv.py --paper paper --out arxiv-submission
"""

from __future__ import annotations

import argparse
import re
import tarfile
from pathlib import Path

ABSTRACT_MAX = 1920  # arXiv abstract character limit


def _braced(text: str, cmd: str) -> str | None:
    """Return the brace-matched argument of ``\\cmd{...}`` (handles nesting / newlines)."""
    m = re.search(rf"\\{cmd}\s*(\[[^\]]*\])?\s*{{", text)
    if not m:
        return None
    i = m.end()
    depth, out = 1, []
    while i < len(text) and depth:
        c = text[i]
        depth += c == "{"
        depth -= c == "}"
        if depth:
            out.append(c)
        i += 1
    return "".join(out).strip()


def _load_macros(paper: Path, main: Path) -> dict[str, str]:
    """Collect no-argument ``\\newcommand{\\name}{value}`` definitions from the source + \\inputs."""
    text = main.read_text(errors="ignore")
    for inp in re.findall(r"\\input\{([^}]+)\}", text):
        f = paper / (inp if inp.endswith(".tex") else inp + ".tex")
        if f.exists():
            text += "\n" + f.read_text(errors="ignore")
    return dict(re.findall(r"\\newcommand\{\\(\w+)\}\{([^{}]*)\}", text))


def _apply_macros(s: str, macros: dict[str, str]) -> str:
    for name, val in sorted(macros.items(), key=lambda kv: -len(kv[0])):
        s = re.sub(rf"\\{name}(?![a-zA-Z])", val, s)
    return s


def _latex_to_text(s: str) -> str:
    s = re.sub(r"\\(emph|textit|textbf|code|texttt)\{([^}]*)\}", r"\2", s)
    s = re.sub(r"\\cite[tp]?\*?(\[[^\]]*\])*\{[^}]*\}", "", s)
    s = re.sub(r"\$([^$]*)\$", r"\1", s)
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    s = s.replace("~", " ").replace("\\&", "&")
    return re.sub(r"\s+", " ", s).strip()


def find_main_tex(paper: Path) -> Path:
    cands = [p for p in paper.glob("*.tex") if "\\documentclass" in p.read_text(errors="ignore")]
    if not cands:
        raise SystemExit(f"no .tex with \\documentclass in {paper}")
    return next((p for p in cands if p.name == "main.tex"), cands[0])


def collect(paper: Path, main: Path) -> tuple[list[Path], list[str]]:
    """Return (files-to-ship, warnings)."""
    text = main.read_text(errors="ignore")
    files, warn = [main], []
    for inp in re.findall(r"\\input\{([^}]+)\}", text):
        f = paper / (inp if inp.endswith(".tex") else inp + ".tex")
        (files.append(f) if f.exists() else warn.append(f"missing \\input: {inp}"))
    figs = re.findall(r"\\(?:includegraphics(?:\[[^\]]*\])?|plotone|plottwo)\{([^}]+)\}", text)
    for fig in figs:
        f = paper / fig
        (files.append(f) if f.exists() else warn.append(f"missing figure: {fig}"))
    bbl = main.with_suffix(".bbl")
    if bbl.exists():
        files.append(bbl)
    else:
        warn.append("no .bbl found — compile with bibtex and ship the .bbl (arXiv may not run your .bst)")
    for bib in re.findall(r"\\bibliography\{([^}]+)\}", text):
        f = paper / (bib + ".bib")
        if f.exists():
            files.append(f)
    return list(dict.fromkeys(files)), warn


def validate(paper: Path, main: Path, abstract: str, files: list[Path]) -> list[str]:
    errs = []
    if len(abstract) > ABSTRACT_MAX:
        errs.append(f"abstract is {len(abstract)} chars (> {ABSTRACT_MAX})")
    text = main.read_text(errors="ignore")
    for inp in re.findall(r"\\(?:input|includegraphics(?:\[[^\]]*\])?|plotone)\{([^}]+)\}", text):
        if inp.startswith("/") or inp.startswith(".."):
            errs.append(f"non-relative path (won't survive on arXiv): {inp}")
    if not any(f.suffix == ".bbl" for f in files):
        errs.append("no .bbl in the package — bibliography will not render on arXiv")
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", default="paper")
    ap.add_argument("--out", default="arxiv-submission")
    a = ap.parse_args()
    paper, out = Path(a.paper), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    main_tex = find_main_tex(paper)
    text = main_tex.read_text(errors="ignore")
    macros = _load_macros(paper, main_tex)
    title = _latex_to_text(_apply_macros(_braced(text, "title") or "TODO: title", macros))
    authors = [_latex_to_text(a) for a in re.findall(r"\\author(?:\[[^\]]*\])?\{([^}]*)\}", text)]
    abm = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", text, re.S)
    abstract = _latex_to_text(_apply_macros(abm.group(1), macros)) if abm else "TODO: abstract"

    files, warns = collect(paper, main_tex)
    errs = validate(paper, main_tex, abstract, files)

    tar = out / "arxiv-source.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        for f in files:
            tf.add(f, arcname=str(f.relative_to(paper)))

    (out / "metadata.yaml").write_text(
        f"""# arXiv submission metadata — fill every TODO, then upload at https://arxiv.org/submit
title: {title!r}
authors:        # full names + affiliations, in order (an AI/LLM is NOT an eligible author)
{chr(10).join('  - ' + repr(x) for x in authors) or '  - TODO'}
abstract: |
  {abstract}
primary_category: TODO   # e.g. astro-ph.IM | astro-ph.GA | astro-ph.HE | astro-ph.CO
cross_lists: []          # 0-3 secondary, e.g. [astro-ph.GA, astro-ph.HE]
comments: TODO           # e.g. "4 pages, 3 figures. Code: github.com/joebarbere/jansky-research"
license: CC BY 4.0       # or CC BY-SA 4.0 | CC BY-NC-SA 4.0 | CC0 | arXiv non-exclusive
report_number: null
journal_ref: null        # add after journal acceptance
doi: null
msc_class: null          # math.* only
acm_class: null          # cs.* only
orcid: TODO              # e.g. 0009-0008-3289-4447 (set on the arXiv account too)
"""
    )

    lines = [
        "# arXiv submission checklist",
        "",
        f"- main file: `{main_tex.name}`",
        f"- package: `{tar}` ({len(files)} files)",
        f"- abstract: {len(abstract)} / {ABSTRACT_MAX} chars",
        "",
        "## Validation",
        *(f"- ❌ {e}" for e in errs),
        *(f"- ⚠️  {w}" for w in warns),
        ("- ✅ no blocking errors" if not errs else "- **fix the ❌ items before submitting**"),
        "",
        "## Upload (manual — arXiv has no submit API)",
        "1. https://arxiv.org/submit → choose license (metadata.yaml).",
        "2. Upload `arxiv-source.tar.gz`; check the AutoTeX preview matches paper/main.pdf.",
        "3. Paste title/authors/abstract/comments; set primary + cross-list categories.",
        "4. Preview & submit (new submitters may need endorsement in the category).",
    ]
    (out / "CHECKLIST.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {tar}, {out/'metadata.yaml'}, {out/'CHECKLIST.md'}")
    print(f"  {len(files)} files | {len(errs)} errors | {len(warns)} warnings")
    return 1 if errs else 0


if __name__ == "__main__":
    raise SystemExit(main())
