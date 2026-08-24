#!/usr/bin/env python
"""Mechanical triage across every paper, for the defects the 2026-08 review round kept finding.

Deep refereeing costs about a session per paper. This is the cheap first pass: it looks only
for failure modes that can be detected without understanding the science, ranks the papers by
what it finds, and says which deserve the expensive treatment.

Every check here corresponds to a real defect found in a real paper this month:

  cited-placeholder   a macro rendering as "--" that the prose actually cites, so the sentence
                      reads as a finished claim with the number gone (four papers shipped this)
  unnamespaced        a mode-dependent macro emitted outside the Syn/Real loop, so an offline
                      rebuild writes synthetic values under a real name -- a WRONG number, not
                      a blank, invisible to the "--" guard (\\tiiNEvents, \\svbNTargets)
  hard-typed          a number in the prose that equals a committed macro's value, i.e. a
                      value that will silently drift when the pipeline is re-run
  bad-doi             a DOI whose Crossref title does not match the .bib title (found two
                      pointing at entirely different papers)
  bare-author         a .bib entry with author = {others} or no author, which renders as
                      "others. 2025, arXiv e-prints"
  no-evidence         a paper whose results file is missing, or marked is_real: false
  overclaim           a verb this repo has had to retract before, near a number

Usage:  uv run python scripts/triage_papers.py [--no-network] [--paper NAME]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "papers"
RESULTS = ROOT / "results"
CACHE = Path("/tmp/triage_crossref_cache.json")

# Slices whose evidence predates the <slug>_metrics.json convention. Each mapping is taken
# from the module's own write call, not guessed: frbstats (the oldest slice) writes the bare
# metrics.json; hi writes rotation_curve.json; spectra's slice was named "uss" internally;
# the torchfdmt paper is backed by the fdmt+singlepulse modules. Aliased here rather than
# renaming the committed files, because the modules and their tests read these paths and a
# rename would churn real evidence for a cosmetic gain.
HARD_TYPED_OK = {
    ("pte2", "0.05"),  # conventional p-threshold; equality with \ptSynFP is coincidence
    ("vgpra", "17.24"),  # Voyager-era literature constant; the injection was chosen to match
    ("vgpra", "16.11"),  # Neptune's Voyager-2 period (Warwick et al. 1989), cited as history
    ("dr20radio", "2.5"),  # "the apparent 2.5x" is a flux RATIO (3.0/1.2), not the 2.5" beam
}

EVIDENCE_ALIAS = {
    "frbstats": ["metrics.json"],
    "hi": ["rotation_curve.json"],
    "spectra": ["uss_metrics.json"],
    "driftsearch": ["drift_metrics.json"],
    "frbperiod": ["period_metrics.json"],
    "torchfdmt": ["singlepulse_metrics.json"],
    # Two legs, two evidence files: the paper's numbers are the synthetic day's, and the real
    # archive day is cited once for the coincidence step's real-data behaviour.
    "ecallisto_pipeline": ["ecallisto_metrics.json", "ecallisto_synthetic_metrics.json"],
}

# Verbs this repo has retracted. "recovers" where the honest word was "responds to"; "pins"
# for a parameter sitting on its prior bound; "confirms" for a non-rejection; "unbiased" for a
# Kaplan-Meier under dependent censoring; "complete and pure" for a test that could not fail.
OVERCLAIM = (
    r"\bpins\b",
    r"\bconfirms\b",
    r"\bunbiased\b",
    r"\bproves\b",
    r"\bdemonstrates that\b",
    r"\bestablishes\b",
    r"\brobust to\b",
    r"\bcomplete and pure\b",
    r"\bexactly reproduc",
    r"\bdefinitive\b",
    r"\bnominal coverage\b",
    r"\bfully compatible\b",
)


def _macros(paper: Path) -> dict[str, str]:
    f = paper / "generated" / "macros.tex"
    if not f.is_file():
        return {}
    out = {}
    text = f.read_text()
    # brace-matched: a value can contain braces (\mathrm{}, exponents), and a [^{}]* pattern
    # silently skips those -- one of the regex traps recorded in CLAUDE.md
    for m in re.finditer(r"\\newcommand\{\\([A-Za-z]+)\}\{", text):
        i = m.end()
        depth, buf = 1, []
        while i < len(text) and depth:
            c = text[i]
            depth += (c == "{") - (c == "}")
            if depth:
                buf.append(c)
            i += 1
        out[m.group(1)] = "".join(buf)
    return out


def _used(tex: str) -> set[str]:
    return set(re.findall(r"\\([A-Za-z]+)(?![A-Za-z])", tex))


def paper_tex_files(paper: Path) -> list[Path]:
    """Every compilable document in a paper directory, not just ``main.tex``.

    RNAAS notes live beside the paper as ``rnaas.tex``. Until 2026-08-21 this script read
    only ``main.tex`` and skipped any directory without one, so **no note in the repo had
    ever been triaged** -- including the one at the head of the submission queue. A gate that
    cannot see a document cannot vouch for it.
    """
    return sorted(
        t for t in paper.glob("*.tex") if "\\documentclass" in t.read_text(errors="replace")
    )


def check_paper(
    paper: Path, *, network: bool, cache: dict, tex_path: Path | None = None
) -> list[tuple[str, str, str]]:
    """Return [(severity, kind, detail)] for one document of one paper."""
    out: list[tuple[str, str, str]] = []
    tex_path = tex_path or (paper / "main.tex")
    if not tex_path.is_file():
        return out
    tex = tex_path.read_text()
    macros = _macros(paper)
    used = _used(tex)

    # 1. a cited macro that renders as a placeholder
    for name, val in macros.items():
        if val.strip() in {"--", ""} and name in used:
            out.append(("HIGH", "cited-placeholder", rf"\{name} renders as '--' and is cited"))

    # 2. mode-dependent macros not namespaced: both Syn and Real variants exist for some stem,
    #    but a bare (un-prefixed) name for the same stem is also emitted
    stems: dict[str, set[str]] = {}
    for name in macros:
        m = re.match(r"^([a-z]{2,4})(Syn|Real)(.+)$", name)
        if m:
            stems.setdefault(m.group(1) + m.group(3), set()).add(m.group(2))
    for stem, seen in stems.items():
        pref = re.match(r"^([a-z]{2,4})", stem).group(1)
        rest = stem[len(pref) :]
        bare = pref + rest
        if len(seen) == 2 and bare in macros and bare in used:
            out.append(
                (
                    "HIGH",
                    "unnamespaced",
                    rf"\{bare} is mode-free but \{pref}Syn{rest}/"
                    rf"\{pref}Real{rest} both exist",
                )
            )

    # 3. hard-typed numbers that duplicate a committed macro value.
    # Adjudicated coincidences are allowlisted: equality of a prose number and a macro is not
    # always provenance. pte2's "p<0.05" is the conventional significance threshold, which the
    # synthetic false-positive rate happens to equal -- substituting \ptSynFP there would
    # claim the threshold was measured. vgpra's "17.24 h" is the Voyager-era literature
    # constant quoted in historical context; \vgSynInjected was CHOSEN to equal it.
    prose = re.sub(r"\\newcommand\{[^}]*\}\{[^}]*\}", "", tex)
    for name, val in macros.items():
        v = val.strip()
        if not re.fullmatch(r"-?\d+\.\d+", v) or abs(float(v)) < 0.01:
            continue  # integers and tiny values collide with section numbers etc.
        if (paper.name, v) in HARD_TYPED_OK:
            continue
        if re.search(rf"(?<![\d.]){re.escape(v)}(?![\d.])", prose) and name in used:
            out.append(("MED", "hard-typed", f"{v} appears in prose and as \\{name}"))

    # 4. bibliography
    bib = paper / "refs.bib"
    if bib.is_file():
        btext = bib.read_text()
        for entry in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", btext, re.S):
            key, body = entry.group(1).strip(), entry.group(2)
            am = re.search(r"author\s*=\s*\{(.*?)\}", body, re.S)
            if am and am.group(1).strip().strip("{}") in {"others", ""}:
                out.append(("HIGH", "bare-author", f"{key}: author = {{others}}"))
            dm = re.search(r"doi\s*=\s*\{([^}]+)\}", body)
            tm = re.search(r"title\s*=\s*\{(.*?)\}\s*,\s*\n", body, re.S)
            if network and dm and tm:
                doi = dm.group(1).strip()
                real = crossref_title(doi, cache)
                if real and _title_mismatch(tm.group(1), real):
                    want_s = re.sub(r"\s+", " ", tm.group(1)).strip()
                    out.append(
                        (
                            "HIGH",
                            "bad-doi",
                            f"{key} ({doi})\n        bib:      {want_s[:78]}\n"
                            f"        crossref: {real[:78]}",
                        )
                    )

    # 5. evidence
    slug = paper.name
    cand = list(RESULTS.glob(f"{slug}*.json"))
    if not cand and slug in EVIDENCE_ALIAS:
        cand = [RESULTS / n for n in EVIDENCE_ALIAS[slug] if (RESULTS / n).is_file()]
    if not cand:
        out.append(("MED", "no-evidence", f"no results/{slug}*.json"))
    else:
        parsed = []
        for c in cand:
            try:
                parsed.append((c, json.loads(c.read_text())))
            except Exception:  # noqa: BLE001 - a malformed evidence file is itself a finding
                out.append(("HIGH", "no-evidence", f"{c.name} is not valid JSON"))
        has_real = any(isinstance(d, dict) and d.get("is_real") is True for _, d in parsed)
        for c, d in parsed:
            if isinstance(d, dict) and d.get("is_real") is False:
                # A deliberately-committed synthetic companion (e.g. vgpra writes its offline
                # run to <slice>_synthetic_metrics.json precisely so it can never clobber the
                # real file) is the DESIGN, not missing evidence -- but only when the real
                # sibling actually exists; a synthetic file standing alone is still a finding.
                if "synthetic" in c.name and has_real:
                    continue
                out.append(("HIGH", "no-evidence", f"{c.name} has is_real: false"))

    # 6. retracted verbs, only where a number is nearby
    for pat in OVERCLAIM:
        for m in re.finditer(pat, prose, re.I):
            window = prose[max(0, m.start() - 120) : m.end() + 120]
            if re.search(r"\d", window) or "\\" in window:
                out.append(("LOW", "overclaim", f"'{m.group(0)}'"))
                break
    return out


_STOP = {"the", "of", "a", "an", "and", "in", "for", "with", "from", "on", "to", "at", "its"}


def _norm(s: str) -> str:
    """Title -> comparable token set.

    Crossref returns markup (``<i>Gaia</i>``, ``H<scp>i</scp>``) and journals differ on
    "H I" vs "HI", so a substring test flags nearly every correct DOI. Strip tags and LaTeX
    braces, then compare tokens.
    """
    s = re.sub(r"<[^>]+>", " ", s)
    s = re.sub(r"[{}\\$]", " ", s)
    return re.sub(r"[^a-z0-9 ]", " ", re.sub(r"\s+", " ", s.lower())).strip()


def _title_mismatch(want: str, got: str) -> bool:
    """True when two titles are too different to be the same paper.

    Containment, not equality: a .bib title is often the subtitle-trimmed form of the
    published one. Flag only when the overlap is poor in BOTH directions -- that is what
    distinguishes "same paper, punctuated differently" from "an entirely different paper",
    which is the failure this is for (a DOI pointing at a colliding-wind binary in a
    paper about radio stars).
    """
    a = {t for t in _norm(want).split() if t not in _STOP and len(t) > 2}
    b = {t for t in _norm(got).split() if t not in _STOP and len(t) > 2}
    if not a or not b:
        return False
    overlap = len(a & b)
    return overlap / min(len(a), len(b)) < 0.6


def crossref_title(doi: str, cache: dict) -> str | None:
    if doi in cache:
        return cache[doi]
    try:
        req = Request(
            f"https://api.crossref.org/works/{doi}",
            headers={"User-Agent": "jansky-research triage (mailto:joe.barbere@gmail.com)"},
        )
        with urlopen(req, timeout=30) as fh:
            title = json.load(fh)["message"]["title"][0]
    except Exception:  # noqa: BLE001 - an unresolvable DOI is reported, not raised
        title = None
    cache[doi] = title
    time.sleep(0.15)  # polite
    return title


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-network", action="store_true", help="skip Crossref DOI checks")
    ap.add_argument("--paper", help="triage a single paper")
    a = ap.parse_args(argv)

    cache = json.loads(CACHE.read_text()) if CACHE.is_file() else {}
    papers = [PAPERS / a.paper] if a.paper else sorted(p for p in PAPERS.iterdir() if p.is_dir())

    rows = []
    for p in papers:
        for tex in paper_tex_files(p):
            found = check_paper(p, network=not a.no_network, cache=cache, tex_path=tex)
            label = p.name if tex.name == "main.tex" else f"{p.name}/{tex.stem}"
            rows.append((label, found))
    CACHE.write_text(json.dumps(cache))

    order = {"HIGH": 0, "MED": 1, "LOW": 2}
    rows.sort(key=lambda r: (-sum(3 - order[s] for s, _, _ in r[1]), r[0]))

    print(f"{'paper':24} {'HIGH':>4} {'MED':>4} {'LOW':>4}  findings")
    print("-" * 92)
    for name, found in rows:
        h = sum(1 for s, _, _ in found if s == "HIGH")
        m = sum(1 for s, _, _ in found if s == "MED")
        low = sum(1 for s, _, _ in found if s == "LOW")
        kinds = sorted({k for _, k, _ in found})
        print(f"{name:24} {h:4} {m:4} {low:4}  {', '.join(kinds) if kinds else 'clean'}")
    print()
    for name, found in rows:
        hi = [f for f in found if f[0] == "HIGH"]
        if hi:
            print(f"== {name}")
            for _s, k, d in hi:
                print(f"   [{k}] {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
