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


# Common LaTeX macros -> plain-text/unicode, so the auto-extracted abstract stays readable
# (arXiv stores the abstract as near-plain text). Applied before the generic macro strip.
_SYMBOLS = {
    r"\approx": "~",
    r"\simeq": "~",
    r"\sim": "~",
    r"\pm": "+/-",
    r"\mp": "-/+",
    r"\times": "x",
    r"\cdot": ".",
    r"\propto": "~prop~",
    r"\equiv": "=",
    r"\leq": "<=",
    r"\geq": ">=",
    r"\ll": "<<",
    r"\gg": ">>",
    r"\lesssim": "<~",
    r"\gtrsim": ">~",
    r"\ell": "l",
    r"\odot": "sun",
    r"\alpha": "alpha",
    r"\beta": "beta",
    r"\gamma": "gamma",
    r"\delta": "delta",
    r"\mu": "mu",
    r"\nu": "nu",
    r"\sigma": "sigma",
    r"\chi": "chi",
    r"\phi": "phi",
    r"\arcdeg": " deg",
    r"\degree": " deg",
    r"\arcsec": "arcsec",
    r"\arcmin": "arcmin",
    r"\sin": "sin ",
    r"\cos": "cos ",
    r"\tan": "tan ",
    r"\log": "log ",
    r"\ln": "ln ",
    r"\bmod": " mod ",
    r"\textbackslash": "\\",
}


def _latex_to_text(s: str) -> str:
    s = re.sub(
        r"\\(emph|textit|textbf|code|texttt|mathrm|mathit|text|textsc)\{([^{}]*)\}", r"\2", s
    )
    s = re.sub(r"\\cite[tp]?\*?(\[[^\]]*\])*\{[^}]*\}", "", s)
    s = re.sub(r"\\(citealt|citeauthor|ref|label)\*?\{[^}]*\}", "", s)
    # sub/superscripts: keep the content inline so e.g. R_0 -> R0, ^{-1/2} -> ^(-1/2)
    s = re.sub(r"_\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\^\{([^{}]*)\}", r"^(\1)", s)
    s = re.sub(r"_([A-Za-z0-9])", r"\1", s)
    for macro, rep in _SYMBOLS.items():
        s = s.replace(macro, rep)
    s = re.sub(r"\\[,;:!> ]", " ", s)  # thin/medium spaces and \>
    s = re.sub(r"\$([^$]*)\$", r"\1", s)  # drop math delimiters, keep contents
    s = re.sub(r"\\[a-zA-Z]+", "", s)  # any remaining macros
    s = s.replace("~", " ").replace("{", "").replace("}", "")
    for esc in ("&", "%", "_", "#", "$"):  # unescape \& \% \_ \# \$
        s = s.replace("\\" + esc, esc)
    s = re.sub(r"\s+([,.;:])", r"\1", s)  # drop spaces left before punctuation by \cite removal
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
        warn.append(
            "no .bbl found — compile with bibtex and ship the .bbl (arXiv may not run your .bst)"
        )
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


# --- arXiv category suggestion ("which groups to submit to") -----------------
# arXiv routes a paper by its *primary* category (one, the best single fit) plus up
# to ~3 *cross-lists* (secondary categories whose readers should also see it). For
# radio astronomy almost everything is under astro-ph; the science domain sets the
# primary and astro-ph.IM rides along whenever the contribution is a tool/method
# (the whole point of this repo's papers). Non-astro cross-lists (eess.SP for pure
# DSP, physics.space-ph for heliophysics) are added only when a *different*
# community should genuinely see it.
ASTRO_CATEGORIES = {
    "astro-ph.GA": "galaxies, Milky Way, ISM, AGN/jets, HI 21cm, Faraday rotation",
    "astro-ph.CO": "cosmology, large-scale structure, isotropy/dipole, source counts, lensing",
    "astro-ph.HE": "FRBs, pulsars/magnetars, transients, accretion/jets, explosive phenomena",
    "astro-ph.SR": "the Sun & solar radio bursts, stars, M/white/brown dwarfs, coherent emission",
    "astro-ph.EP": "planets & the solar system: Jovian/Saturnian/ice-giant radio, magnetospheres",
    "astro-ph.IM": "instrumentation, methods, software, pipelines, data analysis, statistics, RFI",
}

# Non-astro cross-lists a radio paper occasionally warrants (a *different* community).
_NONASTRO_CATEGORIES = {
    "eess.SP": "signal processing — for a pure DSP/algorithm contribution",
    "physics.space-ph": "space physics / heliophysics — solar-wind & magnetosphere radio",
    "physics.ins-det": "instrumentation & detectors — for a hardware/receiver design",
    "gr-qc": "gravitational waves / GR",
    "stat.ML": "machine-learning methodology (e.g. SBI, classifiers)",
}

# Curated (primary, cross_lists) for THIS repo's slices, keyed by paper-dir name.
# Rule of thumb encoded here: primary = the science domain the headline result lives
# in; astro-ph.IM cross-listed because each is a reproducible tool. A few are IM-primary
# because the *tool itself* is the contribution (benchmarks, DSP kernels, pipelines).
SLICE_CATEGORIES: dict[str, tuple[str, list[str]]] = {
    # FRBs & transients
    "frbstats": ("astro-ph.HE", ["astro-ph.IM"]),
    "frbperiod": ("astro-ph.HE", ["astro-ph.IM"]),
    "frbwait": ("astro-ph.HE", ["astro-ph.IM"]),
    "frblens": ("astro-ph.HE", ["astro-ph.CO", "astro-ph.IM"]),
    "driftsearch": ("astro-ph.IM", ["astro-ph.HE"]),
    "lpt": ("astro-ph.HE", ["astro-ph.IM"]),
    "lptv": ("astro-ph.HE", ["astro-ph.IM"]),
    # Pulsars
    "pulsarspec": ("astro-ph.HE", ["astro-ph.IM"]),
    "ppdot": ("astro-ph.HE", ["astro-ph.IM"]),
    "pte2": ("astro-ph.HE", ["astro-ph.IM"]),
    "glitchpop": ("astro-ph.HE", ["astro-ph.IM"]),
    "wdpulsar": ("astro-ph.SR", ["astro-ph.HE", "astro-ph.IM"]),
    # Galaxies / AGN / ISM
    "hi": ("astro-ph.GA", ["astro-ph.IM"]),
    "peaked": ("astro-ph.GA", ["astro-ph.IM"]),
    "southern": ("astro-ph.GA", ["astro-ph.IM"]),
    "offsets": ("astro-ph.GA", ["astro-ph.IM"]),
    "vlbi": ("astro-ph.GA", ["astro-ph.HE", "astro-ph.IM"]),
    "vlass": ("astro-ph.GA", ["astro-ph.HE", "astro-ph.IM"]),
    "stacking": ("astro-ph.GA", ["astro-ph.CO", "astro-ph.IM"]),
    "spectra": ("astro-ph.GA", ["astro-ph.IM"]),
    "fashienv": ("astro-ph.GA", ["astro-ph.CO", "astro-ph.IM"]),
    # Faraday rotation / cosmic magnetism
    "rmsky": ("astro-ph.GA", ["astro-ph.IM"]),
    "rmstructure": ("astro-ph.GA", ["astro-ph.IM"]),
    "rmdipole": ("astro-ph.CO", ["astro-ph.GA", "astro-ph.IM"]),
    # Source counts
    "sourcecounts": ("astro-ph.CO", ["astro-ph.GA", "astro-ph.IM"]),
    # Solar radio
    "solarbursts": ("astro-ph.SR", ["astro-ph.IM"]),
    "windwaves": ("astro-ph.SR", ["physics.space-ph", "astro-ph.IM"]),
    "swaves": ("astro-ph.SR", ["physics.space-ph", "astro-ph.IM"]),
    "triangulate": ("astro-ph.SR", ["astro-ph.IM"]),
    "type3synthesis": ("astro-ph.SR", ["astro-ph.IM"]),
    "typeii": ("astro-ph.SR", ["astro-ph.IM"]),
    "ecallisto_pipeline": ("astro-ph.IM", ["astro-ph.SR"]),
    "ecallisto_census": ("astro-ph.SR", ["astro-ph.IM"]),
    # Stellar coherent emitters (RACS Stokes V)
    "stokesv": ("astro-ph.SR", ["astro-ph.IM"]),
    "stokesv_discovery": ("astro-ph.SR", ["astro-ph.IM"]),
    "svsbi": ("astro-ph.SR", ["astro-ph.IM"]),
    # Planetary radio
    "junodam": ("astro-ph.EP", ["astro-ph.IM"]),
    "skr": ("astro-ph.EP", ["astro-ph.IM"]),
    "vgpra": ("astro-ph.EP", ["astro-ph.IM"]),
    # Instrumentation / DSP
    "torchfdmt": ("astro-ph.IM", ["astro-ph.HE", "eess.SP"]),
    "torchdsp": ("astro-ph.IM", ["astro-ph.HE", "eess.SP"]),
    "rfitrend": ("astro-ph.IM", ["astro-ph.SR"]),
}

# Keyword votes for papers not in the curated table (e.g. a new station paper).
_CATEGORY_KEYWORDS: dict[str, tuple[str, ...]] = {
    "astro-ph.HE": (
        "fast radio burst",
        "frb",
        "pulsar",
        "magnetar",
        "giant pulse",
        "glitch",
        "transient",
        "dispersion measure",
        "dedispersion",
        "repeater",
        "gamma-ray",
        "accretion",
        "blazar",
    ),
    "astro-ph.GA": (
        "galax",
        "milky way",
        "galactic",
        "interstellar",
        " ism ",
        "hi 21",
        "21 cm",
        "21cm",
        "rotation curve",
        "faraday",
        "rotation measure",
        "agn",
        "quasar",
        "radio galaxy",
        "peaked-spectrum",
        "gps",
        "css",
        "h ii",
    ),
    "astro-ph.CO": (
        "cosmolog",
        "large-scale structure",
        "isotropy",
        "dipole",
        "source counts",
        "number counts",
        "gravitational lens",
        "cmb",
    ),
    "astro-ph.SR": (
        "solar",
        " sun ",
        "coronal",
        "type iii",
        "type ii",
        "stellar",
        "m dwarf",
        "white dwarf",
        "brown dwarf",
        "heliospher",
        "flare",
        "coherent emission",
        "chromospher",
    ),
    # NB: no bare "io " — it matches inside "rad-io-" on every radio paper.
    "astro-ph.EP": (
        "planet",
        "jupiter",
        "jovian",
        "saturn",
        "decametric",
        "kilometric",
        "uranus",
        "neptune",
        "magnetospher",
    ),
    "astro-ph.IM": (
        "software",
        "pipeline",
        "reproducible",
        "benchmark",
        "gpu",
        "pytorch",
        "algorithm",
        "detector",
        "rfi",
        "data analysis",
        "statistical",
        "library",
        "framework",
        "tool",
        "method",
    ),
}


def suggest_categories(slice_name: str, title: str, abstract: str) -> tuple[str, list[str], str]:
    """Suggest (primary_category, cross_lists, rationale) for a paper.

    Curated per-slice picks win; otherwise a keyword vote over the title+abstract
    chooses the primary science domain and adds astro-ph.IM whenever the paper reads
    as a tool/method. Always a *suggestion* the human confirms — arXiv moderators can
    reclassify, and the primary is a judgment call at the margins.
    """
    if slice_name in SLICE_CATEGORIES:
        primary, cross = SLICE_CATEGORIES[slice_name]
        return primary, list(cross), f"curated for the '{slice_name}' slice"
    text = f" {title} {abstract} ".lower()
    votes = {c: sum(text.count(k) for k in kws) for c, kws in _CATEGORY_KEYWORDS.items()}
    science = {c: n for c, n in votes.items() if c != "astro-ph.IM" and n > 0}
    is_tool = votes["astro-ph.IM"] >= 2
    if science:
        # Primary = the strongest science domain; then the next domain(s); then
        # astro-ph.IM if the paper reads as a tool. (A science result stays
        # domain-primary with IM as a cross-list — the conventional choice; the
        # curated table above encodes the few genuinely IM-primary papers.)
        ranked = sorted(science, key=lambda c: science[c], reverse=True)
        primary, cross = ranked[0], ranked[1:3]
        if is_tool:
            cross = cross + ["astro-ph.IM"]
    else:
        primary, cross = "astro-ph.IM", []  # pure methods/software, no science domain hit
    return primary, cross[:3], "keyword-inferred — VERIFY against the guide in SKILL.md"


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

    # Auto-fill what we can read from the source: the ORCID in \author[ORCID]{...}, and a
    # comments line seeded from the figure count (page count still needs the human / the PDF).
    om = re.search(r"\\author\[(\d{4}-\d{4}-\d{4}-[\dXx]{4})\]", text)
    orcid = om.group(1) if om else "TODO"
    n_fig = sum(1 for f in files if f.suffix == ".pdf")
    comments = f"TODO pages, {n_fig} figures. Code and data: github.com/joebarbere/jansky-research"

    # Suggest which arXiv groups this paper belongs in (primary + cross-lists).
    primary_cat, cross_cats, cat_why = suggest_categories(paper.name, title, abstract)
    _scopes = {**ASTRO_CATEGORIES, **_NONASTRO_CATEGORIES}
    cat_ref = "\n".join(f"#   {c:<17}{_scopes.get(c, '')}" for c in [primary_cat, *cross_cats])
    cross_yaml = ", ".join(cross_cats)

    tar = out / "arxiv-source.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        for f in files:
            tf.add(f, arcname=str(f.relative_to(paper)))

    (out / "metadata.yaml").write_text(
        f"""# arXiv submission metadata — fill every TODO, then upload at https://arxiv.org/submit
title: {title!r}
authors:        # full names + affiliations, in order (an AI/LLM is NOT an eligible author)
{chr(10).join("  - " + repr(x) for x in authors) or "  - TODO"}
abstract: |
  {abstract}
# categories — SUGGESTED ({cat_why}). Confirm the primary is the single best fit and
# each cross-list genuinely targets a *different* readership (arXiv allows up to ~3):
{cat_ref}
primary_category: {primary_cat}
cross_lists: [{cross_yaml}]
comments: {comments!r}
license: CC BY 4.0       # or CC BY-SA 4.0 | CC BY-NC-SA 4.0 | CC0 | arXiv non-exclusive
report_number: null
journal_ref: null        # add after journal acceptance
doi: null
msc_class: null          # math.* only
acm_class: null          # cs.* only
orcid: {orcid}              # set on the arXiv account/profile too
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
        "3. Paste title/authors/abstract/comments; set the primary + cross-list categories from",
        f"   metadata.yaml (SUGGESTED: {primary_cat}"
        + (f" + [{cross_yaml}]" if cross_yaml else "")
        + ") — verify per 'Choosing the categories' in SKILL.md.",
        "4. Preview & submit (new submitters may need endorsement in the category).",
    ]
    (out / "CHECKLIST.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {tar}, {out / 'metadata.yaml'}, {out / 'CHECKLIST.md'}")
    print(f"  {len(files)} files | {len(errs)} errors | {len(warns)} warnings")
    return 1 if errs else 0


if __name__ == "__main__":
    raise SystemExit(main())
