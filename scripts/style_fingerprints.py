#!/usr/bin/env python3
"""Stage-3a fingerprint runs for the pre-LLM style corpus (plan 89).

Two runs, two committed artifacts:

  --selfscan   fingerprint every papers/<slice>/main.tex (the repo's own prose)
               -> results/stylecorpus_selfscan.json
  --corpus     fingerprint the downloaded Stage-2 sample (LaTeX sources; scanned
               PDFs are handled qualitatively in Stage 3b, not here)
               -> results/stylecorpus_fingerprints.json

Both are offline over local files; --corpus requires the Stage-2 download to have
run (data/style_corpus/arxiv_src/ + selection.json).
"""

from __future__ import annotations

import argparse
import json
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import stylecorpus as sc  # noqa: E402


def _write(out: Path, payload: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(out)
    print(f"wrote {out}")


def run_selfscan(out: Path) -> None:
    per_paper: dict[str, dict[str, float]] = {}
    for main in sorted(Path("papers").glob("*/main.tex")):
        per_paper[main.parent.name] = sc.fingerprint_latex(main.read_text())
    payload = {
        "slice": "stylecorpus",
        "stage": "selfscan",
        "is_real": True,
        "n_papers": len(per_paper),
        "per_paper": per_paper,
        "aggregate": sc.aggregate_fingerprints(list(per_paper.values())),
    }
    _write(out, payload)


def _eprint_tex(path: Path) -> str | None:
    """Concatenated .tex members of an e-print bundle (tar/gz/plain), or None."""
    import gzip

    try:
        if tarfile.is_tarfile(path):
            texts = []
            with tarfile.open(path) as tf:
                for member in tf.getmembers():
                    if member.name.endswith(".tex") and member.size < 5_000_000:
                        fh = tf.extractfile(member)
                        if fh:
                            texts.append(fh.read().decode("utf-8", errors="replace"))
            return "\n".join(texts) or None
        raw = path.read_bytes()
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        if raw[:5] == b"%PDF-":
            return None  # PDF-only submission
        text = raw.decode("utf-8", errors="replace")
        return text if "\\" in text else None
    except Exception as exc:  # noqa: BLE001 - a corrupt bundle is a lost sample point
        print(f"  unreadable {path.name}: {exc}")
        return None


def run_corpus(out: Path) -> None:
    root = sc.corpus_dir()
    selection = json.loads((root / "selection.json").read_text())
    by_id = {s["arxiv_id"]: s for s in selection if s["arxiv_id"]}
    per_era_docs: dict[str, list[dict[str, float]]] = {}
    n_read = n_skipped = 0
    for path in sorted((root / "arxiv_src").glob("*.eprint")):
        arxiv_id = path.stem.replace("_", "/")
        item = by_id.get(arxiv_id)
        if item is None:
            continue
        tex = _eprint_tex(path)
        if tex is None or len(tex.split()) < 500:
            n_skipped += 1
            continue
        per_era_docs.setdefault(item["era"], []).append(sc.fingerprint_latex(tex))
        n_read += 1
    payload = {
        "slice": "stylecorpus",
        "stage": "fingerprints",
        "is_real": True,
        "n_documents": n_read,
        "n_unreadable_or_short": n_skipped,
        "source": "arxiv LaTeX sources from the Stage-2 stratified sample",
        "per_era": {
            era: sc.aggregate_fingerprints(docs) for era, docs in sorted(per_era_docs.items())
        },
        "per_era_n": {era: len(docs) for era, docs in sorted(per_era_docs.items())},
        "all": sc.aggregate_fingerprints([d for v in per_era_docs.values() for d in v]),
    }
    _write(out, payload)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selfscan", action="store_true")
    ap.add_argument("--corpus", action="store_true")
    args = ap.parse_args(argv)
    if args.selfscan:
        run_selfscan(Path("results/stylecorpus_selfscan.json"))
    if args.corpus:
        run_corpus(Path("results/stylecorpus_fingerprints.json"))
    if not (args.selfscan or args.corpus):
        ap.error("pass --selfscan and/or --corpus")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
