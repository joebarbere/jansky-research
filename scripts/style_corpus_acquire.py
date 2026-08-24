#!/usr/bin/env python3
"""Stage-2 acquisition for the pre-LLM style corpus (plan 89).

Four idempotent phases, all resumable (each checkpoints to data/style_corpus/,
which is gitignored; only the manifest under results/ is committed evidence):

  harvest   ADS metadata+abstracts for the FULL corpus, one jsonl.gz per era
  select    seeded stratified pick of the ~2,000-paper full-text sample
  download  arXiv e-print source where an arXiv id exists, else the ADS PDF
  manifest  results/stylecorpus_manifest.json (counts, bytes, checksums)

Usage:
    uv run python scripts/style_corpus_acquire.py            # run all phases
    uv run python scripts/style_corpus_acquire.py --phase harvest
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import stylecorpus as sc  # noqa: E402

MANIFEST_OUT = Path("results/stylecorpus_manifest.json")


def phase_harvest() -> None:
    meta_dir = sc.corpus_dir() / "metadata" / "ads"
    for era in sc.ERAS:
        out = meta_dir / f"{era.label}.jsonl.gz"
        if out.exists():
            print(f"harvest {era.label:>9}: cached ({out.stat().st_size / 1e6:.1f} MB)")
            continue
        q = sc.ads_radio_query(era.lo, era.hi, bibstems=None)
        records = list(sc.ads_harvest(q))
        n_bytes = sc.write_jsonl_gz(records, out)
        print(f"harvest {era.label:>9}: {len(records)} records ({n_bytes / 1e6:.1f} MB)")


def phase_select(seed: int) -> list[dict]:
    sel_path = sc.corpus_dir() / "selection.json"
    if sel_path.exists():
        selection: list[dict] = json.loads(sel_path.read_text())
        print(f"select: cached ({len(selection)} papers)")
        return selection
    meta_dir = sc.corpus_dir() / "metadata" / "ads"
    selection = []
    for era in sc.ERAS:
        records = sc.read_jsonl_gz(meta_dir / f"{era.label}.jsonl.gz")
        picked = sc.stratified_pick(records, era.fulltext_target, seed=seed + era.lo)
        for rec in picked:
            selection.append(
                {
                    "bibcode": rec["bibcode"],
                    "era": era.label,
                    "year": rec.get("year"),
                    "bibstem": (rec.get("bibstem") or ["?"])[0],
                    "doctype": rec.get("doctype"),
                    "arxiv_id": sc.arxiv_id_from_identifiers(rec.get("identifier", [])),
                }
            )
        print(f"select {era.label:>9}: {len(picked)} of {len(records)}")
    tmp = sel_path.with_suffix(".json.part")
    tmp.write_text(json.dumps(selection, indent=1) + "\n")
    tmp.replace(sel_path)
    return selection


def phase_download(selection: list[dict]) -> None:
    done = failed = 0
    for i, item in enumerate(selection):
        path: Path | None
        if item["arxiv_id"]:
            path = sc.fetch_arxiv_source(item["arxiv_id"])
        else:
            path = sc.fetch_ads_pdf(item["bibcode"])
        if path is None:
            failed += 1
        else:
            done += 1
        if (i + 1) % 100 == 0:
            print(f"download: {i + 1}/{len(selection)} (ok={done} failed={failed})", flush=True)
    print(f"download: complete — ok={done} failed={failed}")


def _has_fulltext(item: dict) -> bool:
    root = sc.corpus_dir()
    if item["arxiv_id"]:
        if (root / "arxiv_src" / f"{item['arxiv_id'].replace('/', '_')}.eprint").exists():
            return True
    safe = item["bibcode"].replace("/", "_").replace("&", "+")
    return (root / "ads_pdf" / f"{safe}.pdf").exists()


def phase_topup(selection: list[dict], seed: int) -> list[dict]:
    """Replace selected papers whose full text is unavailable.

    The primary selection round-robins over (journal x doctype) cells, which gives
    the many tiny never-scanned venues equal weight with ApJ/MNRAS — measured
    ~90% 404s in the pre-1992 strata. Replacements are therefore drawn from the
    core-journal subset (where ADS scanned coverage is dense), flagged ``topup``,
    and counted per era in the manifest so the induced drift is visible."""
    meta_dir = sc.corpus_dir() / "metadata" / "ads"
    core = set(sc.BIBSTEMS)
    chosen = {s["bibcode"] for s in selection}
    for era in sc.ERAS:
        have = sum(1 for s in selection if s["era"] == era.label and _has_fulltext(s))
        if have >= era.fulltext_target:
            continue
        pool = [
            r
            for r in sc.read_jsonl_gz(meta_dir / f"{era.label}.jsonl.gz")
            if r["bibcode"] not in chosen and bool(core & set(r.get("bibstem") or []))
        ]
        round_no = 0
        while have < era.fulltext_target and pool:
            round_no += 1
            want = era.fulltext_target - have
            batch = sc.stratified_pick(
                pool, min(2 * want, len(pool)), seed=seed + era.lo + 1000 * round_no
            )
            batch_codes = {r["bibcode"] for r in batch}
            pool = [r for r in pool if r["bibcode"] not in batch_codes]
            for rec in batch:
                if have >= era.fulltext_target:
                    break
                item = {
                    "bibcode": rec["bibcode"],
                    "era": era.label,
                    "year": rec.get("year"),
                    "bibstem": (rec.get("bibstem") or ["?"])[0],
                    "doctype": rec.get("doctype"),
                    "arxiv_id": sc.arxiv_id_from_identifiers(rec.get("identifier", [])),
                    "topup": True,
                }
                ok = (
                    sc.fetch_arxiv_source(item["arxiv_id"])
                    if item["arxiv_id"]
                    else sc.fetch_ads_pdf(item["bibcode"])
                )
                selection.append(item)
                chosen.add(item["bibcode"])
                if ok is not None:
                    have += 1
        print(
            f"topup {era.label:>9}: fulltext {have}/{era.fulltext_target} (pool left {len(pool)})",
            flush=True,
        )
    sel_path = sc.corpus_dir() / "selection.json"
    tmp = sel_path.with_suffix(".json.part")
    tmp.write_text(json.dumps(selection, indent=1) + "\n")
    tmp.replace(sel_path)
    return selection


def _dir_stats(d: Path, suffix: str) -> dict[str, int]:
    files = list(d.glob(f"*{suffix}")) if d.is_dir() else []
    return {"n": len(files), "bytes": sum(f.stat().st_size for f in files)}


def phase_manifest(selection: list[dict], seed: int) -> None:
    root = sc.corpus_dir()
    meta_dir = root / "metadata" / "ads"
    per_era: dict[str, dict] = {}
    for era in sc.ERAS:
        f = meta_dir / f"{era.label}.jsonl.gz"
        records = sc.read_jsonl_gz(f)
        per_era[era.label] = {
            "harvested": len(records),
            "metadata_bytes": f.stat().st_size,
            "metadata_sha256": hashlib.sha256(f.read_bytes()).hexdigest(),
            "selected": sum(1 for s in selection if s["era"] == era.label),
            "selected_with_arxiv_id": sum(
                1 for s in selection if s["era"] == era.label and s["arxiv_id"]
            ),
            "topup_draws": sum(1 for s in selection if s["era"] == era.label and s.get("topup")),
            "with_fulltext": sum(
                1 for s in selection if s["era"] == era.label and _has_fulltext(s)
            ),
        }
    have_src = {p.stem for p in (root / "arxiv_src").glob("*.eprint")}
    have_pdf = {p.stem for p in (root / "ads_pdf").glob("*.pdf")}
    missing = [
        s["bibcode"]
        for s in selection
        if not (
            (s["arxiv_id"] and s["arxiv_id"].replace("/", "_") in have_src)
            or s["bibcode"].replace("/", "_").replace("&", "+") in have_pdf
        )
    ]
    payload = {
        "slice": "stylecorpus",
        "stage": "acquire",
        "is_real": True,
        "selection_seed": seed,
        "selection_sha256": hashlib.sha256(
            json.dumps(selection, sort_keys=True).encode()
        ).hexdigest(),
        "per_era": per_era,
        "arxiv_src": _dir_stats(root / "arxiv_src", ".eprint"),
        "ads_pdf": _dir_stats(root / "ads_pdf", ".pdf"),
        "selected_total": len(selection),
        "selected_missing_fulltext": len(missing),
        "missing_bibcodes": sorted(missing),
    }
    MANIFEST_OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = MANIFEST_OUT.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(MANIFEST_OUT)
    print(f"manifest: {MANIFEST_OUT} (missing fulltext for {len(missing)} of {len(selection)})")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--phase",
        choices=["harvest", "select", "download", "topup", "manifest", "all"],
        default="all",
    )
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    if args.phase in ("harvest", "all"):
        phase_harvest()
    selection: list[dict] = []
    if args.phase in ("select", "download", "topup", "manifest", "all"):
        selection = phase_select(args.seed)
    if args.phase in ("download", "all"):
        phase_download(selection)
    if args.phase in ("topup", "all"):
        selection = phase_topup(selection, args.seed)
    if args.phase in ("manifest", "all"):
        phase_manifest(selection, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
