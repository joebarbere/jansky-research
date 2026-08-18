#!/usr/bin/env python3
"""Stage-1 scoping run for the pre-LLM style corpus (plan 89).

Enumerates the radio-astronomy literature 1933-2021 (ADS per-era counts, needs
ADS_API_TOKEN; arXiv per-era counts, unauthenticated), measures a stratified
per-era PDF size sample, and writes the committed evidence to
results/stylecorpus_scoping.json. Network-heavy and resumable: pass --resume to
reuse counts/sizes already present in the output file.

Usage:
    uv run python scripts/style_corpus_scoping.py                 # full run
    uv run python scripts/style_corpus_scoping.py --skip-ads      # arXiv legs only
    uv run python scripts/style_corpus_scoping.py --resume        # fill in the gaps
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import stylecorpus as sc  # noqa: E402

OUT_DEFAULT = Path("results/stylecorpus_scoping.json")


def _load_partial(out: Path, resume: bool) -> dict:
    if resume and out.exists():
        return json.loads(out.read_text())
    return {}


def _save(out: Path, payload: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--skip-ads", action="store_true", help="no ADS token available")
    ap.add_argument("--skip-sizes", action="store_true", help="counts only")
    ap.add_argument("--size-sample", type=int, default=50, help="size probes per era")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--resume", action="store_true", help="reuse counts/sizes in --out")
    args = ap.parse_args(argv)

    prior = _load_partial(args.out, args.resume)
    ads_counts: dict[str, int] = dict(prior.get("ads_counts_unrestricted", {}))
    ads_restricted: dict[str, int] = dict(prior.get("ads_counts_journal_restricted", {}))
    arxiv_counts: dict[str, int] = dict(prior.get("arxiv_counts", {}))
    sizes: dict[str, list[int]] = {
        k: list(v) for k, v in prior.get("sampled_sizes_bytes", {}).items()
    }
    queries: dict[str, str] = {
        "ads_unrestricted": sc.ads_radio_query(1933, sc.CUTOFF_YEAR - 1, bibstems=None),
        "ads_journal_restricted": sc.ads_radio_query(1933, sc.CUTOFF_YEAR - 1),
    }

    # ---- counts -------------------------------------------------------------
    for era in sc.ERAS:
        if not args.skip_ads and era.label not in ads_counts:
            ads_counts[era.label] = sc.ads_count(sc.ads_radio_query(era.lo, era.hi, bibstems=None))
            ads_restricted[era.label] = sc.ads_count(sc.ads_radio_query(era.lo, era.hi))
            print(f"ADS  {era.label:>9}: {ads_counts[era.label]:>7} "
                  f"({ads_restricted[era.label]} in core journals)")
        if era.hi >= sc.ARXIV_FIRST_YEAR and era.label not in arxiv_counts:
            q = sc.arxiv_radio_query(max(era.lo, sc.ARXIV_FIRST_YEAR), era.hi)
            arxiv_counts[era.label] = sc.parse_arxiv_total(sc.arxiv_get(q))
            print(f"arXiv {era.label:>8}: {arxiv_counts[era.label]:>7}")
        _checkpoint(args.out, ads_counts, ads_restricted, arxiv_counts, sizes, queries)

    # ---- per-era size sample ------------------------------------------------
    if not args.skip_sizes:
        for era in sc.ERAS:
            have = len(sizes.get(era.label, []))
            want = args.size_sample - have
            if want <= 0:
                continue
            era_sizes = sizes.setdefault(era.label, [])
            if era.lo >= sc.ARXIV_FIRST_YEAR:
                q = sc.arxiv_radio_query(era.lo, era.hi)
                for aid in sc.arxiv_sample_ids(q, want, seed=args.seed + era.lo):
                    s = sc.arxiv_pdf_size(aid)
                    if s:
                        era_sizes.append(s)
            elif not args.skip_ads:
                q = sc.ads_radio_query(era.lo, era.hi, bibstems=None)
                for bib in sc.ads_sample_bibcodes(q, want * 2, seed=args.seed + era.lo):
                    if len(era_sizes) - have >= want:
                        break
                    p = sc.fetch_ads_pdf(bib)
                    if p is not None:
                        era_sizes.append(p.stat().st_size)
            print(f"sizes {era.label:>8}: n={len(era_sizes)} "
                  f"mean={sum(era_sizes) / max(len(era_sizes), 1) / 1e6:.2f} MB")
            _checkpoint(args.out, ads_counts, ads_restricted, arxiv_counts, sizes, queries)

    # ---- estimate + final payload ------------------------------------------
    counts_for_estimate = ads_counts if ads_counts else {
        e.label: arxiv_counts.get(e.label, 0) for e in sc.ERAS
    }
    est = sc.size_estimate(counts_for_estimate, sizes, seed=args.seed)
    payload = sc.scoping_payload(
        ads_counts=ads_counts,
        ads_counts_restricted=ads_restricted,
        arxiv_counts=arxiv_counts,
        estimate=est,
        queries=queries,
    )
    payload["sampled_sizes_bytes"] = sizes
    _save(args.out, payload)

    total = sum(counts_for_estimate.values())
    total_gb = float(est["total_gb"])  # type: ignore[arg-type]
    ci = est["ci95_bytes"]
    assert isinstance(ci, list)
    ci_lo, ci_hi = float(ci[0]), float(ci[1])
    print(f"\ncorpus: {total} papers; estimated total PDF size "
          f"{total_gb:.0f} GB (95% CI {ci_lo / 1e9:.0f}-{ci_hi / 1e9:.0f} GB)")
    print(f"wrote {args.out}")
    return 0


def _checkpoint(
    out: Path,
    ads_counts: dict[str, int],
    ads_restricted: dict[str, int],
    arxiv_counts: dict[str, int],
    sizes: dict[str, list[int]],
    queries: dict[str, str],
) -> None:
    """Persist partial progress so --resume can pick up after a network failure."""
    _save(out, {
        "slice": "stylecorpus",
        "stage": "scoping-partial",
        "is_real": True,
        "ads_counts_unrestricted": ads_counts,
        "ads_counts_journal_restricted": ads_restricted,
        "arxiv_counts": arxiv_counts,
        "sampled_sizes_bytes": sizes,
        "queries": queries,
    })


if __name__ == "__main__":
    raise SystemExit(main())
