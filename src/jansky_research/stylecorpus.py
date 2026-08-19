"""Pre-LLM radio-astronomy style corpus — scoping, enumeration, and size estimation.

Stage 1 of plan 89 (``plans/89-traditional-style.md``): enumerate "all radio-astronomy
papers submitted before 2022" across NASA ADS (1933-2021, needs ``ADS_API_TOKEN``) and
arXiv astro-ph (1991-2021, unauthenticated), and turn a stratified per-era size sample
into a ballpark estimate — with a bootstrap CI — of what the *entire* corpus would weigh
as PDFs. The estimate is the committed justification for sampling ~2,000 papers in
Stage 2 rather than mirroring the literature.

Pure logic (query builders, Atom parsing, stratification, the estimator) is tested
offline; everything that touches the network is ``# pragma: no cover`` like every other
slice. The CLI lives in ``scripts/style_corpus_scoping.py``.
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "CUTOFF_YEAR",
    "ERAS",
    "EraStratum",
    "ads_radio_query",
    "arxiv_radio_query",
    "corpus_dir",
    "parse_arxiv_ids",
    "parse_arxiv_total",
    "sample_allocation",
    "scoping_payload",
    "size_estimate",
]

#: First calendar year EXCLUDED from the corpus. ChatGPT shipped 2022-11-30; a
#: 2022-01-01 submission cutoff is safely pre-LLM for published prose.
CUTOFF_YEAR = 2022

#: arXiv holds astro-ph from 1992-04 on; earlier eras are reachable only through ADS.
ARXIV_FIRST_YEAR = 1992


@dataclass(frozen=True)
class EraStratum:
    """One era stratum of the corpus.

    ``fulltext_target`` is the Stage-2 full-text sample size for the stratum (the
    approved plan's allocation, totalling 2,000); Stage-1 size sampling uses the
    smaller uniform ``--size-sample`` count instead.
    """

    label: str
    lo: int
    hi: int  # inclusive
    fulltext_target: int


ERAS: tuple[EraStratum, ...] = (
    EraStratum("1933-1949", 1933, 1949, 25),
    EraStratum("1950s", 1950, 1959, 35),
    EraStratum("1960s", 1960, 1969, 150),
    EraStratum("1970s", 1970, 1979, 150),
    EraStratum("1980s", 1980, 1989, 150),
    EraStratum("1990s", 1990, 1999, 250),
    EraStratum("2000s", 2000, 2009, 400),
    EraStratum("2010-2015", 2010, 2015, 400),
    EraStratum("2016-2021", 2016, 2021, 440),
)

#: Abstract terms that pull in the radio literature without drowning in false hits.
RADIO_ABS_TERMS: tuple[str, ...] = ("radio", "pulsar", "VLBI", "interferometer")

#: ADS keyword-field terms (the curated keyword vocabulary, complements ``abs:``).
RADIO_KEYWORDS: tuple[str, ...] = ("radio continuum", "radio lines")

#: The journals that carried radio astronomy across the eras (ADS bibstems).
BIBSTEMS: tuple[str, ...] = (
    "ApJ",
    "ApJL",
    "ApJS",
    "AJ",
    "MNRAS",
    "A&A",
    "PASP",
    "PASA",
    "AuJPh",
    "Natur",
)


def corpus_dir() -> Path:
    """Root of the (gitignored) style-corpus cache under ``data/``."""
    from jansky_research.data import data_dir

    d = data_dir() / "style_corpus"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _radio_clause() -> str:
    abs_part = " OR ".join(f'abs:"{t}"' for t in RADIO_ABS_TERMS)
    kw_part = " OR ".join(f'keyword:"{k}"' for k in RADIO_KEYWORDS)
    return f"({abs_part} OR {kw_part})"


def ads_radio_query(
    lo: int,
    hi: int,
    *,
    bibstems: tuple[str, ...] | None = BIBSTEMS,
    refereed: bool = True,
) -> str:
    """Build the ADS ``q`` string for radio-astronomy papers in ``[lo, hi]``.

    ``bibstems=None`` gives the unrestricted-journal variant (the corpus total);
    the bibstem-restricted variant scopes the per-journal tables.
    """
    parts = [_radio_clause(), f"year:{lo}-{hi}", "collection:astronomy"]
    if bibstems is not None:
        stems = " OR ".join(f'"{b}"' for b in bibstems)
        parts.insert(0, f"bibstem:({stems})")
    if refereed:
        parts.append("property:refereed")
    return " AND ".join(parts)


def arxiv_radio_query(lo: int, hi: int) -> str:
    """Build the arXiv export-API ``search_query`` for radio papers in ``[lo, hi]``.

    ``cat:astro-ph*`` covers both the pre-2007 flat ``astro-ph`` and the later
    subcategories. The submittedDate window is [Jan 1 ``lo``, Jan 1 ``hi``+1).
    """
    abs_part = " OR ".join(f"abs:{t}" for t in RADIO_ABS_TERMS)
    return f"cat:astro-ph* AND ({abs_part}) AND submittedDate:[{lo}01010000 TO {hi + 1}01010000]"


_ATOM = "{http://www.w3.org/2005/Atom}"
_OPENSEARCH = "{http://a9.com/-/spec/opensearch/1.1/}"


def parse_arxiv_total(atom_xml: str) -> int:
    """Extract ``opensearch:totalResults`` from an arXiv export-API Atom feed."""
    root = ET.fromstring(atom_xml)
    node = root.find(f"{_OPENSEARCH}totalResults")
    if node is None or node.text is None:
        raise ValueError("no opensearch:totalResults in arXiv response")
    return int(node.text)


def parse_arxiv_ids(atom_xml: str) -> list[str]:
    """Extract bare arXiv ids (e.g. ``astro-ph/9601001`` or ``2101.01234``) from a feed."""
    root = ET.fromstring(atom_xml)
    ids = []
    for entry in root.findall(f"{_ATOM}entry"):
        node = entry.find(f"{_ATOM}id")
        if node is None or node.text is None:
            continue
        m = re.search(r"arxiv\.org/abs/(.+?)(v\d+)?$", node.text)
        if m:
            ids.append(m.group(1))
    return ids


def sample_allocation(counts: dict[str, int]) -> dict[str, int]:
    """Cap each era's Stage-2 full-text target by the papers that actually exist.

    No redistribution: a shortfall is reported, not papered over.
    """
    targets = {e.label: e.fulltext_target for e in ERAS}
    return {label: min(targets[label], max(0, counts.get(label, 0))) for label in targets}


def size_estimate(
    counts: dict[str, int],
    sizes: dict[str, list[int]],
    *,
    n_boot: int = 2000,
    seed: int = 0,
) -> dict[str, object]:
    """Estimate the total PDF bytes of the full corpus from a per-era size sample.

    ``counts`` are the per-era paper counts (unrestricted ADS query); ``sizes`` are
    measured PDF byte sizes for the sampled papers of each era. Point estimate is
    sum(count x era mean); the CI bootstraps the sampled sizes within each era
    (resampling papers, NOT refetching), so it carries sampling error only — a
    limitation stated in the payload rather than hidden.

    An era with papers but no measured sizes falls back to the pooled mean of all
    samples and is flagged ``pooled_fallback`` per stratum.
    """
    pooled = [s for v in sizes.values() for s in v]
    if not pooled:
        raise ValueError("no sampled sizes at all — nothing to estimate from")
    rng = np.random.default_rng(seed)
    pooled_arr = np.asarray(pooled, dtype=float)

    per: dict[str, dict[str, object]] = {}
    arrays: dict[str, np.ndarray] = {}
    total_point = 0.0
    for era in ERAS:
        n_papers = int(counts.get(era.label, 0))
        own = sizes.get(era.label) or []
        fallback = n_papers > 0 and not own
        arr = np.asarray(own, dtype=float) if own else pooled_arr
        arrays[era.label] = arr
        mean = float(arr.mean())
        est = n_papers * mean
        total_point += est
        per[era.label] = {
            "count": n_papers,
            "n_sampled": len(own),
            "mean_bytes": mean,
            "est_bytes": est,
            "pooled_fallback": fallback,
        }

    boots = np.empty(n_boot)
    for b in range(n_boot):
        tot = 0.0
        for era in ERAS:
            n_papers = int(counts.get(era.label, 0))
            if n_papers == 0:
                continue
            arr = arrays[era.label]
            tot += n_papers * float(rng.choice(arr, size=arr.size, replace=True).mean())
        boots[b] = tot
    lo, hi = np.percentile(boots, [2.5, 97.5])

    return {
        "per_stratum": per,
        "total_bytes": total_point,
        "total_gb": total_point / 1e9,
        "ci95_bytes": [float(lo), float(hi)],
        "n_boot": n_boot,
        "seed": seed,
        "note": (
            "CI reflects size-sampling error only; corpus enumeration (query recall) "
            "is not bootstrapped."
        ),
    }


def scoping_payload(
    *,
    ads_counts: dict[str, int],
    ads_counts_restricted: dict[str, int],
    arxiv_counts: dict[str, int],
    estimate: dict[str, object],
    queries: dict[str, str],
) -> dict[str, object]:
    """Assemble the committed ``results/stylecorpus_scoping.json`` payload."""
    return {
        "slice": "stylecorpus",
        "stage": "scoping",
        "is_real": True,
        "cutoff_year": CUTOFF_YEAR,
        "eras": [
            {"label": e.label, "lo": e.lo, "hi": e.hi, "fulltext_target": e.fulltext_target}
            for e in ERAS
        ],
        "ads_counts_unrestricted": ads_counts,
        "ads_counts_journal_restricted": ads_counts_restricted,
        "arxiv_counts": arxiv_counts,
        "fulltext_allocation": sample_allocation(ads_counts),
        "size_estimate": estimate,
        "queries": queries,
    }


# --------------------------------------------------------------------------------------
# Network runners (real corpus access; exercised by scripts/style_corpus_scoping.py)
# --------------------------------------------------------------------------------------

ADS_API = "https://api.adsabs.harvard.edu/v1/search/query"
ARXIV_API = "https://export.arxiv.org/api/query"
ARXIV_DELAY_S = 3.5  # export API asks for 3 s between calls; it throttles hard


def _ads_token() -> str:  # pragma: no cover - env/network
    token = os.environ.get("ADS_API_TOKEN", "").strip()
    if not token:
        # The conventional location used by the `ads` python package.
        dev_key = Path.home() / ".ads" / "dev_key"
        if dev_key.exists():
            token = dev_key.read_text().strip()
    if not token:
        raise RuntimeError(
            "No ADS token: set ADS_API_TOKEN or write the token to ~/.ads/dev_key. "
            "Create one free at https://ui.adsabs.harvard.edu/user/settings/token."
        )
    return token


def ads_count(query: str) -> int:  # pragma: no cover - network
    """Number of ADS records matching ``query`` (rows=0, one API call)."""
    import requests

    params: dict[str, str | int] = {"q": query, "rows": 0}
    resp = requests.get(
        ADS_API,
        params=params,
        headers={"Authorization": f"Bearer {_ads_token()}"},
        timeout=60,
    )
    resp.raise_for_status()
    return int(resp.json()["response"]["numFound"])


def ads_sample_bibcodes(query: str, n: int, *, seed: int = 0) -> list[str]:  # pragma: no cover
    """Draw ~``n`` bibcodes spread across the query's result set (network).

    ADS has no server-side random sort; we page at seeded random offsets (sorted by
    bibcode for determinism given a fixed corpus) and draw without replacement.
    """
    import requests

    rng = np.random.default_rng(seed)
    headers = {"Authorization": f"Bearer {_ads_token()}"}
    total = ads_count(query)
    if total == 0:
        return []
    page = 200
    n_pages = max(1, -(-n * 2 // page))  # oversample 2x, then thin
    starts = (
        sorted(int(s) for s in rng.choice(max(total - page, 1), size=n_pages, replace=False))
        if total > page
        else [0]
    )
    pool: list[str] = []
    for start in starts:
        params: dict[str, str | int] = {
            "q": query,
            "rows": page,
            "start": start,
            "fl": "bibcode",
            "sort": "bibcode asc",
        }
        resp = requests.get(ADS_API, params=params, headers=headers, timeout=60)
        resp.raise_for_status()
        pool.extend(d["bibcode"] for d in resp.json()["response"]["docs"])
    pool = sorted(set(pool))
    if len(pool) <= n:
        return pool
    return [str(b) for b in rng.choice(pool, size=n, replace=False)]


def arxiv_get(query: str, *, start: int = 0, max_results: int = 1) -> str:  # pragma: no cover
    """One polite export-API call, returning the raw Atom XML (network).

    ``max_results`` must be >= 1: the API 500s on ``max_results=0`` (measured
    2026-08-17), so counting queries fetch one entry and read ``totalResults``.
    Retries with backoff on 5xx and on the plain-text "Rate exceeded." body.
    """
    import time

    import requests

    params: dict[str, str | int] = {
        "search_query": query,
        "start": start,
        "max_results": max(1, max_results),
    }
    last = ""
    for delay in (ARXIV_DELAY_S, 15.0, 60.0):
        time.sleep(delay)
        resp = requests.get(ARXIV_API, params=params, timeout=60)
        last = resp.text
        if resp.ok and last.strip() != "Rate exceeded." and "api/errors" not in last:
            return last
    raise RuntimeError(f"arXiv export API kept failing for {query!r}: {last[:200]}")


def arxiv_sample_ids(query: str, n: int, *, seed: int = 0) -> list[str]:  # pragma: no cover
    """Draw ~``n`` arXiv ids spread across the query's result set (seeded offsets; network)."""
    rng = np.random.default_rng(seed)
    total = parse_arxiv_total(arxiv_get(query))
    if total == 0:
        return []
    page = 100
    n_pages = max(1, -(-n * 2 // page))
    starts = (
        sorted(int(s) for s in rng.choice(max(total - page, 1), size=n_pages, replace=False))
        if total > page
        else [0]
    )
    pool: list[str] = []
    for start in starts:
        pool.extend(parse_arxiv_ids(arxiv_get(query, start=start, max_results=page)))
    pool = sorted(set(pool))
    if len(pool) <= n:
        return pool
    return [str(i) for i in rng.choice(pool, size=n, replace=False)]


def remote_size(url: str) -> int | None:  # pragma: no cover - network
    """Content-Length of ``url`` without downloading it (HEAD, then Range fallback)."""
    import requests

    resp = requests.head(url, allow_redirects=True, timeout=60)
    cl = resp.headers.get("content-length")
    if resp.ok and cl and int(cl) > 0:
        return int(cl)
    resp = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60)
    cr = resp.headers.get("content-range", "")
    m = re.search(r"/(\d+)$", cr)
    return int(m.group(1)) if m else None


def arxiv_pdf_size(arxiv_id: str) -> int | None:  # pragma: no cover - network
    """Byte size of the arXiv-rendered PDF for ``arxiv_id``."""
    import time

    time.sleep(1.0)  # polite pacing against the PDF CDN
    return remote_size(f"https://arxiv.org/pdf/{arxiv_id}")


def fetch_ads_pdf(bibcode: str, dest_dir: Path | None = None) -> Path | None:  # pragma: no cover
    """Download the ADS-served article PDF for ``bibcode`` (atomic, cached; network).

    Uses the link gateway's ADS_PDF esource (the scanned-literature service for the
    pre-arXiv era). Returns None when ADS holds no full text for the bibcode.
    """
    import requests

    dest_dir = dest_dir or (corpus_dir() / "ads_pdf")
    dest_dir.mkdir(parents=True, exist_ok=True)
    safe = bibcode.replace("/", "_").replace("&", "+")
    target = dest_dir / f"{safe}.pdf"
    if target.exists():
        return target
    url = f"https://ui.adsabs.harvard.edu/link_gateway/{bibcode}/ADS_PDF"
    tmp = target.with_suffix(".part")
    try:
        with requests.get(url, stream=True, allow_redirects=True, timeout=120) as resp:
            # Any failure (404 no fulltext, 5xx on the scanned-article service)
            # means "this sample point is unavailable", not "abort the sweep".
            if not resp.ok or "pdf" not in resp.headers.get("content-type", ""):
                return None
            with open(tmp, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    fh.write(chunk)
    except requests.RequestException:
        return None
    tmp.replace(target)
    return target
