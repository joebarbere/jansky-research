"""Offline tests for the style-corpus scoping logic (plan 89, Stage 1)."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import stylecorpus as sc

ATOM_FIXTURE = """<?xml version="1.0" encoding="UTF-8"?>
<feed xmlns="http://www.w3.org/2005/Atom"
      xmlns:opensearch="http://a9.com/-/spec/opensearch/1.1/">
  <opensearch:totalResults>12345</opensearch:totalResults>
  <entry>
    <id>http://arxiv.org/abs/astro-ph/9601001v2</id>
    <title>An old paper</title>
  </entry>
  <entry>
    <id>http://arxiv.org/abs/2101.01234v1</id>
    <title>A newer paper</title>
  </entry>
  <entry>
    <title>An entry with no id survives parsing</title>
  </entry>
</feed>
"""


def test_eras_cover_1933_to_cutoff_without_gaps_or_overlap() -> None:
    years: list[int] = []
    for era in sc.ERAS:
        assert era.lo <= era.hi
        years.extend(range(era.lo, era.hi + 1))
    assert years == list(range(1933, sc.CUTOFF_YEAR))


def test_fulltext_targets_total_2000() -> None:
    assert sum(e.fulltext_target for e in sc.ERAS) == 2000


def test_ads_query_restricted_and_unrestricted() -> None:
    q = sc.ads_radio_query(1960, 1969)
    assert "year:1960-1969" in q
    assert "property:refereed" in q
    assert 'bibstem:("ApJ"' in q
    assert 'abs:"pulsar"' in q
    assert 'keyword:"radio continuum"' in q
    q_all = sc.ads_radio_query(1933, 2021, bibstems=None, refereed=False)
    assert "bibstem" not in q_all
    assert "refereed" not in q_all
    assert "collection:astronomy" in q_all


def test_arxiv_query_window_is_inclusive_of_hi_year() -> None:
    q = sc.arxiv_radio_query(2000, 2009)
    assert "cat:astro-ph*" in q
    assert "submittedDate:[200001010000 TO 201001010000]" in q
    assert "abs:radio" in q


def test_parse_arxiv_total_and_ids() -> None:
    assert sc.parse_arxiv_total(ATOM_FIXTURE) == 12345
    assert sc.parse_arxiv_ids(ATOM_FIXTURE) == ["astro-ph/9601001", "2101.01234"]


def test_parse_arxiv_total_missing_raises() -> None:
    with pytest.raises(ValueError, match="totalResults"):
        sc.parse_arxiv_total("<feed xmlns='http://www.w3.org/2005/Atom'></feed>")


def test_sample_allocation_caps_at_available_counts() -> None:
    counts = {e.label: 10_000 for e in sc.ERAS}
    counts["1933-1949"] = 7  # scarcer than the target of 25
    alloc = sc.sample_allocation(counts)
    assert alloc["1933-1949"] == 7
    assert alloc["2016-2021"] == 440
    missing = sc.sample_allocation({})
    assert set(missing) == {e.label for e in sc.ERAS}
    assert all(v == 0 for v in missing.values())


def _toy_inputs() -> tuple[dict[str, int], dict[str, list[int]]]:
    counts = {e.label: 0 for e in sc.ERAS}
    sizes: dict[str, list[int]] = {}
    counts["1960s"] = 100
    sizes["1960s"] = [1_000_000, 3_000_000]  # mean 2 MB
    counts["2016-2021"] = 50
    sizes["2016-2021"] = [2_000_000, 2_000_000]
    return counts, sizes


def test_size_estimate_point_value_and_determinism() -> None:
    counts, sizes = _toy_inputs()
    est = sc.size_estimate(counts, sizes, n_boot=200, seed=42)
    # 100 * 2 MB + 50 * 2 MB = 300 MB
    assert est["total_bytes"] == pytest.approx(300e6)
    assert est["total_gb"] == pytest.approx(0.3)
    lo, hi = est["ci95_bytes"]  # type: ignore[misc]
    assert 0 < lo <= hi
    again = sc.size_estimate(counts, sizes, n_boot=200, seed=42)
    assert again["ci95_bytes"] == est["ci95_bytes"]
    per = est["per_stratum"]
    assert isinstance(per, dict)
    assert per["1960s"]["n_sampled"] == 2
    assert not per["1960s"]["pooled_fallback"]


def test_size_estimate_pooled_fallback_is_flagged() -> None:
    counts, sizes = _toy_inputs()
    counts["1970s"] = 10  # papers exist, no sizes measured
    est = sc.size_estimate(counts, sizes, n_boot=50, seed=0)
    per = est["per_stratum"]
    assert isinstance(per, dict)
    assert per["1970s"]["pooled_fallback"]
    assert per["1970s"]["mean_bytes"] == pytest.approx(np.mean([1e6, 3e6, 2e6, 2e6]))


def test_size_estimate_with_no_samples_raises() -> None:
    with pytest.raises(ValueError, match="no sampled sizes"):
        sc.size_estimate({"1960s": 5}, {})


def test_scoping_payload_schema() -> None:
    counts, sizes = _toy_inputs()
    est = sc.size_estimate(counts, sizes, n_boot=50, seed=0)
    payload = sc.scoping_payload(
        ads_counts=counts,
        ads_counts_restricted={k: v // 2 for k, v in counts.items()},
        arxiv_counts={"2016": 1000},
        estimate=est,
        queries={"ads_unrestricted": "q"},
    )
    assert payload["slice"] == "stylecorpus"
    assert payload["is_real"] is True
    assert payload["cutoff_year"] == sc.CUTOFF_YEAR
    assert payload["fulltext_allocation"] == sc.sample_allocation(counts)
    labels = [e["label"] for e in payload["eras"]]  # type: ignore[union-attr]
    assert labels == [e.label for e in sc.ERAS]


def test_corpus_dir_under_data_dir(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("JANSKY_RESEARCH_DATA_DIR", str(tmp_path))
    d = sc.corpus_dir()
    assert d == tmp_path / "style_corpus"
    assert d.is_dir()
