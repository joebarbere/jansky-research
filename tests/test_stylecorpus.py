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


def test_era_of() -> None:
    assert sc.era_of(1933) == "1933-1949"
    assert sc.era_of(1969) == "1960s"
    assert sc.era_of(2021) == "2016-2021"
    assert sc.era_of(1932) is None
    assert sc.era_of(2022) is None


def test_arxiv_id_from_identifiers() -> None:
    assert (
        sc.arxiv_id_from_identifiers(["2021MNRAS.500....1X", "arXiv:2101.01234", "10.1093/x"])
        == "2101.01234"
    )
    assert sc.arxiv_id_from_identifiers(["astro-ph/9601001"]) == "astro-ph/9601001"
    assert sc.arxiv_id_from_identifiers(["1969Ap&SS...4..464G"]) is None
    assert sc.arxiv_id_from_identifiers([]) is None


def test_jsonl_gz_roundtrip(tmp_path) -> None:
    records = [{"bibcode": "b1", "year": 1969}, {"bibcode": "b2", "year": 1970}]
    path = tmp_path / "sub" / "era.jsonl.gz"
    n_bytes = sc.write_jsonl_gz(records, path)
    assert n_bytes == path.stat().st_size > 0
    assert sc.read_jsonl_gz(path) == records
    assert not path.with_suffix(".gz.part").exists()


def _fake_records(n_per_cell: int) -> list[dict]:
    recs = []
    for stem in ("ApJ", "MNRAS", "AJ"):
        for doctype in ("article", "letter"):
            for i in range(n_per_cell):
                recs.append({"bibcode": f"{stem}-{doctype}-{i}", "bibstem": [stem],
                             "doctype": doctype})
    return recs


def test_stratified_pick_spreads_across_cells_and_is_deterministic() -> None:
    records = _fake_records(20)  # 6 cells x 20
    picked = sc.stratified_pick(records, 12, seed=1)
    assert len(picked) == 12
    cells = {(r["bibstem"][0], r["doctype"]) for r in picked}
    assert len(cells) == 6  # round-robin touches every cell before repeating any
    again = sc.stratified_pick(records, 12, seed=1)
    assert [r["bibcode"] for r in picked] == [r["bibcode"] for r in again]
    assert len({r["bibcode"] for r in picked}) == 12


def test_stratified_pick_returns_everything_when_short() -> None:
    records = _fake_records(1)
    assert sc.stratified_pick(records, 100, seed=0) == records


SAMPLE_TEX = r"""
\documentclass{aastex631}
\begin{document}
\begin{abstract}
We measure a thing. Crucially, the result --- obtained with care --- is new.
\end{abstract}
\section{Introduction}
% a comment that should vanish entirely
The flux is $S_\nu = 3$ mJy \citep{smith1970}. \emph{It is a ratio, not a difference.}
First, we do one thing. Second, we do another thing. Finally, we conclude the sweep.
This paper is described by this paper itself, and the reader should note it.
The spectrum was measured by the receiver. I present the census here today.
\begin{itemize}
\item a list item
\end{itemize}
\subsection{The method, and what it is worth}
We delve into the data. Moreover, the fit is shown in Figure 1.
\begin{equation}
x = y
\end{equation}
\end{document}
"""


def test_strip_latex_removes_markup_but_keeps_prose() -> None:
    prose = sc.strip_latex(SAMPLE_TEX)
    assert "comment that should vanish" not in prose
    assert "MATH" in prose  # inline $...$ replaced
    assert "REF" in prose  # \citep replaced
    assert "equation" not in prose.lower() or "x = y" not in prose
    assert "It is a ratio, not a difference." in prose  # \emph unwrapped
    assert "\\" not in prose.replace("\\%", "")


def test_latex_section_titles_and_abstract() -> None:
    titles = sc.latex_section_titles(SAMPLE_TEX)
    assert titles == ["Introduction", "The method, and what it is worth"]
    abstract = sc.latex_abstract(SAMPLE_TEX)
    assert abstract is not None and "We measure a thing" in abstract
    assert sc.latex_abstract("no abstract here") is None


def test_split_sentences_handles_abbreviations() -> None:
    got = sc.split_sentences("The flux (e.g. at 1.4 GHz) is high. It fell later. Fig. 2 shows.")
    assert got[0].startswith("The flux")
    assert len(got) == 3


def test_fingerprint_prose_rates() -> None:
    fp = sc.fingerprint_prose(sc.strip_latex(SAMPLE_TEX))
    assert fp["rule_of_three"] == 1.0
    assert fp["first_singular_per_kw"] > 0  # "I present"
    assert fp["we_per_kw"] > 0
    assert fp["hedge_per_kw"] > 0  # crucially, delve, moreover
    assert fp["self_ref_per_kw"] > 0  # "this paper" twice
    assert fp["reader_addr_per_kw"] > 0
    assert fp["passive_per_sentence"] > 0  # "was measured"
    assert fp["n_sentences"] >= 5


def test_fingerprint_latex_structural_metrics() -> None:
    fp = sc.fingerprint_latex(SAMPLE_TEX)
    assert fp["em_dash_per_kw"] > 0  # --- in the abstract survives via prose? checked on tex
    assert fp["emph_per_kw"] > 0
    assert fp["emph_sentence_start_per_kw"] > 0  # ". \emph{It is..." after REF period
    assert fp["itemize_envs"] == 1.0
    assert fp["n_section_titles"] == 2.0
    assert fp["section_title_mean_words"] == pytest.approx((1 + 7) / 2)
    assert fp["abstract_words"] > 5


def test_aggregate_fingerprints_percentiles() -> None:
    docs = [{"m": float(v)} for v in range(1, 101)]
    agg = sc.aggregate_fingerprints(docs)
    assert agg["m"]["p50"] == pytest.approx(50.5)
    assert agg["m"]["n"] == 100.0
    assert agg["m"]["mean"] == pytest.approx(50.5)
    assert sc.aggregate_fingerprints([]) == {}


def test_lint_paper_flags_only_exceedances() -> None:
    corpus = {
        "em_dash_per_kw": {"p50": 0.0, "p90": 0.5},
        "abstract_words": {"p50": 126.0, "p90": 260.0},
        "rule_of_three": {"p50": 0.0, "p90": 0.0},
    }
    fp = {"em_dash_per_kw": 10.7, "abstract_words": 200.0, "rule_of_three": 1.0}
    found = sc.lint_paper(fp, corpus)
    kinds = {(sev, metric) for sev, metric, _ in found}
    assert ("HIGH", "em_dash_per_kw") in kinds
    assert ("MED", "rule_of_three") in kinds
    assert not any(m == "abstract_words" for _, m, _ in found)  # under p90
    assert sc.lint_paper({"em_dash_per_kw": 0.0}, corpus) == []


def test_macro_names_from_tex() -> None:
    tex = "\\newcommand{\\drFoo}{1.2}\n\\newcommand{\\drBar}{--}\n% \\newcommand{\\x}{y}\n"
    assert sc.macro_names_from_tex(tex) == {"drFoo", "drBar", "x"}


GUARD_TEX = r"""
The fraction is \drFoo\% \citep{smith2020, jones1999} at 5 arcsec. % comment 7.7
We measure \drFoo\ again and 3.14 twice: 3.14. \software{astropy \citep{astropy}}
"""


def test_guard_signature_counts() -> None:
    sig = sc.guard_signature(GUARD_TEX, {"drFoo", "drBar"})
    assert sig["macros"] == {"drFoo": 2}
    assert sig["cites"]["smith2020"] == 1 and sig["cites"]["jones1999"] == 1
    assert sig["numbers"]["3.14"] == 2 and sig["numbers"]["5"] == 1
    assert "7.7" not in sig["numbers"]  # comments stripped
    assert len(sig["software"]) == 1


def test_compare_signatures_prose_only_edit_passes() -> None:
    macros = {"drFoo", "drBar"}
    before = sc.guard_signature(GUARD_TEX, macros)
    reworded = GUARD_TEX.replace("We measure", "The survey yields")
    assert sc.compare_signatures(before, sc.guard_signature(reworded, macros)) == []
    broken = GUARD_TEX.replace("3.14 twice", "3.15 twice")
    problems = sc.compare_signatures(before, sc.guard_signature(broken, macros))
    assert any("3.14" in p for p in problems) and any("3.15" in p for p in problems)
    dropped = GUARD_TEX.replace("\\drFoo\\ again", "it again")
    problems = sc.compare_signatures(before, sc.guard_signature(dropped, macros))
    assert any("drFoo" in p for p in problems)
