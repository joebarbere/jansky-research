# stylecorpus — findings

Plan 89 (`plans/89-traditional-style.md`): scope the pre-LLM radio-astronomy literature as the
empirical basis for the `traditional-style` skill. Stage 1 (scoping, 2026-08-17) is complete;
evidence in `results/stylecorpus_scoping.json`.

## Stage 1 — how big is "all of it"?

**Corpus definition.** Refereed astronomy papers matching a radio clause
(`abs:` radio/pulsar/VLBI/interferometer or `keyword:` "radio continuum"/"radio lines"),
1933–2021, submitted before 2022-01-01 (safely pre-ChatGPT). Enumerated on ADS; arXiv
overlaid for 1992+.

**Counts (ADS, unrestricted journals, refereed):** 113,829 papers total.

| era | papers | core-journal subset | mean PDF (sampled) |
|---|---|---|---|
| 1933–1949 | 207 | 110 | 0.50 MB (n=40) |
| 1950s | 843 | 377 | 0.67 MB (n=50) |
| 1960s | 4,027 | 1,696 | 0.43 MB (n=27) |
| 1970s | 12,434 | 5,283 | 1.11 MB (n=50) |
| 1980s | 16,113 | 7,108 | 0.61 MB (n=50) |
| 1990s | 19,718 | 10,025 | 0.87 MB (n=49) |
| 2000s | 24,203 | 12,911 | 0.45 MB (n=50) |
| 2010–2015 | 16,609 | 9,163 | 1.33 MB (n=50) |
| 2016–2021 | 19,675 | 10,893 | 3.16 MB (n=50) |

arXiv holds 39,207 radio papers 1992–2021 (about half of the 80,205 ADS-refereed papers over
the same span); the arXiv-only slice alone would weigh ~65 GB.

**Ballpark total PDF size of the full corpus: ~138 GB (95% CI 116–164 GB).**
Point estimate is Σ per-era count × sampled mean; the CI bootstraps the per-era size samples.
This is the committed justification for Stage 2's design: a ~2,000-paper stratified sample
(~4 GB at these means) captures the style conventions; mirroring 138 GB buys nothing a style
analysis needs.

**What the estimate does and does not cover.**

- *The CI is size-sampling error only.* Corpus enumeration (query recall) is not bootstrapped:
  the radio clause is a recall proxy whose keyword coverage varies by era, so the *count* — and
  therefore the total — carries an unquantified definitional term. The per-era relative sizes
  are the robust part.
- *Post-1992 sizes are arXiv renders; counts are ADS.* Publisher PDFs for the same papers can
  differ in size from arXiv-rendered ones, so the 2000s+ rows mix a numerator and denominator
  from different services. Direction unknown; likely sub-factor-of-2 at era level.
- *Scanned-era sizes are conditional on ADS serving full text.* The 1960s stratum yielded only
  27/50 sample PDFs; papers ADS cannot serve have unknown sizes and are implicitly assumed to
  match the served ones.
- The rising post-2010 mean (0.45 → 3.16 MB) is real — figure-heavy modern papers — and means
  "all papers as PDFs" is dominated by the *least* traditional era.

**Stage-2 feasibility.** Every stratum's full-text target is satisfiable
(`fulltext_allocation` equals the plan's targets in all nine strata; even 1933–1949 has
207 papers against a target of 25). 266 scanned-era PDFs (186 MB) are already cached under
`data/style_corpus/ads_pdf/` from the size sample and count toward the Stage-2 sample.

**Operational notes.** The arXiv export API 500s on `max_results=0` — count queries must
fetch one entry and read `opensearch:totalResults`; it also rate-limits aggressively
(minutes-long bans), so all calls back off and the CLI checkpoints per era for `--resume`.
The ADS scanned-article service intermittently 504s on individual old articles; a failed
article is a lost sample point, not an abort.
