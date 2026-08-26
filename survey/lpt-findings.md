# Findings — a verified, provenance-carrying LPT population catalogue (plan 35)

`jansky_research.lpt` + `data/lpt_sample.csv` (every value carrying its arXiv provenance — compiled
2026-07 from the discovery papers, cross-checked against the Rea+2026 review arXiv:2601.10393 and
its GitHub data file). GATE-2 note: the review's own Fig. 3 already plots the class on P–Ṗ — the
contribution here is the provenance-typed table + regenerable statistics, not the first diagram.
Companion of GLEAM-X J0704-37 corrected to M3 (not M5).

**v3 update (plan 44, 2026-07-08):** the catalogue was extended in place from 13 → **16** with the
three 2026 discoveries (ASKAP J142431.2−612611; VASTER's J165130.3−450520 and J170036.6−445758),
each coordinate verified against its source-name convention. The stats below are the v3 (N=16)
values; the period-split test moved to Δlog P = 0.176, **p_perm = 0.52 — still not significant**
(the 2 long-period VASTER sources have unreported binary status, weakening the WD-binary/rest
period contrast). The Stokes-V forced-photometry leg lives in the sibling `lptv` slice.

## What the compilation itself found

- **A transcription error in the review's own machine-readable table**: ASKAP J1935+2148 period
  2225.309 s vs the discovery paper's 3225.313±0.002 s (dropped leading digit). Flagged; the
  discovery value is used.
- **Only 2 Ṗ measurements exist in the whole class**: CHIME J0630+25 (+5.2e-12, through a glitch —
  the review conservatively downgrades it to a limit; both framings carried) and CHIME
  J1634+44 (−9.03e-12 — the class's only firm spin-UP, natural for a binary).
- Post-review member included: ASKAP J1745−5051 (accreting WD binary, Rose+2026).

## Population statistics (all regenerate from the CSV)

| quantity | value |
|---|---|
| N (confirmed, v3 2026-07-08) | 16 |
| WD binaries / candidates | 7 |
| X-ray detected | 3 |
| period range | 7.0 min – 6.48 hr (median 73.4 min) |
| Ṗ-constrained objects below the pulsar death line | **9 / 9** |
| WD-binary vs rest period offset | Δlog P = 0.176, **p_perm = 0.52 — NOT significant at N=16** |
| same test on a synthetic real split | p ≈ 0.022 (the test has power; the non-detection is informative) |

The 9/9 below-death-line fact is the class's central puzzle made quantitative; the hinted ~78-min
binary boundary is genuinely open, not yet established.

## Counterpart cross-match (live, 2026-07-02; `lpt.crossmatch_counterparts`)

Per-object VLASS QL2 cone (20″, epoch 1) + LoTSS DR3 forced cutout peak at each LPT position:

- **VLASS (10/13 in Dec>−40 coverage): no persistent 2–4 GHz counterpart above ~0.7 mJy (5σ QL)
  for ANY object.**
- **LoTSS DR3 (3/13 in returned coverage: ASKAP J1935+2148, ILT J1101+5521, CHIME J1634+44): all
  undetected in the mosaic** (peak/rms mJy: 0.30/0.15, 0.01/0.008, −0.02/0.019). ILT J1101's
  non-detection in the *mosaic* is consistent with its burst-only conjunction emission.
- Remaining objects: outside VLASS dec range (3) / LoTSS returned HTTP 500 = no footprint (10).

Reading: the class is burst-only at current survey depths — persistent emission ≲0.7 mJy at
2–4 GHz and ≲0.04–0.76 mJy (5σ) at 144 MHz where covered. Kept in findings (not the paper) since
these numbers don't flow through the macro pipeline yet.

## Honesty rails

- Compilation, not discovery: every number is someone else's measurement with provenance.
- NS dipole quantities (B, τ) NOT assigned to binary members (orbital periods).
- GCRT J1745−3009 flagged as the weakest member; disputed values carry flags in the CSV.
- Selection effects on the period distribution are not modelled (stated).

## Full referee round (2026-08-26): MAJOR REVISION, 18 findings, two BLOCKERs

The underlying work is sound: every metrics field re-derives exactly from the CSV through the
pure functions, no macro renders as `--`, all six checkable DOIs pass Crossref, and the
transcription-bug accusation against Rea+2026's table HOLDS (Crossref settles it: the Caleb et
al. discovery paper is titled "...with a 54-minute period" — 3225.313 s = 53.76 min; the
review's 2225.309 s = 37.09 min contradicts the discovery paper's own title). The null on the
78-min boundary is honest. What blocks acceptance is claim/artifact correspondence.

**BLOCKER 1: the committed figure no longer matches its own data.** figures/lpt_ppdot.pdf was
written 2026-08-04; the CSV was corrected twice since (2026-08-14, 2026-08-21). The figure
plots ASKAP J1832−0911 at Ṗ = 9.0e-12 — a factor of 109 below its true published limit —
under the name the CSV's own flag calls erroneous, and its hard-coded title reads "13
objects, 2026-07" against an abstract saying 16. (The metrics JSON was checked: the
corrections move no rounded metric — median still 73.4 min — so the JSON is stale in
provenance, not value; the figure is stale in value.)

**BLOCKER 2: a catalogue paper with no catalogue.** No table of the 16 objects, no
data-availability statement, no pointer to the CSV, and ~10 of 16 discovery papers uncited —
including the object that sets \lptPmax (6.48 hr), which has no reference anywhere in the
manuscript. The claimed advance is provenance, and none of it is delivered to a reader.

**MAJORs (all verified by exact enumeration/computation):**
- "An EXACT permutation test gives p = 0.5219" — it is 20,000 random shuffles (MC SE ±0.0035);
  the true exact enumeration over all C(16,7) = 11,440 partitions gives p = 0.5258. Enumerate
  (it takes a second) or drop the word.
- Three unknown-status objects (incl. the two longest-period VASTER sources) are counted as
  "no companion", and the labelling moves the headline: Δlog P 0.176 (published) → 0.413 with
  unknowns excluded (exact p 0.526 → 0.167). The bias runs toward the published null, but the
  effect size in the abstract is the smallest of four defensible labellings and no sensitivity
  is reported. Primary should be 7v6 excluding unknowns + the sensitivity table.
- The power demo is easy by construction and misdescribed: "same size" is actually N=13, 4v9
  (defaults), and the injected classes are DISJOINT (Δlog P ≈ 1.03 dex ≈ 6× observed; seed 0
  is a low outlier at 0.349). At the OBSERVED offset the test's power at N=16 is ≈0.058 —
  indistinguishable from its size — and a hard 78-min boundary is already falsified by
  inspection (3 of 7 binaries below it, 3 of 9 non-binaries above; Fisher p = 0.615).
- "9/9 below the death line" is one structural fact, not nine findings: sitting above the line
  would need B = 3.6e16–1.1e20 G (1.5–5 dex above any magnetar; above the NS virial limit for
  the three longest periods), characteristic ages 0.75–2280 yr, drifts up to 1.5e4 s/yr. No
  achievable measurement could have falsified it; and three of the nine are WD binaries the
  paper itself says the criterion doesn't apply to. Report the margin distribution (2–6 dex,
  closest 42×), the field threshold, the death-valley sweep (9/9 across B/P² = 5e10–2e12, min
  margin 2.64×), and 6/6 for the no-companion subset as the claim-carrying count.
- "Verified" overclaims against the table's own history: the CSV's flags record five defects
  found by a downstream audit AFTER the paper was committed (100× Ṗ error, 59″ wrong name,
  three wrong discovery arXiv IDs — one pointing at a CAD-modelling paper, one at a
  dark-matter paper). The manuscript, edited the day after the fixes landed, mentions none —
  while accusing a published review of a transcription error. Disclose and rename to
  "cross-checked".

**MINOR/NIT:** no mechanical pin test protects the CSV that six modules consume (the 100×
error lived 7 weeks; add a pinned-literature cell-by-cell test; lptv.J1839_PERIOD_S duplicates
a CSV cell); lptxray coverage/cones coordinate caches currently consistent (checked, 1e-6 deg)
but unguarded; the metrics JSON's source string says "13 LPTs" above n_lpt=16 (guard keys on
"synthetic" so protection holds; derive the string); deruiter2025 title paraphrased (Crossref:
"Sporadic radio pulses from a white dwarf binary at the orbital period"); the 78-min boundary
hypothesis has no citation; death line cited to Ruderman & Sutherland but implemented as the
B/P² form (Bhattacharya & van den Heuvel 1991) with the 2e11 choice undiscussed; the
accusation sentence's uncertainties and the review-file version/access date live only in the
manuscript (add pdot_err/period_err columns + provenance note); stale gitignored
arxiv-submission tarball holds the 2026-07-02 manuscript and pre-correction figure — delete;
\lptNx emitted, unused; --offline is a no-op flag (and hence no mode-dependent macro hazard —
the one namespacing check this slice passes for free); "no published catalogue carries..." →
"we are not aware".

**Verdict: MAJOR REVISION.** The single change: put the catalogue IN the paper — a
machine-readable 16-row table with per-value provenance and a reference column citing all 16
discovery papers, plus a data-availability statement with a versioned DOI. That one addition
delivers the claimed contribution and exposes/fixes four other findings at once.

**Status: RESOLVED (2026-08-26).** One regeneration from the (already-corrected) CSV plus the
new exact/labelling/margin/power machinery; every referee number reproduced in the committed
pipeline (Δ 0.413/0.176/0.446/0.605 with exact p 0.1667/0.5258/0.0919/0.1905; min margin 42.3;
power 0.062 at the published offset).

1. **BLOCKER 1**: the figure regenerates from the corrected CSV (J1832−0911 at its true Ṗ
   limit under its corrected name) with the title derived from the sample size — "16 objects",
   never hard-coded again. The stale local arxiv-submission tarball was deleted and rebuilt.
2. **BLOCKER 2**: the paper now CONTAINS the catalogue — a generated 16-row deluxetable
   (`generated/table.tex`: period, Ṗ + type incl. m(d), companion status with *unknown shown
   as unknown*, X-ray, and a discovery-reference column citing all 16 papers — 11 new
   Crossref/arXiv-verified bib entries) — plus a data-availability paragraph with the Zenodo
   concept DOI and the one-command rebuild.
3. **"Exact" is exact**: `period_split_stat` enumerates all C(16,7)=11,440 partitions
   (p = 0.5258 where the MC seed said 0.5219) and reports its method string; MC only as a
   fallback above 200k partitions.
4. **The labelling is primary evidence**: `label_sensitivity` commits all four labellings; the
   PRIMARY comparison excludes the three unknown-status objects (7v6: Δ = 0.413, p = 0.167)
   and the paper discloses that the earlier version counted two of the three longest periods
   with the null class silently. No labelling reaches significance — conclusion unchanged,
   effect size honest.
5. **The power statement is measured**: `split_power` (vectorized exact test over simulated
   overlapping classes at the sample's 1.75-dex width) gives power 0.172 at the observed
   primary offset and 0.062 at the published one — essentially the size of the test — quoted
   in place of "the test has power"; the synthetic disjoint demo is rerun at the real
   composition (16, 7v9) and labelled as verifying machinery, not sensitivity. The paper also
   notes a hard 78-min boundary is already disfavoured by inspection (3 binaries below, 3
   non-binaries above).
6. **"9/9" reframed as structural**: `death_line_margins` commits per-object margins (42.3 to
   1.2e6 — no achievable measurement could have falsified the count), the death-valley sweep
   (9/9 across B/P² = 5e10–2e12, min margin 10.6; Chen & Ruderman + Bhattacharya & van den
   Heuvel now cited, the constant named as a choice), and the claim-carrying no-companion
   count (6/6).
7. **"Verified" → "cross-checked"**, with a disclosure paragraph naming the five defects the
   2026-08 downstream audit fixed in this table (100× Ṗ, 59″ name, three wrong arXiv IDs) and
   the accusation's own confirmation route (the discovery paper's "54-minute" title).
8. **The CSV is mechanically pinned**: `test_csv_matches_pinned_literature` asserts every
   period/Ṗ/type/discovery-ID cell (a future edit must change the pin deliberately);
   `test_lptxray_coordinate_caches_match_the_csv` guards the duplicated coordinates;
   `test_lptv_ephemeris_constant_matches_the_csv` guards lptv's hard-coded J1839 period. The
   CSV itself was also repaired: three rows' flags contained unquoted commas that silently
   truncated them for any named-column reader (found when the round-10 rewrite crashed on
   them), and the quoted uncertainties for the accusation-bearing values now live in
   `period_err_s`/`pdot_err_s_s` columns with the review-file comparand recorded in flags.
9. MINORs: the metrics source string derives from the data ("16 LPTs"); deruiter2025's title
   corrected to the published one; the 78-min boundary hypothesis cites rea2026; ruderman1975
   cited for the death-line physics where it belongs; \lptNx used (first LPT X-ray discovery
   sentence, citing wang2025); "no published catalogue" → "we are not aware of".
