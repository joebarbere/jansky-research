# Findings — the pulsar P–Ṗ diagram (ATNF)

`jansky_research.ppdot` turns the ATNF Pulsar Catalogue's periods and spin-down rates into the
P–Ṗ diagram — the H–R diagram of pulsars — deriving surface fields, ages, and the death line, and
classifying the population. Reuses `pulsarspec.is_millisecond` and the project's ATNF VizieR fetch.
Tooling + real-data + recover-a-known done together (VizieR is a single reliable query).

## Data

ATNF Pulsar Catalogue (Manchester et al. 2005) via VizieR `B/psr` (public, no auth): periods `P0` (s)
and period derivatives `P1` (Ṗ). Of the catalogue, **2 052** pulsars have a positive measured Ṗ and
enter the diagram.

## Recover-a-known: the population structure

**Surface field** $B = 3.2\times10^{19}\sqrt{P\dot P}$ G cleanly separates the three populations:

| class | N | median $\log_{10}(B/\mathrm{G})$ | interpretation |
|---|---|---|---|
| millisecond ($P<30$ ms) | 182 | **8.42** | recycled by accretion — low field, fast spin |
| normal | 1 787 | **12.05** | the bulk rotation-powered population |
| magnetar / high-$B$ ($B>10^{13}$ G) | 83 | **13.31** | the top-right of the diagram |

The $\sim$5-orders-of-magnitude ($\sim$10$^{4.9}$, a factor of $\sim$80 000) spread in median field
between millisecond and high-$B$ pulsars is the textbook P–Ṗ structure (Lorimer & Kramer 2004): MSPs
sit at the bottom-left (short $P$, tiny $\dot P$), ordinary pulsars at $B\sim10^{12}$ G, magnetars at
the top-right. We use the standard $B = 3.2\times10^{19}\sqrt{P\dot P}$ G (orthogonal rotator,
$\sin\alpha=1$, $I=10^{45}$ g cm², $R=10$ km); the aligned/equatorial convention gives a coefficient
of $\sim$6.4$\times10^{19}$, i.e. all our fields scale up by $\sqrt2$.

**Death line.** With the constant-$B/P^2$ polar-cap criterion ($B_{12}/P^2=0.2$; Ruderman & Sutherland
1975; Bhattacharya & van den Heuvel 1991), **98.2%** of the catalogue lies *above* the death line — i.e.
almost every catalogued radio pulsar is on the radio-loud side, as expected for a radio-selected sample.
This is **model-dependent** (the line spans a "death valley", Chen & Ruderman 1993), and the $\sim$1.8%
(36 pulsars; this line originally said "~37") below the line are not noise — they are real long-period pulsars (e.g. PSR
J2144$-$3933) that pushed the death-line literature. The Crab pulsar validates the derivations:
$B = 3.8\times10^{12}$ G and $\tau = P/2\dot P \approx 1.26\times10^3$ yr, both textbook.

## Honest assessment & caveats

- **A reproduction, not a discovery.** The tool recovers the well-known P–Ṗ population structure from a
  public catalogue; the contribution is a tested, reproducible pipeline.
- **"Magnetar" is really "high-$B$".** Our $B>10^{13}$ G cut (83 sources) captures both true magnetars
  ($\sim$30 known) and high-field rotation-powered pulsars; it is a position-in-diagram label, not an
  emission-mechanism classification (most true magnetars are X-ray, not radio, sources).
- **The death line is one model.** The constant-$B/P^2$ line is a representative criterion; published
  death lines span a "death valley" (Chen & Ruderman 1993) and depend on the gap model, so 98.2% should
  be read as "almost all, for a standard line," not a precise number.
- **Only pulsars with a positive measured $\dot P$ enter the diagram.** Of the 2536 ATNF entries
  (this bullet originally said "~3 500", a recollection corrected below on 2026-08-24), 484 are
  excluded: 435 have a period but no measured $\dot P$, 35 (mostly globular-cluster pulsars) show a
  *negative* apparent $\dot P$ from cluster acceleration (not real spin-down), 5 record exactly
  zero, and 9 lack a period. An acceleration-corrected treatment would return some of those; here
  they are correctly dropped.
- **Selection, not volume-limited.** The ATNF catalogue is the union of many flux- and
  period-limited surveys (MSPs and faint/long-period pulsars are under-represented), so the *relative*
  population sizes are survey-shaped, not intrinsic; the per-class B-fields are robust to this.
- **$\tau=P/2\dot P$ assumes braking index 3 and $P_0\ll P$** — a characteristic age, not the true age.
- **The synthetic fixture is a round-trip code check.** Its `classify_accuracy` (only reported for the
  offline run) confirms the classifier on injected truth labels; it is *not* a validation of the real
  ATNF results, which have no ground-truth class labels.
- **Reproducible:** `python -m jansky_research.ppdot` regenerates the metrics, the P–Ṗ-diagram figure,
  and the macros from the public VizieR catalogue.

## Referee round on the style conversion (2026-08-24)

**Restyle finding, fixed.** Collapsing an em-dash appositive to commas turned a three-item list
(field, age, luminosity) into a four-item one, so "the standard orthogonal-rotator estimate" read as
a *separate quantity the code computes* rather than a gloss on the $B$ formula. That parenthesis is
the only statement in the paper that the three headline fields are convention-dependent: under the
aligned/equatorial coefficient every quoted value shifts by log10(sqrt2) = 0.15 dex
(8.42 -> 8.57, 12.05 -> 12.20, 13.31 -> 13.46). Re-fenced with a relative clause and semicolons.

## Two pre-existing defects, NOT fixed here (queued for a real-run follow-up)

**1. The Crab "validation" is in the abstract, is hand-typed, and never touches the real data path.**
The abstract says "The Crab pulsar validates the derivations ($B=3.8\times10^{12}$ G,
$\tau\approx1260$ yr)". Neither number is in `results/ppdot_metrics.json` (nine keys, none
Crab-related) or in `generated/macros.tex` (ten macros, none Crab-related). The only Crab check in
the repo is `tests/test_ppdot.py::test_magnetic_field_and_age_crab`, which evaluates the closed forms
at **hand-entered literature values** and asserts a wide window. The real run *cannot* identify the
Crab at all: `fetch_atnf_ppdot` requests the `PSRJ` column and then discards it. Placing the claim in
Results, after the 2052-pulsar census, invites the reader to infer the Crab was found in the analysed
sample. Ask what would make this test fail: only a typo in a two-term formula. It exercises none of
the real path's failure modes (the VizieR `P1` units, row parsing, the positive-Pdot cut).

**Fix (needs a real re-run, VizieR confirmed reachable):** keep `PSRJ`, look up J0534+2200 in the
analysed sample, emit `\ppCrabB`/`\ppCrabTau`. Then "validates" is earned.

**2. "Of the ~3500 entries" is wrong, and the reproducibility universal is false.** Queried directly:
`Vizier B/psr/psr` returns **2536** rows, which is also the number `survey/pulsarspec-findings.md`
records for the same table fetched by the same call. So ppdot's ~3500 is a recollection, not a
measurement, and the two slices disagreed about the size of one table. Neither number was auditable
from `results/`.

Relatedly, the paper says "Every number above is written by the pipeline into the macros this
manuscript `\input`s" and the file's own comment says "none typed by hand". Both are false: the two
Crab numbers, "~3500 entries", "a factor of ~80,000", "~1.8\%" and "~30 true magnetars" are all
typed. In a paper whose claimed contribution *is* reproducibility, a false universal about
reproducibility is the one sentence that has to be right. Note ~80,000 and ~1.8% are arithmetic on
`\ppLogBmagnetar - \ppLogBmsp` and `1 - \ppFracAlive`: correct today, silently stale after a re-run.

**Fix:** emit the fetched row count as a macro, and scope the universal to catalogue-derived numbers.

### Resolved 2026-08-24 (real re-run)

`fetch_atnf_ppdot` now keeps the `PSRJ` column it was already requesting, and
`named_pulsar_derived` looks a pulsar up **in the analysed sample**. The Crab is therefore read from
the fetched catalogue rather than from literature constants, which exercises the `P1` column units,
the row parsing and the positive-Pdot cut --- none of which the old unit test on hand-entered values
could reach.

| quantity | before | after (from the fetched row) |
|---|---|---|
| Crab B | 3.8e12 (typed) | `\ppCrabB` = **3.8e12** |
| Crab tau | ~1260 yr (typed) | `\ppCrabAge` = **1257** yr |
| Crab P | not stated | `\ppCrabPeriod` = 0.0333924 s |
| catalogue size | "~3500 entries" | `\ppNcatalogue` = **2536** |
| field span | "a factor of ~80,000" (typed) | `\ppBSpanFactor` = **78000** (4.89 dex) |
| below the death line | "~1.8%" (typed) | `\ppFracDead` = **1.8** |

The "~3500" was a recollection and wrong; 2536 is the same number `pulsarspec` gets from the
identical call, so the two slices no longer disagree about the size of one table. The hand-typed
arithmetic (~80,000, ~1.8%) is now derived in `run()`, so it cannot go stale after a re-run.

The false reproducibility universal is scoped rather than deleted: "Every **catalogue-derived**
number above is written by the pipeline ... (the ~30 true magnetars and the Crab's ~970-yr
historical age are literature figures, not outputs of this analysis)". The file's own header comment
said the same false thing and was corrected too.

Both re-runs were **purely additive**: every previously published value (n_pulsars 2052,
frac_above_death 0.982, the three median fields) is unchanged, and no macro changed value.

## Full referee round (2026-08-25): MINOR REVISION (upper end), 12 findings

Every catalogue-derived number reproduces: the referee independently re-derived the census
from VizieR `B/psr` TAP (n=2052, 182/1787/83, medians 8.42/12.05/13.31, 36 below the death
line, Crab B=3.794e12 / τ=1256.8 from the fetched row) and the committed figure contains
exactly 2052 marker draws. All five DOIs clean. What holds it up: three claims that outrun
their committed evidence, one undisclosed data property, and the missing per-object table.
All fixes are additive to one real re-run.

**MAJORs:**
1. The "magnetar" median is a property of the cut, not a population: the class is defined by
   the quantity being reported, so 13.31 tracks the threshold directly (12.71 at B>3e12 →
   14.12 at >3e13), and the 4.89-dex "span" moves 4.29–5.70 with it. The measured claim is
   the MSP↔normal separation (3.63 dex, stable to both cuts: 8.37–8.48 / 12.04–12.06).
   Keep the number, change the claim strength, commit the sweep.
2. The catalogue snapshot is ~2017 vintage (publication_date 2017-07-18; no J0901−4046, max
   P0 = 11.79 s) and the paper never says so — which bites the death-line caveat hardest,
   since the sample excludes the entire ultra-long-period class. One sentence in Data +
   `catalogue_version`/`access_date` in the JSON.
3. 98.25% above the death line has no committed threshold sweep though the paper's own caveat
   says the line is "one representative model": measured, "almost all" survives a factor of 2
   in the constant (92.4%) and not a factor of 5 (75.6%). Emit `\ppFracAliveLo/Hi`.
4. "The per-class median fields are robust to this" is asserted with no committed check — and
   is true (S1400-cut excursions ≤ 0.11 dex, the magnetar median moving most because a radio
   flux requirement removes radio-quiet high-B objects, 53/83 survive). Commit it, scoped as
   "robust to a flux cut" (B/psr has no survey-provenance column, so per-survey jackknife is
   unavailable).
5. A 2052-object census committed as 17 scalars — no reader can recompute a single median or
   see which 36 objects are below the line. Commit `results/ppdot_pulsars.csv` (PSRJ, P0, P1,
   B, tau, Edot, class, above_death).

**MINOR/NIT:** "validates" for N=1 with no assertion that can fail on a re-run — widen to
J1939+2134 (B=4.1e8), J1550−5418 (2.2e14), J2144−3933 (2.1e12) with a tolerance assert; the
discard breakdown is asserted not measured (444 null Ṗ / 35 negative / 5 exactly zero, the
last silently dropped by `pd > 0` and mentioned nowhere; 31 of 35 negatives are in globular
clusters); `ppdot.py:200` comment and this file's own "~3 500 entries"/"~37 pulsars" still
assert defects recorded as fixed; `\ppAccuracy` is a dead `--` macro written un-namespaced by
the offline leg only (a placeholder is not defended by `preserve_live_macros` — the
`\tiiNEvents` shape; rename `\ppSynAccuracy` or drop); 0.982-fraction vs 1.8%-percent in one
paper, and `\ppBSpanFactor` = 78000 carries three noise digits (unrounded medians give 77278
— quote 8×10⁴); silent `idx[0]` on name lookup (no ambiguity today, one-line guard); the
gitignored arxiv-submission tarball is stale (2026-06-29, would publish the retracted
values) — regenerate before submission.

**Status: fixes pending.** The single change: commit the per-pulsar table and the two sweeps
— the sweeps mostly *support* the paper; the one claim they do not support is the magnetar
median, which should be described as threshold-set.

**Status: RESOLVED (2026-08-25).** All five MAJORs closed by one real re-run, and every sweep the
referee computed server-side reproduces in the committed pipeline.

1. **Magnetar median reframed as threshold-set.** `magnetar_threshold_sweep` commits the
   referee's table exactly (12.71 at B>3e12 -> 14.12 at >3e13); the abstract and Results now
   present the MSP<->normal separation (`\ppMspNormalSpanDex` = 3.64 dex) as the measured claim
   -- backed by its own committed period-cut sweep (`msp_period_sweep`: MSP median 8.37-8.48
   over 20-50 ms, normal 12.04-12.06) -- and label 13.31 as tracking its own defining cut.
2. **Vintage disclosed.** Data names the frozen snapshot (published 2017-07-18), its 11.8-s
   maximum period, and the absent ultra-long-period class; `catalogue_version` /
   `fetched_utc` / `catalogue_max_p0_s` are in the JSON; the death-line caveat states the
   below-line census is a floor.
3. **Death-valley sweep committed.** `death_line_sweep` over B12/P^2 = 0.05-1.0:
   98.2% at the paper's constant, 92.5% at 2x, 99.6% at half (`\ppFracAlivePctLo/Hi`); the
   Discussion now says "survives a factor of two and not a factor of five" with the numbers.
4. **Flux-cut robustness measured, scoped, committed.** `flux_cut_medians` (S1400 fetched):
   max excursion 0.11 dex, the high-B median moving most (radio-quiet high-B objects removed);
   stated as "under a flux cut", not survey selection, since B/psr has no provenance column.
5. **Per-pulsar table committed.** `results/ppdot_pulsars.csv`: all 2052 analysed rows with
   B, tau, Edot, class, death-line disposition.

MINORs: four named anchors (Crab, B1937+21, J1550-5418, J2144-3933) are read from the fetched
table and *gated* -- `run()` raises if any derived log B leaves its literature window, so a
units regression fails the run (tested with a deliberate 1e6 pdot error); the discard
breakdown is measured (435 null Pdot + 35 negative + 5 exactly zero + 9 lacking P0 = 484,
closing the census arithmetic in Data); the stale `fetch_atnf_ppdot` comment is gone with the
rewrite and this file's "~3 500"/"~37" recollections now carry their corrections inline;
`\ppAccuracy` is dropped from the macro file (offline-only, nothing cited it); the abstract
and Discussion both quote percent (98.2%/1.8%); `\ppBSpanFactor` derives from unrounded
medians at two significant figures (77000, was 78000); `named_pulsar_derived` raises on an
ambiguous name. The stale gitignored arXiv package was regenerated.

Both prior published headline values are unchanged (n 2052, medians 8.42/12.05/13.31,
frac 0.982); the revision added evidence and moved claim strength onto it.
