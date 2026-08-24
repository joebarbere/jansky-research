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
($\sim$37 pulsars) below the line are not noise — they are real long-period pulsars (e.g. PSR
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
- **Only pulsars with a positive measured $\dot P$ enter the diagram.** Of the $\sim$3 500 ATNF entries,
  $\sim$1 400 are excluded: many have no measured $\dot P$, and globular-cluster pulsars often show a
  *negative* apparent $\dot P$ from cluster acceleration (not real spin-down). An acceleration-corrected
  treatment would return some of those; here they are correctly dropped.
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
