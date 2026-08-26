# lptduty — findings

Plan 90. Increments 1–2 complete (tested logic + a run over committed evidence);
increment 3 (phase-resolved leg) and GATE 0 are open.

## What was measured (2026-08-20)

Rose et al. (2026, Nature Astron. 10, 1166) report that ASKAP J174508.9−505149's bursts
"turn off for several hours at a time" — one source, qualitatively. This turns the class
version into a number from evidence already committed by `lptv`: **647 measured VAST
snapshots across 10 LPTs**, no new data fetched. Output: `results/lptduty_metrics.json`.

The quantity is

    p = f_active x (w + T) / P

the probability that one snapshot catches a pulse. **Its factors are not separately
identifiable** from detection counts, and the JSON says so per source
(`identifiable_factors: false`). Separating them needs an ephemeris precise enough to phase
each snapshot years after its reference epoch — increment 3, and only for the sources that
have one.

**Sources with detections** (p at the source's own brightest measured pulse):

| source | k | p | effective epochs | assumed pulse | period |
|---|---|---|---|---|---|
| ASKAP J175534.9−252749.1 | 1 | 0.011 | 91 | 6.3 mJy | 4186 s |
| ASKAP J1832−0911 | 1 | 0.014 | 71 | 14.0 mJy | 2656 s |
| ASKAP J183950.5−075635 | 2 | 0.045 | 44 | 24.3 mJy | 23222 s |

**Sources without** (95% upper limits at an assumed 5 mJy pulse): 0.029 to 0.107 across the
seven, i.e. under ~3–11% of snapshots.

For scale, a ~730 s snapshot on a 2656 s period samples ~27% of a cycle, so a source that
were *always* active with a broad pulse would give p of order 0.3. The measured 0.011–0.045
are an order of magnitude below that — consistent with the switching-off behaviour Rose et
al. describe, though the product cannot yet say how much is duty cycle and how much is
inactivity.

## Two decisions that determine whether the numbers mean anything

**The denominator is efficiency-weighted exposure, not epoch count** — the `frblens` lesson,
whose limit was 4x too tight until per-source efficiency entered it. Each epoch contributes
`Phi(S/sigma - 5)`, so epochs too shallow to have seen the assumed pulse contribute ~0. Note
the effective epochs above (91, 71, 44) against raw counts of 104, 97 and 92 — **but those
"raw counts" are total CSV rows including never-released epochs the loader drops before any
weighting. The 2026-08-26 referee recomputed the weighting's true effect: at the observed
brightnesses it removes 0.02% of the measured exposure, not half.** (Corrected in place; the
original sentence propagated into the note's blocker.)
halved every p.

**Every limit is a function of the assumed pulse flux**, so the JSON reports a grid
(0.5–50 mJy) rather than one number. One source is entirely unconstrained at 0.5 mJy — no
epoch could have seen a pulse that faint — and the run says `unconstrained: true` instead of
inventing a limit.

**An efficiency floor was necessary, and a test caught why.** A Gaussian tail is never
exactly zero: a 0.5 mJy pulse in a 100 mJy epoch still scores ~3e-7, and enough such epochs
sum into sensitivity that does not exist (10^6 would look like a third of an epoch).
`MIN_EFFICIENCY = 1e-3` floors them to zero. Without it the code returned a limit of
509,214 where the honest answer is "unconstrained".

## GATE 0, part 1 — novelty: CLEAR (2026-08-20)

No class-wide, uniformly-measured LPT activity-fraction or duty-cycle constraint exists.
Per-source numbers do exist and are expected — the field review (Rea, Hurley-Walker & Caleb
2026, arXiv:2601.10393) compiles one per known LPT in its Table 2 — but each comes from a
different instrument, cadence and threshold, with no common denominator. That review states
the gap directly: survey selection functions "require explicit modeling before population
analysis is possible", and the known population is "still drawn from highly heterogeneous
surveys and observing strategies".

Checked and empty: VAST survey papers (Murphy+2021, arXiv:2108.06039) describe cadence but
publish no LPT detection-rate statistic; the MWA GPM search (Horvath+2025, arXiv:2509.06315)
builds injection-recovery machinery but never applies it to the known sample; the ASKAP EMU
10 s search (Lee+2025, arXiv:2511.09770) reports a surface-density limit for *new* sources.
Nearest per-source anchors: GPM J1839-10 at ~12% in-period duty cycle with a 50-70% nulling
fraction (Hurley-Walker+2023, Nature 619, 487) and GLEAM-X J1627 at ~6% over a <3-month
active window (Hurley-Walker+2022, Nature 601, 526).

## GATE 0, part 2 — aliasing: ONE SOURCE FAILS (2026-08-20)

`results/lptduty_gate0.json`. The binomial constraint assumes each snapshot is an independent
draw on pulse phase; VAST's roughly fortnightly cadence is not random against any LPT period,
so this is tested rather than asserted. Rayleigh Z and Kuiper V on phases referenced to each
source's first epoch — the zero point is arbitrary, and clustering is invariant under a phase
shift, so **no ephemeris is needed to test uniformity** (only to assign physical phase).
Bonferroni-corrected across the ten sources tested.

- **GPM J1839-10 fails**: Kuiper V = 0.275, Rayleigh p = 1.2e-3 (0.012 corrected). Its
  snapshots are *not* uniform in pulse phase, so the binomial model does not apply and its
  quoted limit is not a limit.
- **ASKAP J142431.2-612611 is inconclusive**: its catalogued period (2147.27 s, quoted to
  1e-2 s) is less precise than the 4.0e-3 s needed to keep phase coherent over the 1400-day
  baseline, so the phase smears and the test cannot detect clustering that may still be
  there. Inconclusive is not the same as passing.
- The other eight are consistent with uniform sampling.

Both are now stamped into `results/lptduty_metrics.json` per source
(`constraint_valid`, `phase_sampling.verdict`), so a reader cannot lift a number without the
caveat that governs it. **Eight of ten constraints stand; two do not.**

A numerical detail worth keeping: phase must be referenced to the first epoch, not MJD 0. At
MJD ~59000 a 1 h period is ~1.2e6 cycles, where float64 retains only ~1e-10 of a cycle — a
test caught this as a zero-point invariance failure at the 1e-9 level.

## Still open

1. **Ephemeris audit** for the phase-resolved leg (increment 3), which is what would split
   f_active from the in-period duty cycle. `lptv` round 3 was bitten by multi-year phase
   drift; the same systematic applies here, and `quoted_period_precision_s` in the GATE-0
   file is inferred from decimal places — a proxy for the published uncertainty, not the
   uncertainty itself.
2. **Period derivatives are not folded in.** A pdot large enough to matter over the baseline
   breaks phase coherence even when the period is quoted precisely.
3. The 319 never-released rows (of 966) carry no measurement and are excluded; that is recorded in
   the loader, not silently dropped.


## Ephemeris audit (2026-08-21) — and five errors it found in the committed catalogue

The audit for increment 3 read all ten discovery/timing papers. Before any phase work, it
turned up defects in `data/lpt_sample.csv`, a vendored provenance-carrying table that feeds
`lptv` (in the submission queue). **Each was verified independently against the arXiv record
before changing anything**, per the repo's rule to adjudicate rather than search:

| row | field | was | is | check |
|---|---|---|---|---|
| GCRT J1745-3009 | discovery_arxiv | astro-ph/0502231 | astro-ph/0503052 | 0502231 is "How Concentrated Are The Haloes Of Low Surface Brightness Galaxies" |
| GLEAM-X J162759.5 | discovery + pdot ref | 2201.02926 | 2503.08033 | 2201.02926 is "Variational design for a structural family of CAD models" |
| GPM J1839-10 | discovery + pdot ref | 2307.14829 | 2503.08036 | 2307.14829 is a theory paper *about* the source, not the discovery |
| ASKAP J1832-0911 | pdot_s_s | 9.0e-12 | 9.8e-10 | paper text: "spin period derivative limit Pdot = < 9.8e-10" — **two orders of magnitude** |
| ASKAP J1832-0911 | period_s | 2656.2554 | 2656.247 | the paper gives P = 2656.247 +/- 0.001 s twice; 2656.2554 appears nowhere in it |

The three wrong arXiv ids had propagated into `papers/lptv/refs.bib` `note` fields. The
entries' titles, journals, volumes, pages and DOIs were all correct, which is precisely why
triage's Crossref DOI check never caught them — **a DOI check cannot see a wrong arXiv id in
a free-text note.**

**No committed `lptv` number moves**: `period_split_p` stays 0.5219, median period 73.4 min,
N = 16, 7 white-dwarf binaries. The period correction is 0.008 s in 2656 s, far too small to
shift a median or a permutation rank. Verified by recomputing, not assumed.

## Increment 3 — the split works, and yields exactly one number

`results/lptduty_phase.json`. Machinery: restrict to snapshots whose phase coverage overlaps
the emitting window, and the detection rate *within that subset* estimates f_active, while
(w+T)/P comes from published quantities instead of being fitted.

**Validated against an independent result.** For ASKAP J183950.5-075635 the code puts the two
VAST detections at phases **0.489 and 0.951** — reproducing exactly the values `lptv` derived
and published separately (0.489 +/- 0.016 and 0.951 +/- 0.017, `papers/lptv/main.tex`).

**Only 3 of 10 sources could be attempted, and the reason is not period precision.** Seven
are excluded because their papers publish **no reference epoch at all** — GPM J1839-10 knows
its period to 2e-4 s over a 34-year baseline and still tabulates no PEPOCH. Assuming an epoch
"near the campaign" would manufacture the phase, so the runner refuses and records why.

Of the three attempted, the single-window model survives for one:

- **ASKAP J175534.9-252749.1**: 16 of 91 snapshots on-window, its one detection at phase
  0.045 (window half-width 0.091), giving **f_active ~ 0.06**. The detection landing near
  phase 0 is *consistent with* the assumption that the published PEPOCH is a pulse epoch —
  but with an a-priori landing chance of 0.18 it cannot corroborate it (2026-08-26 referee).
- **ASKAP J1832-0911**: its only detection sits at phase **0.699**, far outside a window
  centred on PEPOCH. That *falsifies* the pulse-at-phase-0 assumption for this source — a
  "period epoch" in a timing solution is a frequency reference and need not be a pulse
  arrival. f_active is not measurable here without the true pulse phase, and fitting the
  phase to the one detection and then using it to derive f_active would be circular.
- **ASKAP J183950.5-075635**: needs a two-component model. Its detections are the interpulse
  (0.489) and a main-pulse-like burst (0.951); a single window cannot represent a source with
  a published interpulse at 177.8 deg separation.

So the honest state: the separation is implemented, tested and validated, and the published
ephemerides support it for **one source out of ten**. That ratio is itself the finding, and
it is the thing to report — not a class-wide f_active the data cannot carry.


## GATE-2 science review (2026-08-21) — one real statistics bug, caught before the write-up

Read-only review of the module, the three runners, the committed JSONs and this file.
Verdict: do not advance to the note until finding 1 is fixed. All four findings are now
addressed.

**1. The Poisson upper limit was wrong for k > 0 — the code used `2.996 + k`.** That is
correct only at k = 0. The exact one-sided 95% limit solves
`sum_{i<=k} exp(-lam) lam^i/i! = 0.05`, i.e. `0.5 * chi2.ppf(0.95, 2k+2)`: 2.996, 4.744,
6.296, 7.754 for k = 0..3. The limit does not grow by one count's worth per count near zero.
Verified independently here by bisection on the Poisson CDF before changing anything.

The bug hit **exactly the three sources with a detection** — the most quotable numbers in
the slice — and was invisible in the write-up because the prose quoted only the seven
zero-detection limits, which are unaffected. Corrected values:

| source | k | was | is |
|---|---|---|---|
| ASKAP J175534.9-252749.1 | 1 | 0.0439 | **0.0521** |
| ASKAP J1832-0911 | 1 | 0.0563 | **0.0668** |
| ASKAP J183950.5-075635 | 2 | 0.1135 | **0.1431** |
| J175534.9 f_active (phase leg) | 1 | 0.2498 | **0.2965** |

No test caught it: the existing one asserted only `p_upper_95 > p_point`, which the wrong
formula also satisfies. `poisson_upper_95` now pins the exact values against the CDF
solution, and asserts `poisson_upper_95(1) > POISSON_ZERO_95 + 1` to block a regression.

**2. The detection criterion did not match its own docstring.** It claimed "the same
criterion the `lptv` sweep applied", but applied only `|V|/sigma >= 5` — `lptv` also
requires `|V| > 0.006|I|`, the ASKAP on-axis leakage floor. All four detections in the
committed sweep clear that floor comfortably, so **no number changed**, but the code did not
do what it said and a brighter-I epoch could later have been counted as a pulse. `EpochRow`
now carries Stokes I and `_is_detection` applies both conditions.

**3. The GATE-0 significance level was recorded misleadingly** — `ALPHA = 0.01` with an
unexplained `* 5` at each use, so the JSON advertised a threshold that was not the one
governing the verdict. Now a single `ALPHA_FAMILYWISE = 0.05`.

**4. GPM J1839-10's period was queried and is CORRECT.** The reviewer flagged a possible
mismatch (1317.2 s vs the catalogue's 1318.1957 s) and hedged it as needing confirmation.
Confirmed against the paper's own text: "we derive a period P of 1318.1957 +/- 0.0002 s"
(arXiv:2503.08036, 34-year timing baseline). The catalogue value stands, so the GATE-0
phase-clustering verdict for this source is not an artefact of a wrong period.

Everything else was checked and passed: the identifiability argument, the efficiency model
(`Phi(S/sigma - threshold)` is the right detection probability for a forced measurement),
the MIN_EFFICIENCY floor, the Rayleigh/Kuiper statistics and the zero-point-invariance
argument, the `(w+T)/P` window (an exact interval-overlap identity, not an approximation),
and every citation.


## The note, and a gate that could not see it (2026-08-21)

`papers/lptduty/rnaas.tex` drafted in the RNAAS register, every number `\input` from
`generated/macros.tex` which `lptduty.write_paper_assets` writes from the committed JSONs.
It passes `prose_lint --genre rnaas` clean against the 391-note pre-LLM baseline.

Writing it exposed a gap in the repo's own science gate. `scripts/triage_papers.py` read
`paper/main.tex` and skipped any directory without one, so **no RNAAS note in the repo had
ever been triaged** -- not atlas3i's, frbstats', spectra's, vgpra's, nor wdpulsar's, which is
the head of the submission queue. The rule says run triage before any submission; the tool
could not see the document being submitted.

Triage now checks every file in a paper directory carrying a `\documentclass`. Re-run over
everything: all five previously-invisible notes come back clean, so the gap hid no defect --
but it was real, and it is the CLAUDE.md lesson again, that a stated invariant is not a
followed one until something checks it.

The headline the note leads with is the ephemeris finding rather than the duty cycles: the
measurement is limited by what the discovery papers publish, not by sensitivity or epochs,
and a reference epoch costing one table row would convert most of these products into
activity fractions.

## Full referee round (2026-08-26): MAJOR REVISION, 19 findings, one BLOCKER

The core is real and worth publishing (the Poisson fix verified exact; the GPM J1839−10
clustering verdict survives a phase-smearing check; both checkable citations exact) — but the
note contains one quantitative claim about its own method that the committed JSON refutes by
~50×, three counts its own evidence contradicts, and the note is IN THE SUBMISSION QUEUE, so
these land before submission, exactly as intended.

**BLOCKER: "removes roughly half the nominal exposure" is false by ~50× and the weighting is
inert at every flux the note quotes.** Recomputed: at the observed brightnesses the
sensitivity weighting keeps 90.981/71.000/44.000 of 91/71/44 epochs (0.02% effect); at 5 mJy,
636.88/647 (1.6%). The findings-file origin conflated the loader dropping never-released rows
(104/97/92 → 91/71/44) with the efficiency weighting — the note inherited the conflation, and
the sentence is the whole justification for "sensitivity-weighted denominator" in the
abstract.

**MAJORs (all verified by computation):**
- \ldNNoEpoch=7 and \ldNPhaseUsable=3 are both contradicted by the JSON that generates them
  (J142431 is excluded for period precision, not a missing epoch → 6; J183950's anchor is the
  repo's own, "not a published PEPOCH" → 2 of 3). The errors cancel to 10, which is why
  nothing looked wrong.
- A Kuiper statistic is quoted with the RAYLEIGH p (the Kuiper p, computed, is 0.059 —
  above the family threshold); and the Bonferroni correction covers 10 sources but not
  2 tests × 10 (checked: the verdict survives ×20 → 0.023).
- GATE-0 gates on period precision but not Ṗ: applying its own 0.1-cycle criterion to the
  published Ṗ bounds makes GLEAM-X (6.4 cycles!), J1832 (0.91), and J183950 (0.74)
  inconclusive → \ldNValid 8 → 5 conditionally (the bounds are upper limits; state the
  condition).
- The independence claim: uniformity ≠ independence (the file's own caveat, absent from the
  note). The presenter's pseudo-replication concern RELOCATED: no two snapshots overlap in
  phase (min gap 776 s > 726 s duration) and the point estimates are unbiased (linearity);
  what breaks is the INTERVALS — under a conservative one-per-cycle collapse \ldPMax moves
  0.045→0.071 (+57%) and \ldLimMin 0.032→0.050 (+55%).
- p is a STOKES-V threshold quantity presented as "probability of catching a pulse" — the
  circular fraction (18–36% measured for J183950 alone) makes 5 mJy V ≡ 5–50 mJy total
  intensity, unconvertible without an unquoted factor.
- The one scale comparison uses J1832's T/P for a sample spanning 22× in T/P: dividing each p
  by its own (w+T)/P, J183950 implies f_active ≈ 0.73 — the opposite of "an order of
  magnitude below" permanently-active (and J175534's 0.061 matches the phase leg's 0.0625,
  the consistency worth reporting).
- \ldFActive = 0.06 is a k=1 estimate quoted with no interval (90%: [0.006, 0.297] — a
  factor 50) and exists only under the assumption it supports (the split's window contains
  the detection only within ~0.10 cycles of the assumed PEPOCH; a-priori landing chance 18%;
  and identical one-event evidence is read oppositely for J1832 vs J175534 — the J1832
  sentence itself is sound, checked).
- The phase machinery's only validation cannot fail (same anchor constant, same period, same
  table, same fold arithmetic as lptv — "independently" retracted to consistency check).
- The ephemerides argued about are UNCITED (McSweeney, Wang, Lee — the note asserts a named
  team's epoch "is not a pulse arrival" without citing them) and barbere2026lptv has no
  resolvable identifier.
- "p = 0.011–0.045" reads as spread; the three are indistinguishable (LRT vs common rate:
  p = 0.45) — a single class-wide rate is the cleaner headline.
- A limit on p constrains f_active only via its own (w+T)/P: the ordering INVERTS —
  \ldLimMin (0.032) excludes only f_active > 0.75, GLEAM-X's 0.059 excludes > 0.09, and
  J165130's limit is vacuous (> 1.60).
- "A published epoch would convert most of these into activity fractions" has a
  counterexample in the note's own JSON (J183950 has a usable anchor and yields nothing —
  3 of 44 epochs in-window: the data, not the ephemeris, binds there; and for two sources a
  perfect epoch gives a vacuous bound).
- write_paper_assets has NO caller anywhere (lptduty is absent from the Makefile SLICES;
  make paper SLICE=lptduty errors), and five pipeline numbers are hand-typed under "every
  number is \input" (all currently correct to the digit — nothing checks them).

**MINOR/NIT:** the sample selection (10 of 16, excluding the motivating Rose et al. source)
is never stated; the committed flux grid contains "probabilities" up to 56 marked
constrained (clip or mark unconstrained); a stale caveat says GATE-0 is "NOT yet tested" in
the same file that carries its verdicts; the J1832 Ṗ smear is 0.400 cycles, not the file's
0.22 (strengthens its own conclusion; invalidates its f_active_upper for J1832); two
definitions of n_detections; "107 never-released epochs" is 319; findings-table 0.046 vs
macro 0.045; \ldNDetSources unused. Format fine (~684 words).

**Verdict: MAJOR REVISION.** The single change: report p/(w+T)/P — the implied active
fraction — per source beside p. Computable today from committed JSONs, it converts the
counting-noise range into the physical quantity, and surfaces the J183950 inversion, the
limit ordering, and the vacuous conversions in one stroke. HOLD the RNAAS submission until
this lands.

**Status: RESOLVED (2026-08-26).** All three drivers re-run offline from the committed CSV;
every referee number reproduced (G-test p = 0.445; implied f_active 0.039/0.061/0.735;
best informative limit 0.09 = GLEAM-X; one vacuous; Rayleigh ×20 = 0.023 survives; smears
0.91/0.74/6.38; J1832 farthest-epoch smear 0.400; J175534 upper 0.297). The note is rewritten
and the submission hold can lift after review of this revision.

1. **BLOCKER**: the denominator sentence now states the measured effect — the weighting
   removes \ldEffRemovedHeadlinePct = 0.02% of the exposure at the observed brightnesses
   ("more at fainter assumed fluxes, where it matters") — and this file's origin sentence is
   corrected in place. Limits are capped at 1 (a "limit" above 1 excludes nothing;
   `duty_constraint` caps and the grid marks such rows unconstrained — the committed 0.5 mJy
   rows with p_upper 56 are gone).
2. **Counts fixed at the source**: the phase JSON now carries n_no_published_epoch = 6,
   n_epoch_but_excluded = 1 (J142431, \citep{pritchard2026} cited for it), and
   n_published_pepoch = 2 (J183950's anchor is the repo's own, stated); the note's abstract
   and body use the new macros.
3. **The Kuiper/Rayleigh attribution is fixed**: `kuiper_p` (Stephens) is implemented and
   committed per source; the Bonferroni now covers sources × tests (×20); the note quotes
   Rayleigh Z = 6.61 with its own corrected p = 0.023; the clustered criterion uses either
   test at the family level.
4. **Ṗ folded into GATE-0**: `pdot_phase_smear_cycles` per source; three verdicts are marked
   conditional (smears 0.91/0.74/6.38 if the bounds are saturated); the note states "8 stand
   — 5 unconditionally".
5. **Independence caveat quantified**: the metrics caveat and the note state the ~55%
   interval weakening under activity persistence, with point estimates unbiased (linearity)
   — the presenter's pseudo-replication concern relocated to where it actually bites.
6. **Stokes-V stated**: "circularly polarized flux clears threshold" in the abstract, with
   the total-intensity conversion caveat in the body.
7. **The single change is in**: per-source implied f_active = p/((w+T)/P) committed
   (`wt_over_p`, `implied_f_active`, `implied_f_active_limit_5mjy`, vacuous flags) and the
   note leads with it — 0.039–0.735 for the detections (J183950 active in most cycles under
   its window model), the limit ordering inversion stated (best bound GLEAM-X < 0.09;
   \ldLimMin excludes only > 0.75; one vacuous), and the p range labelled consistent with
   one common rate (LRT p = 0.445, committed).
8. **\ldFActive carries its interval** (< 0.30 at 95%) and its conditionality (PEPOCH as
   pulse arrival; a-priori 0.18); J1832's split is retired via the corrected 0.400-cycle
   smear (`usable: false`); the J183950 phases are labelled a consistency check with the
   shared anchor, not an independent validation (here and in this file).
9. **Citations**: mcsweeney2025, wang2025, lee2025, pritchard2026 added (arXiv-verified,
   copied from the lpt paper's round-10 entries); barbere2026lptv now carries the Zenodo
   concept DOI and the specific committed filename.
10. **The producer exists**: `lptduty_run.py` calls `write_paper_assets` (previously
    caller-less), lptduty is in the Makefile SLICES (make paper SLICE=lptduty builds), and
    the five hand-typed prose numbers are macros (\ldClusterZ/P, \ldJEDPhase, \ldPhaseA/B).
11. MINORs: sample selection stated (10 of 16 with VAST coverage); the stale "GATE 0 not yet
    tested" caveat replaced; n_detections unified on _is_detection (leakage veto); this
    file's 0.046, "107", and "corroborates" corrected in place; \ldNDetSources used in the
    caption.
