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
| ASKAP J183950.5−075635 | 2 | 0.046 | 44 | 24.3 mJy | 23222 s |

**Sources without** (95% upper limits at an assumed 5 mJy pulse): 0.029 to 0.107 across the
seven, i.e. under ~3–11% of snapshots.

For scale, a ~730 s snapshot on a 2656 s period samples ~27% of a cycle, so a source that
were *always* active with a broad pulse would give p of order 0.3. The measured 0.011–0.046
are an order of magnitude below that — consistent with the switching-off behaviour Rose et
al. describe, though the product cannot yet say how much is duty cycle and how much is
inactivity.

## Two decisions that determine whether the numbers mean anything

**The denominator is efficiency-weighted exposure, not epoch count** — the `frblens` lesson,
whose limit was 4x too tight until per-source efficiency entered it. Each epoch contributes
`Phi(S/sigma - 5)`, so epochs too shallow to have seen the assumed pulse contribute ~0. Note
the effective epochs above (91, 71, 44) against raw counts of 104, 97 and 92: **roughly half
the exposure is not real at these flux levels**, and an epoch-count denominator would have
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
3. The 107 never-released epochs carry no measurement and are excluded; that is recorded in
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
  phase 0 corroborates the assumption that the published PEPOCH is a pulse epoch.
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
