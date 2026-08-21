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

## Open — do not quote these numbers in a paper yet

1. **GATE 0 novelty** has not been run. A VAST or MWA transient paper may already report
   class-wide activity fractions.
2. **Aliasing is untested.** VAST pointings are not randomly phased against any LPT period.
   If the cadence beats against a period, the binomial model's independence assumption
   fails. This is the single most likely way the numbers above are wrong.
3. **Ephemeris audit** for the phase-resolved leg (`lptv` round 3 was bitten by multi-year
   phase drift; the same systematic applies here).
4. The 107 never-released epochs carry no measurement and are excluded; that is recorded in
   the loader, not silently dropped.
