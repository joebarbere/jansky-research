# Plan 101 — `lvdrift`: the within-burst type III drift law from LV614 LOFAR dynamic spectra

Status: **plan, controls frozen 2026-10-18 before the burst arrays are analysed.** During GATE-0 a
research agent opened 2–3 example arrays to check shape and units only (recorded below); no drift
has been measured from any array. Changes after the first real run go in
`survey/lvdrift-findings.md` with the reason. Idea: `ideas.md` L13.

## Context

Šteinbergs, Kolotkov & Nakariakov (arXiv:2609.26535, HTML v1 read 2026-10-18; accepted ApJS)
catalogue 335 solar radio bursts seen by the Irbene LV614 LOFAR station at 29–60 MHz and fit the
type III frequency-drift law df/dt = A f^a with **A = 0.78 ± 0.91, a = 0.43 ± 0.31** (Fig. 15
caption). Their method (Sect. III.5): **one drift per burst**, the instantaneous rate at the time of
peak flux from a linear regression "in the vicinity of the peak" (window not stated); drifts binned
in 1 MHz bins of peak frequency; least squares "with mean flux … as fitting weights". So the
published law is an **across-burst** relation, and its amplitude is consistent with zero.

The classic laws are steeper and within-burst or population averages over wide bands:
Alvarez & Haddock 1973 (SoPh 29, 197) quoted as df/dt = −0.01 f^1.84 over 75 kHz–550 MHz (as cited
by Šteinbergs et al.); **Zhang et al. 2018** (A&A 618, A165; abstract, checked 2026-10-18):
df/dt = −0.0672 f^1.23 from 1389 type III bursts at **10–80 MHz** (Nançay Decameter Array) —
the same band as LV614.

**Gap.** Nobody has re-analysed the LV614 data (ADS: 0 citations, 2026-10-18). The across-burst
fit cannot constrain the index over a factor-of-two band. The arrays support a **within-burst**
measurement: track each burst's arrival time across its 162 channels and fit the drift law
along the burst. That is a different and better-posed quantity, directly comparable with Zhang
et al. 2018.

## Data (GATE-0, verified)

Zenodo doi:10.5281/zenodo.20157853, one zip, 266,723,235 B, CC-BY-4.0 (MD5 matched by the agent).
`type_{1..5}/npy`: 23, 8, **290**, 6, 5 arrays (the paper lists 293 type III). Arrays are float64,
(time, 162), in SFU, 1 s × 195.3 kHz, subbands 150–311 → 29.3–60.7 MHz; names encode date, start
and stop time and subband range. Per-burst drifts are not in the files; they are in the paper's
Table 3 (drift given for 229 of 293 type III and 7 of 8 type II), joinable on date + start time.

## The within-burst estimator (frozen)

`lvdrift.burst_index(arr, freqs, dt=1.0)`:
1. Background: per-channel median over the array, subtracted; noise σ_ch from the per-channel
   MAD (×1.4826).
2. Smoothing: 3 s × 3 channel boxcar.
3. Ridge: per channel, the time of the smoothed peak, kept if peak/σ_ch ≥ 5 and the peak lies
   ≥ 2 samples from the array edges.
4. Quality: ≥ 40 kept channels spanning a frequency ratio ≥ 1.5, and the kept ridge must be
   monotonic in time to within 2 s (type III bursts drift from high to low frequency).
5. Fit: t_peak(f) = t₀ + C f^(1−a) (the integral of df/dt = −A f^a, a ≠ 1) by non-linear least
   squares with a uniform 1 s timing error; report a and A. Bursts failing step 4 are counted
   and reported, never dropped silently.
6. Population: median of per-burst a with a bootstrap 95% interval (resampling bursts), and a
   joint fit with per-burst t₀, C and one shared a.

## Frozen controls

False-failure rate of every pass rule is simulated in step 1 before the rule is used (CLAUDE.md,
the `hiblend` lesson); if a rule's false-failure rate exceeds 5% it is widened **before** the
real run, with the simulation recorded.

- **P1 — planted bursts (recover-a-known, synthetic).** Bursts with a ∈ {0, 0.43, 1.0, 1.23,
  1.84}, drift amplitudes spanning the paper's type III range, Gaussian time profiles of the
  observed widths, injected at S/N 5–50 into **real noise**: the pre-burst first 20 s of real
  type III arrays (which the drift fit never uses). 50 plants per index. Pass: median recovered
  a within ±0.15 of the planted value for every index.
- **R0 — the paper's numbers (recover-a-known, real).** Their per-burst drift at peak flux,
  reproduced with a linear regression of ridge frequency on time in a ±3 s window around the
  peak: median ratio of ours to Table 3 for the joined bursts in [0.8, 1.25]; type II and III
  means reported against 0.043 and 4.383 MHz s⁻¹. Their across-burst fit is re-run with their
  recipe (1 MHz bins, mean-flux weights) and with an unbinned errors-in-variables fit; both
  reported.
- **N1 — shuffled null.** Each channel of a real burst shifted by an independent random time
  offset (uniform over the array length): the quality cut (step 4) must reject ≥ 95% of shuffled
  bursts.
- **Robustness.** Ridge S/N threshold 4 / 5 / 7; smoothing 1 / 3 / 5; the peak-regression window
  ±2 / ±3 / ±5 s for R0. The headline is the 5 / 3 / ±3 setting; the others are reported.

Interpretation fixed now: the within-burst median a is reported with its interval and compared
with 1.23 (Zhang et al. 2018, same band) and 1.84 (Alvarez & Haddock as cited). The across-burst
index is reported as a separate quantity. If the within-burst index differs from 1.23 by more than
its interval, the note says so and lists the differences in method (peak-time ridge vs Zhang's
estimator, instrument, sample) without choosing between them.

## Deliverables

`src/jansky_research/lvdrift.py` (pure estimator, synthetic plants, join to Table 3, `# pragma`
for download), tests to the 85% floor, `results/lvdrift_*.json` with `source`,
`survey/lvdrift-findings.md`, `papers/lvdrift/` RNAAS note (software-cites `jansky-research`).

## Order of work

1. Estimator + synthetic tests + false-failure simulation of every rule (no real arrays).
2. P1 on real noise windows. 3. N1, R0, then the within-burst run. 4. Findings, GATE-2 review,
note.

## Risks

- **Timing resolution.** At 4 MHz s⁻¹, 1 s spans ~20 channels; slow bursts are well sampled,
  fast ones poorly. P1 measures where the estimator stops working; the result is stated for the
  bursts where it works.
- **Selection.** The sample is NOAA-list driven and the brightest of grouped bursts (Sect. II.2);
  stated as such.
- **Table 3 parsing.** Read from the paper's HTML table; joined counts reported.
