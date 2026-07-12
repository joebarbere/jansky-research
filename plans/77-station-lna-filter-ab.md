# 77 — LNA-first vs filter-first: measured noise-figure/overload A/B for the H-line front end

Status: 📋 planned (hardware-gated) — needs a **dedicated discrete parts buy** (~$85: SPF5189Z
wideband LNA board ~$5, standalone 1420 MHz bandpass filter ~$35, Nooelec SAWbird+ H1 ~$45) plus
bench instruments (LiteVNA 64, tinySA Ultra+, HackRF, BG7TBL noise source) delivered and
health-checked; bench-only except the optional on-sky feed leg — scheduled during station
commissioning. Updated 2026-07-11 for the integrated-feed architecture (see Context).

## Context

Every amateur H-line build guide repeats the folk trade-off: LNA-first minimizes system noise
figure (Friis), filter-first protects the LNA from out-of-band overload by urban transmitters —
but the station-science sweep found nobody has *published measured numbers* for both orderings on
the same parts. The amateur literature stops at qualitative advice; rigor (a real measured A/B
with stated uncertainties) is the novelty axis, matching this repo's house style.

**2026-07-11 architecture update.** The station's front end is no longer a discrete LNA + filter
chain: it is the sealed **KrakenRF Discovery H-line active feed** (2× QPL9547 LNA around 2× SAW —
an LNA-first design whose first stage sees the unfiltered sky, relying on the QPL9547's high IP3;
see `station/hydrogen-line-receiver.md`, which notes the integrated feed *fixes* the ordering a
discrete build has to choose). Two consequences for this slice:

1. **The A/B no longer chooses the station's configuration** — that decision is baked into the
   feed. The framing shifts to a pure methods contribution: publish the measured trade nobody
   has, on cheap discrete parts anyone can buy, and anchor it against the commercial integrated
   units people actually fly.
2. **The A/B needs its own small parts buy** (the build dropped the discrete LNA/filter): the
   SPF5189Z board and a standalone 1420 MHz bandpass filter, plus the SAWbird+ H1 as the
   bench-measurable integrated reference. All three are inventoried as optional/spare items in
   the build note ([[hydrogen_line_telescope_build]] § Research program); this slice makes them
   load-bearing.

Measurement arms:

- **A — LNA-first (discrete):** SPF5189Z → 1420 MHz filter.
- **B — filter-first (discrete):** 1420 MHz filter → SPF5189Z.
- **C — integrated commercial reference (bench):** Nooelec SAWbird+ H1 — SMA in/out, so it runs
  the identical cabled bench sequence; shows where a popular integrated unit lands on the same
  axes as A/B.
- **D — the station's own feed (on-sky, system-level; optional but valuable):** the Discovery
  feed's input is the antenna element inside the radome — it **cannot be bench-injected**. Its
  characterization is a ground/cold-sky Y-factor system-temperature measurement through the full
  chain during commissioning, reported alongside (not as) the controlled A/B. No published
  system-Tsys measurement of this feed exists either — it dates the note to the Discovery Dish
  era and grounds the A/B in flown hardware.

## Deliverables

- `src/jansky_research/lnabench.py`: `cascade_noise_figure` (Friis prediction for either
  ordering), `y_factor_nf` (hot/cold Y-factor reduction with error propagation),
  `overload_margin` (1 dB compression / desense vs injected interferer power),
  `sky_ground_tsys` (on-sky Y-factor reduction for arm D),
  `synthetic_bench` fixture (known NF/gain stages → all reductions recover them),
  `run`/`_figure`/`_write_macros`/`_main`.
- Tests to the 85% floor (synthetic bench fixtures — no hardware needed for tests);
  `papers/lnabench/`; `survey/lnabench-findings.md`; wiring (Makefile/Snakefile).

## Approach

0. GATE 0: A/B parts (SPF5189Z, standalone 1420 filter, SAWbird+ H1) + bench gear on hand and
   triaged per `station/test-equipment.md` (VNA S21 of the filter, LNA gain check, HackRF health
   sequence, noise-source sanity check).
1. Tooling + synthetic recover-a-known: build the reductions offline; a simulated two-stage
   chain with known NF/gain must round-trip through `y_factor_nf` and `cascade_noise_figure`;
   a simulated hot/cold sky pair must round-trip through `sky_ground_tsys`.
2. Bench leg (`# pragma: no cover`): S21 sweeps of arms A/B/C; Y-factor NF per arm (ambient vs
   50 Ω cold reference, BG7TBL stimulus — uncalibrated ENR, so deltas not absolutes); desense
   curves with a HackRF out-of-band tone stepped in power at representative urban interferer
   frequencies (cell/FM bands). Identical sequence for all three cabled arms.
3. On-sky leg (arm D, during commissioning): ground vs cold-sky Y-factor through the Discovery
   feed + full chain → system temperature with propagated uncertainty; compare against the
   feed's published/estimated NF and the A/B-derived expectations for an LNA-first architecture.
4. Compare measured NF deltas to the Friis prediction; quantify the overload margin each
   ordering buys; place the integrated units (C bench, D on-sky) on the same axes; state what
   the numbers imply for urban-rooftop builders choosing between discrete orderings and
   integrated feeds.
5. GATE-2 science review → RNAAS-style methods note in `papers/lnabench/`.

## Verification

Synthetic bench round-trips recover injected NF/gain/compression (and synthetic Tsys) to stated
tolerance; measured Friis prediction vs Y-factor agreement is the real-data anchor; checks
green; GATE-2 sign-off.

## Risks & mitigations

- **Hardware delivery slips** → bench-only slice at core; all software + synthetic legs ship
  first, the bench leg lands whenever parts do; arm D rides the commissioning schedule and is
  optional to the paper.
- **Hobby-grade instrument accuracy** (tinySA/LiteVNA absolute cal, BG7TBL uncalibrated ENR) →
  report deltas between orderings, not absolute NF claims; propagate instrument spec
  uncertainties into every number.
- **Urban RFI contaminating the Y-factor cold reference** → use a 50 Ω terminated reference on
  the bench, not sky; the overload test injects a controlled tone, so ambient RFI is excluded.
- **Discovery feed is not bench-injectable** (antenna-integrated input) → structural: arm D is
  system-level on-sky only and framed as an anchor measurement, never as part of the controlled
  A/B; the discrete arms carry the comparison.
- **Cold-sky leg weather/RFI sensitivity (arm D)** → repeat on ≥2 nights; report spread as the
  uncertainty, not a single-shot number.
