# 77 — LNA-first vs filter-first: measured noise-figure/overload A/B for the H-line front end

Status: 📋 planned (hardware-gated) — needs the H-line front-end parts (LNA + 1420 MHz filter) and
bench instruments (LiteVNA 64, tinySA Ultra+, HackRF) delivered and health-checked; bench-only, no
sky required — scheduled during station commissioning

## Context

Every amateur H-line build guide repeats the folk trade-off: LNA-first minimizes system noise
figure (Friis), filter-first protects the LNA from out-of-band overload by urban transmitters —
but the station-science sweep found nobody has *published measured numbers* for both orderings on
the same parts. The amateur literature stops at qualitative advice; rigor (a real measured A/B
with stated uncertainties) is the novelty axis, matching this repo's house style. The station's
own design doc (`station/hydrogen-line-receiver.md`) already flags "filter order is itself an
experiment" and both configurations will be built; `station/test-equipment.md` inventories the
bench. This is the first station slice: it needs no sky, only the commissioning bench.

## Deliverables

- `src/jansky_research/lnabench.py`: `cascade_noise_figure` (Friis prediction for either
  ordering), `y_factor_nf` (hot/cold Y-factor reduction with error propagation),
  `overload_margin` (1 dB compression / desense vs injected interferer power),
  `synthetic_bench` fixture (known NF/gain stages → both reductions recover them),
  `run`/`_figure`/`_write_macros`/`_main`.
- Tests to the 85% floor (synthetic bench fixtures — no hardware needed for tests);
  `papers/lnabench/`; `survey/lnabench-findings.md`; wiring (Makefile/Snakefile).

## Approach

0. GATE 0: all front-end parts + bench gear on hand and triaged per
   `station/test-equipment.md` (VNA S21 of the filter, LNA gain check, HackRF health sequence).
1. Tooling + synthetic recover-a-known: build both reductions offline; a simulated two-stage
   chain with known NF/gain must round-trip through `y_factor_nf` and `cascade_noise_figure`.
2. Bench leg (`# pragma: no cover`): S21 sweeps of both orderings; Y-factor NF per ordering
   (ambient vs cold reference); desense curves with a HackRF out-of-band tone stepped in power
   at representative urban interferer frequencies (cell/FM bands).
3. Compare measured NF deltas to the Friis prediction; quantify the overload margin each
   ordering buys; state which configuration the Philadelphia rooftop should fly and why.
4. GATE-2 science review → RNAAS-style methods note in `papers/lnabench/`.

## Verification

Synthetic bench round-trips recover injected NF/gain/compression to stated tolerance; measured
Friis prediction vs Y-factor agreement is the real-data anchor; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Hardware delivery slips** → bench-only slice; all software + synthetic legs ship first, the
  bench leg lands whenever parts do.
- **Hobby-grade instrument accuracy** (tinySA/LiteVNA absolute cal) → report deltas between
  orderings, not absolute NF claims; propagate instrument spec uncertainties into every number.
- **Urban RFI contaminating the Y-factor cold reference** → use a 50 Ω terminated reference on
  the bench, not sky; the overload test injects a controlled tone, so ambient RFI is excluded.
