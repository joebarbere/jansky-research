# 77 — LNA-first vs filter-first: measured noise-figure/overload A/B for the H-line front end

Status: 🧊 **deferred to backlog (2026-07-11)** — see `survey/candidate-gaps.md` § Station
hardware methods. Restructured after the integrated-feed architecture change: the station's
sealed Discovery feed fixes its own ordering, so no measurement outcome here changes any station
decision, and the slice requires a dedicated ~$85 discrete-parts buy purely for the paper. The
one station-relevant piece — the on-sky ground/cold-sky Tsys measurement of the Discovery feed
(formerly arm D) — **moved to plan 78's commissioning leg**, where it feeds the calibration
substrate and plan 79's error budget. Pick this plan back up only if the fast-RNAAS credit is
wanted; it remains the quickest citable result in the station track and the bench never expires.

## Context

Every amateur H-line build guide repeats the folk trade-off: LNA-first minimizes system noise
figure (Friis), filter-first protects the LNA from out-of-band overload by urban transmitters —
but the station-science sweep found nobody has *published measured numbers* for both orderings on
the same parts. The amateur literature stops at qualitative advice; rigor (a real measured A/B
with stated uncertainties) is the novelty axis, matching this repo's house style.

**2026-07-11 architecture update (and why this is now backlog).** The station's front end is the
sealed **KrakenRF Discovery H-line active feed** (2× QPL9547 LNA around 2× SAW — an LNA-first
design whose first stage sees the unfiltered sky; see `station/hydrogen-line-receiver.md`, which
notes the integrated feed *fixes* the ordering a discrete build has to choose). The A/B therefore
no longer informs the station: it survives only as a pure methods contribution on parts bought
for the purpose.

Measurement arms (if picked up):

- **A — LNA-first (discrete):** SPF5189Z → 1420 MHz filter.
- **B — filter-first (discrete):** 1420 MHz filter → SPF5189Z.
- **C — integrated commercial reference (bench):** Nooelec SAWbird+ H1 — SMA in/out, so it runs
  the identical cabled bench sequence; shows where a popular integrated unit lands on the same
  axes as A/B.

Parts: SPF5189Z board (~$5), standalone 1420 MHz bandpass filter (~$35), SAWbird+ H1 (~$45);
bench: LiteVNA 64, tinySA Ultra+, HackRF, BG7TBL noise source.

## Deliverables (if picked up)

- `src/jansky_research/lnabench.py`: `cascade_noise_figure` (Friis prediction for either
  ordering), `y_factor_nf` (hot/cold Y-factor reduction with error propagation),
  `overload_margin` (1 dB compression / desense vs injected interferer power),
  `synthetic_bench` fixture (known NF/gain stages → both reductions recover them),
  `run`/`_figure`/`_write_macros`/`_main`.
- Tests to the 85% floor (synthetic bench fixtures — no hardware needed for tests);
  `papers/lnabench/`; `survey/lnabench-findings.md`; wiring (Makefile/Snakefile).
- (The on-sky `sky_ground_tsys` reduction that used to live here ships with plan 78.)

## Approach (if picked up)

0. GATE 0: A/B parts (SPF5189Z, standalone 1420 filter, SAWbird+ H1) + bench gear on hand and
   triaged per `station/test-equipment.md` (VNA S21 of the filter, LNA gain check, HackRF health
   sequence, noise-source sanity check).
1. Tooling + synthetic recover-a-known: build the reductions offline; a simulated two-stage
   chain with known NF/gain must round-trip through `y_factor_nf` and `cascade_noise_figure`.
2. Bench leg (`# pragma: no cover`): S21 sweeps of arms A/B/C; Y-factor NF per arm (ambient vs
   50 Ω cold reference, BG7TBL stimulus — uncalibrated ENR, so deltas not absolutes); desense
   curves with a HackRF out-of-band tone stepped in power at representative urban interferer
   frequencies (cell/FM bands). Identical sequence for all three cabled arms.
3. Compare measured NF deltas to the Friis prediction; quantify the overload margin each
   ordering buys; place the integrated unit on the same axes; state what the numbers imply for
   urban-rooftop builders choosing between discrete orderings and integrated feeds. Cross-link
   plan 78's measured Discovery-feed Tsys as the flown-hardware anchor.
4. GATE-2 science review → RNAAS-style methods note in `papers/lnabench/`.

## Verification

Synthetic bench round-trips recover injected NF/gain/compression to stated tolerance; measured
Friis prediction vs Y-factor agreement is the real-data anchor; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Hobby-grade instrument accuracy** (tinySA/LiteVNA absolute cal, BG7TBL uncalibrated ENR) →
  report deltas between orderings, not absolute NF claims; propagate instrument spec
  uncertainties into every number.
- **Urban RFI contaminating the Y-factor cold reference** → use a 50 Ω terminated reference on
  the bench, not sky; the overload test injects a controlled tone, so ambient RFI is excluded.
- **Staleness while deferred** → none that matters: the folk question is decades old and the
  parts are commodity; the gap only closes if someone else publishes the measurement first
  (re-run the novelty check before picking up).
