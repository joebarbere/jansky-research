# 76 — EDGES raw-data averaging/flagging robustness check (calibration held fixed)

Status: 📋 planned (not started, **parked — high effort; pick up only if keen**) — GATE 0
pending: full-text novelty pass + raw-data-URL verification (the fable-ideas scan ran
egress-blocked; see the standing caveat there)

## Context

The EDGES 78-MHz absorption feature (Bowman et al. 2018) has been contested chiefly on
calibration and foreground-model grounds; SARAS 3 disputes it observationally. A narrower,
rarely-asked question is open (fable-ideas F39): holding the EDGES team's *own calibration
fixed*, how robust is the averaged spectrum to day-selection and RFI-flagging choices alone?
This is emphatically NOT a recalibration — no receiver modelling, no new foreground physics; it
is a data-selection sensitivity analysis on the public raw/intermediate products, using their
released calibration solutions verbatim. Recover-a-known is built in: with their selections,
their averaged spectrum must be reproduced first. High effort for a note-scale product, hence
parked unless keen.

## Deliverables

- `src/jansky_research/edgescheck.py`: `fetch_edges_raw` (public release, network, `# pragma`),
  `apply_released_calibration` (their solutions applied verbatim — no refitting),
  `average_spectrum` (day-selection + flagging-mask parameterised), `selection_sweep` (jackknife
  over days, flagging-threshold sweep, season splits; feature amplitude/width per variant),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/edgescheck/`; `survey/edgescheck-findings.md`; wiring.

## Approach

0. **GATE 0:** verify the EDGES raw/intermediate data and calibration solutions are publicly
   downloadable at file level and sufficient to reproduce the averaged spectrum without
   re-deriving calibration; full-text pass on the reanalysis literature (Hills+, Sims & Pober,
   SARAS 3) to confirm the selection-robustness axis is genuinely unclaimed. Either failure →
   stays parked.
1. Tooling + synthetic recover-a-known: a mock day-set with an injected absorption feature and
   known day-to-day systematics round-trips through the averaging/flagging machinery.
2. Anchor leg: reproduce their published averaged spectrum with their day selection and flags —
   the mandatory anchor before any variation.
3. Sweep leg: jackknife day-selection, flagging-threshold, and season-split variants; report the
   feature's amplitude/width distribution across variants, whatever it shows.
4. GATE-2 (scope discipline: any drift into calibration or foreground-model claims is out of
   bounds; multiplicity across variants) → paper (note-scale robustness report).

## Verification

Their averaged spectrum reproduced from their selections first — no variant is run before that
anchor passes; the variant distribution is reported plainly, stable or not; checks green; GATE-2
sign-off.

## Risks & mitigations

- **NOT a recalibration — say it loudly** → calibration solutions applied verbatim and frozen;
  the paper's scope statement leads with this; GATE-2 enforces it.
- **High effort, note-scale ceiling** → parked by default; the GATE-0 access check is cheap and
  runs first so effort is only spent if the data actually support the scope.
- **A stable result is unexciting; an unstable one invites overreach** → both outcomes are
  pre-framed as selection-sensitivity measurements, nothing more.
