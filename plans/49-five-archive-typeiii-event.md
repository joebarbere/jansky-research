# 49 — Five-archive simultaneous type III event at cycle-25 max (completes `type3synthesis`)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm all five
archives (e-Callisto, OVRO-LWA, Wind/WAVES, PSP/RFS, SolO/RPW) are simultaneously live

## Context

Exactly one five-instrument simultaneous type III analysis exists (the 2019-04-15 event:
PSP+STEREO-A+Wind+e-Callisto+EOVSA, arXiv:2306.01910) — pre-OVRO-LWA and well off solar maximum.
All five modern archives are now simultaneously public: the e-Callisto qkl tree, OVRO-LWA
spectrograms, Wind/WAVES, PSP/RFS, and SolO/RPW, and cycle-25 maximum supplies bursts daily.
This slice completes the merged `type3synthesis` arc: find **one** burst crossing all five bands
(a time-boxed hunt, 3 weeks of search effort, pre-committed), then run the existing
drift-to-distance ladder end-to-end on it. Deliverable: mutual consistency of five independent
distance/speed estimates on one real event — or an honest bounded null ("no five-archive event
found in N searched days"). Near-total reuse of four merged slices (`type3synthesis`,
`solarbursts`, `ecallisto_census`, the Wind/WAVES reader); the only new code is an OVRO-LWA
spectrogram reader.

## Deliverables

- `src/jansky_research/type3event.py`: `read_ovro_lwa_spectrogram` (the one new reader),
  `fetch_archives` (five-archive windowed download, `# pragma`), `find_crossband_candidates`
  (e-Callisto/OVRO-LWA burst lists × spacecraft-band activity windows → coincidence shortlist),
  `run_distance_ladder` (existing drift-to-distance machinery applied per archive),
  `consistency_table` (five estimates + errors), `synthetic_crossband_burst` (one injected
  drifting burst rendered into all five formats → ladder round-trip),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/type3event/`; `survey/type3event-findings.md`; wiring.

## Approach

0. **GATE 0:** verify live access + current formats for all five archives (e-Callisto qkl tree,
   OVRO-LWA spectrograms, Wind/WAVES, PSP/RFS CDFs, SolO/RPW CDFs); full-text pass on
   arXiv:2306.01910 and its citing papers to confirm no cycle-25 multi-archive event study has
   landed; **pre-commit the 3-week search time box in writing here before starting**.
1. Tooling + synthetic recover-a-known: one synthetic drifting burst injected into all five
   formats; the ladder must return five mutually consistent distance/speed estimates.
2. The hunt (time-boxed): shortlist candidate days from e-Callisto/OVRO-LWA burst activity ×
   PSP/SolO data-availability windows; confirm cross-band crossing on the shortlist.
3. Real leg: run the full ladder on the found event; the five-way consistency table. If the box
   expires empty: the bounded null (days searched, coverage fractions) is the paper.
4. GATE-2 science review: per-archive frequency-to-distance model assumptions stated uniformly,
   spacecraft-position light-travel corrections, honest framing of the null branch.
5. Paper: one-event five-archive consistency test (or the bounded null).

## Verification

The recover-a-known is reuse of the existing drift-to-distance ladder on the found event — each
archive's estimate must land within its stated errors of the synthetic round-trip behaviour;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **The hunt can fail** → the 3-week time box is pre-committed above; the bounded null (search
  coverage + non-detection) is a stated acceptable deliverable, not a failure.
- **Archive-format drift** (five independent providers) → GATE-0 pins formats; readers tested on
  vendored samples per archive.
