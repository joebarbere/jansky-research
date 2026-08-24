# Findings — e-Callisto streaming ingest + cross-station coincidence QC

`jansky_research.ecallisto_catalog` is the worker behind the e-Callisto Airflow-on-Podman ingest
pipeline (`plans/31`). It scans each station's 15-minute dynamic spectrum for a drifting type III ridge
(reusing the `solarbursts` tools) and — the step that makes the output trustworthy — vets the
per-station candidates by **cross-station coincidence**.

## Why coincidence is the key QC

e-Callisto is 150+ heterogeneous ground stations, uncalibrated and RFI-heavy, so *any single station*
produces spurious ridges. The physical discriminant: a real solar burst radiates to the whole sunlit
hemisphere and is recorded at the **same universal time** by many stations, while RFI/local artefacts are
single-station. `coincident_events` clusters the day's candidate rows in peak time (single-linkage,
60 s tolerance for clock offsets + burst duration) and confirms a cluster spanning ≥2 distinct stations.
For the real path, each candidate's peak time is converted to **UT-of-day** (file-start + local peak) so
stations whose 15-minute files begin at different UTs are compared on one clock.

## Recover-a-known (synthetic): coincidence confirms the real burst, rejects RFI

A synthetic day — a real type III injected at **4 stations** at a common UT, plus **3 single-station**
interference events at distinct times, plus 3 quiet stations:

| quantity | value |
|---|---|
| stations scanned | 10 |
| burst candidates (per-station) | 7 (4 real + 3 RFI) |
| **coincidence-confirmed events** | **1** (the real burst, 4 stations, drift −6.9 MHz/s) |
| single-station candidates rejected | 3 |

The coincidence step recovers exactly the injected event and rejects every single-station spurious
candidate — the QC works.

## Real-data reality check (2011-09-14): coincidence is coverage- and detection-limited

A quick real run — the 11:45 UT window of the `solarbursts` recover-a-known day (2011-09-14, whose
11:50 type III was cleanly fit at **BIR**), scanned across the six stations that happened to have an
11:45 file (BLEN7M, BLENSW, DARO, HUMAIN, INPE, MRO) — produced **no coincidence-confirmed event**: not
one of those six registered a clean type III ridge in that window (drifts present but `r2` below the
0.5 threshold, or the wrong sign). This is the honest limitation, not a bug: (i) BIR, the station that
caught the event, was not in that six-station subset; (ii) whether any *other* station saw a given burst
depends on sunlit coverage and pointing; and (iii) individual-station detection on uncalibrated,
RFI-heavy e-Callisto data is threshold-sensitive. So a real multi-station coincidence needs the *full*
active-station set for the window and event-tuned detection — which is exactly the census work the
coincidence QC enables, not something a six-station snapshot with the synthetic-tuned thresholds
delivers. The **coincidence logic** is validated on the synthetic day above; its **real yield** is a
coverage-limited lower bound.

## Honest assessment & caveats

- **Candidates → events, not yet a census.** The coincidence promotes candidates to confirmed events;
  a full multi-cycle *occurrence census* (with completeness corrections, station-coverage weighting, and
  calibration caveats) is the natural next step — the pipeline's per-day reduce scales for it.
- **Type III only.** The drift-based detector targets fast negative-drift ridges (type III); type II
  (slower, shock-driven) would need a second template.
- **Coincidence depends on station coverage.** A burst seen by only one active station cannot be
  confirmed; the confirmed rate is a lower bound set by how many stations observed the event.
- **Reproducible:** `make ecallisto-day DATE=...` runs a day's scan + coincidence without Airflow; the
  DAG's `reduce_day` writes both the per-day candidate CSV and the confirmed-events CSV idempotently.

## The paper's only result was drawn from the wrong run (found + fixed 2026-08-23)

Found by the referee round on the traditional-style conversion, not by the conversion itself.
This is the un-namespaced mode-dependent macro failure in its purest form.

`_write_macros` emitted seven names (`\ecNevents`, `\ecMaxEventStations`, `\ecNbursts`,
`\ecNscanned`, `\ecMedDrift`, `\ecBurstFrac`, `\ecNrfiRejected`) from **both** legs, and `run()`
wrote both legs' metrics to one `results/ecallisto_metrics.json`. The real 2011-09-14 archive day
ran last, so the committed macros held real values --- while **all twelve macro uses in the paper
describe the synthetic day**. The abstract typeset as:

> On a synthetic day (a real burst at **0** stations plus single-station interference), the
> coincidence step confirms exactly **0** event and rejects **8** single-station candidates.

and the Method attributed the real day's `8 of 1512` and `-5.941 MHz/s` to a synthetic run of ten
stations. `preserve_live_macros` could not help: both legs wrote real values under one name, so
there was nothing to arbitrate. Correct synthetic values are 4 stations, 1 event, 3 rejected,
7 of 10 flagged, -6.911 MHz/s.

**Fixed:** macros namespaced `\ecSyn*`/`\ecReal*`; each leg writes its own results file
(`ecallisto_synthetic_metrics.json` + `_catalog.csv`, allowlisted like `vgpra_synthetic_*`); the
paper repointed to the synthetic namespace.

Two things this taught that are not in CLAUDE.md yet:

1. **`preserve_live_macros` does not accumulate names, only values.** It rewrites the lines the new
   run emits and silently drops any existing macro the new text never mentions. The two-namespace
   design works only because each writer emits *both* namespaces, filling its own and leaving the
   other as `--`. A writer that emits only its own namespace **deletes** the other leg's numbers.
2. **The un-namespaced test was the lock.** `test_run_offline_writes_artifacts` asserted
   `\ecNbursts in macros` --- it passed precisely because the defect was present, and would have
   failed on the fix. Third instance of this pattern in the repo.

**Still open:** the real archive day is now cited once (it flagged `\ecRealNbursts` candidates among
`\ecRealNscanned` station-days and confirmed `\ecRealNevents`), but it remains a single day and does
not characterise the coincidence step's real-data performance.
