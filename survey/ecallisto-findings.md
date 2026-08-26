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

## Full referee round (2026-08-26) — pipeline paper: MAJOR REVISION, 15 findings, three BLOCKERs

The framing is right, the macro namespacing fix has held (offline rerun cannot clobber the real
JSON or macros — verified live), and the synthetic validation is MORE robust than claimed (the
referee ran 60 seeds: 60/60 recover, 0 of 180 pure-noise spectra false-flagged). But the real
leg is not measuring what its sentences say.

**BLOCKER 1: the committed figure is the real leg's output under the synthetic leg's caption.**
figures/ecallisto.pdf extracts as "Coincidence QC: 0 confirmed / OOTY … BIR" while the caption
describes "\ecSynMaxEventStations = 4 detections at a common time (red, dashed) confirmed as one
event". `_figure` writes one path from both legs with no guard, and `make reproduce` runs the
real leg last. The paper's only figure contradicts its caption and abstract.

**BLOCKER 2: \ecRealNscanned = 1512 double-counts every file (the index regex matches each
filename twice: 763 distinct listed → 1526 returned) and only 559 distinct spectra were actually
fetched — the count is 2.7× the spectra analysed, and "station-days" is wrong by a further ~96×
(files are 15-minute spectra). burst_fraction inherits it.**

**BLOCKER 3: 64% of scanned rows analyse a DIFFERENT file than the one they are labelled with.**
ingest_day lists a filename, extracts HHMM, then fetch_ecallisto re-resolves by
closest-preceding-start — files starting at HHMM59 resolve to the PREVIOUS file (15 min early),
and focus-code siblings collapse onto the first match. Referee simulation against the live
index: 970/1516 rows (64.0%) fetch a different file; the mechanism reproduces the committed 1512
to 0.3%. The UT peak time is computed from the LABEL's start second while the data came from
another file — so up to 64% of t_peak_s values fed to the 60 s coincidence are wrong by up to
~15 minutes. The same defect is in the production DAG.

**MAJORs:** the committed CSV omits t_peak_s and file — the two fields coincident_events
consumes (the innerrc lesson), which is also how the referee could show the "8 candidates" are
4+4 byte-identical duplicate rows = TWO real detections (\ecRealNbursts wrong by 4×); no
attempted denominator + the bare-except silent drop is undisclosed here while the census sibling
now discloses it (14 of 1526 entries silently dropped on this very day); "the DAG and make
ecallisto-day produce identical rows" is false (the DAG reimplements the worker inline with
opposite error semantics, NaN→None, and MAX_FILES=12 vs the uncapped run); the synthetic
validation cannot fail on the axis it claims (RFI and bursts are the SAME synthetic_burst
function differing in seed; nearest injected gap to the 60 s boundary is 200 s); single-linkage
clustering has no span cap (8 stations at 50 s spacing — true span 350 s — confirm as ONE event:
"a 60-second tolerance" overstates); no false-positive rate anywhere (two chance-coincident RFI
stations promote 20/20 seeds at dt ≤ 59 s — the QC's reliability IS the chance-coincidence rate,
computable in closed form, unquoted); the real-data null (the QC has NEVER confirmed a real
event) is buried in the Discussion and the findings-file explanation is stale — BIR is present
in the committed run with a candidate (13 of 15 stations produced NO candidate: the null
measures the single-station threshold, not station coverage).

**MINOR/NIT:** single-seed "recovers exactly" should quote the referee's ensemble (60/60,
0/180); the figure is the one artifact with no merge guard (verified live); the census sibling's
120 s tolerance is unremarked here; the same three spectra are "injected bursts" in §3 and
"interference" in §4; "60-second" and "two distinct stations" are hand-typed against defaults.
Citations and macros all verify clean.

**Verdict: MAJOR REVISION.** The single change: make ingest_day scan the file it listed
(de-duplicate the listing, fetch by filename), then re-run 2011-09-14 and re-derive every
\ecReal* macro — retiring all three blockers and turning the buried null into the paper's best
sentence: of 15 stations across a full day, only two produced a type III candidate, and they
were not coincident.

**Status: RESOLVED (2026-08-26) — pipeline paper.** ingest fixed at the root and the real day
re-run end to end; the referee's predictions confirmed by the committed pipeline: 763 files
listed (their 763), 761 scanned (2 fetch failures now counted, never silent), **2 candidates**
(the phantom 8 were duplicate rows), BIR's at t_peak = 42666 s = 11:51:06 UT — consistent with
the known 11:50 UT type III (the referee's "plausible" identification, now auditable from the
committed t_peak_s) — and OOTY's 3.2 h earlier; 0 coincident events.

1. **BLOCKER 1 (figure/caption)**: each leg writes its own figure (ecallisto_syn.pdf /
   ecallisto_real.pdf); the old shared path is deleted; the synthetic caption describes the
   synthetic panel again and the real day gets its own figure with the honest reading.
2. **BLOCKERS 2+3 (double-count + wrong-file analysis)**: `list_day_files` de-duplicates (the
   index HTML renders each filename twice); `scan_file` fetches EXACTLY the listed file
   (`fetch_ecallisto` gained a filename= path; the HHMM re-resolution and its 64% mismatch are
   documented in its docstring as the lesson). n_scanned now counts spectra; "station-days"
   language corrected.
3. **The committed catalogue can audit the headline**: t_peak_s and file columns added — the
   two fields coincident_events consumes.
4. **The denominator exists**: n_files_listed / n_stations_listed / n_fetch_failed / station
   scope committed and quoted; the bare-except is counted, not silent.
5. **DAG = CLI**: both call the shared `scan_file`; the abstract's identical-rows claim is now
   scoped and true (the DAG stamps the date, JSON-sanitises NaNs, and caps its fan-out — all
   stated).
6. **The validation can fail**: the fixture gained a constant narrowband carrier, a broadband
   zero-drift impulse (30 s from the carrier — INSIDE the tolerance, so a morphology failure
   would be promoted), and a reverse-drift ridge; all three rejected on morphology. The
   single-seed sentence is replaced by the committed 60-seed ensemble (60/60 exact recovery,
   0 quiet false flags, 0 contaminant false flags, drift sd 0.015).
7. **Single-linkage honesty**: each confirmed event carries its span (a chained cluster can
   exceed the tolerance — stated); `chance_coincidence_rate` commits the closed-form expected
   chance events per day (0.0292 synthetic; 0.0014 for the real day's counts) — the QC's
   actual reliability number, previously never quoted; the census sibling's 120 s tolerance is
   cross-referenced in print.
8. **The real-data null is the abstract's second sentence**, with the corrected explanation:
   13 of 15 stations produced no candidate on a day with a known burst — the bottleneck is the
   single-station detection threshold, not station coverage (this file's earlier six-station /
   BIR-excluded explanation is superseded by the committed evidence).
9. Hand-typed 60 s / 2 stations → \ecTolS/\ecMinStations macros; per-leg \ecSynSource /
   \ecRealSource markers (fixing a downgrade-guard see-saw that had blocked the synthetic leg
   from updating its own namespace); "injected bursts" vs "interference" vocabulary unified.
