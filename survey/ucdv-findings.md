# ucdv findings — RACS Stokes-V two-epoch census of the Gaia UCD sample (plan 59)

Running log, newest section last. Plan: `plans/59-racs-v-ucd-census.md` (controls frozen 2026-10-15,
PR #344). Every number below is in a committed `results/ucdv_*.json` unless marked otherwise.

## GATE 0 (2026-10-15, scout, PR #343)

Novelty: Perger, Seli & Vida (arXiv:2506.21169, full text) is VLASS Stokes I only; RACS-low2
(arXiv:2606.16182, full text) measures V at Stokes-I catalogue positions and names no UCDs. ADS
keyword queries found no targeted forced-V UCD census. Data: Reylé et al. (2018) via VizieR
`J/A+A/619/L8` — **14,915 rows** confirmed on download (Gaia DR2 ids, positions at J2015.5).

## Pre-data amendments to the frozen controls (both before any real pixel was evaluated)

1. **C1 null rule.** The first frozen version (95% sign interval, one null per target-epoch,
   n ≈ 80) fails a *perfect* pipeline 6.6% of the time (17% at n = 40; 20,000 simulations). The
   synthetic offline run tripped it on pure noise (std 0.83). Now 5 nulls per target-epoch
   (n ≥ 400) and a 99% interval: 0.8% false-fail. Cutouts enlarged to 0.06° so a 2′ null plus its
   90″ annulus, and the 3′ bright-source guard, fit. (CLAUDE.md: the `hiblend` gate lesson.)
2. **Provenance split.** The first version split set A/B by Gaia `source_id` parity. **All
   14,915 Reylé ids are even**, so every target would have been in calibration set A and nothing
   could certify. Now the first bit of SHA-256 of the decimal id: 50.2% B overall.

Both were caught by building the machinery against the plan before opening data — the order the
slice pattern asks for.

## Footprint (CASDA ObsCore metadata only)

RACS V products: 1,177 low2 (λ 0.338 m, MJD 59600–59899) and 1,702 low3 (λ 0.318 m, MJD
60200–60499). RACS-low1 has none. **12,751** UCDs lie in both a low2 and a low3 V footprint
(low2 stops at Dec +48°); 6,374 in set B. Median distance 87 pc; 20 within 20 pc, 1,181 within
50 pc. (A first footprint test passed the antipodal region too — every target "covered" — caught
by the count and fixed by orienting each edge from the polygon centroid.)

## S0 — throughput and scope (results/ucdv_s0_throughput.json)

40 random set-B targets × 2 epochs with a single SODA job per I+V pair (astroquery's default is
a job per product, polling every 20 s; here 3 s). **80/80 pairs ok, 0 failures, 22.0 s per
pair.** The full both-epoch sample projects to 156 h > 120 h, so the frozen rule sets the census
to **d < 107.8 pc: 9,805 targets, 119.99 h**, chosen from parallax and timing only. Option (b),
batched multi-position jobs, was not tried: (a) already reaches 77% of the sample, and the 2,946
dropped are the most distant, whose limits are weakest. Recorded as a deviation from the plan's
"(a) and (b)".

## C1, C2 on real cutouts, set B (results/ucdv_controls_c1c2.json)

First evaluation of any pixel value; target V itself not evaluated.

- **C1 null, V (400 positions): pass.** 51.0% positive (99% interval 43.5–56.5%), std(V/σ) 0.96,
  0 at ≥ 5σ. The same positions in Stokes I: 51.0% positive, std 1.01 — I is genuinely forced too.
- **C2 injection (480 per level, target pixels + null sites): pass.** Completeness 2.9 / 47.7 /
  96.5 / 100 / 100% at 3 / 5 / 7 / 10 / 20 σ; median recovered/injected 0.982 at ≥ 10σ.
  Target pixels alone (80 per level): 95.0% at 7σ, ratio 0.972. The 2–3% low bias is the
  nearest-pixel read of a beam centred off the pixel centre; inside the frozen 0.95–1.05.
- **Image noise:** median V σ 0.149 mJy/beam (0.116–0.269) over the 80 cutouts, so 3σ limits
  near 0.45 mJy — deeper than the 0.25 mJy assumed in early estimates.

Deviation: C2 asked for 200 injections per level "at the target pixel"; 80 target cutouts give 80,
so null sites in the same cutouts were added to reach 480. Both are reported.

## C3 recover-a-known (results/ucdv_controls_c3.json)

- **C3(a): pass.** All 67 `source_type == 'star'` rows of the RACS-low2 V catalogue (CASDA
  `AS110.racs_low2_v_v01`, 221 rows) fetched in their own SBID (stars flare, so the most-central-
  tile rule is overridden by `obs_id`), forced V at the catalogue's fitted position: median
  forced/catalogue **1.021** (IQR 1.00–1.09; frozen range 0.85–1.15), **0 handedness flips in 67**,
  local σ_V / catalogue noise 0.98. Largest outlier: RACS-LOW2 J125545.3−445152, 0.60 → 0.81 mJy
  with σ 0.12 (1.7σ apart).
- **C3(b) WISE J0623−0456: not run yet.** It is a T8 dwarf, too faint for Gaia and so not in the
  Reylé sample; the plan makes it "reported, not pass/fail" and requires the full text of
  arXiv:2306.15219 to be read before any flux from it is quoted. Pending, and not a gate for the
  census.

All three frozen controls pass on real data before the census starts.

## Census launch (2026-10-15)

50-slot trial first (part of the census, kept): 50/50 ok, 21.0 s per slot, the 250 per-slot nulls
52% positive with std(V/σ) 1.10 and none at ≥ 5σ, median σ_V 0.159 mJy/beam, ~11 MB of float32
cutouts per 50 slots (≈ 4.4 GB for the census). Target V not evaluated. The full run then started
as `ucdv-census.service` (systemd user unit, resumable) from a detached worktree pinned to
`38426ad`, 9,805 targets × 2 epochs = 19,610 slots, ~115–120 h expected.

## Census stall and restart (2026-10-16)

The census ran at ~165 slots/h until 07:54 EDT on Oct 16, then every CASDA login failed for
~11.5 h; the runner's fixed 60 s retry made ~285 failed OPAL logins before it was noticed and
stopped at 19:28 (1,694 of 19,610 slots done, all `ok`). A manual login then succeeded — no lockout.
Retries now back off exponentially (60 s doubling, cap 30 min; ~29 logins for an outage that
long) and log timestamped errors. Failed slots were never written, so no data was affected.
After the restart (19:29) logins succeeded but no slot was written: CASDA's SODA queue was leaving
cutout jobs PENDING indefinitely (a manual probe sat PENDING > 4 min; a second was aborted at
60 s, still PENDING), and astroquery's job poll has no timeout. `_cutout_urls` now runs the same
SODA steps with a 10-minute deadline, aborts the stuck job and raises, so the backoff applies and
the census waits out CASDA outages on its own (~2 attempts/h while down).

## CASDA data-access outage (2026-10-16 → at least 2026-10-17 09:15 EDT)

Since ~08:00 EDT Oct 16 no cutout job has left PENDING (23 timed-out attempts by Oct 17 08:51,
one every ~40 min under the backoff). A plain file-staging job fails immediately with
`errorSummary: "Error: A problem occured obtaining access to the requested item(s)"` (job
e44f96ee-b313-4b8f-b548-10f86d203037), while login and TAP work. This is CASDA's archive back end,
not this code or the account. No notice on casda.csiro.au or research.csiro.au/casda. The census
service waits it out and resumes on its own; 1,694 / 19,610 slots are banked.

## C3(b) WISE J0623−0456 — literature read (arXiv:2306.15219, full text, 2026-10-17)

Detected in **RACS-mid** (SBID 21060, MJD 59216): Stokes I 4.17 ± 0.41 mJy beam⁻¹, |V|/I
66.3 ± 9.0% (Sec. II.1). **Not detected in RACS-low** (SBID 8592), 5σ limit 1.81 mJy beam⁻¹
(Sec. II.1; Table 1 note). RACS-mid position 06ʰ23ᵐ09.28ˢ −04°56′22.8″, proper motion
(−0.93 ± 0.01, +0.17 ± 0.02) ″ yr⁻¹, distance 11.44 ± 0.38 pc (Sec. II.4); period 1.912 h,
B > 0.71 kG, emission cut off near 1.9–2.0 GHz (abstract, Sec. III). It is a T8 dwarf, absent
from Gaia and so from the Reylé sample. For this census it is a rotationally pulsed source with
no prior ~0.9 GHz detection, seen in two 15-min snapshots against a 1.9 h period — a
non-detection is the likely outcome and is not a pipeline failure. `run_c3b` measures it once
CASDA's archive recovers; reported, not pass/fail, as the plan states.

## Census complete (2026-10-25) and frozen rules applied (2026-10-26)

19,600 of 19,610 target-epochs fetched (99.95%); the 10 missing are low3 images CASDA refuses
("You do not have access to any of the requested data files", in every retry pass). Frozen rules
(`analyze_census`, read-only dry run, 2026-10-26):

- C1 null monitor on all census nulls: n = 97,967, V > 0 in 50.09% (99% interval 49.59–50.41%),
  std V/σ 1.009, one value ≥ 5σ (Gaussian expectation 0.06; allowed ≤ 3.06). **Passes.**
- 19,576 of 19,600 slots sane; median 3σ |V| limit 0.57 mJy/beam.
- One slot with |V|/σ_V ≥ 5; 0 candidates by the photometry alone, 1 **held** (Stokes I detected,
  cutout leakage floor undetermined: < 3 field sources above 50σ). No target differs between
  epochs under the frozen two-epoch rule (9,771 pairs).

## The held slot: 2MASS J11231964−0509045 (vetting, 2026-10-26)

*Everything in this section after the census row is post hoc follow-up, not part of the frozen
census.*

- **Census row:** Gaia DR2 3786550031235785344, low3, SBID 56518 (MJD 60310.89): V = +0.773 ±
  0.132 mJy/beam (5.8σ), I = 1.14 ± 0.19 (5.9σ), V/I ≈ 0.68 (`classify_emitter`:
  `highly_circular`; far above ASKAP leakage of order 1%, so the undetermined floor cannot
  explain it). low2 (SBID 38679): V −0.1σ, I −0.2σ. Image inspection: compact, beam-shaped I and V
  at the propagated position (peak within one 2.5″ pixel); low2 noise only.
- **Identification** (research agent, VizieR/SIMBAD, locators in its report; key items):
  SIMBAD LM*, M6V (Davoudi et al. 2025, AJ 170, 213, Table 1, SpeX S/N 69, ±0.5); M7 photometric
  (Reylé 2018, J/A+A/619/L8 `sptphot`) and M7/dM7 from four LAMOST catalogues (S/N 14–17). So it
  is at the M6/M7 boundary: a UCD only by the loose definition. Rotation P = 0.129 d (Newton et
  al. 2018 J/AJ/156/217; Fetherolf et al. 2023 J/ApJS/268/4); TESS flares (Ren et al. 2026,
  J/ApJS/282/29: 4 flares; Seli et al. 2025); Hα emission EW 13–17 Å (LAMOST). No radio-star
  catalogue lists it (Driessen et al. 2024 SRSC, Wang et al. 2024, V-LoTSS); no paper on it.
  **Its RACS-low3 Stokes I source is catalogued but unidentified:** Thomson et al. 2026 (PASA 43,
  89) RACS-LOW3 J112319−05090, Fpeak 1.13 mJy, 3.1″ from the propagated position.
- **Every other ASKAP observation of the position** (`followup_lightcurve`,
  `results/ucdv_followup_j1123.json`): 36 observations with I+V images; 6 refused by CASDA, 1
  (SBID 56437, low3 the day before) NaN at the target pixel (tile edge), **29 measured**.
  Significant epochs, all **positive** V:

  | SBID (project) | MJD | MHz | V (mJy/beam) | I/σ |
  |---|---|---|---|---|
  | 56518 (RACS-low3, the census row) | 60310.9 | 943 | +0.773 ± 0.132 (5.8σ) | +5.9 |
  | 61512 (VAST) | 60422.5 | 887 | +0.485 ± 0.150 (3.2σ) | +2.8 |
  | 64973 (VAST) | 60544.2 | 887 | +0.449 ± 0.154 (2.9σ) | +2.9 |
  | 82408 (EMU, deep) | 61091.6 | 943 | +0.174 ± 0.030 (**5.9σ**) | +5.1 |
  | 83189 (VAST) | 61125.7 | 887 | +0.502 ± 0.139 (3.6σ) | −0.8 |

  The other 24 measured epochs (2020–2026, 855–1655 MHz) span −1.8σ to +2.2σ; none ≤ −2.5σ.
- **Statistics (post hoc):** excluding the census epoch, the best epoch (EMU 82408, 5.9σ,
  one-sided p = 1.9e-9) survives 28 looks at p ≈ 5e-8; the 28 combine to Stouffer z = +3.3. The
  EMU image shows a compact V source within 2″ of the propagated position (one 2″ pixel) and an
  I counterpart at 5σ. Two independent observations 2 yr apart detect the star at > 5σ in V with
  the same sign.
- **Reading:** a new detection of strongly circularly polarized, intermittent radio emission from
  a 3.1-h-rotating M6/M7 flare star at 21 pc, consistent in sign over 2024–2026 — the behaviour
  expected of electron-cyclotron maser emission (interpretation, not measured: no time
  resolution within the 15-min snapshots has been examined).
- **Still open:** GATE-2; whether the RACS VIII/Pritchard V tables (not on VizieR) list it; a
  time-resolved look within SBIDs 56518 and 82408 (needs visibilities); the census paper's UCD
  definition given the M6/M7 type.

## GATE-2 (2026-10-26): PASS-WITH-CHANGES, changes applied

Independent read-only science review of the census and the detection; worktree confirmed clean
afterwards. Findings and what was done:

1. **Census null valid** (NOTE). Forced photometry reads one fixed pixel; nulls V > 0 50.09%,
   std 1.009. The nulls show a mild heavy tail at 3σ (155 with |z| ≥ 3 vs ~132 Gaussian, +17%);
   the frozen C1 criterion tests only 5σ. Stated, not a failure.
2. **Snapshot-fraction limits (MUST-FIX, fixed).** The denominator used each slot's measured
   |V| + 3σ_V, which is data-dependent (drops positive fluctuations and the detection's own slot)
   and not the detectability the 5σ rule implies; held detections never entered the numerator;
   the L = 1e14 row (n_eff = 2, "limit" 1.50) constrained nothing. Now: sensitivity = 5σ_V (and
   7σ_V, C2's 96.5% completeness point) at the slot's distance, independent of V; held and
   candidate slots enter the numerator at their measured L; limits > 1 are not reported; quoted
   L = 3e14, 1e15, 1e16. Results (5σ / 7σ): 3e14: n_eff 17 / 4, one detection, limit 0.28 / —;
   1e15: 458 / 225, limit **0.0065** / 0.013 (was 0.0035 under the old definition); 1e16:
   17,886 / 13,152, limit **0.00017** / 0.00023. The detection's L = 4.23e14 erg s⁻¹ Hz⁻¹ is a
   macro. The two superseded 1e14 keys were deleted from `results/ucdv_metrics.json` by hand
   (the merge guard keeps keys by design).
3. **Trial factors (SHOULD-FIX, reworded).** The best follow-up epoch (EMU 82408, p = 1.9e-9 one-
   sided) survives 28 looks at p ≈ 5e-8. The Stouffer +3.3 over 28 epochs is carried by the same
   four selected epochs: the other 24 sub-2.5σ epochs give mean z = +0.08 (Stouffer +0.37), so
   there is **no sub-threshold excess** and the Stouffer number is a restatement, not independent
   support. The census epoch itself carries a look-elsewhere factor of 19,600.
4. **EMU noise tail (SHOULD-FIX, done, post hoc).** Forced V at blank positions ≥ 16″ apart in
   the saved cutouts, avoiding the target and I > 5σ: EMU 82408, 465 positions, V > 0 49%, std
   0.95, max |z| 3.57, 2 at ≥ 3σ (1.3 expected); RACS 56518, 391 positions, 50%, 1.03, 2.97, 0
   (1.1 expected). The two 5.8–5.9σ values sit well outside their own images' noise.
5. **Chance coincidence (SHOULD-FIX).** Census nulls: ≤ 1 event ≥ 5σ in 97,967, i.e. ≲ 0.2
   expected among 19,600 target slots (an upper bound; consistent with zero). Catalogue density
   (~221 RACS-low2 V sources over the survey, reviewer's area assumption) gives ~1e-4 within 3″
   over all 9,805 targets. Proper motion moves the star only ~0.5″ over 2024–26, so variability,
   not motion, is what links the emission to the star.
6. **Astrometry (SHOULD-FIX, reworded).** I peak ~0.6″ and V peak one 2.5″ pixel from the
   propagated position; the catalogued RACS-low3 I source is 3.1″ away, ~2σ for a 6σ source in a
   13″ beam. Wording: "consistent within the positional uncertainty", not a sub-pixel claim.
7. **V/I (NOTE).** ~0.7 ± 0.2 (I is only 5.9σ); ≳ 30× typical ASKAP leakage, so leakage is
   excluded even though this cutout's floor is undetermined.
8. **Framing (SHOULD-FIX, adopted for the note).** "An M6–M7 dwarf in the Reylé UCD sample
   (spectroscopic type M6V)", not "an ultracool dwarf". Not "new"/"first" until the Pritchard
   2021/2024 and RACS VIII V tables are checked; until then, "not in the radio-star catalogues
   checked (SRSC, Wang et al. 2024, V-LoTSS)". "Forced-V measurements show", not "discovered".
   "Intermittent" is inferred from the non-detections; electron-cyclotron maser is
   interpretation only.
9. **Two-epoch rule (NOTE).** This star's low2/low3 pair differs by ~5.9σ in V but is not flagged
   because the frozen comparison includes the 10% scale term and the spectral bound; one sentence
   in the note.
10. **Accounting (SHOULD-FIX).** The note states 19,610 slots → 10 refused → 19,600 fetched → 24
    failed the per-cutout sanity check → 19,576 with limits. The follow-up count is 29 measured.

## C3(b) WISE J0623−0456 measured (2026-10-26; reported, not a gate)

`run_c3b` (`results/ucdv_controls_c3b.json`): low2 SBID 39068, V = −2.43 ± 0.15 mJy/beam
(V/σ = −16.3), I/σ = 8.7, |V|/I ≈ 0.7; low3 SBID 56420, V = −0.96 ± 0.15 (−6.6σ), I/σ = 6.0,
|V|/I ≈ 0.5. The pipeline recovers this known emitter in both epochs, with one sign. The 10-17
expectation above ("a non-detection is the likely outcome") was wrong; \citet{rose2023} reported
no RACS-low1 detection, and its RACS-mid circular fraction was 66%. Sign convention versus Rose et
al. not compared (their sign not read).

## RNAAS note drafted (2026-10-26)

`papers/ucdv/rnaas.tex`: abstract ~120 words, body ~890; every number a macro (detection,
follow-up, C3(b), accounting and fraction limits added to `_MACROS`); figure = (a) fraction limit
vs L, (b) the star's V/σ in every ASKAP epoch. Prior-report check before writing (research agent,
full text/tables): not in Pritchard et al. 2021 (Table 2), Pritchard et al. 2024 (Table 4, all
76 stars), or the RACS VIII V catalogue (CASDA TAP `AS110.racs_low2_v_v01`, 221 rows, none within
6′); SRSC, Wang et al. 2024 and V-LoTSS negative earlier. Wording therefore "absent from the
catalogues of …", not "first". Four reference titles reconstructed from truncated ADS output were
wrong (Davoudi, Thomson, Duchesne, Pritchard 2024) and were replaced with the Crossref titles.
Style lint (RNAAS genre) and triage clean; compiles in the tectonic container, 3 pages.
