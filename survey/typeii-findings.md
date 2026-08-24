# Findings — OVRO-LWA type II burst detector + census (plan 50, fable-ideas F13)

`jansky_research.typeii`: a slow-drift + fundamental/harmonic type II solar radio burst detector
for the OVRO-LWA archive, with a synthetic recover-a-known and a CME cross-match. The published
OVRO-LWA real-time detector (arXiv:2603.25446) is type-III-only; this fills the type II gap.

## GATE 0 (2026-07-09)

- **Data ACCESSIBLE as FITS arrays**: OVRO-LWA Level-1 beamforming dynamic spectra, FITS,
  13.4–86.9 MHz, 768 ch, 256 ms, ~0.6 GB/day, 2024-04→present; pattern
  `ovro-lwa.lev1_bmf_256ms_96kHz.YYYY-MM-DD.dspec_I.fits` (freq list + time list + Stokes-I dspec).
- **Novelty PASS**: arXiv:2603.25446 is a YOLO type III detector (type-III-only). No OVRO-LWA type
  II census exists; RSTN Cycle-24 catalogue (Lawrance+2024, 429 events) is a different instrument
  + cycle; arXiv:2512.21846 is an N=10 case study. Genuinely unclaimed.
- **Type II signature pinned**: slow drift <~0.1 MHz/s (to ~1) vs type III tens–hundreds;
  fundamental+harmonic ratio ~2; band-split ~1.23; duration minutes (61% <15 min); recover-a-known
  = Gopalswamy+2005 fast(>900 km/s)-and-wide(>60°) CME driver bias.
- Cross-match catalogues open: LASCO CME v2 (CDAW), GOES flares (NCEI NetCDF/CSV), SILSO.

## Detector

Three physical discriminators on a background-subtracted dynamic spectrum: (1) slow Theil-Sen drift
rate (0.01–2 MHz/s, the 2 a guard band below type III); (2) duration ≥90 s; (3) ridge coherence
|corr(t,f)| ≥ 0.55, which rejects scattered noise (low corr) and horizontal RFI (near-zero drift).
Harmonic score (co-emission at 2×/0.5× the ridge freq) reported for purity, NOT required
(single-lane type IIs exist). Narrowband RFI masked by channel-variance outlier. Reuses the
`solarbursts` coronal-density drift model (driven by a slow shock, not a fast beam) + its
`background_subtract`.

## Synthetic recover-a-known (offline, in CI)

Mixed-difficulty census — 24 injected type II (strong+harmonic, single-lane, near-threshold-SNR),
16 type III contaminants, 8 RFI-only:

- **Purity 1.0** — no type III or RFI misclassified as type II.
- **Completeness is honestly SNR-dependent** (NOT a single saturated number):

  | injected SNR | 2 | 2.5 | 3 | 4 |
  |---|---|---|---|---|
  | completeness | 0.33 | 0.63 | 0.92 | 1.0 |

  Strong bursts fully recovered; near-threshold genuinely missed — the honest performance floor.
- CME cross-match recovers the injected fast-and-wide bias (0.82 fast, 0.86 wide, median 1116
  km/s) — **a WIRING CHECK**: these echo the injected CME distribution and validate that the
  temporal match picks the right CME among decoys, NOT an independent Gopalswamy reproduction
  (that needs real events).

## The real census (ran 2026-07-09/10): a FALSE-POSITIVE-DOMINATED NULL

**The data is on AWS Open Data** (public bucket `ovro-lwa-solar`, `spec_fits/<YYYY>/<YYYYMMDD>.fits`;
no login/CAPTCHA — the Turnstile only gates the query UI). Each daily file is ~1.7 GB (4D I/V
dynamic spectrum, ~15–85 MHz, ~0.26 s), so the real leg **streams**: `stream_dspec` opens the S3
FITS lazily (astropy `use_fsspec`), range-reads ONLY the Stokes-I plane in frequency chunks, and
block-averages to ~4 s bins on the fly — a day is processed **entirely in memory, nothing on disk**
(~30 s/day at the measured 30 MB/s), then freed. Ran **all 765 observing days 2024-04→2026-07,
0 failures**.

**Result: 331 candidates — and the honest purity test says they are false-positive dominated, NOT
a type II census:**

- **Rate 0.43/day** — ~4× the RSTN Cycle-24 metric type II rate (~0.11/day). First red flag.
- **CME association is background-like**: matched-CME median **478 km/s** ≈ the background CME
  median (**379 km/s**, only 6% fast) — the general CME population, not the fast (≥900 km/s) wide
  drivers real type II require. Observed match rate (0.55) is BELOW the chance rate (0.64: at solar
  max, any random time has a CME within ±2h). The ±2h window is physically right (metric type II
  precedes the CME's LASCO first-appearance) but near-powerless at solar max — the null rests on
  the SPEED distribution, not the match-rate deficit.
- **drift ⊥ CME speed** (Pearson r = **0.085**): the burst drift encodes shock kinematics, so a
  real sample would correlate with CME speed. It doesn't → the matches are coincidence.
- **83% window-saturated**: 0.831 of "bursts" fill the entire 15-min window (half at exactly the
  window length). Real type II are varied-duration transients; a window-filling ridge is persistent
  RFI/background tracked as a slow drift — the classic false-positive signature.
- **Harmonic cut makes it WORSE** (harm≥0.5 → 12% fast, 352 km/s): the "harmonic" candidates lack
  real fundamental+harmonic structure.
- **The flare-gated subset is a CONFOUND, not a rescue**: X-flare-coincident candidates look
  textbook (13/13 CME-matched, ~967 km/s, 69% fast-wide) — but that's the flare-size↔CME-speed
  relation (matched speed tracks flare CLASS: X 967 / M 544 / C 381 / none 488), not the detector
  finding type II (drift still ⊥ CME speed). The May-2024 Gannon-storm-week pilot (which looked
  like a clean 83%-fast recovery) is the same confound at its extreme — that week was wall-to-wall
  fast-wide CMEs, so even false positives coincided with real fast drivers. Does not generalize.

**Honest bottom line:** a blind slow-drift spectral detector does NOT yield a usable type II census
in OVRO-LWA's RFI-heavy 15–85 MHz band — exactly why the archive's published detector is
type-III-only. Not a detector defect (synthetic completeness/purity are as reported) but a
real-data PURITY limit from RFI/background; a reliable census needs interferometric imaging (source
at the eruption site) or human vetting. A handful (~order a few) of the 331 may be real type II but
are not separable with the available spectral features. **We claim no census, no occurrence rate,
no individual detection.** `purity_diagnostics` makes this conclusion pipeline-generated
(`association_is_background_like`, `drift_cme_speed_corr`, `window_saturation_frac`).

- **Real FITS format confirmed on first contact** (4D PRIMARY array, axes time/freq/1/Stokes, Jy);
  `parse_lwa_dspec` handles it (+ a fallback to the older table layout), unit-tested on both.
- The streamed real census runs via `scripts/typeii_real.py --start 2024-04-01 --end 2026-07-09`
  (needs the `typeii` extra: fsspec + aiohttp). CDAW LASCO CMEs + HEK GOES flares + SILSO fetched
  automatically; CDAW only publishes through ~2026-02 (2026-03→07 lack CME data — noted).

## GATE-2 — detector build (PASS with required fixes, all applied)

- **R1**: dropped the false claim that `ecallisto_census.coverage_corrected_rate` is an active
  reuse (it's now genuinely used in the real leg's occurrence piece).
- **R2**: the 1.0/1.0 was an easy-synthetic ceiling → added single-lane + near-threshold cases;
  completeness is now 0.92 with the reported SNR curve.
- **R3**: labelled the synthetic CME association a wiring check (echoes the injection).
- Suggested, applied: in-memory-FITS test for `parse_lwa_dspec`; single-lane + weak-SNR assertions.

## GATE-2 — real census (PASS with required fixes, all applied)

The reviewer independently re-ran the confound tests on the 331-event list and **confirmed the
null is correct** (not a missed real sub-population). Required fixes, all applied:

- The **flare-gated enrichment is a confound** (X-flare subset looks textbook fast-wide) — proven
  by drift⊥CME-speed (r=0.09) and matched-speed-tracks-flare-class; must be in the writeup, is.
- The **window-saturation FP signature** (83% fill the window) — must be reported, is.
- The **CME match is near-powerless at solar max** (chance 0.64) — the null rests on the speed
  distribution + confound, not the match-rate deficit; stated.
- **Must NOT claim** a census/rate/detection, or that the flare-gated enrichment is a real
  sub-population, or that occ r=0.07 is a meaningful cycle-phase null (low leverage) — none claimed.
- The confound evidence (`drift_cme_speed_corr`, `window_saturation_frac`) is now
  pipeline-generated in `purity_diagnostics`, so the artifact carries it, not just the prose.

## Reproduce

Offline (detector + synthetic + tests): `uv run python -m jansky_research.typeii --offline --out "$(mktemp -d)"` (never `--out .` — that is the invocation that once gutted the real metrics)
Real (streamed, in memory, no disk): `uv run --extra typeii python scripts/typeii_real.py --dates
2024-05-14 2024-05-15 --cme data/typeii/lasco_cme.csv`.


## Regeneration note (2026-08-04)

The real census was re-run end-to-end (same pipeline, same date range) after the 2026-07-30
offline-build incident clobbered the original metrics: the archive now serves **768** observing
days (3 backfilled by OVRO since the original 765-day run) yielding **332** candidates (was 331).
All conclusions unchanged — the extra day's candidate joins the false-positive-dominated
population; the purity confounds re-derive to the same verdict. The paper macros now carry the
regenerated numbers; this file's prose above describes the original run.

## Full referee round (2026-08-24): MAJOR REVISION, 15 findings, one BLOCKER

The null is real and verified (window saturation 0.828 recomputed exactly; 167/332 at the window
length; the flare-class trend confirmed). What blocks: the apparatus around it.

**BLOCKER: the tracked arxiv.yaml still says "48 days, zero failures"** -- the synthetic event
count in place of 768 real days, the exact number from the original \tiiNEvents clobber, in the
hand-authored submission abstract that sits outside the macro system's protection, already baked
into the built arxiv-submission/metadata.yaml.

**MAJORs:** the entire synthetic leg (purity 1.0, the SNR curve) has NO committed results file --
structurally, since preserve_live_results correctly refuses to let the offline leg write into the
real metrics file, and the slice has nowhere else to put it; meanwhile the paper claims "every
synthetic number above is pipeline-generated (results/typeii_metrics.json -> macros)", false for
every one of them. Purity 1.0 is one seed and 24 negatives (95% bound: 11.7%), quoted with no
denominator or interval. The purity test cannot fail the way the real data fails: the injected
contaminants (fast type III, zero-drift RFI) sit orders of magnitude from the acceptance band,
while the class that dominates the census -- coherent slow-drift window-filling background, |df/dt|
0.012-0.064 -- was never injected; and the 24-negative test is ~25x too small to see the census's
~0.5%/window false rate, so "this is not a detector defect" is unsupported. **The match-rate
deficit is an artifact of the CDAW coverage gap**: coverage ends 2026-03-31, 57 detections (17%)
lie beyond it with no possible match, and restricted to the covered span the comparison REVERSES
(0.676 > chance 0.619) -- the committed association_is_background_like boolean is wrong as
computed. The harmonic-cut "confirming" sentence rests on 8 CMEs and reverses at >=0.7 (median
469, frac_fast 0.25). The SNR-completeness curve was measured at 300 s / 0.5 s bins / Gaussian
noise vs the deployed 900 s / 4.1 s / real band, with every injection window-saturating and one
duration/speed -- so the duration and drift cuts have no measured completeness. The "~4x RSTN"
compares a solar-max-only rate to a whole-cycle average, and n_windows_swept is unrecorded. The
matched sample is 3.3x enriched in fast CMEs (z=7.8) and the one macro showing it
(\tiiRealFracFast=0.194) is the one macro left uncited while \tiiRealBgFracFast IS cited in the
"essentially the background" sentence.

**MINOR/NIT:** 14 candidates are the same structure in adjacent 50%-overlap windows (dedup: 318);
the Method says a ~1-hr asymmetric match window where the code uses +/-2 h symmetric; r=0.058
lacks its CI (-0.09, +0.20); the committed figure has an INVERTED frequency axis (drift appears
upward in a slow-drift paper) and renders nulls as zero-height bars, and the paper includes no
figure at all; lawrance2024 pages=75 is an article number; the findings file itself tells readers
to run --out . (the gutting command).

**Status: fixes pending** (all recomputations possible from the committed event list; purity
ensemble offline).

### Resolved 2026-08-24 (revision): the contaminant that can fail was injected, and it failed

All 15 findings addressed. The centrepiece: `synthetic_slow_background` injects the class the
real band actually contains (coherent, slowly drifting, window-filling ridges with |df/dt| drawn
from the real candidates' 0.01-0.07 MHz/s range), and a 30-seed ensemble measures the detector
against every class:

| class | trials | false positives | 95% upper bound |
|---|---|---|---|
| fast type III | 480 | 0 | 0.8% |
| narrowband RFI | 240 | 0 | 1.6% |
| **slow-drift background** | 240 | **91** | **44%** |

The old "this is not a detector defect" is retracted in the paper: the detector is clean against
contaminants that cannot enter the acceptance band and false-positives at ~38% against the one
that can, so the unusable census is a JOINT property of the band and the detector -- which is the
stronger, mechanistic form of the paper's null (the census's 83% window saturation and the
ensemble's failure mode now corroborate each other).

**The match-rate deficit is retracted as a coverage artifact.** purity_diagnostics is now
coverage-aware: with a fresh CDAW fetch (through 2026-05; 4407 CMEs), 303 of 332 detections fall
inside coverage and the covered match rate is **0.68, above the 0.612 chance rate** -- the old
deficit came from counting post-coverage candidates as unmatched. The committed
`association_is_background_like` now keys on the speed comparison and window saturation, not the
match rate. The speed test itself is sharpened: matched median 476.5 vs background 372.0, matched
fast-fraction 0.184 vs 0.057 (KS p ~ 0) -- a real but modest enrichment, exactly the flare-CME
confound's signature, nowhere near the >=900 km/s driver population; the previously-uncited
`\tiiRealMatchFracFast` is now in the same sentence as the background value.

Also fixed: the arxiv.yaml BLOCKER ("48 days, zero failures" -> 768 -- the original clobber's
number surviving outside the macro system); the synthetic leg now writes its own
`results/typeii_synthetic_metrics.json` (allowlisted; the Reproducibility sentence had claimed a
provenance that was false for every synthetic number); the harmonic-cut sweep is committed with
its denominators (8 matched CMEs at >=0.5, trend reversing at >=0.7 -- "confirming" -> "consistent
with, on counts too small for the word confirm"); adjacent-window duplicates counted (332
detections = 320 distinct structures, stated); the "4x RSTN" comparison scoped as cycle-phase
inflated (whole-cycle average vs solar-max-only census); the +/-2 h symmetric match window stated
correctly in the Method; r = 0.044 quoted with its pair count; the committed figure's inverted
frequency axis fixed and null-as-zero bars removed; lawrance2024's article number moved out of
pages; and this file's own reproduce command no longer tells the reader to run `--out .`.

Values that moved with the fresh CME catalogue (three more months of coverage): matched median
472.5 -> 476.5, background 378.0 -> 372.0, r 0.058 -> 0.044, chance 0.619 -> 0.612. The verdict is
unchanged; every diagnostic is now recomputable from the committed event list plus the stated
CDAW months.
