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

## Status: no real census run yet — but the data is NOT access-blocked (AWS Open Data)

**Corrected 2026-07-09.** Initially I thought the portal (`ovsa.njit.edu/lwadata-query`) was a hard
block: it's a JS SPA behind a Cloudflare Turnstile bot challenge. But reading the app bundle
revealed the actual data is on **AWS Open Data** — public bucket `ovro-lwa-solar`, path
`spec_fits/<YYYY>/<YYYYMMDD>.fits`, directly downloadable with **no login and no CAPTCHA** (the
Turnstile only gates the interactive query UI). Verified: the bucket lists and serves the FITS over
plain HTTPS (`s3_dspec_url`).

- The daily files are **large (~1.7 GB; a 4D I/V dynamic spectrum, ~15–85 MHz, ~0.26 s)**, so the
  real leg **streams** rather than downloads: `stream_dspec` opens the S3 FITS lazily (astropy
  `use_fsspec`) and range-reads ONLY the Stokes-I plane in time-chunks, block-averaging each to
  ~4 s bins on the fly (`_downsample_time`), so a day is processed **entirely in memory with
  nothing written to disk** (peak ≈ one reduced chunk, never the 1.7 GB file). `sweep_day` slides
  a 15-min window across the day; `real_census` frees each day before streaming the next.
- **Real FITS format confirmed on first contact** (`20240514.fits`): a 4D PRIMARY array, axes
  (time, freq, 1, Stokes), Jy, freq from FREQMIN/FREQMAX (GHz). `parse_lwa_dspec` was rewritten to
  this real layout (with a fallback to the older table layout), unit-tested on both.
- The deliverable remains the validated detector + method; the streamed real census (and the SILSO
  occurrence-vs-cycle-phase piece) runs via `scripts/typeii_real.py --dates ... --cme ...` (needs
  the `typeii` extra: fsspec + aiohttp for the lazy cloud reads). A multi-year census is a compute
  follow-on, NOT an access or disk-space block.

## GATE-2 (PASS with required fixes, all applied)

- **R1**: dropped the false claim that `ecallisto_census.coverage_corrected_rate` is an active
  reuse (it's deferred with the real leg — disclosed).
- **R2**: the 1.0/1.0 was an easy-synthetic ceiling → added single-lane + near-threshold cases;
  completeness is now 0.92 with the reported SNR curve.
- **R3**: labelled the CME association a wiring check that echoes the injection, not a Gopalswamy
  reproduction (docstring + source + paper).
- **R4**: carried all disclosures (Turnstile blocker, no real census, SNR-curve caveat, injection-
  echo) into this findings file and the paper.
- Suggested, applied: in-memory-FITS test for `parse_lwa_dspec` (the riskiest untested code);
  single-lane + weak-SNR classification assertions; guard-band comment on DRIFT_SLOW_MAX.

## Reproduce

Offline (detector + synthetic + tests): `uv run python -m jansky_research.typeii --offline --out .`
Real (streamed, in memory, no disk): `uv run --extra typeii python scripts/typeii_real.py --dates
2024-05-14 2024-05-15 --cme data/typeii/lasco_cme.csv`.
