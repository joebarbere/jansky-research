# dr20radio — findings (in progress)

Slice plan: `plans/88-dr20-radio-census.md`. Target: the first radio-counterpart census of
the SDSS-V DR20 Black Hole Mapper quasars (VLASS north / RACS south).

## Increment 1 + northern leg — done 2026-08-07

Tooling (7 offline tests: synthetic fraction recovery to ±0.03 after measured false-match
correction; carton-circularity exclusion demonstrated to matter) + the real VLASS leg:
DR20 `spAll-lite` (506,550 rows) → **202,691 clean quasars** (`CLASS=QSO`, `ZWARNING=0`)
north of −40°, 280 radio-carton objects excluded from fractions.

| quantity | E2 | E3 | any epoch |
|---|---|---|---|
| matched (2.5″) | 8,465 | 8,900 | 9,466 |
| raw fraction | 4.18% | 4.39% | **4.67%** |
| measured false-match | 0.009% | 0.009% | — |
| corrected | 4.17% | 4.38% | — |

- **Redshift trend:** ~5.1% (z<0.5) → ~4.1% (z≈2–2.5) → 6.9–7.2% (z>3) — the classic
  flux-limited-selection rise at high z; reported per selection class in the paper, never
  as "the quasar population" (BHM targeting is heterogeneous).
- **Observatory nuance:** 74,462 of the *northern-sky* census quasars were observed from
  LCO (du Pont) — telescope hemisphere ≠ sky hemisphere; the southern-first claim is about
  sky south of −40°, where only RACS reaches.
- **An expectation corrected (kept honestly):** the radio-carton objects matched VLASS at
  only 31–35%, not the ~100% the plan anticipated — and that is physics, not failure. The
  cartons were selected from RACS (888 MHz) / LOFAR (144 MHz); steep-spectrum sources fade
  below VLASS's 3 GHz depth. The ~100% validation expectation transfers to the RACS leg
  (increment 2), and the 31–35% VLASS rate becomes a free cross-frequency measurement. The
  synthetic test masked this by planting carton counterparts in the same catalog — the
  synthetic will gain a two-survey variant in increment 2.
- E3 interim lists (VLASS Memo 22) use a simplified schema (binary `Flag`, no
  Duplicate/Quality flags, no empty islands) — reader handles each epoch's own convention.

Committed evidence: `results/dr20radio_north.json`.

## Increment 2 — the southern leg (the categorical first) — done 2026-08-07

RACS-low DR1 fetched strip-wise from CASDA TAP (115 × 1° Dec strips, resumable, 2,123,638
sources cached locally), matched at 5″ (RACS's 25″ beam motivates the wider radius; false
match measured as always).

| block | n quasars | matched | corrected fraction |
|---|---|---|---|
| **δ ≤ −40° (SDSS × RACS: first ever possible)** | **73,074** | 2,930 | **3.95%** |
| −40° < δ ≤ +30° (VLASS cross-check band) | 174,171 | 6,626 | 3.76% |

- The deep-south census is **pure LCO** (73,074/73,074) — the du Pont exclusivity the
  categorical-first claim rests on, confirmed in the data itself.
- Overlap-band consistency: RACS 3.76% vs VLASS 4.67% in comparable sky — same order, with
  the difference carrying the frequency/depth information the luminosity-matched contrast
  (increment 3) will quantify properly.
- Deep-south z-trend: elevated at z<0.5 (6.6% — host/low-z AGN effects), flat ≈3.4–3.9%
  through the bulk, rising at z>4 on small numbers (11%, n=73).
- **Carton validation, split by selecting survey (the increment-1 lesson paying off):**
  racsradio cartons match RACS at **88.9%** (56/63 — the near-expected pipeline validation;
  the residual is plausibly catalog-version/Galactic-cut differences between the targeting
  product and the DR1 galactic-cut source list, noted for the paper), while lofarradio
  cartons match at **17.2%** (5/29) — 144 MHz-selected steep-spectrum sources mostly fade
  below 888 MHz, the same cross-frequency physics as the 31–35% VLASS figure. The pooled
  66% that first came out was two clean populations averaged; the two-survey synthetic now
  regression-tests this whole class of confusion.

Committed evidence: `results/dr20radio_south.json`.

## Increment 3 — the luminosity-matched north/south contrast — 2026-08-07

Fluxes wired through the whole chain (`crossmatch` now returns the counterpart index; the
VLASS loaders carry `Peak_flux`; RACS fluxes were already cached). Both surveys' flux limits
(VLASS 1.0 mJy — CIRADA reliability threshold; RACS 3.0 mJy — Hale et al. 2021 95%
completeness) convert per-quasar to rest-1.4 GHz luminosity limits (α = −0.7 K-correction,
Planck18); a quasar counts only if matched AND its counterpart clears the COMMON (max) limit.

| leg | raw corrected | above common L limit |
|---|---|---|
| North (VLASS any-epoch) | 4.67% (raw) | **4.06%** |
| Deep south (RACS, δ≤−40°) | 3.95% | **2.82%** |
| Overlap band (RACS) | 3.76% | 3.12% |

The gap survives luminosity matching — but the overlap band (same-ish sky as VLASS's 4.67%)
shows RACS recovering less at fixed L, so a survey-side component is real: candidate causes
are the α = −0.7 assumption (flat-spectrum cores are relatively brighter at 3 GHz, so the
common limit filters them asymmetrically), RACS's 25″ beam (blending/position shifts vs the
5″ radius), and hemisphere targeting-mix differences (northern BHM cartons lean
eROSITA-selected). The paper reports the contrast with these stated as caveats — it is a
measured survey comparison, not a hemisphere-physics claim.

Figure + macros: `paper_assets` (`--paper`) renders `dr20radio_fractions.pdf` (raw and
luminosity-matched fractions vs z, Wilson bars) and `generated/macros.tex` from the
committed evidence only.

## GATE-2 + paper — 2026-08-07

`papers/dr20radio/main.tex` (AASTeX, 4 pp) drafted from macros generated off the committed
evidence; science-reviewer pass returned **no blockers**, four should-fixes, all applied:

- **DR20 byline**: the data-release paper's author is the SDSS Collaboration, not
  "Almeida et al." — refs.bib and the plan corrected.
- **RACS completeness is two numbers** (Hale et al.: ~3 mJy source-count based, ~5 mJy
  simulation based): both stated; the luminosity-matched contrast repeated at the
  conservative 5 mJy — north 3.45% vs south 2.39% (primary: 4.06% vs 2.82%).
  **Superseded 2026-08-12 (referee round):** this was reported as showing the gap robust to
  the choice. It shows nothing of the kind — the common luminosity limit is the RACS one in
  *both* legs, so raising the flux floor rescales north and south identically and the ratio
  cannot move (1.4391 → 1.4434). The parameter the contrast actually depends on is the
  assumed spectral index: sweeping α over 0…−1 moves the gap 0.23 → 1.66 pp. See the
  `luminosity_matched_alpha` block in `results/dr20radio_{north,south}.json`.
- **North carton validation split by selecting survey** (same lesson, third appearance):
  at 3 GHz, RACS-selected cartons recover at 49%, LOFAR-selected at 27% — the previously
  quoted pooled ~31% was again an average over two populations. The cross-frequency fading
  table is now complete: selecting-survey 89%, 888 MHz→3 GHz 49%, 144 MHz→888 MHz 17%,
  144 MHz→3 GHz 27%*. (*The 144 MHz class recovering slightly better at 3 GHz than at
  888 MHz reflects the different comparison samples/depths, not spectra turning over —
  noted, not interpreted.)
- **Overlap-band luminosity-matched number surfaced** (3.12% vs north's 4.06% on shared
  sky) — the survey-side-effects claim now cites its evidence.
- Nits: RACS footprint wording (−80° not the TAP query bound), Wilson 68% stated, dropped
  out-of-z-range quasars (<0.15%) noted in Limitations, FIRSTCARTON-only exclusion caveat
  added, north z<0.5 elevation acknowledged, obs breakdown now census-only, plan status
  updated.

## Next

- Increment 2: RACS southern leg (the categorical first — SDSS quasar spectra south of
  −40° × RACS), racsradio-carton validation against its selecting survey, two-survey
  synthetic variant, north/south contrast at matched luminosity limits.

## Referee round 2 (2026-08-12) — after the α demotion

The revised draft went back to the referee. Verdict: **major revision**, 16 findings. All the
new α macros recomputed correctly, but the paragraph written to *repair* the first-round
finding contained a fresh error, which is the finding worth remembering:

1. **BLOCKER (my own regression).** The new sentence gave the northern α-sweep range as
   `\drLumNorthFlatPct--\drLumNorthConsPct` = 3.06–3.45%. `\drLumNorthConsPct` is the 5 mJy
   conservative variant — a *different axis*, and the very one the same paragraph had just
   declared incapable of testing this. The real range is 3.06–**4.37%**. The range endpoint
   had no macro, so I reached for the nearest-looking name. Fixed by emitting
   `\drLumNorthSteepPct`/`\drLumSouthSteepPct`: **if a range needs an endpoint, the endpoint
   gets its own macro.**
2. **"Unchanged to three decimal places" was false** — 1.4391 vs 1.4434 differ *in* the third
   decimal. Now quoted as the two ratios via `\drRatioFid`/`\drRatioCons`, a 0.3% shift.
3. **The stated mechanism was false at the end of the sweep that sets the headline.** The
   common limit is the larger of the two K-corrected limits; RACS binds for α ≳ −0.9, but
   between −0.7 and −1 the VLASS limit overtakes it (3 mJy at 888 MHz → 1.90 mJy at 1.4 GHz,
   against 1 mJy at 3 GHz → 2.14 mJy). So at α = −1 about a quarter of the gap comes from the
   *south falling*, not the north rising. Now stated with the crossover.
4. **The southern denominator was never intersected with the RACS footprint.** `deep_south` is
   a pure declination cut, so quasars in RACS-low DR1's |b| ≲ 5° hole entered the denominator
   as guaranteed non-detections. Measured rather than assumed: **52 objects, 0.071%**, none
   below the −85° floor. Now computed by `run_south` into `racs_footprint` and reported in
   Limitations as a bound. The bias is ~50× smaller than the fraction itself.
5. **The carton "same-objects" claim was wrong** — the 888 MHz leg is a small, more southerly
   subset of the 3 GHz leg's sky, and the rates differ by ~1σ on samples of tens. Downgraded
   from "does *not* track the frequency jump / depth dominates" to a consistency statement.
6. **A sign error in the K-correction would have inverted the paper and no test would see it.**
   `tests/test_dr20radio.py::test_k_correction_sign_and_alpha_monotonicity` now pins the sign,
   the RACS→VLASS crossover, and monotonicity in α.

Also fixed: Limitations/Introduction/figure caption did not carry the demotion (the paper had
it in the abstract and summary only); "nearly flat at ~4%" described a decline significant at
many σ; `\drNorthAnyPct` is raw, not chance-corrected.
