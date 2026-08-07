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

## Next

- Increment 2: RACS southern leg (the categorical first — SDSS quasar spectra south of
  −40° × RACS), racsradio-carton validation against its selecting survey, two-survey
  synthetic variant, north/south contrast at matched luminosity limits.
