# 30 — Synthesis paper: a type III electron beam from the corona to 0.4 AU, geometrically validated

Status: ✅ done (module + figure + paper + GATE-2; real reproduce gives corr 0.989, ratio 2.18, reach 0.49 AU)

## Context

Four slices already track solar type III radio bursts with the same **drift-to-distance** idea — the
emission sits near the local plasma frequency, so the frequency drift, inverted through a density model,
gives the electron beam's heliocentric distance and speed:

| slice | instrument | band | reach | method / event |
|---|---|---|---|---|
| `solarbursts` | e-Callisto (ground) | 20–90 MHz | corona, ~1.5–3 R⊙ | Newkirk; 2011-09-14, 0.14 c, R²=0.90 |
| `windwaves` | Wind/WAVES RAD2 | 1–14 MHz | 2.4→10.2 R⊙ (Alfvén surface) | Leblanc; 2003-10-28 X17, 0.083 c |
| `swaves` | STEREO/WAVES HFR | 0.125–16 MHz | 2.3→82.6 R⊙ = **0.38 AU** | Leblanc; **2013-05-15**, 0.150 c |
| `triangulate` | STEREO-A+B direction-finding | — | 15→106 R⊙ | **independent 3D geometry**; **2013-05-15** |

The whole drift-to-distance ladder rests on an *assumed* density model (Newkirk in the corona, Leblanc
in the heliosphere). The fresh contribution is that **one of these distances is independently validated
by geometry**: `swaves` and `triangulate` analyse the **same 2013-05-15 event**, so for that burst we
have both a STEREO/WAVES density-model distance (to 82 R⊙) and a STEREO-A+B triangulated geometric
distance (15–106 R⊙) — and they agree, corr(r_geom, r_plasma) = **0.989**. No single slice makes this
point; together they do. This paper unifies the four into one reproducible framework spanning the
corona to 0.4 AU and uses the triangulation as the geometric check on the density-model distance the
whole approach depends on.

This is a synthesis/methods paper for **astro-ph.SR**, genuinely fresh (the geometric validation of the
density-model distance ladder + a single reproducible corona→IP framework), and honest: three of the
four use *different events*, so it is "one framework validated across four instruments and regimes,"
not literally one beam tracked end-to-end (see the optional GATE 0 upgrade).

## Deliverables

- `src/jansky_research/type3synthesis.py` — a thin orchestration module (no new physics): imports the
  four slices, runs each (offline-synthetic for CI, real for `reproduce`), and assembles
  - a **unified distance ladder figure**: heliocentric distance vs emission frequency for all four,
    on one log–log panel from ~100 MHz (corona) to 0.125 MHz (0.4 AU), with the Newkirk and Leblanc
    model curves and the four instruments' coverage shaded;
  - the **2013-05-15 cross-check panel**: `swaves` density-model r(f) vs `triangulate` geometric r,
    with the correlation and ratio;
  - a small **systematics table** common to all four (harmonic vs fundamental factor-2, peak-time vs
    onset speed bias, average-density-model error, the R²/independent-time-sample caveat).
  - `run`/`_figure`/`_write_macros` (every headline number from the four slices' metrics) / `_main`.
- `tests/test_type3synthesis.py` — offline, composes the four synthetic fixtures; 85% floor.
- `papers/type3synthesis/` (AASTeX, astro-ph.SR); `survey/type3synthesis-findings.md`.

## Approach

0. **GATE 0 (optional, headline upgrade) — the single-event hunt.** Search for one type III caught
   simultaneously by e-Callisto **and** Wind/WAVES **and** STEREO/WAVES (a strong 2011–2014 event near a
   STEREO-favourable longitude). If found, the paper becomes "one beam tracked from 100 MHz to 0.125 MHz
   across four instruments with an independent geometric distance" — much stronger. If not, ship the
   framework-synthesis (below); record the search either way.
1. **Orchestration + unified figure.** Build `type3synthesis` over the four slices' existing `run()`
   outputs; produce the distance-ladder figure and the 2013-05-15 cross-check; emit macros.
2. **The argument.** (i) The same drift-to-distance method, with a corona (Newkirk) → heliosphere
   (Leblanc) density-model handoff, spans ~1.5 R⊙ to 0.4 AU across four public instruments; (ii) for
   2013-05-15 the density-model distance is **independently confirmed geometrically** (triangulation,
   r=0.989), validating the model the ladder depends on; (iii) the shared systematics (harmonic,
   peak-time, average density) bound all four consistently.
3. **GATE-2 science review** before write-up: the cross-event framing must not over-claim a single
   beam; the geometric validation must be stated as event-specific; deceleration / speed differences
   between events explained honestly (different events, peak-time bias, density excursions).
4. **Write-up** `papers/type3synthesis/` — a synthesis/reproducibility contribution: a unified,
   geometrically-validated drift-to-distance framework, corona → 0.4 AU. This is the **lead arXiv
   submission** (replacing frbstats in the shortlist).

## Verification

- `make test` / `cov` green on the four synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` produces the unified ladder figure + the 2013-05-15 cross-check + macros.
- (Real-data) `make reproduce` regenerates the four slices and the synthesis from public archives; the
  2013-05-15 geometric/plasma correlation reproduces; GATE-2 sign-off.

## Risks & mitigations

- **Over-claiming "one beam" (highest) →** lead with "one framework across four regimes"; only the
  GATE-0 single-event version may claim a single beam; the geometric validation is explicitly
  event-specific (2013-05-15).
- **Different events / inconsistent speeds →** discuss honestly (peak-time bias, deceleration, density
  excursions); the robust cross-instrument claim is the *distance ladder* and the *geometric validation*,
  not a single speed.
- **Scope (it touches four slices) →** the module only orchestrates existing tested code; no new physics,
  so the risk surface is the figure/argument, gated by GATE-2.
