# Findings — a type III electron beam from the corona to 0.4 AU, geometrically validated

`jansky_research.type3synthesis` unifies four separately-validated slices into one reproducible
**drift-to-distance** framework for solar type III bursts, and adds the element no single slice provides:
an **independent geometric check** on the density-model distance the whole method rests on.

## The unified ladder (each leg a recover-a-known on public data)

| regime | slice | instrument | density model | reach (recover-a-known) | speed |
|---|---|---|---|---|---|
| corona | `solarbursts` | e-Callisto (ground) | Newkirk | ~1.5–3 R⊙ | 0.14 c |
| inner heliosphere | `windwaves` | Wind/WAVES RAD2 | Leblanc | 2.4 → 10.2 R⊙ (Alfvén surface) | 0.083 c |
| interplanetary | `swaves` | STEREO/WAVES HFR | Leblanc | 2.3 → 82.6 R⊙ = **0.38 AU** | 0.150 c |
| geometric | `triangulate` | STEREO-A+B direction-finding | **none** | 15 → 106 R⊙ | — |

The same drift-to-distance method, with a Newkirk (corona) → Leblanc (heliosphere) density-model handoff,
spans more than three decades in frequency (~100 MHz → 0.125 MHz) and tracks the beam from the low corona
to ~0.4 AU.

## The centerpiece: a model-free check on the model distances

Every drift-to-distance result is only as good as the assumed `n_e(r)`, and that assumption is almost
never checked independently. Here it is: **`swaves` and `triangulate` analyse the same 2013-05-15 event**,
so the STEREO/WAVES plasma-frequency distance `r_plasma(f)` and the STEREO-A+B geometric (triangulated)
distance `r_geom` are two *independent* estimates of the same quantity. Across the triangulated channels
they correlate at **r = 0.989**. Both estimators decrease with frequency by construction, so *some*
correlation is trivial (a linear ramp already gives r ≈ 0.75 vs Leblanc); what the measured value adds
is agreement in the **log–log curvature** of the density-model track over two decades, not just a shared
monotone trend. The absolute geometric scale runs ~2× high (ratio ≈ 2.18): direction-finding outward
bias on a degrees-wide source, a likely active-region density enhancement, **and radio-wave scattering
shifting the apparent interplanetary source outward of the plasma level** (Krupar et al. 2015) all push
the same way. The robust statement is that the two independent estimators **track each other in shape**,
not that they agree in absolute calibration.

## Honest assessment & caveats

- **A synthesis/reproduction, not a discovery.** Each leg is an existing recover-a-known; the new element
  is the unified framework + the geometric validation of the model distance.
- **Different events.** The corona (2011-09-14), inner-heliosphere (2003-10-28 X17), and interplanetary
  (2013-05-15) legs are *different bursts*, so this is one method validated across four regimes, **not one
  beam tracked end to end**. The stronger demonstration — a single event caught simultaneously by
  e-Callisto + Wind + STEREO — is the natural next step (an optional GATE-0 event hunt).
- **The cross-check is event-specific and shape-only.** It validates the 2013-05-15 interplanetary
  distances; the absolute geometric scale is direction-finding–biased, so it checks the distance *shape*,
  not the calibration.
- **Inherited per-leg systematics.** Harmonic vs fundamental (factor 2 in distance), peak-time vs onset
  speed bias, the average density model, and the coarse-cadence fit statistics are all bounded in the
  component slices and carried through here.
- **Two stitched density models, not one.** The Newkirk (corona) and Leblanc (heliosphere) models are
  not designed to join continuously; at the ~15–20 MHz handoff they differ in heliocentric radius by
  ~50% (white-light active corona vs Wind type III drifts). The ladder is two stitched regimes.
- **Not the first geometric check.** Prior STEREO direction-finding work (e.g. Krupar et al.) has
  compared triangulated positions to density-model predictions; the new element here is the *same-event,
  reproducible, open-pipeline* comparison, not direction-finding as a density-model check per se.
- **Reproducible:** `python -m jansky_research.type3synthesis` regenerates the unified ladder + the
  2013-05-15 cross-check figure and the macros; `make reproduce` runs the four recover-a-known events on
  the public archives (e-Callisto, SPDF Wind/WAVES + STEREO/WAVES + STEREO L3 DF).
