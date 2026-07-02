# Findings — Jovian DAM occurrence census from Juno/Waves (plan 37)

`jansky_research.junodam`: an occurrence census over the (CML, Io-phase) plane from the public
Juno/Waves calibrated dataset (doi:10.25935/6jg4-mk86; daily CDFs, 110 chan to 40.5 MHz with
shipped per-channel Background/Sigma), using Horizons for the **sub-Juno** System III CML (the
naive IAU W_III formula errs by up to ~40° at Juno — GATE-0 verified) and Lieske (1987) for Io.

## Recover-a-known (synthetic orbit)

Real CML/Io rates + Bernoulli activity in the canonical Io-A/B/C/D boxes: injected contrast 8.75,
recovered **7.2** (cell-edge smearing accounts for the gap); 324/324 plane cells covered in a
synthetic month.

## Multi-orbit real leg (7 months / 210 days, 2016–2019, v01+v02, 8.6 GB local)

| quantity | value |
|---|---|
| 15-s bins / active | 1,209,600 / 44,294 (3.7%) |
| **duty cycle by range quartile** | **12.9% → 1.5% → 0.20% → 0.07% — proximity dominates by ~180×** |
| Io-box contrast, aggregate | **1.12** |
| per-month contrasts (7 months) | 1.56, **2.22**, 0.87, 0.35, 0.93, 0.70, 0.84 — median 0.87, both sides of 1 |
| distance-resolved contrast (near→far) | 1.24 / 0.36 / 1.37 / 0.71 (far quartiles activity-starved) |
| Io-box contrast (full month, CORRECT Φ_Io convention) | **2.22** |
| far-half Io contrast | **1.55** |

**GATE-2 caught a real convention blocker**: Io phase was computed as Λ_Io−CML instead of the
standard Φ_Io = CML+180°−Λ_Io — displacing every canonical box (except Io-B by coincidence) and
faking a weak 1.38 contrast. With the correct convention the single Feb-2017 month gave contrast 2.22 — but the multi-orbit
extension shows that was the high tail of month-to-month scatter: **across 7 months the
canonical boxes do not coherently organise Juno-frame occurrence** (median 0.87, aggregate
1.12, no strengthening with distance). One orbit can mislead; the census needed the spread.
Lesson recorded: **a synthetic round-trip cannot validate frame conventions** (it injected and
recovered in the same wrong frame). A beaming-model fit over the full v02 archive is the
definitive follow-on.

## Caveats

- One orbit; uneven (CML, phase) exposure (per-cell exposure carried; <3-visit cells masked).
- 10%-of-band activity criterion is one choice (suppresses narrowband RFI, dilutes narrowband DAM).
- No Io/non-Io separation; no beaming/proximity modelling — demonstrated, not modelled.
- Units are V²m⁻²Hz⁻¹ (occurrence needs only detection vs the shipped background).
- Reproduce: download a month of CDFs to `data/junodam/` (URL pattern in module), then
  `uv run python -m jansky_research.junodam --out .` (Horizons access needed for CML).
