# Findings — Jovian DAM occurrence census from Juno/Waves (plan 37)

`jansky_research.junodam`: an occurrence census over the (CML, Io-phase) plane from the public
Juno/Waves calibrated dataset (doi:10.25935/6jg4-mk86; daily CDFs, 110 chan to 40.5 MHz with
shipped per-channel Background/Sigma), using Horizons for the **sub-Juno** System III CML (the
naive IAU W_III formula errs by up to ~40° at Juno — GATE-0 verified) and Lieske (1987) for Io.

## Recover-a-known (synthetic orbit)

Real CML/Io rates + Bernoulli activity in the canonical Io-A/B/C/D boxes: injected contrast 8.75,
recovered **7.2** (cell-edge smearing accounts for the gap); 324/324 plane cells covered in a
synthetic month.

## Real month (2017 February, 27 daily CDFs ≈ 1 GB; one FAIL day of 28)

| quantity | value |
|---|---|
| 15-s bins / active | 155,520 / 6,244 (4.0%) |
| **near-half vs far-half duty cycle** | **7.72% vs 0.31% — proximity dominates by ~25×** |
| Io-box contrast (full month, CORRECT Φ_Io convention) | **2.22** |
| far-half Io contrast | **1.55** |

**GATE-2 caught a real convention blocker**: Io phase was computed as Λ_Io−CML instead of the
standard Φ_Io = CML+180°−Λ_Io — displacing every canonical box (except Io-B by coincidence) and
faking a weak 1.38 contrast. With the correct convention (+ the Io-C CML wrap through 0°):
contrast **2.22** full-month, **1.55** far-half. Reading: Io organisation genuinely persists
from orbit, reduced vs ground-based — consistent with part of the canonical structure encoding
observer beaming. Lesson recorded: **a synthetic round-trip cannot validate frame conventions**
(it injected and recovered in the same wrong frame). Multi-orbit v02 census is the follow-on.

## Caveats

- One orbit; uneven (CML, phase) exposure (per-cell exposure carried; <3-visit cells masked).
- 10%-of-band activity criterion is one choice (suppresses narrowband RFI, dilutes narrowband DAM).
- No Io/non-Io separation; no beaming/proximity modelling — demonstrated, not modelled.
- Units are V²m⁻²Hz⁻¹ (occurrence needs only detection vs the shipped background).
- Reproduce: download a month of CDFs to `data/junodam/` (URL pattern in module), then
  `uv run python -m jansky_research.junodam --out .` (Horizons access needed for CML).
