# Findings — RACS Stokes-V discovery: two-epoch forced photometry of the nearest M dwarfs

`jansky_research.stokesv_discovery` + `scripts/stokesv_discovery_real.py` deliver plan 33: forced
I+V photometry of the nearest CNS5 M dwarfs at TWO RACS-mid epochs (MJD 59233 / 60769 — a 4.2-yr
same-band pair at 1367.5 MHz, found at GATE 0 after RACS-low1 V images turned out not to exist),
with Gaia PM propagation to each epoch, per-epoch leakage floors, and signed two-epoch variability.

## Recover-a-known (synthetic epoch pair)

Selection completeness/purity 1.0/1.0; flare-flagging completeness/purity 1.0/1.0; and the
headline rationale: **22.9% of injected emitters are invisible to any single epoch**
(`single_epoch_miss_frac`) — the quantified version of the stokesv paper's variability limit.

## Real sky (60 nearest CNS5 M dwarfs, 115-min CASDA run)

| quantity | value |
|---|---|
| targets measured (≥1 epoch) | 39 (CASDA staging failures account for the rest) |
| complete epoch pairs | 19 |
| ≥5σ V detections | **1 system: GJ 65 (BL+UV Ceti)** — V=9.26±0.15 (2021) → 7.11±0.15 mJy (2025) |
| inter-epoch ΔV significance | **10σ** at 1367.5 MHz over 4.2 yr |
| everything else | quiescent; median 5σ V limit **0.83 mJy** |

GJ 65 is a *recovery* (it is in the RACS-low2 Paper VIII blind V catalogue —
RACS-LOW2 J013906.5−175647, 4″ away — and in SRSC): the pipeline finds the prototype coherent
emitter and adds the mid-band two-epoch V change that blind single-epoch catalogues don't provide.
No candidate survived the novelty bar (checked against Paper VIII/SRSC/SIMBAD before "new" could
be used). PM propagation validated live: Barnard's star moved 39.8″ between the epochs.

## Honest caveats

- Two epochs bound variability, not timescale (flare vs secular).
- Staging availability limits the census: 35% of targets got no data, ~half of measured targets
  lack the second epoch (~52% of target-epoch slots failed; rows kept with notes; resumable).
- Both GJ 65 snapshots may be burst states (bursts last minutes-hours vs 15-min integrations);
  its quiescent V floor is unconstrained here.
- Leakage floor is per-field; a beam-position-dependent model would sharpen faint candidates.
- Reproduce: `uv run python scripts/stokesv_discovery_real.py` (needs `CASDA_USERNAME` +
  `~/.casda_pw`; ~2 h; resumable) then `uv run python -m jansky_research.stokesv_discovery --out .`.
