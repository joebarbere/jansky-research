# Findings — Voyager 2 PRA re-derivation of the Uranus & Neptune radio rotation periods (plan 46, F9)

`jansky_research.vgpra`: a blind, modern periodogram re-derivation of the ice-giant radio rotation
periods from the open Voyager 2 PRA encounter data, reusing the merged `frbperiod` Rayleigh-Z²
machinery. Lamy+2025 (Nat. Astron. 9, 658) moved Uranus's period by 28 s using HST UV aurora; the
1986/1989 **radio** values (Uranus 17.24±0.01 h, Neptune 16.11 h) were never reanalysed with modern
statistics, and the Cecconi+2017 PRA refurbishing covers only Jupiter/Saturn.

## GATE 0 (2026-07-10)

- **Dataset IDs pinned + egress verified** (curl 200): `VG2-U-PRA-3-RDR-LOWBAND-6SEC-V1.0`
  (`DATA/VG2_URN_PRA_6SEC.TAB`, 49 MB) and `VG2-N-PRA-3-RDR-LOWBAND-6SEC-V1.0`
  (`DATA/VG2_NEP_PRA_6SEC.TAB`, 79 MB) on PDS-PPI, direct HTTP.
- **Format** (from the .LBL + real bytes): fixed-width ASCII, **line-based** (2284 data bytes + LF;
  NOT the 2286 `RECORD_BYTES` stride — a trap). Each row = one 48-s major frame: `DATE` (YYMMDD) +
  `SECOND` (of day) + 8 sweeps × (1 status word + 70 channels) in **millibell** (=0.01 dB); sweep k
  starts at `SECOND + 6k` s. 70 low-band channels `f_i = 1326.0 − 19.2·i` kHz. Neptune 36028 rows ≈
  21 d; Uranus ≈ 13 d — tens of rotations. First 8 Neptune rows vendored at
  `tests/data/vg2_nep_pra_6sec_sample.tab` for the parser test.
- **Novelty PASS**: Lamy+2025 is HST-UV (not radio); arXiv:2604.19863 is radio *occultations*
  (atmospheric geometry); Cecconi+2017 = Jupiter/Saturn only. No modern radio PRA re-derivation of
  the U/N periods exists.

## Method (and a real-data method upgrade)

- **Flux series**: sum linear power (`10**(mB/1000)`) over 100–1000 kHz per 6-s sweep.
- **First attempt — burst-epoch Rayleigh** (`detect_bursts` + `frbperiod` fold): a NULL on the real
  data (Z² ≈ 0.5–3, i.e. noise level). Not a real null — a *method* limit: discrete burst epochs
  throw away the continuous modulation. (Synthetic recover-a-known passed, so the machinery is
  sound; the real emission just isn't clean phase-clustered bursts.)
- **Upgrade — Lomb-Scargle on the continuous flux** (`flux_period_posterior`, PRIMARY): the rotation
  modulates a continuous, red-noise-dominated emission, so LS on the (log, 0.1-h-binned,
  quadratically-detrended) flux is far more sensitive. The analytic LS FAP is **meaningless** here
  (autocorrelated samples + huge N → absurdly small), so uncertainty = a **rotation-block
  bootstrap** (resample whole rotation-length blocks). Burst-Rayleigh kept as a secondary
  cross-check.
- **Band-stability** (`band_stability`): an ACHROMATICITY check — a coherent signal (rotation OR a
  geometry/beaming/window modulation) is achromatic, so a small sub-band spread confirms the peak is
  a real coherent signal, NOT noise. It does **not** separate right from wrong (GATE-2 catch: on the
  real data Uranus's *wrong* peak is *more* band-stable than the truth); the right/wrong verdict
  rests on comparison to the historical prior.
- **Synthetic recover-a-known** (offline, CI): searched in the SAME wide 14–20 h window used on the
  real data, both methods recover an injected 17.24 h period (LS 17.22 h, Rayleigh 17.24 h). This is
  the control that makes the real null interpretable.

## The real reanalysis (ran 2026-07-10): a CONTROLLED NULL — NEITHER period recovered

Blind LS over an un-tuned wide 14–20 h window on both encounter volumes (Uranus 179,688 spectra /
13 d; Neptune 288,224 / 21 d):

| planet | blind LS peak | combined σ | sub-band spread | offset vs hist | recovers? |
|---|---|---|---|---|---|
| **Uranus**  | 18.44 h | 1.94 h | **1.76 h** (wanders) | +1.20 h | **NO** |
| **Neptune** | **20.0 h** (railed to bound) | 1.94 h | 0.70 h | +3.89 h | **NO** |

- **The synthetic recovers 17.22 h in the same wide window; the real data recover neither period.**
  So the failure is a **data limitation**, not a pipeline bug (the typeii control pattern).
- **Uranus** peaks at 18.44 h, offset +1.2 h from 17.24 h and NOT band-robust (sub-band peaks span
  18.4/18.0/**14.3** h). Consistent with 17.24 h only within the vacuous ±1.9 h band.
- **Neptune** rails to the 20 h search bound — the signature of a periodogram dominated by the
  red-noise / flyby-envelope continuum, not a rotational line. A *narrow* window can place the peak
  near 16 h, but that only shows the result is window-dependent; with an un-tuned wide window there
  is no stable rotational peak. (The earlier narrow-window "16.04 h ≈ 16.11 h" was a window
  artifact — corrected after the widen.)
- **Recovery criterion** (pipeline-generated `recovers_hist`): within 1 combined-σ AND band-robust
  (spread < 0.6 h) AND not railed. Neither planet qualifies.
- **Precision headline** (planned from day one): the blind uncertainties (~2 h) are **hundreds of
  times coarser** than the 28-s HST shift — even a success could never approach modern precision.

**Honest bottom line**: a blind Lomb-Scargle of the Voyager PRA total-power flux recovers a clean
injected rotation in synthetic tests but **recovers neither real ice-giant period** — the auroral
total-power flux is not a clean rotational sinusoid over these short flybys, so the historical
determinations' beaming/magnetic-longitude modelling was essential, and no naive periodogram on this
data can approach the modern HST precision. A reproducible modern-statistical demonstration of *why*
the original methods were needed.

## GATE-2 (2026-07-10) — PASS on honest framing, with required fixes (all applied)

The review confirmed the coarseness/precision message is delivered well, and caught real problems,
all fixed:

- **The band-stability narrative was inverted** (the key catch). I had called achromaticity "the
  load-bearing honesty gate" that tells a rotation measurement from a wrong-period peak — but on the
  real data Uranus's *wrong* peak (spread 0.17 h) is *more* band-stable than Neptune's; achromaticity
  never even bound `total_unc`. Reframed: `band_stability` is a coherence check (signal vs noise),
  NOT a right/wrong discriminator; the verdict rests on the historical prior.
- **Window-dependence exposed**: widening the search window (as the reviewer's edge-railing point
  prompted) moved the Neptune peak from 16.04 h to 14.9→20 h. The earlier "Neptune recovered" was a
  **narrow-window artifact**. Switched to one un-tuned wide 14–20 h window for both planets; result
  is now correctly a null for both. `recovers_hist` now requires proximity AND band-robustness AND
  not-railed.
- **"First independent radio confirmation" overclaim** → it is the same (only) Voyager dataset;
  changed to "modern blind re-analysis," and the whole framing to a controlled null.
- **Uranus 18.4 h "beaming" asserted too strongly** → softened to a coherent non-rotation modulation
  (geometry/beaming and/or a spectral-window sidelobe).
- **Citations**: Lecacheux+1993 DOI `10.1029/93GL02958`→`10.1029/93GL03117`; Neptune uncertainty
  ±0.02→±0.05 h (Warwick); Lamy pages 658–666→665; Cecconi PRE7→PRE8.
- **Verified good**: millibell→`10**(mB/1000)` conversion; time reconstruction (the same code
  recovers the synthetic and would recover Neptune if the signal were clean, ruling out a time-axis
  bug); dropping the analytic LS FAP under red noise; the rotation-block bootstrap.

## Reproduce

Offline (metric + synthetic recover-a-known + tests + real-sample parser test):
`uv run python -m jansky_research.vgpra --offline --out .`
Real (downloads both PDS-PPI volumes): `uv run python scripts/vgpra_real.py --cache /tmp/vgpra`
(writes `results/vgpra_metrics.json`, `is_real=True`).
