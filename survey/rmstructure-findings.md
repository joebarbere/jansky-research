# Findings — RM structure functions from SPICE-RACS (plan 36, first leg)

`jansky_research.rmstructure` extends `rmsky` (Taylor+2009, 37k RMs) to the SPICE-RACS grids with
the analysis the DR2 release paper (arXiv:2605.16917) did not do: **noise-debiased second-order RM
structure functions per Galactic-latitude bin** (SF = ⟨ΔRM²⟩ − ⟨σ₁²+σ₂²⟩, Haverkorn+2004
convention; source-level bootstrap errors; recorded pair subsampling).

## GATE 0 (verified live 2026-07-02)

- **DR2 IS public**: `spice-racs.dr2.fits.gz` (4.97 GB) on CSIRO DAP collection csiro:64891, no
  auth (48-h presigned S3 URLs). Not on CASDA TAP or VizieR yet.
- **DR1 on CASDA TAP**: `AS110.spice_racs_dr1_corrected_cut_v02` = **24,758 spectral-component rows** (live count; the
  polarimetric RM detections are a subset — ~5.8k in the DR1 paper, 7,707 at our S/N≥8 cut), columns incl. `l, b, rm, rm_err, snr_polint`, no auth.

## Recover-a-known (synthetic screen)

- Pure-noise sky debiases to SF≈0 (not the 2σ² floor) — the debiasing works.
- Injected 2° coherence + 5× plane amplitude boost: enhancement ratio 4.64±0.35 recovered;
  low-|b| SF plateau ≫ high-|b| (6526 vs 511 rad²m⁻⁴); half-plateau break 3.7° ≈ injected 2° ×
  2√(ln2)≈1.67 (theory 3.33°; the log-binned estimator snaps to the nearest bin, ±30%).

## First real leg (DR1, 7,707 RMs at S/N≥8)

| quantity | value |
|---|---|
| high-|b| SF plateau | **5,755 rad²m⁻⁴** → RM dispersion ~54 rad m⁻² |
| half-plateau scale | **0.54°** (bin-limited; may reflect S/N-cut clustering or tile geometry — null test needed) |
| low-|b| column | **empty — DR1's corrected-cut excludes the Galactic plane** |

The plane exclusion is the honest headline limitation: the disc–halo SF contrast the method is
built for is validated on synthetics only until the DR2 run (a 5 GB download + the same code).
The plateau mixes Galactic power with intrinsic+extragalactic RM scatter → **upper bound only**:
literature high-|b| Galactic dispersion is ~9–15 rad/m² (Mao+2010, Taylor+2009), so the 54 rad/m²
is intrinsic-scatter-dominated. Prior art: Stil+2011 did NVSS per-region SFs; ours is the first
from SPICE-RACS, not the first ever.

## DR2 full-sky leg (DONE 2026-07-02)

337,548 RMs (S/N≥8) from the local DAP file (`data/spice-racs.dr2.fits`, gunzipped for memmap):

| quantity | DR2 value |
|---|---|
| plane enhancement ratio (|b|<10 / >60) | **11.11 ± 0.09** (Taylor 2009: ~5.4 from 37k RMs) |
| SF plateau, disc (|b|<10°) | **60,846 rad²m⁻⁴** |
| SF plateau, halo (|b|>10°) | **2,612 rad²m⁻⁴** (→ σ≈36 rad/m², still intrinsic-dominated) |
| disc–halo fluctuation-power contrast | **~23×** |
| break scales | 1.41° (disc) vs 0.87° (halo) |

The intrinsic-scatter floor is latitude-independent, so the disc–halo DIFFERENCE (~58,000
rad²m⁻⁴) is dominated by the Galactic magneto-ionic medium — the measurement the method paper
staked out. Caveat added: disc sightlines also gain depolarisation-selected populations; a
DEFROST-style separation would tighten the difference argument. Pairs at this scale are drawn by
random sampling (unbiased; fraction recorded). Reproduce: download csiro:64891, gunzip, then
`uv run python -m jansky_research.rmstructure --dr2 --out .`.
