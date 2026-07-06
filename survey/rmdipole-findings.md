# Findings — the first RM dipole/isotropy test (plan 38)

`jansky_research.rmdipole` asks a question no published paper has asked: **is the extragalactic
rotation-measure sky statistically isotropic at dipole order?** Every cosmic-dipole-anomaly test
to date is source-count/flux based (Böhme+ arXiv:2509.16732; arXiv:2509.18689; RMP colloquium
arXiv:2505.23526; Mittal & Lewis arXiv:2605.27520); SPICE-RACS DR2 (arXiv:2605.16917) is the
first RM catalogue big enough to change that.

## GATE 0 (full-text pass, 2026-07-05 — the fable-ideas scan ran egress-blocked)

- All five load-bearing arXiv IDs verified against their abstract/full-text pages; none does an
  RM dipole. The DR2 release paper itself contains no dipole/isotropy analysis ("dipole",
  "anisotropy", "isotropy" absent from its full text).
- Kill-condition sweep: **no published or preprint catalogue-scale RM dipole/isotropy test
  exists.** Closest prior art (cite, different question): uniform-cosmological-B-field RM dipole
  — Kronberg & Simard-Normandin 1976, Vallée 1990, Kolatt 1998, all at ~10² RMs; Mao+2010
  (Galactic vertical field toward the poles); Mtchedlidze+ arXiv:2511.19508 (simulation only).
- **Kinematic expectation for RM** (GATE-2 correction of the GATE-0 claim "none exists"): an
  observer boost Doppler-rescales RM_obs = RM(1+β·n̂)² — a signed dipole of amplitude
  2β ≈ 2.5×10⁻³, two orders below this test's sensitivity and largely absorbed by local
  subtraction. No *paper* derives it, but the physics is trivial. Framed as an isotropy test;
  a detection at measurable amplitude could not be kinematic.
- Data verified on disk: `data/spice-racs.dr2.fits` = 9,294,225 components × 125 columns incl.
  `rm, rm_err, l, b, goodRM_flag, nn_rm_med, nn_rm_count` (DAP csiro:64891).

## Design (deltas from plan 38)

- **healpy dropped**: the dipole is fit per-source (LSQ on the monopole+dipole design matrix),
  no pixelisation, no new dependency. Binned maps appear only in the figure.
- **Residuals**: DR2's own `nn_rm_med` nearest-neighbour GRM subtraction (the Malik+2026
  arXiv:2605.16924 convention) is the primary path; 5°-|b|-band median subtraction is the
  cross-check. Local subtraction absorbs any *signed* mean-RM dipole by construction → the
  tested statistics are dipoles in residual **power** (rm²−σ², noise-debiased, primary) and
  **|residual|** (robust companion). A Vallée-style uniform-field limit needs model-based
  subtraction — deferred, noted in the paper.
- Significance from **footprint-preserving scrambles** (permute residuals among real positions
  within 5° Dec bands): preserves the Dec≤+49° footprint and Dec-dependent systematics exactly,
  destroys all RA structure. 999 scrambles → p-floor 10⁻³.

## Recover-a-known

- Synthetic (CI, offline): injected power dipole A=0.4 on a DR2-like footprint recovered in
  amplitude and direction (<15° apex error); no-dipole control stays null; pure-noise power
  debiases to σ₀² exactly; a σ-map dipole is found by the noise stat and NOT by the debiased
  power stat; piled-up outliers are caught by the clip diagnostic.
- **On the real footprint**: injected A=0.3 at the 102,830 real |b|≥45 positions → recovered
  0.281±0.009, apex 5.7° from injected, p=0.001. The pipeline sees a real dipole if one is there.

## Result (999 scrambles; full table in `results/rmdipole_metrics.json` + the paper)

| leg | amp | apex (RA,Dec) | sep. CMB | p |
|---|---|---|---|---|
| power, nn, \|b\|≥45 (primary) | 0.529±0.065 | (72°, +9°) | 97° | 0.001 |
| **power, nn, \|b\|≥45, clip top-1% \|r\|** | 0.311±0.022 | (86°, −56°) | 80° | **0.933** |
| abs, nn, \|b\|≥30 (widest, robust) | 0.142±0.006 | (56°, −45°) | 100° | **0.956** |
| noise (σ²) dipole, \|b\|≥45 | 0.037±0.006 | (232°, −37°) | 65° | 0.001 |
| latitude-band variants | 0.26–0.77 | various | 100–153° | 0.001 |

**Headline (honest):** the extragalactic RM sky is **isotropic at dipole order in its
distributional core**. The formally significant full-sample power dipole is carried *entirely*
by the top 1% of |residual| — clipping 1,029 of 102,830 sources moves p from 0.001 to 0.93.
The one marginal robust-stat entry (abs, |b|≥45 nn, p=0.023) is disclosed in the paper: no
trials correction across 8 variants (~15–20% chance of one such p), apex 82° from the CMB, no
persistence at |b|≥30 — read as an echo of the same tail. No sky-statistic variant's apex is
closer than 80° to the CMB dipole apex: **no RM counterpart to the source-count dipole
excess.** The latitude-band subtraction yields inflated, unstable
dipoles — the expected signature of longitude-dependent Galactic residue, and the reason the
nn subtraction is primary. The noise map has its own small (3.7%) dipole pointing elsewhere —
survey depth non-uniformity, excluded as the power-dipole carrier.

## Interpretation limits (GATE-2 material)

- A positive RM-power dipole would have had no unique cosmological reading (no kinematic
  expectation); the value here is the isotropy bound + the honest null machinery.
- The ~10³ tail sources that carry the full-sample signal are unidentified: plausible carriers
  are imperfect local GRM subtraction, leakage survivors, and nπ-ambiguity outliers. Identifying
  them (tile clustering? bright-source proximity? `rm_width` morphology?) is the natural
  follow-on and would double as a DR2 quality diagnostic.
- The scramble preserves Dec-band means, so it tests the **RA-projected** dipole only:
  conservative for a general dipole, blind to a polar-axis-aligned one (sensitivity to a
  CMB-aligned dipole is nearly full because the CMB apex sits at Dec −6.9°). Stated in the
  paper (GATE-2 required fix).
- The permutation null also assumes within-band exchangeability; the measured degree-scale RM
  correlations (rmstructure) violate it, so small p-values strictly read "RA structure exists,"
  not "a dipole exists" — consistent with (and supporting) the systematics reading.
- An RA-dependent systematic aligned with a dipole could in principle masquerade. The tail-clip
  result makes this moot for the core (nothing survives to explain), but it caps how strongly
  the tail anisotropy itself can be interpreted.

## Reproduce

`uv run python -m jansky_research.rmdipole --n-scramble 999 --out .` (needs the local 9.3 GB
DR2 FITS; ~10 min CPU). Offline CI leg: `--offline`. Everything in the paper flows from
`results/rmdipole_metrics.json` → `generated/macros.tex` + `generated/legs_table.tex`.
