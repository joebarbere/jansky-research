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

## Full referee round (2026-08-24): MAJOR REVISION, 16 findings

The presenter/referee round (paper-presenter + paper-referee agents). All macros resolve, the
legs table matches the JSON, and the arithmetic checks --- the problems are what the evidence
does not contain.

**The five MAJORs:**
1. **The null has no sensitivity.** p=0.933 on the clipped leg means the scramble null's own
   median amplitude exceeds the observed 0.31 --- and the number needed to state an excluded
   amplitude (`null_amps`) is computed and stripped before commit (`rmdipole.py:406`). The
   frblens lesson: "we saw nothing" is a constraint only after dividing by what could be seen.
   Fix: commit null percentiles, restate the headline as "RA-projected dipoles with |p|/m > X
   excluded at 95%".
2. **Table 1's +/- and p contradict each other** (0.311 +/- 0.022 next to p=0.933 is a "14-sigma"
   amplitude 93% consistent with null): the bootstrap SE is within-realization fit precision, the
   rmstructure shape. State which column carries significance.
3. **The tail clip was never tested against a genuine dipole** and 0.99 is a single unswept
   choice; a variance dipole preferentially populates the clipped tail by construction. The one
   unit test injects 40-sigma contaminants and cannot fail. Fix: clip the injection leg; sweep
   0.95/0.98/0.99/0.995.
4. **The injection control is single-seed** (0.2812 +/- 0.0087 vs injected 0.3 is 2.2 SE low,
   explained by assertion), Gaussian where the real field is heavy-tailed --- amp 0.28 on the
   Gaussian injection gives p=0.001 while amp 0.31 on the real clipped field gives p=0.933 at the
   same positions, so the injection says nothing about detectability in the real residuals.
5. **"Isotropic at dipole order" is unscoped in the abstract** though the test is blind to the
   Dec-projected component (the Methods say so; the headline does not).

Eleven MINOR/NIT: the ~6-deg direction-recovery claim traces to a function default, not committed
evidence; the nn path's own |b| amplitude rise goes unremarked while the same signature indicts
the latitude path; the tail's apex is quoted but never fit; n_scramble is not in the JSON; four
bib entries cite preprints now published (incl. the DR2 data citation) and boehme2025/secrest2025/
malik2026 lack DOIs; the kinematic "two orders of magnitude" is a factor 2.4 in the fitted
statistic; the committed figure plots the rejected leg and is not in the paper.

**Status: fixes pending** (data local: data/spice-racs.dr2.fits).

### Resolved 2026-08-24 (revision): the null becomes a limit, and every diagnostic came back clean

All 16 findings addressed; the real leg re-run with 999 scrambles and eight new legs. Purely
additive on the original eleven legs.

**The headline is now an exclusion limit.** The scramble null's amplitude percentiles are
committed per leg (the old code stripped them): the clipped-core null's 95th percentile is
**0.354**, so the paper states "RA-projected power-dipole amplitudes |p|/m > 0.354 excluded;
nothing below that constrained" instead of "isotropic". The RA-projection scope is now in the
abstract, not only the Methods.

**Every referee-demanded diagnostic resolved in the paper's favour, and is now evidence:**
- *Clip circularity*: the same 0.3 injection pushed through the 0.99 clip is attenuated to 0.223
  but still detected at p=0.001 -- the clip cannot manufacture the isotropic core from a real
  dipole. The quantile is swept 0.95-0.995: p = 0.78-0.93 throughout.
- *Single-seed injection*: twenty realizations give **0.3006 +/- 0.0091 (bias +0.0006)** -- the
  injection is unbiased, the seed scatter equals the bootstrap SE, and the old "partial-sky
  bias" explanation of the one low draw (0.2812) was wrong; the ensemble refutes it and the
  paper now says so.
- *Real-field detectability*: a 0.3 dipole painted onto the real heavy-tailed residuals
  (sign-randomised first) is detected at p=0.001, though its recovered amplitude (0.556) shows
  the amplitude scale is uncalibrated on the real field -- conclusions are therefore stated in
  detection/exclusion terms.
- *The tail's own apex*: fitted alone, the 1029 clipped sources have amp 0.52, p=0.002, apex
  **(25 deg, +65 deg)** -- beyond the survey's Dec +49 edge, 65 deg from the full-sample apex the
  paper had implicitly assigned them, supporting the leakage/edge-systematics reading.
- *nn-path leakage*: the neighbour-count requirement is inert (N identical at nn_min 3/5/10),
  now stated as a property of the catalogue rather than a cut (finding 14 resolved: the >=5
  mask does nothing).

Bib: the four now-published preprints completed (incl. the DR2 data citation, PASA 43 e089) and
three missing DOIs added. The kinematic expectation restated in the fitted statistic (4beta).
The apex range corrected to 80-153 deg; the look-elsewhere clause now says the 17% is
conservative in the null's favour; n_scramble and seed committed per leg; the figure plots the
conclusion-carrying clipped leg and its null. One mechanical LaTeX lesson: macro names cannot
carry digits (\rmdRealClipNull95 silently broke the preamble; renamed NullNinetyFive).
