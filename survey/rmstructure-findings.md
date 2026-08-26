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

**GATE-2 delta upgrade**: the survey's own `goodRM_flag` is now applied — **333,173 components**
(matches the DR2 paper's pre-dedup goodRM count to ONE row, an independent filter validation;
field-overlap duplicates sit below the smallest SF bin). The unflagged 337,548 sample was
leakage-contaminated; corrected numbers below (from the local DAP file):

| quantity | DR2 value |
|---|---|
| plane enhancement ratio (|b|<10 / >60) | **11.17 ± 0.10** (statistical only; Taylor09 ~5.4 — NOT apples-to-apples: NVSS nπ caps plane |RM|, southern sky adds inner-Galaxy sightlines; both push SPICE-RACS higher) |
| SF plateau, disc (|b|<10°) | **62,065 rad²m⁻⁴** (upper bound; secular plane gradients mixed in) |
| SF plateau, halo (|b|>10°) | **2,284 rad²m⁻⁴** |
| disc–halo fluctuation-power contrast | **~27×** |
| break scales | 2.29° (disc) vs 3.7° (halo) — the unflagged sample's 0.5° halo break was a LEAKAGE ARTIFACT |

### Latitude ladder (six |b| bins, the resolved profile)

| |b| | n | σ_RM (rad/m²) | σ_Gal (floor-subtracted) |
|---|---|---|---|
| 0–5° | 22,220 | **214.7** | 214.4 |
| 5–10° | 27,089 | 104.9 | 104.2 |
| 10–20° | 56,801 | 66.3 | 65.2 |
| 20–30° | 53,936 | 40.0 | 38.2 |
| 30–50° | 91,197 | 20.2 | 16.3 |
| 50–90° | 81,930 | 11.9 | (floor bin: ≡0; true polar Galactic term → lower bound) |

Monotonic ×18 fall; the polar endpoint sits at the literature intrinsic+extragalactic floor
(Mao+2010 ~9–15), licensing the quadrature floor subtraction (DEFROST-lite; the polar bin as
floor estimate). The intrinsic-scatter floor is latitude-independent, so the disc–halo DIFFERENCE
(~58,000 rad²m⁻⁴) is dominated by the Galactic magneto-ionic medium — the measurement the method paper
staked out. Caveat added: disc sightlines also gain depolarisation-selected populations; a
DEFROST-style separation would tighten the difference argument. Pairs at this scale are drawn by
random sampling (unbiased; fraction recorded). Reproduce: download csiro:64891, gunzip, then
`uv run python -m jansky_research.rmstructure --dr2 --out .`.

## Full referee round (2026-08-24): MAJOR REVISION, 18 findings, two BLOCKERs

The honesty discipline held (verbs, upper-bound framing, macros all resolve; every abstract number
traces). The problems are what the evidence does not contain.

**BLOCKERs:** (1) The abstract's closing claim -- the quality flags "proved material", moving the
high-|b| break from 0.5 to 3.70 deg -- rests on a hard-coded 0.5 whose run was never committed; a
1.3% cut credited with a three-bin shift, asserted with one number and an unevidenced leakage
mechanism. (2) The sample definition contradicts the release's own published counts: DR2's 8-sigma
post-dedup count is 2.5e5 and its 6-sigma count 3.4e5, while this paper's "S/N >= 8" sample is
333,173 -- and the code resolves the S/N column by unrecorded fallback chain, so nobody can say
which column the cut used or whether it bit.

**MAJORs:** the headline 11.17 +/- 0.10 is an i.i.d. source bootstrap on a field the paper itself
shows is spatially correlated (coherence 2.29 deg -> effective N ~patches, error understated ~6x)
and ~25% duplicated -- the slice's own recorded lesson, on the real leg; fix is a spatial block
jackknife. The enhancement ratio is never defined in the paper, and the recover-a-known validates
a DIFFERENT statistic (fixture pole cut |b|>15 vs the science |b|>60). The floor
latitude-independence premise is challenged in the findings file (depolarisation-selected disc
populations) and the caveat never reached the manuscript. The floor value licensing the
subtraction is never quoted and moves sigma_Gal 13.5-18.1 across its 9-15 range. plateau_err,
n_pairs, pair_fraction are computed and discarded, so the six-figure plateaus and the "monotonic"
ladder carry no uncertainties. The 30-seed ensemble exists only as mean/std in macros (per-seed
ratios uncommitted; though 30 seeds IS enough for "threefold": 3.17 +/- 0.42). The SD 1.11 is
quoted where the response claim needs the SEM (10.6 sigma); 3.15-vs-5 is unexplained in the paper
(band-average over a 5-deg profile cannot equal the peak boost). Duplicates are defended for the
SF and undefended for the two statistics that carry the paper.

**MINOR/NIT:** committed metrics predate the current code (missing keys); \rmsSource says
"synthetic" in the file carrying real values (the guard's marker inverted -- namespace it); DR1
numbers evidenced nowhere and thomson2023 uncited at the sentence attributing them; "quadrant
signs" overstates rmsky's own docstring ("conflates" quadrants); the stale arXiv package has a
dangling-citation abstract ("(per-region NVSS structure functions exist; )"); 11.17 +/- 0.1
precision mismatch; README rows stale (says "awaits the public DR2 file").

**Status: RESOLVED (2026-08-24).** The full DR2 leg re-ran with every fix in place.

**Blocker 2 (sample definition), dissolved by one dedup.** The 333,173 goodRM rows contain
tile-overlap repeats: 246,508 unique `cat_id` — exactly the release's published post-dedup
8-sigma count. `load_spice_racs_dr2` now dedups (one row per `cat_id`, keeping the observation
nearest its tile centre, `snr_polint` tie-break) and commits the whole cascade with the S/N
column named: 9,294,225 raw → 338,313 (snr_polint ≥ 8) → 337,548 (finite) → 333,173 (goodRM) →
246,508. The paper's old defence ("duplicates sit below the smallest SF bin") was wrong twice
over: same-source pairs at ~0 separation inject pure noise power into the smallest bins, and
26% of sources were double-weighted in every median.

**The headline error was understated 11x, not ~6x.** On the deduped sample the enhancement
ratio is 10.97 with a leave-one-block-out jackknife (601 blocks of 10 deg) SE of **1.1**; the
i.i.d. source bootstrap says 0.10. Quoted as 11.0 ± 1.1 (display precision matched). The new
test proves the jackknife exceeds the bootstrap on the correlated fixture, so the failure mode
is pinned in CI, not just fixed once.

**Blocker 1 (quality-flag claim), now measured on committed runs.** Unflagged + deduped:
high-|b| half-plateau scale 0.87 deg; flagged: 2.29 deg. The claim survives in direction; the
old numbers (0.5 vs 3.70) were pre-dedup and the flagged one changed materially — which is
exactly why an uncommitted variant run was a blocker.

**The rest:** per-seed ensemble committed (results/rmstructure_synthetic.json, 30 ratios,
allowlisted as synthetic-by-design); ratio defined in the abstract; fixture scope stated (pole
cut 15 vs 60 deg; why the band median must sit below the injected peak boost of 5); ladder rows
carry plateau_err, sigma_rm_err, n_pairs, pair_fraction; the floor is quoted (12.3 ± 0.21) with
its leverage made explicit — plane bins insensitive (201.9 vs 201.1 under a floor nearly twice
as large), the 30–50 deg bin dominated (18.3, spanning 16.1–20.1 across literature floors
9–15); provenance macros namespaced (\rmsRealSource / \rmsSynSource under the live \rmsSource
the guard reads); thomson2023 cited at the DR1 sentence; "quadrant signs" → "coarse net RM sign
per quadrant"; abstract cut to arXiv's 1920-char limit and the arXiv package rebuilt clean
(the dangling-citation abstract is gone); README row updated.

## Full referee round (2026-08-26): MAJOR REVISION, 21 findings, one BLOCKER

The core is sound and unusually honest: the sample cascade is fully committed and consistent
(ladder bins sum to 246,508 exactly), all 29 macros resolve, every ladder row reproduces from
its committed inputs to the last digit, the 30-seed ensemble is committed per-seed, "responds
to" is the right verb, and all five journal citations Crossref-verify. What blocks it: one
abstract claim rests on an uncommitted noise-dominated statistic, and the paper applies its
hardest-won lesson to the headline while leaving the condemned estimator under every other
error bar.

**BLOCKER: the quality-flag claim ("unflagged, the break falls 2.29° → 0.87°, leakage
masquerading as small-scale structure") rests on an 11-valued threshold-crossing statistic**
(bin centres of a 1.62×-per-bin grid), whose underlying curves are neither committed (run()
discards sf/sf_err/n_pairs for this comparison) nor uncertainty-bearing, evaluated on
UNMATCHED random pair sets, on a curve the paper's own figure shows dipping non-monotonically
by 3× at exactly those separations (~160 pairs in the 0.87° bin). The flag removes 1.17% of
rows yet the mechanism claimed is a separation-independent contaminant that moves the plateau
only 4%. Commit both curves on a matched pair set with an uncertainty on the break, or the
claim leaves the abstract.

**MAJORs:**
- EVERY error except the headline uses the bootstrap the paper condemns as 11× too small:
  plateau_err (hence sigma_rm_err for all six bins, incl. the abstract's floor ±0.21) is the
  same i.i.d. source bootstrap. Jackknife the ladder or label the errors shot-noise-only.
- The block jackknife's own tuning is never varied (10° is the only value run) while Honest
  Limits itself says plane gradients span "tens of degrees" — and the (l,b) blocks are not
  equal-area (polar blocks ~1/10 the solid angle, and the poles are the denominator). Sweep
  block size (5–30°) or use HEALPix; quote the largest/stabilised SE.
- THE CI TEST LOCKS THE OUTLIER IN: `assert 3.0 < ratio < 7.0` passes only because seed 0
  (4.64) is the 6th-highest of 30 — 17 of the 30 committed seeds fall OUTSIDE the window (min
  0.95). Assert on the ensemble mean/SD instead.
- The depolarization/goodRM selection-vs-|b| caveat (recorded as unresolved in this file after
  the previous round) is STILL absent from the manuscript, and the data to bound it is one line
  away (run() already loads goodrm=False; commit the unflagged count per ladder bin).
- The statistic's definitional sensitivity is measured in the SIBLING slice and unmeasured
  here: rmsky commits alt-bins 5°/70° moving the same ratio +39% (5.4→7.5), cut300 and e_RM
  variants — none run on DR2, where the quoted error is 10%.
- The only figure mislabels its axis (real leg plots RA vs Galactic b under "lon"), saturates
  at ±40 against a 202 rad/m² plane σ (the headline contrast is invisible), and never plots
  the six-bin ladder that is the paper's result.
- A direct offline run in the repo root silently swaps the real figure for the synthetic one
  (JSON and macros are guarded; the figure is not).

**MINOR/NIT:** the floor bin's σ_Gal=0.0 is a construction recorded as a measurement (mark
null/by-construction; the literature floor range allows a polar Galactic term up to 8.4 rad/m²
— which the quoted 16.1–20.1 span happens to bracket, and the headline ratio has NO floor
subtraction at all — both worth saying); 11.0 vs 202.3/12.3=16.4 sit unreconciled in the
abstract (different statistics, different bins); the Taylor comparison drops rmsky's committed
comparable jackknife (±0.59) and "recover from Taylor" should be "measure from"; stil2011 is
in refs.bib and the findings file ("not the first ever") but uncited in the paper; ~330 of
"601 blocks" contribute zero (report the effective count); prose cascade omits the committed
n_after_finite step and "26% duplicates" uses the rows denominator (35.2% of unique); two
Honest-Limits macros + the DR1 counts (7,707/5,818) are in no committed JSON; the "release's
own post-dedup count" validation cites no table (and this file previously attached the same
rhetoric to a DIFFERENT number); the fixture plateau ratio (12.8) is referenced but never
printed; true_coherence_deg is a bare NaN literal (RFC 8259); thomson2023 pages=e040 is an
article-number (the lawrance2024 shape); stale arxiv tarball predates the artifacts; the
figure caption doesn't say which leg is printed.

**Positive findings recorded (do not re-run):** the SF coherence-scale test is NOT vacuous;
preserve guards are wired for JSON+macros; the relative jackknife error matching rmsky's
(10% vs 11% on a 6.6× smaller sample) is the signature of field-dominated error — say it; and
the paper's single best unused asset is already committed: the fixture block-jackknife SE
(1.1) reproduces the true 30-seed field-to-field scatter (1.114) to 1% while the bootstrap
understates it 3.7× — the genuine recover-a-known for the ERROR MODEL that licenses ±1.1 on
the real leg. Lead with those three numbers.

**Verdict: MAJOR REVISION.**
