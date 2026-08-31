# Findings — Milky Way HI rotation curve (tangent-point method)

`jansky_research.hi.run` builds the inner-Galaxy rotation curve by the tangent-point method from the
Leiden/Argentine/Bonn (LAB) HI 21 cm survey (Kalberla et al. 2005), fetching the $(b, v)$ slices at
$\ell = 10°$–$80°$ from VizieR (VIII/76), reading the $b=0$ spectrum, extracting the terminal
velocity, and applying $R = R_0\sin\ell$, $V(R) = v_\mathrm{term} + V_0\sin\ell$. We use the modern
Reid et al. (2019) Galactic constants $R_0 = 8.15$ kpc and $V_0 = 236$ km/s (their best-fit circular
speed at the Sun, not the older IAU 220 km/s). LAB velocities are already in the LSR frame, so no
solar-motion correction is applied.

## Result: an approximately flat (non-Keplerian) rotation curve

| $\ell$ (°) | $R$ (kpc) | $V$ (km/s) | note |
|----|------|------|------|
| 10 | 1.42 | 228 | bar region (R<4): unreliable |
| 20 | 2.79 | 220 | bar region |
| 30 | 4.07 | 247 | |
| 40 | 5.24 | 254 | |
| 50 | 6.24 | 264 | |
| 60 | 7.06 | 260 | |
| 70 | 7.66 | 252 | |
| 80 | 8.03 | 263 | |

Beyond the bar-dominated inner region, from $R \approx 4$ to $8$ kpc the rotation speed stays at
**$V \approx 257 \pm 6$ km/s** — it does **not** fall off as the Keplerian $V \propto R^{-1/2}$ that
the visible mass alone would give. (There is a gentle rise from $R\sim3$ to $\sim6$ kpc, likely
non-circular streaming, on top of a broadly flat curve.) Recovering this **approximately flat,
non-Keplerian** curve from public HI data with a ~30-line tool is the **validation** — a known,
foundational result, not a new measurement (`paper/figures/rotation_curve.pdf`).

## Honest limitations

- **Absolute level is ~9% high** ($257$ vs $V_0 = 236$ km/s). The main cause is the
  **terminal-velocity estimator**: a fixed 2 K brightness threshold *overestimates* $v_\mathrm{term}$
  relative to inflection-point spectral fitting (McClure-Griffiths & Dickey 2016 find threshold
  crossings $\sim7$ km/s higher), which inflates $V$; the rest is the simple estimator and
  non-circular motions. The **shape (flat, non-Keplerian) is robust; the absolute level is the soft
  part.**
- **Circular-orbit assumption.** $V = v_\mathrm{term} + V_0\sin\ell$ assumes the tangent-point gas is
  on a circular orbit; spiral-arm / bar-driven streaming biases individual points (not just adds
  scatter) and cannot be removed without a model.
- **Bar region excluded.** Points at $R \lesssim 4$ kpc ($\ell = 10°, 20°$) are dominated by the
  Galactic bar's non-circular motions and are flagged unreliable and dropped from the mean.
- **Not, by itself, a dark-matter detection.** An inner-Galaxy flat curve is *consistent* with a dark
  halo but can also arise from the disk's own mass; the unambiguous dark-matter case needs the curve
  **beyond $R \approx R_0$**, where the visible-mass contribution falls — a regime requiring a
  different method and not covered here.
- **First quadrant, $b=0$, inner Galaxy only.** One quadrant and the mid-plane spectrum; the
  tangent-point method works only for $R < R_0$.

## Bottom line

A clean, honest **validation**: a small, tested, offline-capable tool recovers the approximately flat
inner Milky Way rotation curve from public LAB HI data, with the absolute normalisation ~9% high (a
known terminal-velocity-estimator overestimate) and the non-Keplerian shape robust. A tooling +
validation contribution — and explicitly *not* a standalone dark-matter detection.

## Full referee round (2026-08-24): MAJOR REVISION, 13 findings

Evidence chain verified consistent (JSON -> macros -> prose; DOIs all clean; the make-figures
clobber is structurally closed). The problems are in what the numbers mean.

**MAJORs:** (1) The load-bearing citation is misattributed: McClure-Griffiths & Dickey 2016's
~7 km/s is the offset between the Clemens 1985 *CO* terminal curve and their fitted *HI* curve --
a cross-tracer, cross-survey comparison attributed only partly to the threshold method; they
never measure threshold-vs-fit on the same HI data. (2) Even at face value ~7 km/s explains only
about a third of the measured excess (11-28 km/s per longitude, mean ~20.6); the abstract's
single-cause "because" overstates, and the honest split in this findings file was dropped from
the paper; MG&D publish a uniform Q1 terminal curve the attribution could be tested against, and
their error-function edge fit could be run on the cached slices -- neither was done. (3) The
estimator is a global velocity maximum above 2 K with no contiguity requirement; the l=70/80
points (implied v_term 29.9/30.9 km/s vs flat-curve expectations ~14/~4) look like a
local-emission floor, at longitudes where V0 sin(l) is already 88% of the plotted value
(dropping both moves the mean only 256.6 -> 256.2 -- say so). (4) "Flat" is never quantified:
from the committed JSON, slope = 2.8 +/- 1.8 km/s/kpc and a Keplerian normalised at the
innermost point predicts 176 vs 263 measured at R=8 -- both one-liners that strengthen the
paper. (5) threshold_k=2.0 and flat_radius_min_kpc=4.0 are unswept and uncited (the l=30 point
sits 0.07 kpc above the cut; excluding it moves the mean to 258.6 -- no committed evidence
shows robustness).

**MINOR/NIT:** the +/-6 is np.std ddof=0 (ddof=1: 7; SEM: 2.9) and is scatter, not an
uncertainty, orthogonal to the correlated threshold systematic; the synthetic edge is 0.6 km/s
wide so threshold ~= inflection by construction (a realistic 5-10 km/s edge would measure the
bias in-house); the LSR-frame sentence glosses the standard-vs-measured solar motion mismatch;
N=6 never stated, no (l, R, V) table, and Results opens with a dangling "Table-style," artifact;
figure lacks error bars/Keplerian curve and draws the mean through excluded bar points;
arxiv-submission stale (pre style pass); "sub-megabyte survey" and "~30-line tool" both
overstate; \hiVflat is still un-namespaced mode-dependent (230 offline / 257 real), protected
today only by the source-marker rule.

**Status: fixes pending** (LAB slices cached; findings 2, 4, 5 fixable from the committed JSON
plus one edge-fit re-run).

**Status: RESOLVED (2026-08-25).** The misread citation is replaced by an in-house measurement
that strengthens the paper. An error-function edge estimator (`terminal_velocity_edge`) now runs
on the same spectra as the threshold estimator: the measured threshold-vs-edge offset is
**13.8 km/s** on the flat sample -- twice the misattributed ~7 -- and the edge-fit flat mean is
**242.9 +/- 8.0 km/s, within 2.9% of the Reid 2019 V0**, so the 9% excess is demonstrated to be
almost entirely the estimator, on this data, with MG&D 2016 now cited only for the edge-fit
method and the Clemens-CO provenance of their ~7 km/s stated. Flatness is quantified (slope
2.8 +/- 1.8 km/s/kpc; Keplerian at the outermost point 176 vs 263 measured). The synthetic
fixture's edge is widened to a realistic 5 km/s, so it now reproduces the wing overshoot
offline (13.6 km/s) instead of hiding it behind a 0.6 km/s edge. Sensitivities committed:
threshold sweep 1.5-5 K moves the mean 258 -> 252; drop-innermost 259; drop-outermost-two 256.
The +/-6 is relabelled as point-to-point scatter (ddof=1: 7 km/s) with the SEM (2.9) quoted and
its orthogonality to the shared estimator systematic stated. A per-longitude (l, R, V_thr,
V_edge) table is pipeline-written into the paper; both bar-excluded points are marked. Every
spectrum used is verified to be a single contiguous run above 2 K (the l=70/80 wing worry is
addressed by the edge fit, and the drop-outer-two variant is committed). Macros are namespaced
(\hiSyn*/\hiReal*), the figure draws both estimators + the Keplerian curve with bar points open,
the frame-mismatch sentence is honest, the size claims are scoped, the "Table-style," artifact
is gone, and the arXiv package is rebuilt clean.

## Dense resampling + VGPS cross-validation (2026-08-31)

The 8-longitude sample was the referee's finding 5 left half-fixed: the parameters were swept,
but the *sampling* was never questioned. Resampled to \hiNlong = 71 longitudes
($\ell = 10$--$80\arcdeg$, every $1\arcdeg$ -- LAB is all-sky and each slice is ~370 kB, so the
cost is one fetch loop) and compared against the tabulated VGPS terminal velocities of
McClure-Griffiths & Dickey 2016 (VizieR `J/ApJ/831/124/table1`, 748 rows over
$18.4 < \ell < 67.0\arcdeg$ at $0.065\arcdeg$). The comparison is made on `v_term` itself, so no
$R_0$/$V_0$ choice enters it.

**The headline result is a validation the paper did not previously have.** Over the 49
overlapping longitudes the edge-fit terminal velocities agree with VGPS to
**+1.04 +/- 0.15 km/s** (sd 1.05), and the slopes agree too: LAB edge
**+2.52 +/- 0.55** against MG&D's **+3.82 +/- 0.18** on their own 579 points beyond the same cut
(restricted to their longitude range and their sampling reduced to our 1-degree grid, the two
give +3.59 +/- 0.76 and +3.57 +/- 0.71 -- indistinguishable). A 0.5-degree-beam all-sky single-dish
survey reproduces an interferometric survey's terminal-velocity curve at the km/s level, in
normalisation and in shape. Robust to the matching window: +/-0.25/0.5/1.0 deg, mean or median,
all give +0.97 to +1.17; nearest-point +1.17.

**And it retracts the flatness claim.** The published slope was 2.8 +/- 1.8 km/s/kpc at N=6 --
1.6 sigma, reported as "consistent with flat". At N=51 it is **+1.79 +/- 0.57 (3.1 sigma)** on the
threshold estimator and **+2.52 +/- 0.55 (4.6 sigma)** on the edge fit. The curve rises. Nothing
about the data changed; the six-point sample simply could not resolve a rise this gentle, so the
old sentence was a statement about the sample's power dressed as a statement about the Galaxy.
MG&D's own curve rises at the same rate, so the *correct* reading is that the dense LAB curve
reproduces the reference in shape -- a better result than the one being retracted, arrived at by
giving the test enough points to fail. This is the mirror of the fashienv lesson: run the test even
when you expect it to confirm the claim. The title carried "Flat" and has been changed.

The Keplerian contrast (176 vs ~257 at the outermost point) is untouched, so the
dark-matter-relevant statement stands. The flat *level* also survives: 256.6 (scatter 5.3, SEM 0.7)
and 243.0 (scatter 5.5, SEM 0.8) at N=51, against 257 and 243 at N=6 -- the old point estimates
were right, only their error bars and the slope were underpowered.

**The threshold estimator is confirmed biased, and the bias is not a constant.** vs VGPS it sits
**+15.26 +/- 0.56 km/s** high with sd **3.92** (range +7.5 to +29.9), so no single offset would
correct it. The gap closes monotonically as the threshold climbs the profile edge --
1.5 K +17.2, 2 K +15.3, 3 K +13.2, 5 K +10.6, 10 K +7.5, **20 K +4.2**, 40 K +0.6 -- and 20 K is
the threshold MG&D seed their own fit from, so the two pipelines converge where they should. That
is a prediction the data could have refused and did not. The overshoot correlates with the fitted
edge width (r = 0.44, slope 1.02 over all 71; r = 0.74, slope 1.86 over the 49 overlap), against an
erfc-model prediction of ~1.5: the relation is real but loose, because the fitted width is noisy.

**Honest caveat, stated in the paper.** MG&D's estimator is the same family as ours -- a sum of two
error functions seeded from a 20 K threshold, on continuum-masked latitude-averaged spectra. So the
+1.04 bounds survey-to-survey and pipeline-to-pipeline differences (beam, sampling, absorption
treatment) and does *not* test the estimator family; both would move together under a different
definition of the terminal velocity. Claiming an independent-method confirmation here would be the
`dr20radio` vacuous-robustness-check error in a new costume.

**Two defects the dense sample exposed that the 8-point sample hid.** (1) 5 of the 71 spectra have
non-contiguous emission above 2 K -- the paper's flat assertion that "in every LAB spectrum used
here the emission above 2 K is a single contiguous velocity run" was true of the 8 and is not true
of the 71; it now reports the count. (2) The edge fit fails to converge on 1 sightline
($\ell = 12\arcdeg$), which is dropped and stated.

**Fixture change.** The offline fixture used a single 5 km/s edge width, so it could not exercise
the width-bias relation the paper now reports; it varies 3--8 km/s and reproduces the scaling
(slope 1.80, r = 0.96). The offline leg also now runs the same comparison machinery against the
injected curve in place of the VGPS one, recovering it to 0.09 km/s -- so the comparison code is
covered without a network, which is where the previous version had no coverage at all.

**Macro hygiene.** A non-finite metric used to reach the macro file as the literal string `nan`
(the 40 K sweep point has no crossing on a 30 K synthetic fixture); every numeric macro now goes
through a formatter that emits `--` instead, which the arXiv assembler already blocks on. The
`\hiSyn*` namespace was refilled by running the offline leg with `out=<tmpdir>` and calling
`_write_macros` on the real path, per CLAUDE.md; `preserve_live_macros` held every real value.
