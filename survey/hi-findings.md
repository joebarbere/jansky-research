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
