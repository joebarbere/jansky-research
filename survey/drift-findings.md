# Findings — CPU-only SETI drift-search injection-recovery benchmark (+ honest Voyager check)

`jansky_research.driftsearch` builds a reproducible, offline injection-recovery benchmark on the
pure-NumPy Doppler-drift search in `jansky.seti`, and includes an honest real-data check against the
Voyager-1 file. Run with `python -m jansky_research.driftsearch` (100 trials/cell by default).

## 1. Injection-recovery benchmark (the deliverable)

Inject synthetic drifting tones over a grid of S/N × drift rate, run the brute-force de-drift
search, and measure the recovered fraction. With 100 trials/cell and a detection threshold of 10
(noise-only best S/N peaks near ~5):

- **50% completeness at injected S/N $\approx 1.3$; 90% at $\approx 1.5$** (drift-averaged; the
  binomial scatter at the crossing is $\sim\pm0.1$ at 100 trials/cell). Recovery rises from 0 below
  S/N $\sim$0.75 and saturates at 100% by S/N 2 (heatmap: `paper/figures/drift_recovery.pdf`).
- **False-positive rate $< 0.9\%$** (0/400 noise-only realisations; 95% Clopper-Pearson upper
  limit) at this threshold.
- Recovery is flat across the tested drift rates (0–0.6 chan/sample), as expected for a brute-force
  search whose grid includes the true drift.

**What's genuinely new.** `setigen` (Brzycki et al. 2022) already provides synthetic signal
injection compatible with BL data, and turboSETI efficiencies have been characterised with it. The
contribution here is narrower but real: a **self-contained, hardware-independent benchmark cell that
fixes both the injection model and the detector** in one reproducible unit — a portable reference,
not a new injector. `results/drift_metrics.json` holds the full curve.

## 2. Voyager-1 real-data check — an honest NEGATIVE result

Pointing the same detector at the Breakthrough Listen GBT open-data file of Voyager 1
(`validate_voyager`) is a cautionary result, not a success:

- The brightest channel is **band-centre (channel N/2)** with a value $\sim10^{3}$× any real tone —
  the spectrometer **DC-spike artifact**, not the spacecraft. A naive brightest-channel "detection"
  reports this and is wrong.
- At the **documented Voyager-1 carrier (8420.216 MHz; Estévez 2021)** — searched with a wide drift
  grid, since the carrier drifts $\sim-0.69$ Hz/s (several channels/sample here) — the detector
  returns only **S/N $\approx 4.85$**, indistinguishable from a blank window (**4.59**). It does
  **not** recover the carrier (`recovered = False`).
- So the `jansky.seti` detector, validated on injected tones in clean synthetic noise, **fails on the
  real Voyager-1 data**: it is fooled by the DC spike, and the real drifting carrier amid the BL
  data's structure is beyond this teaching-grade tool. Recovering Voyager needs proper SETI tooling
  (blimpy/turboSETI). This honestly **bounds where the tool works.**

(Earlier in development this check appeared to "detect Voyager at S/N $2\times10^{5}$" — that was the
DC spike. Catching and reporting it is the point.)

## Honest limitations

- **Relative, not absolute.** The benchmark's "S/N" and "drift rate (chan/sample)" are the internal
  units of the `jansky.seti` synthetic model, not calibrated janskys or Hz s$^{-1}$. It compares
  detectors on a common reference; it is not an on-sky sensitivity.
- **No cadence / RFI realism.** A single tone in Gaussian noise — no ON/OFF cadence, no realistic
  RFI; it measures raw recovery, not RFI rejection.
- **Real-data gap (Section 2).** The detector does not generalise to real BL data as-is.
- **No discovery, and none expected** (the only real "signal" near the data is Voyager, which this
  tool does not even recover).

## Bottom line

A reproducible CPU-only injection-recovery **benchmark** for the `jansky.seti` drift detector (50%
completeness at S/N $\approx 1.3$, false-positive rate $<0.9\%$), plus an **honest negative
real-data check**: the same detector is fooled by the band-centre DC spike and does not recover the
true Voyager-1 carrier. A tooling/benchmark contribution with its limits stated plainly — no
overclaimed validation.

## Referee round on the style conversion (2026-08-23)

Verdict *minor revision*. All three limitation paragraphs survived the dissolution of
`\section{Limitations}` intact (none carried a magnitude). Two MAJORs, both fixed:

1. **A one-detector, one-file result had been restated as a class-level claim.** Dropping the
   inferential "So" left the indefinite subject `A detector validated on injected tones in clean
   synthetic noise fails on the real Voyager-1 data` --- i.e. *any* such detector. Made definite
   ("This detector, validated on ..., fails ...").
2. **A limitation softened in the direction the evidence does not support.** `does not generalize
   ... as-is` had become `as configured`, implying a settings choice might succeed --- which the
   same paragraph rules out (the drift grid *was* widened for the -0.69 Hz/s carrier, and the
   conclusion is that purpose-built SETI software is required). Restored, and the clause order
   restored so the failure leads rather than trailing a scope note.

Also fixed: the abstract's Voyager test had become additive ("We also point ...") rather than a
test of whether the benchmark transfers, with the non-detection demoted to a participle; the
"relative, not absolute" caveat no longer sat next to any number, so the Discussion's headline S/N
now says "in the model's internal units"; the section stated no reason a null is worth reporting.

**Title:** "Honest" dropped. The conversion removed all six body instances, leaving it an orphan
the paper no longer sustained; in journal register it is a claim about the authors, not the result,
and `\shorttitle` already omitted it. `arxiv-submission/metadata.yaml` is gitignored and
regenerates from the .tex, so nothing needed syncing.

**Pre-existing, still open:** the Voyager numbers (4.85, 4.59, 8420.216 MHz, -0.69 Hz/s) and the
95% Clopper-Pearson limit ~0.9% are typed into `main.tex` despite the file's own comment that no
number is hand-typed. The 0.9% needs N = 400, which lives only in the `fpr_trials=400` default in
`src/jansky_research/driftsearch.py`, not in `results/drift_metrics.json`. Rule-of-three on 0/400
gives 0.75%, so 0.9% is conservative rather than wrong.

## Full referee round (2026-08-26): MAJOR REVISION, 13 findings, two BLOCKERs

The synthetic benchmark leg is sound and byte-identically reproducible (the referee re-ran it:
committed JSON and figure reproduce exactly at seed 0 / 30 trials; all four DOIs Crossref-clean;
macros all match). The Voyager leg — half the title — is not a null.

**BLOCKER 1: the "Voyager-1 null" is a targeting error; the detector RECOVERS the carrier.**
The hard-coded VOYAGER_CARRIER_MHZ = 8420.216 maps to channel 419016 of the cached file — blank
sky (peak MAD-z 2.76). The actual Voyager signal is at channel 747929 = **8419.29703 MHz**,
unambiguous: MAD-z 205, telemetry-subcarrier sideband doublets at ±22.50 kHz, and a peak-channel
walk of 2.444 channels/sample = **−0.3741 Hz/s** (not the paper's hand-typed −0.69). Running the
module's own search there: **S/N = 997.6 at best_drift 2.45** (one grid step from the measured
slope) vs blank 4.59. Every null sentence — abstract, Discussion, title, two README rows — is
wrong, not overstated. The referee reproduced the paper's published 4.85/4.59 first (the
committed narrative is faithful to the code as written; the code searched the wrong place).
Offered hypothesis, not verified: the 0.9197 MHz offset ≈ 32.7 km/s at X-band — a
barycentric/topocentric frame mismatch. Do not construct a replacement constant; locate the
carrier in the data.

**BLOCKER 2: two hand-typed physical constants are wrong for this file** (8420.216 MHz;
−0.69 Hz/s), and the refs.bib estevez2021 note propagates the frequency into a citation
annotation that is untrue of this file.

**MAJORs:** the Voyager leg has NO committed evidence, no test, no path in make reproduce (the
"REAL public data" target runs the synthetic leg only) — under a "no number is typed by hand"
header carrying eight hand-typed numbers; the benchmark's configuration (64×512 waterfall,
41-point drift grid, 1.5-channel line width, 400 FPR trials, seed) is the deliverable and
appears nowhere in paper or JSON — and the FPR is entirely a function of it; the per-drift
completeness matrix is discarded (only the drift-average is committed) so "flat across drift" is
not auditable — referee recomputed it and the claim HOLDS (1.250/1.297/1.297); the recovery
figure's imshow y-axis is quantitatively wrong (non-uniform S/N rows on a linear extent: the
50% crossing reads ≈1.85 off the figure vs the caption's 1.3); the FPR test could not have
failed (noise-only best-S/N distribution: mean 4.12, σ 0.33, max 5.28 in 300 draws vs threshold
10 ≈ 18σ — "cluster near S/N~5" overstates; 0/400 distinguishes nothing).

**MINOR/NIT:** the DC-spike bug this paper exists to report has NO regression test (the failure
mode is live: the module's own _snr on the DC window returns 1.985e5 — the referee reproduced
the original spurious detection); trials mismatch (JSON 30, run() default 100, findings doc
"100" — re-run at both: crossings move ≤0.006, published rounding unaffected; findings-doc
scatter estimate ±0.1 overstates the true ≈0.03 SE); no source/is_real marker on JSON or macros
(a bare run(".") silently rewrites the evidence with different-depth numbers, no guard fires);
the "~0.9%" is the two-sided CP endpoint labelled one-sided (one-sided = 0.746%) and N=400
exists only as a code default; on-grid injection checked and NOT inflating completeness
(off-grid at half-step: C50 1.279–1.292 vs 1.250 — within scatter; worth one stated row);
brzycki2022 title truncated.

**Verdict: MAJOR REVISION.** The single change: point the Voyager search at the Voyager signal
(locate the carrier in the file — 8419.29703 MHz, ±22.5 kHz sidebands, −0.374 Hz/s), commit the
resulting dict, and rewrite Section 4 and the title around the recovery at S/N ≈ 10³ — keeping
the genuinely valuable caution that the band-centre DC spike is 2400× the carrier and a
brightest-channel search reports the artifact. The corrected paper is stronger than the
submitted one.
