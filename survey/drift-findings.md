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
