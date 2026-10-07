# Plan 98 — `hiblend` v2: a calibration that must pass a held-out test before C2 and β are rerun

Status: **plan, controls frozen 2026-10-06 before any v2 code ran on real data.** Follows plan 97,
whose frozen outcome was **ambiguous** (C2 failed) and stays so whatever this plan finds. Changes
after the first v2 run go in `survey/hiblend-findings.md` with the reason.

## Why, and what is already seen

Plan 97's post-hoc diagnostics (`survey/hiblend-findings.md` step 4,
`results/hiblend_diagnostics.json`) found that the C3 calibration does not fit isolated targets
it was not trained on. The held-out median residual by S/N quintile reached +0.039 dex (8.9σ) in
the top quintile and −0.025 dex in the second, about the size of the effect the test looks for.
The null slope also fell from 0.256 ± 0.087 to 0.117 ± 0.085 when the null sample was reweighted
to the primary sample's S/N.

**This plan is not blind.** Before any of it was written, the following had been seen:
- the plan-97 calibration residuals of every isolated target, binned by S/N (D0);
- β_null, overall and by S/N tercile (step 3, D1);
- the primary β = 0.385 ± 0.136 (step 3).

A pass here is therefore weaker evidence than a pass on first contact would have been. Two
choices answer that:
1. The new calibration is certified on held-out data, under a rule fixed here (C3′).
2. C2 is made **stricter** than in plan 97, so that a calibration change that merely moves the
   null slope cannot pass.

Provenance (CLAUDE.md, *an oracle that fed a fit may never certify the result*):
- **Production side.** The calibration's functional form was chosen after seeing D0 on all
  isolated targets, and the density term after seeing D3.
- **Certification side.** C3′ uses an RA-strip split that the form was not tuned on. The form
  itself was still motivated by those targets, so C3′ certifies the fit, not the choice of form.

## What changes (frozen)

1. **Calibration form.**
   - Log S/N enters as a **linear spline**, with knots at the 10th, 20th, …, 90th percentiles of
     log S/N in the calibration sample. This replaces plan 97's quadratic, which cannot follow
     the rise in the top quintile.
   - The other terms stay as in plan 97: linear in log W50, log S_g and dec/90.
   - **New term:** linear in log(1 + N₁₅). D3 measured −0.013 ± 0.004 for it in isolated
     targets.
   - The covariates stay as plan 97 built them: S_g not blend-subtracted, S/N from the flux
     error.
2. **C3′, held-out calibration (new gate; it runs first).**
   - Split the isolated sample by RA strip, with even and odd strips defined by `floor(ra / 2°)`.
     Fit on even strips and evaluate on odd, then the reverse.
   - In each direction, take the held-out median residual in each of 10 S/N deciles, with its
     bootstrap SE.
   - **Pass** requires both directions to have every decile median |Δ| < 0.010 dex (under a
     quarter of the ≈0.045 dex equal-flux signal), and χ² of the 10 medians against 0 with
     p > 0.01.
   - **Fail** stops the run, with the outcome *"the survey-to-survey flux scale cannot be
     calibrated to the precision this test needs"*. That is a reportable negative about the
     method.
3. **C0 and C1 rerun under the new calibration**, with the plan-97 criteria unchanged. C0 needs
   ≥ 80% detection of β = 1 at ≥ 3σ. C1 needs β = 1.0 ± 0.3 for a planted signal, and β
   consistent with 0 when nothing is planted. A more flexible calibration can absorb signal, and
   C1 is what measures that.
4. **C4 is unchanged** (same match), and is carried over from plan 97.
5. **C2, strengthened.** β_null must satisfy |β_null| < 2σ:
   - overall, **and**
   - **separately in each of the three null-sample S/N terciles** (the D1 split).

   Plan 97's failure was concentrated in the lowest tercile (0.38 ± 0.12). A calibration that
   only diluted it overall would fail here.
6. **Seed 98**, set before running. The bootstrap counts are unchanged.

## Outcomes, stated in advance

- **C3′ fails:** stop. Report that calibration is the limiting systematic; no β is quoted.
- **C0 or C1 fails:** ambiguous. Report which control failed; no claim.
- **C2 fails (overall or any tercile):** ambiguous. The environment-linked systematic is not
  calibration alone; no claim.
- **All pass, and β > 0 at ≥ 3σ:** *blending supported, after recalibration.* The write-up must
  say that the analysis is a second attempt designed after plan 97's samples were seen.
- **All pass, and β consistent with 0:** *blending not supported.* Report the 95% upper limit on
  β.

**Secondary (reported, not part of any claim):** β of the primary sample by S/N tercile.

## Deliverables

- `hiblend.calibration_design_v2`, `heldout_calibration_check` (C3′), and `run_gated_v2`, with
  tests on the synthetic sky. One test plants an S/N nonlinearity of the D0 shape and requires
  plan 97's form to fail C3′ and v2's form to pass.
- `scripts/hiblend_real.py --v2`, writing `results/hiblend_v2_metrics.json` with `source` set.
  Plan 97's results file is not touched.
- `survey/hiblend-findings.md` step 5.
