# Findings — ultra-steep-spectrum source hunt (TGSS × NVSS)

Run of `jansky_research.spectra.run` over a high-latitude extragalactic field near the North
Galactic Pole (centre RA $180.0°$, Dec $+30.0°$; $3°$ cone), cross-matching **TGSS ADR1**
(147.5 MHz; Intema et al. 2017) with **NVSS** (1.4 GHz; Condon et al. 1998) on VizieR, 15" match
radius. This is the honest assessment, revised after the GATE-2 science review. The full matched
catalogue is committed at [`survey/uss_candidates.csv`](uss_candidates.csv) (456 sources).

## Sanity checks (the method works)

- **456 matched sources**; **median $\alpha_{150}^{1400} = -0.73$** — consistent with the typical
  $-0.7$ to $-0.8$ range for flux-limited extragalactic radio samples. The closest comparator,
  de Gasperin, Intema & Frail (2018; the 1.4M-source TGSS×NVSS index catalogue), finds a weighted
  mean $-0.79$; our median sitting slightly flatter is consistent with TGSS flux-scale inflation
  making some steep sources look flatter. The distribution peaks near $-0.8$ as expected.

## The finding: 6 ultra-steep-spectrum candidates, none a known high-$z$ radio galaxy

Six sources have $\alpha < -1.3$ (the classic HzRG selection; De Breuck et al. 2000). **None has a
measured redshift or a high-$z$ radio-galaxy classification in either NED or SIMBAD** — all six are
HzRG *candidates*.

| # | RA | Dec | $\alpha$ | $S_{150}$ (mJy) | match sep | NED | SIMBAD | note |
|---|------|------|------|------|------|------|------|------|
| 1 | 181.0497 | +29.3853 | $-1.65$ | 126 | 5.8" | — | — | no counterpart in either DB |
| 2 | 179.5148 | +31.6903 | $-1.37$ | 100 | **9.4"** | — | NVSS radio src | loose match → lower confidence |
| 3 | 181.0344 | +30.2634 | $-1.34$ | 131 | 5.4" | WISE IR (IrS) | — | faint IR counterpart |
| 4 | 177.8709 | +28.8189 | $-1.32$ | 191 | 3.5" | WISE IR (`*`) | NVSS radio src | NED's nearest is a star (likely chance IR align) |
| 5 | 180.2691 | +29.1648 | $-1.31$ | 75 | 1.0" | — | — | tight match, no counterpart |
| 6 | 177.0305 | +29.3276 | $-1.31$ | 135 | 0.3" | — | NVSS radio src | tight match |

Where SIMBAD "matches", it knows the source *only* as a bare `NVSS J…` radio entry — no
classification, no redshift. So none of the six is catalogued as a high-$z$ radio galaxy.

## Manual cross-check against de Gasperin (2018) — the candidates mostly do NOT survive

The recommended follow-up (GATE-2 review, S-3) was a manual cross-check against published TGSS×NVSS
catalogues. The authoritative one is **de Gasperin, Intema & Frail (2018)**, the 1.4M-source
flux-scale-**corrected** TGSS×NVSS spectral-index catalogue (`J/MNRAS/474/5008/spidxcat`). All six
candidates are already in it; comparing their catalogued index to ours:

| my $\alpha$ | de Gasperin SpIndex | sep | USS in de Gasperin? |
|------|------|------|------|
| $-1.65$ | $-1.52\pm0.12$ | 1.9" | **yes** |
| $-1.37$ | $-1.21\pm0.10$ | 4.6" | no |
| $-1.34$ | $-1.16\pm0.07$ | 2.5" | no |
| $-1.32$ | $-1.27\pm0.05$ | 0.7" | no |
| $-1.31$ | $-0.92\pm0.00$ | 3.2" | no |
| $-1.31$ | $-1.30\pm0.07$ | 0.9" | borderline |

**Two conclusions, both deflating:**

1. **Our indices are systematically steeper** than the corrected catalogue by a mean
   $\Delta\alpha = -0.15$ (range $-0.01$ to $-0.39$). This is exactly the TGSS ADR1 flux-scale
   inflation the review warned about (B-1), magnified by this field straddling the Dec $+30°$ edge
   of the rescaled coverage. Using the corrected fluxes, **only 1 of the 6 is genuinely USS**
   ($\alpha<-1.3$); a 2nd is borderline ($-1.30$) and the other 4 are ordinary steep sources.
2. **All six are already in a published 1.4M-source catalogue**, so **none is a new source.**

The automated Saxena (2018) check could not be validated (its VizieR table did not resolve), but it
is now moot: de Gasperin already settles the question.

## Revised conclusion — no novel discovery

The honest outcome of the manual cross-check is that **there is no novel ultra-steep-spectrum
discovery here.** The one source de Gasperin also calls USS ($\alpha=-1.52$) is a *known* entry in
their catalogue with no measured redshift — at best a long-standing HzRG candidate, not something
this search newly revealed. The other five were USS only in our *uncorrected* TGSS measurement and
are not USS once the flux scale is handled properly.

This is a genuinely useful (if humbling) result: it is a **cautionary, fully reproducible
demonstration that raw TGSS ADR1 × NVSS USS selection is dominated by the TGSS flux-scale
systematic**, and it validates directly against the authoritative catalogue. The deliverable is the
open, tested, CPU-only tool **plus this honest cross-check** — emphatically not a discovery.

## Honest limitations (this is a candidate list, not a discovery)

- **The candidates do not survive the manual cross-check** (see the de Gasperin section above):
  only 1 of 6 is genuinely USS once the TGSS flux scale is corrected, and all 6 are already in a
  published catalogue. TGSS × NVSS USS selection is an established technique (de Gasperin et al.
  2018; Saxena et al. 2018), and this run did not improve on it.
- **Chance-coincidence rate:** at 15" with the NVSS source density in this field, $\sim$1–2 false
  matches are expected across the full 456-source sample. The tight USS matches ($\le 5.8$") are
  therefore robust; candidate 2's lone 9.4" association is the most likely to be spurious.
- **TGSS flux-scale systematic** (Hurley-Walker 2017, arXiv:1703.06635): $\sim$15% typical, up to
  $\sim$40–50% in places, biasing $\alpha$ by $\sim$0.1–0.2. The rescaled TGSS-RSADR1 covers only
  Dec $\le +30°$, so the **upper half of this cone uses the uncorrected ADR1 scale** and may carry a
  larger, unquantified bias — relevant to candidate 2 (Dec $+31.7°$) and to 5, 6 sitting at the
  $-1.31$ threshold.
- **Two-point index only.** The 4 "inverted" ($\alpha>0$) sources are *not* peaked-spectrum/GPS
  claims — at 150 MHz vs 1.4 GHz a positive index is more likely a calibration/variability/
  resolution artifact; confirming a turnover needs a third frequency (e.g. VLASS 3 GHz).
- NED's nearest object to candidate 4 is classified as a star — most likely a chance infrared
  alignment, not the radio source's host.

## Bottom line

- **No novel discovery.** The manual cross-check against the authoritative de Gasperin (2018)
  catalogue reduced the six raw USS candidates to a single confirmed-USS source — itself already
  catalogued and without a redshift. The other five were flux-scale artifacts.
- **What is real and worth keeping:** an open, tested, CPU-only TGSS×NVSS USS-selection tool; the
  textbook population median ($-0.73$) as a validation; and a clean, reproducible demonstration that
  *uncorrected* TGSS ADR1 fluxes inflate $\alpha$ by $\sim0.15$ and manufacture spurious USS
  candidates — which is exactly why a cross-check against the corrected de Gasperin catalogue is
  mandatory. Reported as a negative/cautionary result, honestly, rather than dressed up as a find.

## Full referee round (2026-08-25): MAJOR REVISION, 16 findings, four BLOCKERs

The tool is sound, the committed evidence is internally consistent, and the six-source manual
lookup reproduces to the digit on a fresh VizieR query. But the paper's causal claim is
refuted by its own field, and the honest paper hiding underneath is stronger.

**BLOCKER 1: the reference catalogue is not flux-scale corrected, and its authors say so.**
main.tex/rnaas.tex call de Gasperin (2018) "the authoritative, flux-scale-CORRECTED TGSS×NVSS
catalogue"; §4 of their own paper says the Hurley-Walker 15% TGSS systematic "can
systematically bias our spectral index values by ~0.06" — the catalogue is built from the same
uncorrected ADR1, so a raw-minus-de-Gasperin difference is structurally incapable of measuring
the TGSS flux scale (which is common to both), and the quoted 0.15 would require a 41% flux
inflation nothing claims.

**BLOCKER 2: the claimed systematic offset is absent from the population.** Referee-side
cross-match of all 456 committed positions against the non-limit de Gasperin entries (387
pairs): mean Δα = +0.004 ± 0.005 — no offset; per-source scatter 0.093 = 1.2× the combined
formal errors. The −0.15 exists only in the six sources selected on the noisy variable (a
no-free-parameter selection model predicts −0.05 ± 0.035 of the observed −0.11) — the
rmstructure lesson in a new costume: the number measures selection on the extreme.

**BLOCKER 3: the declination-edge mechanism is refuted by the same 387 pairs** (Dec ≤ +30:
+0.006 ± 0.006; Dec > +30: +0.002 ± 0.007 — indistinguishable; the largest of the six offsets
is INSIDE the rescaled coverage; and there is no mechanism — only the correction's
availability changes at +30°, and neither side applies it).

**BLOCKER 4: the two headline numbers (Δα ≈ −0.15; "only 1 of 6 survives") are hand-typed**
under "no number is typed by hand": no macro, no JSON key, no query log; reference_spindex()
is pragma-no-cover and called from NOWHERE — the cross-check was a one-time manual lookup.
(To the author's credit: the referee re-ran it and all six values are correct. The provenance
is what is broken.)

**MAJORs:** one of the six "reference" indices is a NON-DETECTION LIMIT (Scode=L, e=0.000 —
"±0.00 should have stopped the analysis"), supplying the largest term (−0.392; drop it: −0.11)
and breaking "all 6 already appear in it" (that entry is a TGSS non-detection DISAGREEING with
our 9.2σ ADR1 detection); the selection function is reported at half strength — the reference
catalogue holds 11 USS sources in the same cone and the naive cut recovers 1, so the method is
17% PURE and 9% COMPLETE, and "manufactures spurious candidates" is wrong at the population
level (6 found where the reference has 11); the citable mechanism is in de Gasperin's own
Fig. 3: NVSS and TGSS are equally sensitive at exactly α = −1.3, so the joint sample is
truncated for everything steeper (to be flagged at −1.65 needs S₁₅₀ ≥ 102 mJy — 51% of the
sample); the published reproduction command (`python -m jansky_research.spectra`, no
coordinates) runs the SYNTHETIC leg into the repo root and clobbers the unguarded
uss_candidates.csv + both paper figures; the field centre/radius appear in neither manuscript
(\usSource unused and would mis-render); "validates the method" via the median is a test that
cannot fail (dG cone median −0.749 vs ours −0.73 — real but powerless; the same run is 9%
complete); the offline fixture injects USS 3σ from the cut.

**MINOR/NIT:** "all six already catalogued so none is new" is a tautology (the real check —
SIMBAD/NED, no redshifts — holds on re-query but rests on another uncommitted manual lookup,
annotate_known also having zero callers); limitation (i) mis-states Hurley-Walker (a zero-mean
15–17% SCATTER → 0.06–0.07 in α, not an inflation of 0.1–0.2; the docstring's "40–50%" not
found in the source); the RNAAS's only figure does not display the note's own conclusion (the
corrected values are not plotted); the ~1–2 chance matches are right but computed nowhere
(1.48); survey/ and results/ carry duplicate candidate CSVs, only one guarded; the printed
equation says 150 MHz where the code uses 147.5; this file's own "flatter" sign error
(inflating S₁₅₀ makes α STEEPER); stale arXiv package.

**Status: fixes pending.** The single change: run reference_spindex() over all 456 sources,
commit the confusion matrix, and let the paper say what it is entitled to: a TGSS×NVSS
catalogue USS cut is ~17% pure and ~9% complete because the two surveys are equally sensitive
at exactly the threshold — a stronger cautionary result than the one currently claimed, one
function call away.
