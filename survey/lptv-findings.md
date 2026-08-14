# Findings — LPT catalogue v3 + Stokes-V forced photometry (plan 44)

`jansky_research.lptv` (+ 3 new rows in `data/lpt_sample.csv`, + `scripts/lptv_real.py`): extends
the merged `lpt` catalogue to v3 and runs the first systematic multi-epoch forced Stokes-V
photometry at all LPT positions.

## GATE 0 (2026-07-08)

- **3 new rows verified and transcribed** from the discovery papers, coordinates decoded from each
  source name (all agree with the CSV to <1″ — the provenance discipline that caught the Rea+2026
  review's 2225-vs-3225 s typo):
  - ASKAP J142431.2−612611 (arXiv:2603.07857): 216.130, −61.4364, P=2147.27 s (36 min), possible
    WD binary unconfirmed, XMM non-detection; already has per-source pol (circular ~8%, RM −222).
  - ASKAP J165130.3−450520 (VASTER, arXiv:2606.20067): 252.8761, −45.0888, P=23317.9 s (6.48 hr).
  - ASKAP J170036.6−445758 (VASTER, arXiv:2606.20067): 255.1525, −44.9661, P=16895.9 s (4.69 hr).
  - The GATE-0 agent flagged ASKAP J1745-5051 as a possible 4th addition — but it is ALREADY in
    the v2 CSV (row 13). So v3 = 13 + 3 = **16**.
- **Stokes-V novelty PASS**: no targeted multi-epoch LPT V survey exists. RACS-low2 Paper VIII
  (arXiv:2606.16182) is a BLIND V catalogue (did not target LPTs) — a fence. Per-source pol
  published for a handful (GLEAM-X J1627 ~90% linear; ASKAP J1935 >70% circular weak state; J1424
  circular ~8%; CHIME J1634 ~100% circular) — fences to cite, not the systematic survey.
- Rea+2026 review lists 12 confirmed + 2 WD pulsars, no population synthesis, no systematic V.
- CASDA RACS-low1/low2/mid V products cover the (southern, Dec −45 to −62) LPT positions.

## Catalogue leg (DONE)

Population stats at **N=16** (all regenerate from the CSV):

| quantity | value |
|---|---|
| N (confirmed, v3) | 16 |
| WD binaries / candidates | 7 |
| X-ray detected | 3 |
| Ṗ measurements | 2 |
| period range | 7.0 min – 6.48 hr |
| median period | 73.4 min |
| below death line / Ṗ-constrained | 9 / 9 |

**Period-split test (WD-binary vs rest): Δlog median 0.176, permutation p = 0.52.** The hinted
~78-min binary-boundary is **still not significant at N=16** — the plan's headline question,
answered honestly. The two long-period VASTER additions (binary status not reported) and the
36-min J1424 do not sharpen the split; reported, not spun.

## Recover-a-known (offline, in CI)

- Injected |V|/I=0.6 circularly-polarized point source into a synthetic RACS cutout →
  `measure_circular_pol` recovers V within noise, classifies `highly_circular`, correct handedness
  (LCP for V<0); blank field stays below 5σ. `summarize_v_sweep` unit-tested: detections, deepest
  limits, leakage vetting (a formally-significant V below 0.6% of I rejected), inter-epoch
  handedness flips.

## V leg (DONE — CASDA sweep, 229 min, 191 rows, 154 good, 1 failed)

`scripts/lptv_real.py` reused the plan-41 `wdpulsar_real.py` CASDA machinery (obscore query,
complete I+V grouping, SODA cutout, retry-with-relogin, resume-by-CSV) at all 16 LPT positions.
Forced I+V per (LPT, obs_id) across RACS low+mid; 15/16 LPTs covered (ILT J1101+5521 at Dec +55.5
is outside the RACS southern footprint → uncovered). Each LPT had 5–22 epochs.

**Not all-limits — 1 secure + 1 candidate single-epoch circular detection + 1 confusion vetoed:**

| target | epoch | I (mJy) | V (mJy) | \|V\|/I | V sig | offset | verdict |
|---|---|---|---|---|---|---|---|
| ASKAP J174508.9-505149 (accreting CV) | mid | 21.6 | −3.25 | 15% | 21.6σ | 0.73″ | **secure detection** (on-source, 25× above leakage floor, known bright CV) |
| ASKAP J165130.3-450520 (2026 VASTER) | mid | 4.37 | +2.56 | 59% | 12.5σ | 3.2″ | **candidate** — real 12.5σ V signal but 3.2″ off (3–6× the astrometric budget); association not certain |
| ASKAP J183950.5-075635 (longest P) | low | 240.8 | +96.1 | 40% | 294σ | 5.3″ | **confused** — 260× the source's median I, off-centre → nearby source, vetoed |

- Both signals are **single-epoch** — each source is a limit in all its OTHER epochs. So these are
  burst states caught in snapshots, NOT persistent circularly-polarized counterparts. Framed as
  such: persistent circular pol is not a class property, but RACS occasionally catches a burst.
- **Confusion veto** (offset > 4″ AND detection-epoch I > 10× the source's median I): flags
  J183950's 240 mJy / 5.3″ peak as a confusing source, not the LPT. The leakage veto alone
  (|V| > 0.6% I) does NOT catch confusion. **Secure/candidate split** (offset < 2″ = secure): the
  on-centre CV is secure; the 3.2″ VASTER peak is a candidate. Median 3σ V limit **0.474 mJy**
  across the 15 covered LPTs (0.36–0.77); 0 handedness flips.

## GATE-2 (PASS with required fixes, all applied)

Reviewer independently re-derived every detection from the raw CSV. Verdict: J1745 a clean secure
detection; J1651 a real 12.5σ signal but overclaimed. Fixes applied:
- **R1 — J1651 downgraded to candidate** (was "more strikingly" in the abstract — backwards; it's
  the *weaker* claim). The 3.2″ offset (~3–6× the astrometric budget for an 18σ source) is now
  disclosed; a new `secure`/`candidate` split (offset < 2″) encodes it, with macros
  lvRealNVSecure=1, lvRealNVCandidate=1.
- **R2 — confusion-veto thresholds disclosed** as heuristics set with the sweep in view (4″/10×
  confusion, 2″ secure); the paper notes the confused peak is far from the boundary in both
  dimensions while the candidate sits just inside — which is *why* it's a candidate.
- Suggested, applied: off-axis-leakage note (J1745 clears even a pessimistic several-% off-axis
  leakage; leakage wouldn't switch on for one epoch); a near-boundary confusion unit test (3.2″/
  8.6× → candidate); the period-split p-direction (0.27→0.52, *away* from significance); per-source
  pol citations incl. CHIME J1634 ~100% circular (arXiv:2507.05139).

## Referee round (2026-08-13) — peak search + a polarization claim that is a flux limit

**BLOCKER 1 confirmed: 152/154 epochs have I > 0**, median offset 10.15″ (more edge-weighted
than uniform-in-disc) — the same peak-search defect as wdpulsar, through the same helper,
under a title claiming "Forced Photometry". The forced re-run (search_arcsec=0.0) is queued
behind wdpulsar's on the same CASDA session.

**BLOCKER 2: "persistent circular polarization is not a class property" is a statement about
polarization drawn from sources undetected in Stokes I in 149/154 epochs.** The implied
fractional-polarization limits are 49–114% — the census does not exclude 100% circular
polarization for most targets. Sharpest case: CHIME J1634+44, published at ~100% circular
polarization (the class's natural recover-a-known), appears as "<0.41 mJy" — a flux
non-detection reported without noting the recovery failure. The claim must become "no LPT
shows persistent circularly polarized emission above ~0.5 mJy in RACS snapshots".

Also for the rewrite: duty cycle appears nowhere (the four longest-period sources have
27–43% phase coverage — their non-detections constrain nothing); there are **no RACS-low1
data** in the census (the abstract cites mcconnell2020 for data that never enter; 50 of 81
low-band epochs are RACS-low3, unnamed and uncited); the headline 0.474 mJy is the median
*best-epoch* limit — the per-epoch median is 0.726 (90th pct 2.09), and for an intermittent
source the per-epoch number is the relevant one; the confusion veto's denominator is a
noise-max floor and its offset clause passes 92% of all rows; "an exact permutation test" is
a 20,000-draw Monte Carlo whose quoted p = 0.5219 differs from the true exact p = 0.5258
(11,440 partitions, both computed); and the N=13 → 0.27 prose number is unbacked (it
reproduces: exact p = 0.2692). "First systematic multi-epoch" needs narrowing to "first
targeted, multi-epoch, uniform V compilation" against vaster2026 (itself a multi-epoch ASKAP
variability pipeline that discovered one of these sources) and the published racslow2 blind
V catalogue covering 15/16 positions.

## Forced re-measurement complete (2026-08-14) — and the vetoed source is a pulse

The re-sweep (191 epochs, `search_arcsec=0.0`) inverted the paper's most interesting call.
Under forced photometry at the catalogue positions (median offset 0.85″; blank-sky epochs now
48% positive in I, as they must be):

- **ASKAP J1745−5051** confirms: 15% circular, 21.6σ, at-position.
- **ASKAP J165130.3−450520** *strengthens*: 67% circular at 11.9σ, at-position — the 3.2″
  positional ambiguity that made it a "candidate" is gone. Now a detection.
- **ASKAP J183950.5−075635** — vetoed by the old analysis as "a nearby confusing source" —
  is a **single-night pulse recovery**: I = 164 mJy, V = +65 mJy (40% circular, 199σ) at
  MJD 60335, with all nine other epochs consistent with zero including the following night.
  The two adjudication queries the referee prescribed settle it: the survey's own blind
  leakage-corrected V catalogue has **no source within 60″** (different epoch — the source
  was off), and the I catalogue has **no source brighter than 4 mJy within 2 arcmin** — there
  is no neighbour to confuse with. The 5.3″ "offset" that triggered the veto was a fifth of
  the RACS-low beam. This is the source doing what its published pulses do.

The paper now reports **three** single-epoch burst detections and rescopes its negative
honestly: the non-detections are limits on circularly polarized *flux* (best-epoch median
0.47 mJy, per-epoch median 0.73), not on fractional polarization — implied |V|/I limits are
49–114%, and CHIME J1634+44 (published ~100% circular) is among the flux non-recoveries.
Duty cycle is now arithmetic in the text (~15-min integrations; the four longest-period
sources have 27–43% phase coverage, so their non-detections constrain little; J183950 is the
existence proof — ten epochs, one pulse). Epoch composition corrected (no RACS-low1 anywhere;
50/81 low-band epochs are the 2024 pass); "exact permutation" corrected to what it is; the
novelty claim narrowed to "first targeted, multi-epoch, uniform".

## Phase-fold check of the J183950 recovery (2026-08-14) — SUPERSEDED same day, see below

Folding all ten epoch MJDs (committed in results/lptv_realtargets.csv) on the catalogued
23,221.74-s period (data/lpt_sample.csv), relative to the detection epoch:

| MJD | phase rel. detection | I (mJy) |
|---|---|---|
| 59243.100 | 0.744 | 0.17 |
| 59615.100 | 0.930 | −0.07 |
| 59675.900 | 0.163 | −0.04 |
| 60312.200 | 0.791 | −0.15 |
| 60313.200 | 0.512 | 0.11 |
| **60335.100** | **0.000** | **164.03** |
| 60336.100 | 0.721 | −0.43 |
| 60352.100 | 0.256 | 0.37 |
| 60628.300 | 0.977 | −0.32 |

A ~15-min RACS integration is ±0.019 in phase. **No non-detection epoch samples the on-pulse
window** — the nearest (phase 0.977) misses it by ~9 minutes; the epoch one night after the
detection lands at phase 0.72. The single-pulse light curve is therefore fully consistent
with a periodic emitter sampled off-pulse nine times and on-pulse once, which is exactly what
a 6.45-h-period source with narrow pulses looks like in survey snapshots. (The epoch at
phase 0.512 — the interpulse phase if the published interpulse sits at 0.5 — shows nothing,
consistent with the interpulse being weaker/intermittent.) Caveat: this fold uses the
catalogued period with no published zero-point; a phase-connected comparison against the
discovery ephemeris (arXiv:2501.09133) is the decisive follow-up and is listed below.

## Research opportunities opened by the recovery

1. **Phase-connect against the published ephemeris** (discovery paper's T0 + P): if MJD
   60335.1 lands in a predicted pulse window, the recovery is conclusive beyond argument.
2. **VAST epochs**: the VAST survey carries many more ASKAP I+V epochs of this field than
   RACS. Extending the sweep to VAST obscore products could yield more pulses, a pulse-rate
   estimate, and a duty-cycle measurement — the natural lptv v2, using the existing forced
   machinery unchanged.
3. **Dilution arithmetic**: 164 mJy averaged over a ~15-min integration implies an intrinsic
   pulse flux of ~0.5–1.6 Jy for published minutes-scale pulse widths — directly comparable
   with the discovery paper's pulse fluxes, and a check that the recovery is energetically
   consistent.
4. **Handedness**: our V is +65 mJy (one sign, 40%); compare against the published
   polarization behaviour of pulse vs interpulse.


## The identity check supersedes the fold (2026-08-14, later)

Fetching the discovery paper's source (arXiv:2501.09133, Lee et al.) settles the recovery
beyond any statistical argument: its observation table lists **SBID 57929 — 2024-01-26
03:12:35, RACS-low, 887.5 MHz — as the discovery observation itself**, and that is the same
scheduling block as our detection epoch. Our census blindly re-found the discovery pulse.

Quantitative cross-checks, all consistent:
- Their published pulse: **0.70 Jy decaying to 0.03 Jy** across the 15-min block. Averaged
  over the integration that is ~0.2 Jy — against our block-integrated forced **I = 164 mJy**.
- Their published circular fraction: **37%** — against our **|V|/I = 40%**.
- Their period 23,221.740 ± 0.332 s matches the catalogue value our fold used, and with the
  paper's exact start time our epoch sits ~12 min from the predicted pulse centre, within
  their reported ToA jitter.

Two corrections to the earlier entry. First, the phase-fold table used the CSV's `epoch_mjd`,
which is quantized to 0.1 day (±72 min) — at a 6.45-h period that quantization swamps
minute-level phase statements, so the fold was decoration, not evidence; the identity makes
it unnecessary. Second, the framing "caught an unknown pulse" is wrong in the other
direction: this is a **blind end-to-end recover-a-known** — the pipeline, pointed at a
catalogue position with no knowledge of the epoch, re-found a published Jy-level pulse at
199σ with polarization intact, after its earlier peak-search configuration had wrongly
vetoed the same signal as confusion. That is the stronger and honest claim, and the paper
now makes it with the discovery paper cited (lee2025 added to refs.bib).

Consequence for the VAST follow-up idea: still live — the discovery paper's own archival
search (their Methods) found no additional archival pulses, so a VAST-epoch extension should
be framed as extending their archival null with V-sensitive forced photometry, not as
virgin territory.


## Round-2 referee (2026-08-14) — major revision, applied same day

Seventeen findings; verdict major. The campaign's lesson held for the second time running: the
rehabilitation paragraph written the same morning contained three of the errors ("blind" for a
test that cannot fail; numbers living only in a .bib comment; the beam-attenuated 164 mJy used
where the peak-search 241 mJy is the right comparator). The two blockers were mine from the
first rewrite: the abstract's rescoped-negative numbers (median I/sigma = 1.9; |V|/I limits
49–114%) still traced to the superseded peak-search CSV — on the forced CSV the median signed
I/sigma is −0.14 and the fractional-polarization constraint is vacuous (undefined for half the
sample) — and the offset argument for J165130.3's promotion was circular (under locked
photometry offset_arcsec is pixel rounding, bounded at half a pixel diagonal; 0.85″ was also
the survey median misattributed to one source, the row value being 0.44″).

Everything applied:
- Rescoped negative recomputed from the forced CSV and stated honestly: Stokes I consistent
  with zero at every non-detection position → the census constrains polarized *flux* only and
  places no |V|/I constraint on any undetected source (stronger and simpler than 49–114%).
- J165130.3 promotion re-based on recovered flux at the locked pixel (83% of peak-search I,
  95% of V) instead of the inert offset.
- J183950: "blind end-to-end validation" → guaranteed-not-blind recovery of the published
  pulse in its own discovery block (SBID 57929 identity), explicitly excluded from any
  pulse-rate reading; decay-average arithmetic shown (213 mJy exp. average vs 241 mJy
  peak-search, 13%); forced 164 mJy identified as beam-attenuated (0.68 ratio in both I and V
  — itself the single-point-source signature); beam fraction corrected (a third of ~15″, not
  a fifth); "matches" → "consistent with"; 199σ qualified as statistics-only. Evidence
  committed: results/lptv_j183950_adjudication.json.
- Confusion veto described as superseded (structurally inert under locked photometry);
  secure/candidate and handedness-flip language dropped as vacuous.
- Stale two-detection fragments fixed (Results opening, contribution sentence, figure caption
  — which promised detection markers the figure never had).
- Epoch composition: no RACS-low1 anywhere; abstract now names RACS-low2 + RACS-mid + the
  2023–24 low pass (the findings file had recorded this as fixed in round 1; the paper hadn't
  been).
- "Thirteen sets of limits" → twelve + one uncovered. Phase coverage 27–75% cumulative (the
  27–43% range had silently excluded J170036.6 at 74.6%), independence assumption stated,
  "single-visit" corrected (~4% per visit).
- Validation paragraph honest about being a smoke test far from decision boundaries; parent
  stokesv recovery flagged as 12″ peak-search mode. Reproducibility claim scoped to what is
  actually pipeline-generated. Instrument parameters cited at point of use.
- **Real catalogue error found by the referee's name-vs-coordinate sweep**: the ASKAP
  J1832−0911 row carried an erroneous name (J183244.5−091121, 59″ from its coordinate). The
  coordinate matches the published VLBA position (18:32:48.4589 −09:11:15.297, wang2025) to
  2″, so the nine-epoch limit is on the real source; the *name* was the transcription error.
  Corrected in data/lpt_sample.csv, results CSV/JSON, v_table.tex; disclosed in the paper.
  Catalogue coordinate quantization (3-decimal degrees = 3.6″ for 13/16 rows) now stated as a
  systematic.
- wang2025 and chime1634 now cited; wang2025 coordinates verified via Crossref (Nature 642,
  583–586).
