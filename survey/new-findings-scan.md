# Survey — slices with genuine *new-findings* potential (2026-06)

A focused scout (three cited `radio-research-assistant` syntheses) of candidate slices that could
yield **new findings** (not reproductions) while reusing the existing `spectra`/`peaked`/`vlass`
tooling. Honest verdicts; the realistic amateur unit of contribution is a vetted candidate list or a
value-added catalogue, not a guaranteed discovery.

## A. RACS Stokes-V search for coherent radio emitters (radio stars / M-dwarfs / UCDs / pulsars)
**Verdict: conditional GO — highest genuine-discovery ceiling.** Circular polarization is a near
-unambiguous flag of coherent emission; each new V-detected star is a real, citable find.

- **State of the art (crowded but not closed):** Pritchard+2021 (33 stars, RACS-low; arXiv:2102.01801),
  Pritchard+2024 (multi-epoch VAST; arXiv:2312.11031), Driessen+2024 Sydney Radio Star Catalogue
  (839 stars; arXiv:2404.07418), RACS-low2 Paper VIII just published 2026-06 (formal V catalogue, 61
  stars/85 pulsars; arXiv:2606.16182). ASKAP teams own the data + calibration.
- **Genuine gaps:** (1) **forced-photometry V survey of a curated M-dwarf/UCD input list** (LSPM,
  CARMENES, Reylé 10-pc) across RACS-low2 + RACS-mid V images — pushes below the blind 5σ threshold;
  (2) **RACS-low1 vs RACS-low2 V variability** (no published 2-epoch V transient comparison);
  (3) multi-band V spectral characterisation of known V stars.
- **Data:** CASDA (`astroquery.casda`, SODA cutouts; free OPAL account). All RACS V products public,
  no embargo. **Killer systematic:** off-axis Stokes-I→V **leakage** (≤1% on-axis, ~12% at beam
  edges; RACS-high holography lost). Must impose a per-beam leakage floor (7× median |V/I|).
- **Reuse:** forced peak photometry (CADC→CASDA, near-direct), cross-match, SIMBAD/NED vetting,
  two-point index — all reused. New: CASDA auth wrapper, leakage-floor estimator (~100 lines),
  proper-motion (Gaia DR3) confirmation.
- **Amateur unit:** ~5–20 new circularly-polarized stars/UCDs + upper-limit burst-rate stats.
- **Risk:** ASKAP teams publish RACS-mid/high standalone V catalogues within 1–2 yr; pick the
  forced-photometry-of-a-target-list or 2-epoch-V-variability angle to complement, not duplicate.

## B. Southern multi-band spectral-curvature catalogue (GLEAM-X + 3 RACS bands)  ← natural `peaked` sequel
**Verdict: GO (scoped) — highest tooling reuse; fixes our northern slice's core limitation.**

- **The structural advantage:** GLEAM-X DR2 gives **20 in-band sub-bands over 72–231 MHz**
  (Ross+2024, arXiv:2406.06921) → a **real, measured turnover**, not the TGSS upper-limit our northern
  `peaked` slice was forced into. Combined with RACS-low/mid/high (887.5/1367.5/1655.5 MHz; Duchesne+
  2023–2026) → 5–7 flux points over a ×23 frequency span, both sides of the SED for many sources.
- **Genuine gap:** no published **all-sky (Dec −80°→+30°) homogeneous GLEAM-X + all-3-RACS curvature
  catalogue.** Callingham+2017 (arXiv:1701.02771) predates RACS; GLEAM-X DR2 flags 18,869 curved
  sources but without RACS GHz anchor; RadioSED II (Kerrison+2025, arXiv:2509.21233) does it well but
  only over 300 deg² of Stripe 82; RACS-high paper defers curvature fitting to future work.
- **Data:** all public — GLEAM-X DR2 & GLEAM & SUMSS on VizieR; RACS via CASDA. No embargo.
- **Reuse:** `spectral_index`, `crossmatch`, `peak_frequency` (parabolic log-log fit), `classify_sed`
  reuse nearly verbatim; the ≥5-point SEDs let the parabolic fitter measure a real turnover and even
  support an SSA/FFA model. Adds a southern USS (α<−1.2) high-z-radio-galaxy candidate list for free.
- **Pitfalls:** resolution mismatch (GLEAM-X ~45″ vs RACS ~10–25″ → compactness/single-component cut),
  ~10% cross-survey flux-scale ties, low-ν ionospheric position errors (use per-source errors).
- **Amateur unit:** a value-added **scout catalogue** of southern GPS/CSS/HFP + USS candidates
  (~10⁵ multi-band SEDs) — interim all-sky-south list until the GLEAM-X team's version.
- **Risk:** Kerrison/Ross/Murphy group may publish a full-sky version in 6–12 months; an independent
  reproducible pipeline still has cross-check/citation value, and it directly upgrades our own slice.

## C. Multiwavelength radio stacking / forced photometry (Gaia/eROSITA/WISE × VLASS/RACS)
**Verdict: NO-GO on VLASS Quick-Look; only conditional GO with VLASS *Single-Epoch* or RACS + an
injection-recovery simulation.** Lowest reuse, highest competition risk.

- **Killer systematic:** for sub-threshold sources VLASS QL flux sits in the **dirty residual**, not
  CLEAN-restored — a 20–40% residual/restored PSF bias on top of the 10–15% QL scale error, and the
  closest precedent (Perger+2023, arXiv:2311.01128) doesn't model it. A referee will require an
  injection-recovery calibration. (Note: the "Bonnassieux 2023 ~1.4× VLASS bias" from the older scan
  is **unverified**; the verified anchors are White+2007 for FIRST and the VLASS QL Users' Guide.)
- **Better substrate:** VLASS Single-Epoch images (self-cal, deeper CLEAN, public via CADC SODA) or
  RACS-mid (clean PAF imaging) largely remove the QL pathology.
- **Reuse:** cutout fetch + forced photometry reused; but needs a real image-plane median-stack,
  injection-recovery sim (~500 lines, the gating item), confusion/rms estimation, bright-neighbour mask.
- **Risk:** eROSITA/Gaia/WISE × radio stacking is the obvious next paper for big collaborations
  (MPE/CSIRO/CIRADA) and likely already in prep — competitive window measured in months.

## Recommendation
- **For the most genuine *new-findings* (discovery) ceiling → A (RACS Stokes-V), forced-photometry of
  a curated M-dwarf/UCD target list or the 2-epoch V-variability angle.** Each detection is a real find.
- **For the cleanest continuation of our work + highest reuse + a real measured turnover → B (southern
  GLEAM-X×RACS curvature).** It literally fixes the upper-limit limitation `peaked` documented.
- **C is parked** unless we move to VLASS-SE/RACS and commit to an injection-recovery harness.

Runners-up from the earlier `data-source-scan.md` not re-scouted here: Astrogeo VLBI multi-decade AGN
variability (HIGH upside, new toolchain) and Wind/STEREO WAVES burst re-classification (solar).
