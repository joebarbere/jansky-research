# Before/after worked conversions (plan 89)

Three conversions from this repo's own papers, applying the strip (S1–S11) and add (A1–A13)
rules of `style-guide.md`. Every macro invocation, `\cite` key, and numeric literal is
preserved exactly; only prose changes. Target register: modern corpus (1990s–2021,
AASTeX/MNRAS). Abstract length target: the corpus median (fingerprints p50 = 163 words,
mean 154) against our selfscan p50 of 247.

---

## Conversion 1 — dr20radio abstract

### BEFORE (verbatim, `papers/dr20radio/main.tex` lines 15–44; 374 words)

```latex
SDSS-V Data Release 20 delivers the Black Hole Mapper's spectroscopically confirmed quasar
sample, including the first optical SDSS spectra ever taken from the southern hemisphere
(BOSS on the du~Pont at Las Campanas). Because VLASS stops at $\delta = -40^\circ$ and
southern SDSS quasar spectra did not exist before DR20, the pairing ``SDSS spectroscopic
quasars $\times$ RACS south of $-40^\circ$'' was impossible until this release; no radio
cross-match of any SDSS-V quasar catalog has been published. I present the first such
census. North of $-40^\circ$, \drNorthCensus\ clean-redshift DR20 quasars matched against
the VLASS Quick Look epoch catalogs at \drAlphaVlassRad\arcsec\ give per-epoch radio-detection fractions
of \drNorthEtwoPct\% (E2) and \drNorthEthreePct\% (E3) after subtracting a
\emph{measured} chance-coincidence rate (\drNorthFmPct\% from position-shift trials), or
\drNorthAnyPct\% matched in any epoch (raw, uncorrected). South of $-40^\circ$ --- sky where these are the
first SDSS quasar spectra --- \drSouthCensus\ quasars (all observed from Las Campanas)
matched against RACS-low DR1 at 5\arcsec\ give \drSouthPct\% (chance rate \drSouthFmPct\%).
Radio-targeted open-fiber cartons are excluded from all fractions as circular; used instead
as a positive control, RACS-selected cartons recover at \drCartonRacsPct\%
($n=\drCartonRacsN$) against their selecting survey, while LOFAR-selected cartons recover
at only \drCartonLofarPct\% ($n=\drCartonLofarN$) at 888\,MHz, and at 3\,GHz VLASS
recovers RACS-selected cartons at \drCartonVlassRacsPct\% and LOFAR-selected at
\drCartonVlassLofarPct\% --- recovery rates for the same targeting classes across
surveys differing in both frequency and depth. Comparing the two hemispheres above a common
rest-frame 1.4\,GHz luminosity limit requires a radio spectral index, so I measure it from
quasars the surveys share rather than assuming one: $\alpha \simeq -0.6$ to $-0.9$, with the
spread set by which survey selects the sample rather than by measurement error. The
comparison reduces algebraically to a pair of flux cuts, so a difference in percentage points
inherits the overall normalisation and falls whenever either limit is raised; the ratio does
not, and is \drRatioMeas --- stable to \drRatioLimSpreadPct\% across every flux-limit variant,
and uncertain over \drRatioAlphaLo--\drRatioAlphaHi\ from the index alone. It is a property
of this catalogue--survey pairing, consistent with survey-side effects (limit definitions,
beam-size mismatch, hemisphere targeting mix), and I attribute none of it to astrophysics.
```

### AFTER (~185 words; every `\drXxx` macro and every number retained)

```latex
We present the first radio cross-match of an SDSS-V Black Hole Mapper quasar catalogue,
including the first SDSS quasar spectra taken from the southern hemisphere (Data Release 20).
North of $\delta = -40^\circ$, \drNorthCensus\ clean-redshift quasars matched against the
VLASS Quick Look epoch catalogs at \drAlphaVlassRad\arcsec\ give radio-detection fractions of
\drNorthEtwoPct\% (E2) and \drNorthEthreePct\% (E3) after subtraction of a measured
chance-coincidence rate (\drNorthFmPct\% from position-shift trials), or \drNorthAnyPct\%
matched in any epoch (raw, uncorrected). South of $-40^\circ$, where no SDSS spectra
previously existed, \drSouthCensus\ quasars matched against RACS-low DR1 at 5\arcsec\ give
\drSouthPct\% (chance rate \drSouthFmPct\%). Radio-targeted open-fiber cartons are excluded
from all fractions and used as a positive control: RACS-selected cartons recover at
\drCartonRacsPct\% ($n=\drCartonRacsN$) against their selecting survey and at
\drCartonVlassRacsPct\% in VLASS at 3\,GHz, while LOFAR-selected cartons recover at
\drCartonLofarPct\% ($n=\drCartonLofarN$) at 888\,MHz and \drCartonVlassLofarPct\% at
3\,GHz. The radio spectral index measured from quasars the two surveys share is
$\alpha \simeq -0.6$ to $-0.9$, with the spread dominated by sample selection. Above a
common rest-frame 1.4\,GHz luminosity limit the north/south detection ratio is \drRatioMeas,
stable to \drRatioLimSpreadPct\% against the flux-limit choices but spanning
\drRatioAlphaLo--\drRatioAlphaHi\ over the index range. We attribute this ratio to
survey-side effects (limit definitions, beam-size mismatch, targeting mix) and claim no
hemispheric difference in quasar radio properties.
```

### Rationale, line by line

1. Opening two sentences of framing ("delivers...", "Because VLASS stops at... was
   impossible until this release; no radio cross-match... has been published. I present the
   first such census.") → one "We present the first radio cross-match of..." sentence.
   **S6** (result-first abstract; framing moves to the Introduction, which already carries
   it) and **A1** (the "We present <data> of <object>" opener — cf. ApJS 229: "We present a
   catalog of hard X-ray sources..."). The "first" claim survives — factual firsts are
   corpus-legal (1995MNRAS.272L...1G) — but stated once, not three times.
2. "I present" → "We present"; closing "I attribute" → "We attribute". **§4 voice rules**:
   corpus single authors use authorial "we"; dr20radio's body-"I" (2.63/kw) sits above the
   corpus p90 (1.44). (If the author prefers to keep "I" as a signature, that is defensible —
   but it is the above-p90 choice, so the conversion shows the median one.)
3. "a \emph{measured} chance-coincidence rate" → "a measured chance-coincidence rate".
   **S2**: corpus emph-for-stress is zero in 2016–21; the word "measured" carries its own
   weight next to "(from position-shift trials)".
4. "South of $-40^\circ$ --- sky where these are the first SDSS quasar spectra ---" →
   "South of $-40^\circ$, where no SDSS spectra previously existed,". **S1** (em-dash
   parenthetical → comma clause) and **A3**; also removes the third repetition of "first".
5. The carton sentence: the trailing em-dash coda ("--- recovery rates for the same
   targeting classes across surveys differing in both frequency and depth") is deleted;
   the four rates are regrouped by targeting class so the sentence needs no gloss. **S1**,
   **S10** (let the numbers be the punch). All four macros and both $n$ macros kept.
6. "Comparing the two hemispheres ... requires a radio spectral index, so I measure it from
   quasars the surveys share rather than assuming one" (method narration) → "The radio
   spectral index measured from quasars the two surveys share is ...". **S6**: method
   rationale out of the abstract; the *result* ($\alpha \simeq -0.6$ to $-0.9$) stays.
   "with the spread set by which survey selects the sample rather than by measurement
   error" → "with the spread dominated by sample selection" — same claim, half the words.
7. "The comparison reduces algebraically to a pair of flux cuts, so a difference in
   percentage points inherits the overall normalisation and falls whenever either limit is
   raised; the ratio does not, and is \drRatioMeas ---" → "the north/south detection ratio
   is \drRatioMeas, stable to ... but spanning ...". **S6** again: the algebraic argument is
   Section 3.4's job; the abstract keeps the number, its two stability figures, and the
   attribution. The em-dash before "stable" becomes a comma (**S1**).
8. "It is a property of this catalogue--survey pairing, consistent with..., and I attribute
   none of it to astrophysics." → "We attribute this ratio to survey-side effects (...) and
   claim no hemispheric difference in quasar radio properties." Same honest negative, in the
   corpus's committed register (**S8**: one clear commitment, no re-hedging; the parenthesis
   holds the list, **A3**). "claim no X" is the corpus's own form ("we do not claim the
   detection ... as strongly", 1995MNRAS.274..701G).
9. Length: 374 → ~185 words. Still above the corpus p50 (163) because all 19 macro-borne
   numbers are contractually retained; it now sits inside the corpus p50–p75 band (163–226)
   instead of at our selfscan p90.

---

## Conversion 2 — dr20radio Results passage (§4.3)

### BEFORE (verbatim, lines 260–291)

```latex
\subsection{The cross-hemisphere comparison, and what it is worth}
\label{sec:contrast}

Above the common limit (Section~\ref{sec:method-lum}) the fractions are \drLumNorthMeasPct\%
in the north and \drLumSouthMeasPct\% in the south at the measured index. That difference is
the weakest number in this paper, and three properties bound how it should be read.

\emph{It is a ratio, not a difference.} Since the comparison is a pair of flux cuts, raising
either survey's limit deepens both and shrinks both fractions; a difference in percentage
points inherits that normalisation and falls, while the ratio does not. Across every limit
variant tried --- the conservative \drSlimRacsCons\,mJy RACS limit and VLASS limits of
\drSlimVlassConsA\ and \drSlimVlassConsB\,mJy --- the ratio moves only
\drRatioLimLo--\drRatioLimHi\ (\drRatioLimSpreadPct\%) while the percentage-point gap moves
by a quarter. The contrast is \drRatioMeas, stable against both surveys' flux-limit choices.

\emph{Its uncertainty is the spectral index, and that index is selection-dependent.} Over the
index's range (Section~\ref{sec:method-alpha}) the ratio spans
\drRatioAlphaLo--\drRatioAlphaHi; at $\alpha = 0$ it would be \drRatioAlphaFlat. Nothing else
in the analysis moves it comparably: the population scatter in $\alpha$, tested by giving
every quasar its own index, shifts each fraction by $\sim$0.3 percentage points and their
ratio hardly at all.

\emph{What remains is not attributable to the sky.} The northern leg is counted down to
$1.2\times$ its catalogue's own threshold while the southern is counted at a 95\%
completeness limit; RACS's 25\arcsec\ beam blends and shifts faint sources against a
5\arcsec\ radius, the same effect visible in the \drCartonRacsPct\% positive control; the BHM
targeting mix differs by hemisphere; and \drNorthOutsideRacsPct\% of the northern census lies
outside the RACS footprint entirely, so the two legs are not matched-sky samples. I therefore
report \drRatioMeas\ as a property of this catalogue--survey pairing, and claim no
hemispheric difference in quasar radio properties. Readers wanting a cross-hemisphere
measurement rather than a cross-survey one should wait for RACS-mid and VLASS at comparable
depth and resolution.
```

### AFTER

```latex
\subsection{Cross-hemisphere comparison above the common limit}
\label{sec:contrast}

Above the common limit (Section~\ref{sec:method-lum}) the fractions are \drLumNorthMeasPct\%
in the north and \drLumSouthMeasPct\% in the south at the measured index. Because the
comparison is a pair of flux cuts, raising either survey's limit deepens both cuts and
shrinks both fractions, so a difference in percentage points inherits the overall
normalisation and falls, while the ratio does not. We therefore quote the contrast as a
ratio. Across every limit variant tried (the conservative \drSlimRacsCons\,mJy RACS limit
and VLASS limits of \drSlimVlassConsA\ and \drSlimVlassConsB\,mJy) the ratio moves only
\drRatioLimLo--\drRatioLimHi\ (\drRatioLimSpreadPct\%), while the percentage-point gap moves
by a quarter. The contrast is \drRatioMeas, stable against both surveys' flux-limit choices.

The dominant uncertainty is the spectral index, which is itself selection-dependent
(Section~\ref{sec:method-alpha}): over the index's range the ratio spans
\drRatioAlphaLo--\drRatioAlphaHi, and at $\alpha = 0$ it would be \drRatioAlphaFlat.
Nothing else in the analysis moves it comparably. The population scatter in $\alpha$,
tested by giving every quasar its own index, shifts each fraction by $\sim$0.3 percentage
points and their ratio hardly at all.

Several survey-side effects remain and are not separable here. The northern leg is counted
down to $1.2\times$ its catalogue's own threshold while the southern is counted at a 95\%
completeness limit; RACS's 25\arcsec\ beam blends and shifts faint sources against a
5\arcsec\ radius, the same effect visible in the \drCartonRacsPct\% positive control; the
BHM targeting mix differs by hemisphere; and \drNorthOutsideRacsPct\% of the northern
census lies outside the RACS footprint entirely, so the two legs are not matched-sky
samples. We therefore report \drRatioMeas\ as a property of this catalogue--survey pairing,
and claim no hemispheric difference in quasar radio properties. A cross-hemisphere
measurement, as opposed to a cross-survey one, will require RACS-mid and VLASS at
comparable depth and resolution.
```

### Rationale

1. Title "The cross-hemisphere comparison, and what it is worth" → "Cross-hemisphere
   comparison above the common limit". **S4**: corpus titles are noun phrases naming the
   operation (median 2.9 words; "Zero witty/conclusory section titles", 2016–21).
2. "That difference is the weakest number in this paper, and three properties bound how it
   should be read." — deleted. **S3** (self-referential epistemics, our 2.24/kw vs corpus
   0.55/kw). Its content is not lost: the weakness is *shown* by the two paragraphs that
   attach the cause (index selection-dependence) and the magnitudes
   (\drRatioAlphaLo--\drRatioAlphaHi; \drRatioLimSpreadPct\%) — the corpus form ("should be
   treated only as lower limits because...", arXiv:2106.16211).
3. The three "\emph{...}." run-in mini-headings → plain topic sentences ("We therefore quote
   the contrast as a ratio." / "The dominant uncertainty is the spectral index..." /
   "Several survey-side effects remain..."). **S2**: sentence-initial \emph has corpus
   p90 = 0; the corpus carries paragraph logic in topic sentences, not typography.
4. Both em-dash asides ("--- the conservative ... ---") → parentheses. **S1/A3**: the aside
   holds data, so it is a parenthesis (corpus asides are "units, acronym definitions, refs
   and values").
5. "I therefore report" → "We therefore report" (**§4**), and "Readers wanting a
   cross-hemisphere measurement ... should wait for" → "A cross-hemisphere measurement ...
   will require RACS-mid and VLASS at comparable depth and resolution." **S9** (reader
   address, corpus p90 = 0) + **A10** (end on what data would settle it: "Higher-resolution
   observations are desired to distinguish...", arXiv:2004.09369).
6. Unchanged: every macro (\drLumNorthMeasPct, \drLumSouthMeasPct, \drSlimRacsCons,
   \drSlimVlassConsA, \drSlimVlassConsB, \drRatioLimLo, \drRatioLimHi, \drRatioLimSpreadPct,
   \drRatioMeas, \drRatioAlphaLo, \drRatioAlphaHi, \drRatioAlphaFlat, \drCartonRacsPct,
   \drNorthOutsideRacsPct), every number ($1.2\times$, 95\%, 25\arcsec, 5\arcsec,
   $\sim$0.3, $\alpha = 0$), both `\ref`s, the label, and every claim including the
   honest negative ("claim no hemispheric difference").

---

## Conversion 3 — lptv abstract opening

### BEFORE (verbatim, `papers/lptv/main.tex` lines 13–27)

```latex
Long-period radio transients (LPTs) --- sources pulsing on minute-to-hour periods
\citep{hurleywalker2022,caleb2024} --- are a young class whose progenitors span neutron stars and
white-dwarf binaries \citep{rea2026}. We extend a provenance-carrying LPT catalogue to \lvNLpt\
confirmed members with the three 2026 discoveries \citep{deruiter2026,vaster2026}, each coordinate
verified against its source-name convention. At \lvNLpt\ members the hinted $\sim$78-minute period
boundary between white-dwarf-binary and other LPTs is \textbf{still not significant} (permutation
$p = \lvPeriodSplitP$). We then perform what has not been done before: systematic multi-epoch forced
Stokes-V (circular-polarization) photometry at every LPT position across the RACS
\citep{mcconnell2020} epochs that cover it: RACS-low2 \citep{racslow2}, RACS-mid, and the
2023--24 low-band pass.
```

### AFTER

```latex
Long-period radio transients (LPTs) are sources pulsing on minute-to-hour periods
\citep{hurleywalker2022,caleb2024}, a young class whose proposed progenitors span neutron
stars and white-dwarf binaries \citep{rea2026}. We extend a provenance-carrying LPT
catalogue to \lvNLpt\ confirmed members with the three 2026 discoveries
\citep{deruiter2026,vaster2026}, each coordinate verified against its source-name
convention. At \lvNLpt\ members the suggested $\sim$78-minute period boundary between the
white-dwarf-binary members and the rest remains not significant (permutation
$p = \lvPeriodSplitP$). We also present the first systematic multi-epoch forced Stokes-V
(circular-polarization) photometry at every LPT position, using the RACS
\citep{mcconnell2020} epochs that cover each source: RACS-low2 \citep{racslow2}, RACS-mid,
and the 2023--24 low-band pass.
```

### Rationale

1. "LPTs) --- sources pulsing ... --- are a young class" → "LPTs) are sources pulsing on
   minute-to-hour periods ..., a young class ...". **S1**: the definitional em-dash pair
   becomes an ordinary appositive; the sentence is now the corpus's plain-background opener
   (**A1** alternative form — "Wolf-Rayet (WR) stars are evolved massive stars, presumably
   on their way to becoming supernova", arXiv:2008.03725). All three `\citep` keys keep
   their positions relative to the claims they support.
2. "\textbf{still not significant}" → "remains not significant". **S2**: bold never appears
   in corpus prose in any era (1960s notes: "Bold face never appears in prose"); the verb
   "remains" carries "still". The number ($p = \lvPeriodSplitP$, macro exact) is the
   emphasis.
3. "the hinted ... boundary between white-dwarf-binary and other LPTs" → "the suggested ...
   boundary between the white-dwarf-binary members and the rest". "Suggested" is the corpus
   hedge verb (**A5** family); grammatical smoothing only, claim identical.
4. "We then perform what has not been done before: systematic multi-epoch forced Stokes-V
   ... photometry" → "We also present the first systematic multi-epoch forced Stokes-V ...
   photometry". **S7**: the novelty drumbeat ("what has not been done before") becomes the
   corpus's factual-first form ("the first study of polarized submillimetre emission from
   the Sgr B2 giant molecular cloud", 1995MNRAS.272L...1G). One "first", stated as a
   property of the data set, not of the authors' daring. The colon-headline construction
   ("We then perform X: Y") also falls under **S10**.
5. "across the RACS epochs that cover it" → "using the RACS epochs that cover each source"
   — repairs the dangling "it" while the double-colon chain ("...that cover it: RACS-low2
   ..., RACS-mid, and...") is reduced to one colon introducing the list, the corpus's
   grammatical-list form (**A3**; lists "introduced grammatically", 1980s notes).
6. Unchanged: \lvNLpt (twice), \lvPeriodSplitP, all six citation keys, $\sim$78-minute,
   2026, 2023--24 (an en-dash range, which is corpus-correct and untouched), and the
   honest negative itself — the non-significance stays in the abstract, per §5 of the
   style guide ("validations, limits, and negatives stay; only their voice changes").
```
