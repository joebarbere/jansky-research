# Traditional-style conversion guide (plan 89)

Empirical basis: nine era note files (`data/style_corpus/stylenotes/1933-1949.md` …
`2016-2021.md`, ~70 papers read in full with verbatim quotes), quantitative fingerprints from
1,137 pre-LLM LaTeX papers (`results/stylecorpus_fingerprints.json`), and a self-scan of this
repo's 45 papers (`results/stylecorpus_selfscan.json`). Every rule below carries a corpus
number or a verbatim quote with its bibcode/arXiv id. Per-1000-word rates are written "/kw".

---

## 1. The target register

**Era-invariant principles (hold from 1934 to 2021):**

- **Openings state facts about nature, instruments, or the field — never the paper's
  importance.** 1947: "Radio-frequency radiation originating outside the earth's atmosphere
  was first discovered by Jansky at a frequency of 18 megacycles per second" (1947ApJ...105..235).
  2020: "Wolf-Rayet (WR) stars are evolved massive stars, presumably on their way to becoming
  supernova" (arXiv:2008.03725, ApJL 900 L3). No era has an "In recent years X has attracted
  considerable attention" opener, and "novel" was never observed ("'new' appears... 'novel'
  never" — 1990s notes, 9 papers).
- **Measurements are stated flat; hedging is reserved for interpretation and graded per
  claim.** "The mean rms deviation in position in one hour ... is found equal to
  0.64 ± 0.20 minutes of arc" vs "it seems as if the magnitude of the positional fluctuations
  is larger when the intensity is relatively small" (1958ApNr....6...75M). 2016–21: "hedges
  attach to the interpretation, never to the number." The ladder is graded — "clearly show" >
  "strong indications" > "possibly" in one abstract sentence (1995MNRAS.274..701G) — and
  asymmetric: pick "very likely" or "we can say nothing", and say which.
- **Limitations are woven in-line, with magnitudes, at the point of use.** "Together all these
  biases may lead to underestimate the estimated parameter up to a factor of $\sim$10"
  (arXiv:2106.16211, MNRAS 507). Zero papers in any era have a "Limitations" section.
- **Plain connectives repeat unashamedly.** MNRAS 507 uses "However" 18 times; "nobody varied
  it." The workhorses are However / Thus / Hence / Also / In addition / Note that.
- **Negatives and failures are first-class results.** A table titled "Negative results of H₂O
  search" (1973A+A....26..487L); "Unfortunately, despite this adopted strategy, the first
  mini-survey ... was highly contaminated by ghost rays" (arXiv:1703.00021, ApJS 229).
- **The paper ends on the result or on what data would settle it, not on significance.**
  "Higher-resolution observations are desired to distinguish..." (arXiv:2004.09369, ApJ 895).

**The modern target specifically (1990s–2021, AASTeX/MNRAS register — what conversions aim at):**

- Abstract: one paragraph, corpus median 163 words (mean 154; p75 = 226; 2016–21 median ~220
  and 7–11 sentences). Opens "We present/report <data> of <object> with <instrument>";
  numbers with uncertainties by sentence 3; one hedged interpretive claim at the end.
- Sections: canonical Introduction / Observations (or Data) / Results / Discussion /
  Conclusions; titles are noun phrases, median 2.9 words (fingerprints
  `section_title_mean_words` p50 = 2.86).
- Voice: "we" for decisions and claims, passive for procedure — and both at HIGHER rates than
  our papers (see §4). Sentences long and additive (corpus median sentence 24–26 words), with
  high variance: 50-word evidence sentences beside short verdicts.
- Punctuation: the aside is a parenthesis holding data — "(hereafter Paper I)", "($1''$
  corresponds to 85 pc)" — not an em-dash. Em-dashes ≈ zero (see strip list #1).
- Texture: roadmap paragraph at the end of the intro (7 of 8 papers in 2016–21); "Note that"
  as the caveat marker; "(private communication)" with a name as a citable source;
  conclusions as numbered findings each restating its own number ± error.

---

## 2. Strip list — LLM tells, ordered by diagnostic strength

Each entry: (a) corpus evidence, (b) replacement, (c) verbatim corpus counter-example.

### S1. Em-dash parentheticals and pivots  [strongest single tell]
(a) Corpus `em_dash_per_kw`: p50 = 0, p75 = 0, p90 = 0.66 across 1,137 papers — **the median
pre-LLM paper contains zero em-dashes**. 2016–21 notes: "effectively zero across 68,000 words
of prose." Our papers: mean 10.7/kw, p50 10.7 — 16× the corpus p90; dr20radio 10.2, frblens 14.5.
(b) Replace with a parenthesis (if the aside is data), a comma clause, or two sentences. The
only period-attested appositive break is a *spaced double hyphen*, used a handful of times per
paper at most.
(c) "the so-called ``dirty'' image -- contains artifacts" (arXiv:1912.04970, PASP 132);
"compact stellar remnants $-$ white dwarfs, neutron stars, and black holes $-$" (2016–21 corpus).

### S2. \emph for rhetorical stress; \emph-opened sentences
(a) Corpus `emph_per_kw`: p50 = 0, p90 = 0.92; `emph_sentence_start_per_kw` p90 = 0 — a
sentence *opening* with \emph is effectively never observed. 2016–21: "`\emph`: literally zero
in all eight papers"; italics exist only for mission names (*NuSTAR*), coined terms, variables.
Our papers: emph 7.9/kw mean; emph-sentence-starts 1.72/kw mean (p50 1.62).
(b) Delete the emphasis, or restructure so the stressed word lands by position; use a plain
topic sentence in place of an italic run-in mini-heading; for a caveat, "Note that ...".
(c) The rare legitimate use is semantic, on a load-bearing physics word: "\emph{volume-averaged}",
"\emph{internal} depolarization" (arXiv:1210.4237, ApJ 766). Bold in prose is worse: "Bold face
never appears in prose" (1960s notes) — lptv's `\textbf{still not significant}` has no corpus
precedent in any era.

### S3. Self-referential epistemics (the paper rating its own numbers)
(a) Corpus `self_ref_per_kw` p50 = 0.55; ours p50 = 2.24 — 4× the corpus median. No corpus
paper meta-narrates its own evidence ("That difference is the weakest number in this paper",
dr20radio §4.3). Corpus limits name a **cause and a magnitude**, not a ranking.
(b) Replace the meta-sentence with the cause-attached caveat itself, at the point of use.
(c) "the outflow parameters estimated here should be treated only as lower limits because it
might include several observational biases" (arXiv:2106.16211); "we use the cataloged RCS
values as an order of magnitude guide only" (arXiv:2006.04327, PASA 37). Even abstract-level
self-criticism is about the *fits*, not the paper: "the fits are not particularly good or
unique and are very model dependent" (1995MNRAS.273..812H).

### S4. Editorializing section titles
(a) 2010–15: "Titles are nouns naming the object, data, or operation. Editorializing titles
are essentially absent." 2016–21: "Zero witty/conclusory section titles." Corpus title median
2.86 words. Our titles like "The cross-hemisphere comparison, and what it is worth" and
"Southern census (RACS): the first of its kind" (dr20radio) have no corpus analogue.
(b) Noun phrase naming the object, quantity, or operation. The maximum attested editorial
content is a hedge or a genuine question.
(c) "Determination of the period of the binary producing the pinwheel" (arXiv:2008.03725);
"Hints of superorbital variability of LS I +61 303 in hard X-ray" (arXiv:1402.6159); question
headings exist but are real open questions: "What is a dwarf galaxy?" (astro-ph/0006290),
"Is HS 240 an interstellar bubble?" (1989A+A...221..311W).

### S5. Dedicated Limitations sections
(a) "No paper has a 'Limitations' or 'Caveats' section" — stated independently in the 2000s
("zero instances"), 2010–15, and 2016–21 notes; nothing like it in 1933–1995 either.
dr20radio has `\section{Limitations}`.
(b) Dissolve into the sections where each limit bites, each with its magnitude; at most one
closing scope paragraph (2–3 sentences) in the Conclusions.
(c) "Although only two sources are studied in this work ... Further comparison of embedded
multiple protostellar systems is needed to confirm these results." (arXiv:1805.05205, A&A 617);
"There are some limiting factors, however." as a mid-paper paragraph opener (arXiv:1001.4731).

### S6. Methodology-narrating 250+-word abstracts
(a) Corpus `abstract_words` p50 = 163, p75 = 226, p90 = 273. Ours: p50 = 247, mean 280, p90 =
381; dr20radio 374, lptv ~490. Half our abstracts sit above the corpus 90th percentile.
(b) Result-first: "We present <data> of <object> with <instrument>"; a number with its error by
sentence 3; framing and method rationale move to the Introduction; end on the one interpretive
claim, hedged, or on the stated limit of what can be concluded.
(c) "We present the results of the first study of polarized submillimetre emission from the
Sgr B2 giant molecular cloud. ... Detections were obtained for six out of 10 positions studied,
with percentage polarizations of 0.8–2.6 per cent." (1995MNRAS.272L...1G).

### S7. Marketing adverbs and novelty claims
(a) Absent in every era file: no "novel", "state-of-the-art", "robust" (as praise),
"crucially/notably/importantly/interestingly," as sentence adverbs, no "underscores/highlights
the importance of", no "paves the way". 1970s: the strongest rhetorical figure in nine papers
is one "not only exists, but dominates" (1973A+A....25..303S).
(b) State the fact that makes it new; guard superlatives with scope; attach "crucial"-class
words only to a concrete operation.
(c) "To our knowledge, these wings correspond to the fastest molecular outflow ever observed
in the Galaxy" (1989A+A...222L...1C); "Here we report the most secure identification of a
central image" — not "the first", because a competitor exists and is discussed by name
(astro-ph/0312136, Nature 427).

### S8. Symmetric blanket hedging
(a) Every era file independently: "No symmetric both-sides hedging ('While X, it is important
to note that Y')". 2016–21: "LLM hedges are evenly distributed politeness; era hedges are
lumpy — zero on measurements, doubled ... on exactly the claims the authors doubt." Note our
papers actually *under*-hedge by modal count (hedge_per_kw ours 0.29 vs corpus 0.45): we state
everything flat and then meta-comment (S3) instead of hedging the interpretive sentences.
(b) Zero hedges on numbers; one graded, cause-attached hedge per interpretive claim; commit
when the evidence allows.
(c) "There is no doubt, however, that, on average, the spectrum steepens outwards along the
jet" beside "we doubt its reality and label it in X" in the same literature
(1989A+A...217...44F; 1989A+A...216...31T); "This however is not conclusive, since the
observations were not performed simultaneously." (arXiv:1004.3058, AJ 139).

### S9. Reader address and rhetorical steering
(a) Corpus `reader_addr_per_kw` p90 = 0 (mean 0.019). dr20radio: "Readers wanting a
cross-hemisphere measurement ... should wait for RACS-mid" (0.33/kw, above the corpus p90).
(b) Recast as a statement about what data would settle it, or — in calibration/survey papers
only — as advice to *users of the data*, which is period-attested.
(c) "Observers wishing to use these data for calibration of observations taken outside of the
2001–2008 epoch should keep in mind the variability and related geometric issues described
above." (arXiv:1001.4731, ApJS 192).

### S10. Rule-of-three constructions and punchy verdict rhetoric
(a) 2016–21: "no rhetorical triads; no punchy short summary sentences ('The result is
striking.'); no colon-headline sentences ('The implication is clear: ...')". Corpus
`rule_of_three` p90 = 0. (Short sentences exist — "This is quite likely." (1995MNRAS.275..309E)
— but they carry a *judgment about the science*, not applause for the paper's own result.)
(b) Let the number be the punch; join genuinely serial facts with "and" or number them.
(c) "Low mass dwarfs exists, but are rare." (astro-ph/0006290 — a four-word verdict *after*
the evidence, typo shipped and all).

### S11. Elegant variation of connectives
(a) 2016–21: "These papers repeat 'However'/'Thus'/'Also' verbatim dozens of times; an LLM
rotates Moreover/Furthermore/Additionally/Nevertheless." 1970s: "Almost no signpost adverbs...
never 'Notably,' 'Importantly,' 'Crucially,'".
(b) Reuse However / Thus / Hence / Also / In contrast / Note that without shame.
(c) MNRAS 507 (arXiv:2106.16211): "However" ×18 in one paper.

---

## 3. Add list — traditional habits an LLM omits

- **A1. The abstract opener + numbers by sentence 3.** "We present a catalog of hard X-ray
  sources in a square-degree region surveyed by \textit{NuSTAR} in the direction of the Norma
  spiral arm. ... Twenty-eight sources are firmly detected and ten are detected with low
  significance" (arXiv:1703.00021). Alternative attested opener: 1–2 sentences of plain
  background fact, then the "we" sentence (arXiv:2008.03725).
- **A2. The roadmap paragraph.** Present in **7 of the 8** 2016–21 papers (absent only in the
  4-page ApJL) and 6 of 8 in 2010–15. Verbatim: "The paper is organized in the following
  manner. In Section~2, we present a brief introduction to the data. ... A conclusion of this
  study is presented in Section~5." (arXiv:2106.16211). Passive form equally attested: "In
  Section 2 the observations and their calibration are described." (arXiv:1010.3790).
- **A3. Parenthetical asides instead of dashes.** Corpus parentheticals run 9.5–25.5/kw
  (2016–21 table) and hold *data*: "(hereafter Paper I)", "(FOV)", "($1''$ corresponds to
  85 pc)", "(no fit involved)". Move every em-dash aside into one of these or into its own
  sentence.
- **A4. Agentive belief verbs — used once, then committed.** "While this revised limit is
  larger than previously reported, we believe it to be more robust and still represents the
  best current constraint" (arXiv:1301.5906, MNRAS 433); "we feel that our final catalog ...
  should be extremely reliable" (astro-ph/0401133, AJ 127); "is thought to be fortuitous"
  (arXiv:0709.3873, PASA 24). One belief verb where evidence is genuinely sub-decisive; no
  re-hedging afterwards.
- **A5. Causal hedges.** The hedge names its cause: "This however is not conclusive, since the
  observations were not performed simultaneously." (arXiv:1004.3058); "The brightening of C3
  at the epoch 1984.30 may not be real, because of the contamination of sidelobes on the map."
  (1989A+A...216...31T).
- **A6. The we/passive voice split** — "we" for decisions and claims, passive for procedure:
  "we performed simulations of stray light contamination and focused our observations on
  three areas" vs "The data were calibrated in the standard manner using the CASA ... package"
  (arXiv:2008.03725). See §4 for the rates.
- **A7. Woven-in limitations with magnitudes.** "Together all these biases may lead to
  underestimate the estimated parameter up to a factor of $\sim$10" (arXiv:2106.16211); "The
  large fractional uncertainty of ∼±31% in the parallax dominates the fractional uncertainty
  in the radius... the luminosity ... is of little value and has not been listed."
  (arXiv:0709.3873).
- **A8. Plain repeated connectives and "Note that".** "Note that" as a load-bearing caveat
  marker in text *and* captions (1.36/kw in arXiv:2106.16211): "Note that although a satellite
  can appear in two consecutive observation IDs, it appears in the above plot as a single
  datum." (arXiv:2006.04327).
- **A9. Numbered conclusions restating each number ± error.** "1. In the six years covered by
  our maps, two jet components ... emerged from the core in succession. ... Their angular
  velocities were 0.17±0.05 and 0.11±0.02 mas/yr, respectively." (1989A+A...216...31T);
  hand-set `$\bullet$` lists with each bullet re-hedged (arXiv:2106.16211).
- **A10. Ending on the observation that would settle it.** "Higher-resolution observations are
  desired to distinguish..." (arXiv:2004.09369); "Further VLBI observations are needed to
  confirm its existence." (1989A+A...216...31T).
- **A11. Rhetorical questions only as genuine open questions**, answered "unknown", never by
  the paper's own result: "Where does the matter come from? How can the matter lose its
  angular momentum to accrete?" (arXiv:2004.09369); "What drives these jets and winds at
  speeds of a few hundreds of km/s is not known." (arXiv:1208.3351).
- **A12. Community scaffolding.** "(hereafter X)" contractions; "Paper I/II" series
  continuity; "in prep"/"submitted" forward references; "(Brian Grefenstette, personal
  communication, May 7, 2014)" with a name and date (arXiv:1703.00021); `\facility{}`,
  `\software{}` liturgy; the referee thanked ("We thank an anonymous referee for a careful
  revision of our paper that improved its clarity.").
- **A13. Era flavor (optional, use sparingly).** The spaced `--` as the *only* appositive
  break; "$\sim$" glued to order-of-magnitude numbers and "of order 1,300 years" phrasing
  (arXiv:2008.03725); quoted-term introductions "the so-called ``dirty'' image"; self-contained
  captions that decode every symbol in words and may carry a caveat ("The blue, black and red
  dashed lines are drawn `by-eye' for visualization (no fit involved)", arXiv:2106.16211).

---

## 4. Voice rules

Fingerprint comparison (corpus "all", N=1137, vs selfscan aggregate, N=45):

| metric | corpus p50 | our p50 | direction |
|---|---|---|---|
| `we_per_kw` | 7.53 (2016–21: 8.27) | 5.30 | corpus uses **more "we"** |
| `passive_per_sentence` | 0.333 | 0.205 | corpus uses **more passive** |
| `self_ref_per_kw` ("this paper/census...") | 0.55 | 2.24 | we self-refer 4× more |
| `first_singular_per_kw` | 0.098 (p90 1.44) | 0 (mean 0.45; dr20radio 2.63) | comparable, ours lumpy |
| `mean_sentence_words` | 26.2 | 31.3 | our sentences longer |

- The conversion direction is counter-intuitive: **add** "we" and **add** passive
  simultaneously, by deleting the third channel our papers overuse — inanimate self-reference
  ("This census measures...", "This paper performs...") and meta-commentary. Decisions become
  "we" sentences; procedure becomes passive.
- The split is functional, in every era: "instrument/observation/reduction prose is
  impersonal-passive; analysis and judgement are 'we'" (1960s counts: 17/18 passive clauses in
  a Methods page, 14 "we" verbs in 20 analysis sentences, 1968ApJ...154....3D).
- Single authors: the modern norm is authorial "we" even alone (1995MNRAS.275..217M: 9 "we", 0
  "I" in a single-author Discussion; 1970s: "'I' occurs ONLY in acknowledgments"). "I" in the
  body is attested but era-bound (1989A+A...209L...1W "I therefore propose"; 1955ApJ...121..367T
  "I have used") — dr20radio's body-"I" at 2.63/kw sits above the corpus p90 (1.44); keeping it
  is a defensible signature, but the corpus-median choice is "we"/passive with "I" reserved for
  the acknowledgments.
- "The paper/This paper" as agent is fine **once** ("This paper describes an automated masking
  algorithm...", arXiv:1912.04970) — as a roadmap or genre marker, not as a recurring narrator.
- Never "we" for the community; third parties are "the authors" or named.

---

## 5. What never changes in a conversion

- **Macros, verbatim.** Every `\drFoo`, `\lvFoo`, `\svbFoo`... invocation is preserved
  character-for-character, including trailing `\ ` spacing and `\%` context. Macros are the
  papers' evidence channel (`report.preserve_live_macros`, the namespacing rules in CLAUDE.md);
  a style pass that touches one can silently blank or swap a published number.
- **`\cite`/`\citep`/`\citet` keys and their grammatical role.** `\citet` is the subject of its
  sentence and is not deletable (CLAUDE.md lesson). Keys never edited.
- **Numeric literals, units, signs.** Every hard-coded number ($-40^\circ$, 888\,MHz, 5\arcsec,
  $p$-values) survives exactly. "Write the number after you run it" cuts both ways: a style
  pass writes no numbers at all.
- **`\software{}` / `\facility` blocks, keywords, `\input{generated/macros}`,** the
  `janskyresearch` + `jansky` citations, ORCID/affiliation front matter.
- **Section count and order may change; content claims may not.** Dissolving a Limitations
  section into Results/Discussion (S5) moves sentences — every caveat, magnitude, and negative
  must land somewhere. Deleting a hedge that gates a claim is a science edit, not a style edit.
- **The repo's honesty conventions are compatible with the target register — restate, don't
  remove.** The corpus states limits plainly, in-line, with magnitudes: "We failed to detect
  the high velocity flow in the 2-1 line of ¹³CO down to a limit of 0.06 K (3 σ)"
  (1989A+A...222L...1C); "our best upper limits ... may have been incorrectly reported as low
  as (50 mK)$^2$" (arXiv:1301.5906) — a corpus paper retracting its own number in the abstract.
  Validations, limits, and negatives stay; only their voice changes (meta-commentary → cause +
  magnitude at point of use).
- After any conversion: rebuild and run `uv run python scripts/triage_papers.py`; diff the
  PDF-extracted numbers against the pre-conversion build if in doubt.

---

## 7. Second-tier tells (found by the pilot's blind A/B, 2026-08-18)

The fashienv pilot passed the mechanical lint (all fingerprint metrics under corpus p90) and
still carried tells a blind judge caught. These are register-level, invisible to rate
metrics — check for them by reading, after the linter is clean:

- **Aphoristic X-not-Y antitheses.** "the limit is bias, not variance, and resampling cannot
  see bias". Corpus form: state which, cite why, no epigram.
- **"Headline result" and other blog/PR register.** Corpus papers rank results with "our main
  result" or simply order them; "headline", "takeaway", "the punchline" never appear.
- **Colon-pivot mini-sentences.** "The follow-on step is mechanical:", "We also flag a real
  sensitivity:". Corpus form: "We note that ...", or just say the thing.
- **Anthropomorphized methods as sentence subjects.** "The estimator has a known weakness
  that this measurement runs into"; "what this measurement needs and does not yet have".
  Corpus form: "The estimator is biased when ...", "a footprint mask is not available".
- **Triads with a "not merely" flourish.** "from a different telescope, sky, and pipeline,
  not merely a consistency check". Two items or an unadorned list.
- **Frame-setting opener sentences.** "Two facts set the scope of this analysis." Open with
  the first fact instead.
- **Dev-speak in prose.** "a one-line change", "the leak", "pipeline" as agent. Say what the
  operation is in astronomy terms.
- **Litotes as understatement.** "so the leak is not small" — give the magnitude.

The strip rules S1–S11 remove what the metrics see; this list is what a referee still smells.
A conversion is done when a blind reader picks it as pre-2022 with at most a couple of these
remaining.

---

## 8. RNAAS notes (genre addendum, 2026-08-18)

Empirical basis: `results/stylecorpus_rnaas.json` (391 pre-LLM notes, 2017–2021, the
arXiv-deposited subset of all 1,035 — baseline conditional on deposit) and
`data/style_corpus/stylenotes/rnaas.md` (13 notes read in full; sample is 2017–2018).
Lint with `prose_lint.py <paper> --file rnaas.tex --genre rnaas`.

The genre is not a compressed paper; it is a different shape:

- **Median 830 words (p90 1,480). The abstract convention BROKE mid-corpus:** the
  journal made abstracts *required* on 2020-05-01 (verified against the RNAAS
  instructions, 2026-08-18; current word limit 1,500). Measured in our sample: 7% of
  pre-2020-05 notes carry an abstract vs 88% after, median ~103 words when present
  (p90 147). For any new or converted note: **abstract required, keep it ≲100–130
  words**. Do not trust the aggregate `abstract_words` p50 = 0 in
  `results/stylecorpus_rnaas.json` for this metric — it averages across the policy
  break; all other metrics are unaffected by it.
- **The title carries the verdict** ("TRAPPIST-1e HAS A LARGE IRON CORE"); the opener
  states one plain fact and the result arrives within the first paragraph. No
  pre-summary opening, no closing recap — notes stop on content.
- **Structure**: often a single blank `\section{}` or 2 short sections (p50 = 2,
  p90 = 5), 3–5 paragraphs, exactly one exhibit that carries the note, then
  acknowledgments.
- **Voice runs hotter than papers**: "we" p50 8.77/kw (papers 7.53), first-person
  singular p90 2.73/kw — signed opinion and even exclamation marks ship. Hedges stay
  local and graded; criticism is named with the arithmetic shown.
- **Same hard zeros as papers**: em-dash p90 1.52/kw, `\emph` p90 2.46/kw, hedge
  vocabulary p90 1.55/kw. No roadmaps, no signposting, no intensifiers.
- Recurring note species, useful when drafting: correction, self-update, utility
  list, proof-of-concept, small measurement.
