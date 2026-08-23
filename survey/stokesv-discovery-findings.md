# Findings — RACS Stokes-V discovery: two-epoch forced photometry of the nearest M dwarfs

`jansky_research.stokesv_discovery` + `scripts/stokesv_discovery_real.py` deliver plan 33: forced
I+V photometry of the nearest CNS5 M dwarfs at TWO RACS-mid epochs (MJD 59233 / 60769 — a 4.2-yr
same-band pair at 1367.5 MHz, found at GATE 0 after RACS-low1 V images turned out not to exist),
with Gaia PM propagation to each epoch, per-epoch leakage floors, and signed two-epoch variability.

## Recover-a-known (synthetic epoch pair)

Selection completeness/purity 1.0/1.0; flare-flagging completeness/purity 1.0/1.0; and the
headline rationale: **22.9% of injected emitters are invisible to any single epoch**
(`single_epoch_miss_frac`) — the quantified version of the stokesv paper's variability limit.

## Real sky (60 nearest CNS5 M dwarfs, 115-min CASDA run)

| quantity | value |
|---|---|
| targets measured (≥1 epoch) | 39 (CASDA staging failures account for the rest) |
| complete epoch pairs | 19 |
| ≥5σ V detections | **1 system: GJ 65 (BL+UV Ceti)** — V=9.26±0.15 (2021) → 7.11±0.15 mJy (2025) |
| inter-epoch ΔV significance | **10σ** at 1367.5 MHz over 4.2 yr |
| everything else | quiescent; median 5σ V limit **0.83 mJy** |

GJ 65 is a *recovery* (it is in the RACS-low2 Paper VIII blind V catalogue —
RACS-LOW2 J013906.5−175647, 4″ away — and in SRSC): the pipeline finds the prototype coherent
emitter and adds the mid-band two-epoch V change that blind single-epoch catalogues don't provide.
No candidate survived the novelty bar (checked against Paper VIII/SRSC/SIMBAD before "new" could
be used). PM propagation validated live: Barnard's star moved 39.8″ between the epochs.

## Honest caveats

- Two epochs bound variability, not timescale (flare vs secular).
- Staging availability limits the census: 35% of targets got no data, ~half of measured targets
  lack the second epoch (~52% of target-epoch slots failed; rows kept with notes; resumable).
- Both GJ 65 snapshots may be burst states (bursts last minutes-hours vs 15-min integrations);
  its quiescent V floor is unconstrained here.
- Leakage floor is per-field; a beam-position-dependent model would sharpen faint candidates.
- Reproduce: `uv run python scripts/stokesv_discovery_real.py` (needs `CASDA_USERNAME` +
  `~/.casda_pw`; ~2 h; resumable) then `uv run python -m jansky_research.stokesv_discovery --out .`.

## Referee round (2026-08-12) — three blockers, all confirmed by recomputation

**The photometry was not forced.** The Method said "forced peak flux at the propagated pixel";
`measure_circular_pol` took the brightest Stokes-I pixel within 12″ and read V there. On blank
sky that is a noise-maximum statistic, and the census showed it: **I > 0 for 54 of 54**
quiescent targets (p = 2⁻⁵⁴ for a genuine fixed-pixel measurement) at median I/σ = 2.2. So each
quiescent "limit" was measured up to about one synthesized beam from the star. A genuinely
forced mode (`search_arcsec <= 0`, reading the pixel containing the position) is now
implemented, the docstring — which claimed "forced ... at a locked (ra, dec)" and then
described a peak search in the next sentence — is corrected, and a test pins the distinction:
forced photometry on noise goes negative about half the time, the 12″ search never does.
Re-measuring the census in forced mode needs a CASDA re-fetch of all 120 cutouts (none are
cached) and is the outstanding item.

**The 10σ GJ 65 decline is image noise only.** I fell 15.99 → 11.87 mJy (−26%) and V 9.26 →
7.11 (−23%), leaving V/I constant to 3.4% — the signature of a flux-scale difference between
two independently calibrated observations, not of a change in emission state. With a per-epoch
scale term: **5.2σ at 3%, 3.5σ at 5%, 1.8σ at 10%**. The committed table records no tile or
beam position, so the systematic cannot be bounded from shipped evidence. Reported now as a
marginal decline, and the title's "a GJ 65 Variability Recovery" is now "a GJ 65 Recovery".

**One system, four rows.** CNS5 424 and CNS5 425 (BL and UV Ceti, 2″ apart, one beam) carry
**bit-identical** photometry at both epochs. `\svdRealDet` counts systems, but the prose said
"target rows", inflating one measurement into two in the place a skimming reader looks.

**Two arXiv-assembler bugs, fixed in the assembler rather than papered over.** The generated
abstract contained `9.267.11 mJy` — a flux that does not exist — because `\rightarrow` had no
symbol mapping and the generic `\[a-zA-Z]+ → ''` sweep fused the operands of
`$V=9.26\rightarrow7.11$`. And `\citealt`, which is natbib's *textual* form, was deleted
alongside the parenthetical `\citep`, leaving `(it appears in the blind V catalogue, )`. Both
now handled; arrows map to `->` and `\citealt`/`\citeauthor` resolve from refs.bib like
`\citet`. This affected every paper, not just this one.

Also: RACS-mid was cited to `mcconnell2020`/`hale2021`, which are both RACS-**low** — the
design paper and the low-band source catalogue. Duchesne et al. 2023 (PASA 40, e034,
Crossref-verified) is the RACS-mid data release, and is where the beam, astrometry and
flux-scale accuracy this paper needs are characterised. The abstract's claim that a leakage
floor "guards" the killer systematic is corrected: it is exercised in simulation only and
cannot be computed for this sample (the floor is a multiple of the *bright*-source median
|V|/I, and this census's median I/σ is 2.2). And the synthetic validation's 1.0/1.0 is
restated as what it is — injected emitters sit 3.5–14× above the floor, contaminants at about
a seventh of it, so it establishes the selection arithmetic, not sensitivity near the boundary.

Open: the forced re-measurement; per-target limits as a machine-readable table; recording
tile/beam position so the flux-scale term can be bounded; and a real-data figure.

## Correctness pass (2026-08-22) — the peak search is now stated, and a macro stopped lying

**`\svdRealDet` counted neither systems nor rows.** The paper read "This is one independent
measurement, not two, and \svdRealDet\ counts systems rather than rows", and the macro is 2.
`n_v_detections_5sig` is `len(set(det))` over catalogue *names*: the four >=5-sigma rows are
(CNS5 425, mid1/mid2) and (CNS5 424, mid1/mid2), the two entries of one binary 2 arcsec apart
in the same beam with bit-identical photometry. So the chain is **4 rows -> 2 CNS5 entries ->
1 system**, and the macro sits in the middle of it. The sentence now says exactly that, and
the estimator carries a comment recording what the number counts.

**The census is a 12-arcsec peak search and the Method said otherwise.** `stokesv.py` defaults
`search_arcsec = 12.0`, which "finds the brightest Stokes-I pixel within that radius"; only
`search_arcsec <= 0` is genuinely forced. The Method claimed measurement "at the propagated
position", and two occurrences of "forced" survived the 2026-08-12 round that removed it from
the title and abstract.

The signature is in the committed evidence and was re-verified here: **all 54 quiescent rows
have I > 0**, at a median I/sigma_I of 2.19, where a genuine fixed-pixel measurement on noise
would go negative about half the time (p = 2^-54 = 5.6e-17).

Fixed by stating it rather than by re-measuring: the Method now describes the 12-arcsec search,
"forced" is gone from the two places it was a method claim, and Honest limits (v) says plainly
that the quiescent limits are **upper limits on a beam-scale maximum, not on the star**, quotes
the 54/54 result, and records that the V value is read wherever that maximum fell. The forced
re-measurement (~2 h, CASDA credentials, a re-fetch of all 120 cutouts) remains the outstanding
work on this census.
