# `dmzcal` findings (plan 99)

## Step 0 — provenance (2026-10-10, before any posterior)

The decisions, counts and caveats are recorded in `plans/99-dmzcal-coverage.md`, in the section
"Step 0 outcome and pre-posterior refinements". The evidence is
`results/dmzcal_provenance.json`, produced by `scripts/dmzcal_provenance.py` with FRBs/FRB
pinned at 996fcda and zdm at e0f985b.

Headline: the certification sample is 36 bursts (CHIME 19, DSA 11, ASKAP 4, MeerKAT 2). Four
facts found along the way are worth keeping:

- **zdm's out-of-the-box redshift state cannot be audited for provenance.** `load_state()`
  defaults to `HoffmannHalo25`, which is unpublished. Its documented p(z|DM) script imports
  `imp`, so it does not run on Python ≥ 3.12. Anyone who calls zdm today without choosing a
  state gets parameters whose fit sample they cannot know.
- **The published fit sample is wider than "CRAFT".** Hoffmann et al. 2025 also used the DMs of
  25 DSA-110 commissioning bursts (redshifts for only 3). The DSA bursts that remain available
  for certification come from the later Sharma+2024 sample.
- **The FRBs/FRB host table under-flags spectroscopic redshifts.** The `Spectrum` column is blank
  for all 19 CHIME/KKO gold hosts of Leung et al. 2025, which do have spectra. A naive
  "flagged-spectrum only" rule would have left 20 certification bursts and an UNDERPOWERED
  verdict caused by a metadata gap.
- **One redshift is transcribed wrongly.** FRB 20231201A is 0.119 in FRBs/FRB and 0.1119 in the
  primary table. This could be reported upstream.
