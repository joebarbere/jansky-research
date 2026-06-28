---
name: radio-source-lookup
description: Gather what's known about a radio source or sky position across catalogues and archives — SIMBAD/NED for identity, NVSS/FIRST/VLASS/RACS for radio flux and spectral index, and the VO via astroquery/pyvo. Use when the user names a source (e.g. Cygnus A, 3C 273, an SRSC radio star) or a position and wants its radio properties.
---

# Look up a radio source

Assemble a concise profile of a named source or sky position from open catalogues. (Ported from the
`jansky` course's skill; adapted to reuse this repo's own helpers.)

## Procedure

1. **Resolve identity & position** with `astroquery.simbad` (and NED for extragalactic objects):
   coordinates, object type, redshift. Use `astropy.coordinates.SkyCoord` for the position.
2. **Radio measurements** via `astroquery.vizier` cone searches — reuse
   `jansky_research.spectra.fetch_survey` (TGSS/NVSS) and the patterns in `stokesv.fetch_racs_i` /
   `fetch_radio_star_measurements` (RACS, SRSC) rather than re-querying by hand:
   - **NVSS** (`VIII/65`) — 1.4 GHz flux over most of the sky.
   - **FIRST** (`VIII/92`) — 1.4 GHz high-resolution.
   - **RACS** (`J/other/PASA/38.58`, RACS-low Stokes I) for the southern sky; **TGSS** 150 MHz.
   - Then estimate a **spectral index** from two or more frequencies with
     `jansky_research.spectra.spectral_index` (our own tested helper).
3. **Which instruments can observe it** — if `../jansky` is checked out, cross-reference its
   `docs/telescopes.md` (declination & frequency coverage) and `docs/resources.md` for archives.

**Always wrap network calls in try/except** and fall back gracefully; if offline, say so and report
what would be queried.

## Report

A short profile: identity/position/type, the radio flux densities found (with frequency and survey),
an estimated spectral index if ≥2 points exist, and pointers to archives/instruments. Cite the
surveys (NVSS = Condon et al. 1998, FIRST = Becker et al. 1995, RACS = McConnell/Hale et al. 2020/21,
TGSS = Intema et al. 2017).

For heavy or multi-archive data hunts (finding actual observation datasets to download), hand off to
the **archive-scout** agent.
