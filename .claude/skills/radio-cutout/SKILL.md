---
name: radio-cutout
description: Fetch a radio-survey FITS image cutout at a sky position, routing to the archive that serves it (CADC SODA for VLASS — working/default; Data Central for GLEAM/GLEAM-X and RACS Stokes I; CASDA for RACS Stokes V). Use when a slice needs an image cutout at an RA/Dec.
---

# Radio image-cutout fetcher

One place that knows **which archive serves which survey** and how to pull a cutout — distilled from
the research slices, which kept re-deriving this. The **CADC SODA** path is the validated workhorse;
other surveys are routed to their archive with honest notes.

## Archive map (verified during the slices)

| Survey | Archive | Status |
|---|---|---|
| **VLASS** (2–4 GHz Quick-Look) | **CADC SODA** | **working — this script** (verified: 3 epoch cutouts at M87) |
| any other **CADC** collection | CADC SODA | working — pass `--collection <NAME>` |
| **GLEAM / GLEAM-X** (72–231 MHz) | **Data Central** POST cutout API | Stokes I; `fits:true` + band PK (see below) |
| **RACS-low/mid/high** Stokes I | Data Central / CASDA | Stokes I total intensity |
| **RACS Stokes V** | **CASDA only** | use the `casda-cutout-fetch` skill (CASDA often down) |

## Usage

```
# VLASS (CADC SODA) — the working default
uv run python .claude/skills/radio-cutout/fetch_cutout.py \
    --survey vlass --ra 187.7059 --dec 12.3911 --size-arcmin 1 --out cutouts

# any other CADC collection
uv run python .claude/skills/radio-cutout/fetch_cutout.py \
    --collection VLASS --ra <ra> --dec <dec> --size-arcmin 1.5 --out cutouts
```

Each saved file is validated to parse as a FITS before keeping. Exit codes: `0` ok · `2` bad usage ·
`5` no image at that position · `6` survey routed elsewhere (GLEAM/RACS — see notes the script prints)
· `7` download/validation failure · `8` unexpected error.

## The non-CADC routes (notes, not yet wired into the script)

- **GLEAM / GLEAM-X / RACS-I — Data Central.** POST to
  `https://datacentral.org.au/api/services/cutout/` with JSON: `ra`/`dec` as strings **with a unit**
  (`"30.0d"`), `radius` in **arcsec** (≥1), `data_releases` and `greyscale_bands` as **integer PKs**
  (enumerate them from the endpoint's `OPTIONS` → `actions.POST.greyscale_bands.grouped_choices`),
  and `fits: true`. RACS on Data Central is **Stokes I only** (no V). For *catalogue* GLEAM-X/RACS
  fluxes (not images) prefer VizieR — `jansky_research.southern.fetch_gleamx` / `fetch_racs_bands`.
- **RACS Stokes V — CASDA.** The only public source; use the `casda-cutout-fetch` skill (CADC does not
  hold RACS). As of 2026-06 CASDA's query backend was down (TAP/SIA2 + web cutout all 0/errors).

## Reuse

The CADC path is the same one `jansky_research.vlass.fetch_vlass_cutouts` uses for forced photometry;
this skill generalises it to arbitrary CADC collections + a name-filter, and adds FITS validation +
clear exit codes. For Stokes-V forced photometry on a cutout, feed the array + WCS to
`jansky_research.stokesv.measure_circular_pol`.
