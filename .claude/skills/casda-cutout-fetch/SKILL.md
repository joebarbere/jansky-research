---
name: casda-cutout-fetch
description: Download RACS Stokes-V (or I) FITS cutouts from the CASDA web Cutout Service by driving the browser with Playwright and an OPAL login. Use when the CASDA VO services (TAP/SIA2) are down or awkward and you need a RACS cutout at a position. Best-effort browser automation; aggressively guarded.
---

# CASDA Cutout Service fetcher (Playwright)

**Reality first (state this to the user):** CASDA's machine-friendly VO services (TAP catalogue,
ObsCore, SIA2) have been unreliable — returning `relation "..." does not exist` and SIA2 500 NPEs —
which blocks the normal `astroquery.casda` `query_region → cutout` flow. This skill routes around that
by driving the **web** CASDA Cutout Service (https://data.csiro.au/domain/casdaCutoutService) with a
real browser: OPAL login → position+survey+radius cutout → download the FITS. It is **best-effort**:
when CASDA's discovery backend is itself down, the *web* service also returns zero results for every
position (verified against the service's own example target), and no automation can download what the
backend won't surface — the skill detects this and exits cleanly with code 5.

## Prefer the API first (when to use this skill)

**Try `astroquery.casda` (`query_region` → `cutout`) before this skill.** When CASDA's backend is
healthy, the programmatic API is the proper, robust tool and needs no browser. This skill exists as a
**fallback** for the case we actually hit: the VO services (TAP/SIA2) misbehaving while the web portal
still authenticates. Honest caveat: its **download step has never been exercised end-to-end** — the
CASDA outage has so far prevented any cutout from being produced — so login + navigation + parsing are
verified but the final download/validation is not. It is kept because (a) it encodes the full working
CASDA web-login + cutout-form flow, (b) it's a genuinely independent access path, and (c) it costs
nothing when unused (not imported by the package; no CI/test impact). Remove it if a future CASDA where
the API is reliable makes the redundancy pointless.

## Prerequisites

- **Playwright + Chromium** (vendored deps; install once):
  `uv pip install playwright && uv run playwright install chromium`
  (`requirements.txt` in this folder pins playwright; on a bare server you may also need
  `uv run playwright install-deps`, which wants root.)
- **OPAL account** (free, self-service at https://opal.atnf.csiro.au/register). Username = email.
- **Password file** at `~/.casda_pw`, mode `0600` — the script reads it and never logs or
  screenshots it (the login screenshot is taken *before* the password is typed). Never pass the
  password on the command line.

## Usage

```
uv run python .claude/skills/casda-cutout-fetch/fetch_cutout.py \
    --ra 232.8191 --dec 34.6908 \
    --surveys RACS-Mid,RACS-High \
    --size-arcmin 2 \
    --username joe.barbere@gmail.com \
    --out casda_cutouts
```

- `--surveys` defaults to `RACS-Mid,RACS-High` — the bands that carry **Stokes V** (RACS-Low DR1 has
  none). `--pol-regex` (default matches Stokes V) picks which result rows to download.
- RA/Dec are **decimal degrees**. `--size-arcmin` is the cutout radius.
- `--headed` shows the browser for debugging; `--budget-seconds` / `--step-timeout` bound the run.

## What it does (the flow it automates)

1. Load the DAP, dismiss the "Acknowledgement of Country" modal.
2. **Sign In → Sign in with OPAL Account**, fill Username/Password, submit; **confirm** a `Sign out`
   marker appears (else exit 4).
3. Navigate to the results URL the cutout form builds:
   `…/casdaCutoutService/results?surveys=RACS-Mid&surveys=RACS-High&ra=…&dec=…&size=…`.
4. Parse the per-survey `RACS-* DR1 N` tallies. If all zero → exit 5 (likely the backend outage).
5. For each result row matching `--pol-regex`, trigger the download, **validate the FITS magic bytes**
   (`SIMPLE  =`) and a minimum size, and save to `--out`.

## Exit codes

| code | meaning |
|---|---|
| 0 | success — at least one FITS downloaded and validated |
| 2 | bad usage / missing or empty password file |
| 3 | browser could not launch (missing system libs) |
| 4 | OPAL login not confirmed |
| 5 | **no results** — usually the CASDA discovery backend is down (web service returns 0 for every position) |
| 6 | results found but none matched the Stokes/pol filter / had no download link |
| 7 | a download failed FITS validation (truncated / not a FITS) |
| 8 | unexpected error (page structure changed, timeout) — see the `FAIL_*.png/.html` diagnostics |

## Guards (aggressively cautious by design)

- Verifies the password file exists and is non-empty *before* launching anything.
- Per-step timeouts **and** a global wall-clock budget — it cannot hang.
- Explicit login-success check; explicit "no results" check; FITS-magic + size validation on every
  download.
- On any failure it writes `FAIL_<stage>.png` + `FAIL_<stage>.html` to `--out` for inspection.
- Read-only on the site: it only navigates and downloads, never edits or deletes.
- Distinct exit codes so a caller can tell "outage" (5) from "broken automation" (8).

## Status (recorded 2026-06)

Verified live: **login, navigation, form/results parsing, and the guards all work**; the script
reaches the results page and parses it correctly. **It could not complete a download** only because
CASDA's discovery backend was returning **0 results for every position — including the cutout
service's own example target (PSR B1919+21)** — i.e. the same outage that takes down TAP/SIA2,
manifesting in the web service. The skill is kept ready: when CASDA recovers, re-run it and it should
download. See `survey/stokesv-findings.md` for the full data-access investigation.
