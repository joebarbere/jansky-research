#!/usr/bin/env python3
"""Download RACS Stokes-V (or I) FITS cutouts from the CASDA web Cutout Service via Playwright.

Best-effort browser automation of the CSIRO Data Access Portal CASDA Cutout Service
(https://data.csiro.au/domain/casdaCutoutService), used when the CASDA VO services (TAP/SIA2) are
down or awkward. It logs in with an OPAL account, submits a position+survey+radius cutout, waits for
the result, and downloads the matching FITS file(s).

Design stance: **aggressively cautious**. Every step has a bounded timeout and an explicit
success/failure check; a failure screenshots + dumps the page and exits with a distinct, documented
code rather than hanging or silently "succeeding". It only reads and downloads — it never mutates
anything on the site.

Credentials: never passed on the command line. Username via --username (or $CASDA_USERNAME); the
password is read from a 0600 file (--pw-file, default ~/.casda_pw) and is never logged or screenshotted
(the login screenshot is taken *before* the password is typed).

Exit codes:
  0  success: at least one FITS downloaded and validated
  2  bad usage / missing credentials file
  3  could not launch the browser (missing system libs etc.)
  4  login failed (OPAL auth not confirmed)
  5  no results found — typically the CASDA discovery backend is down (the web service returns 0
     results for *every* position, including its own example, when the backend is out)
  6  results found but no FITS matched the requested Stokes/pol filter
  7  a download started but the file failed validation (not a FITS, or truncated)
  8  unexpected error (page structure changed, navigation timeout, ...)

This is the automation half of the `casda-cutout-fetch` skill; see SKILL.md.
"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

PORTAL = "https://data.csiro.au/domain/casda"
CUTOUT = "https://data.csiro.au/domain/casdaCutoutService"
RESULTS = "https://data.csiro.au/domain/casdaCutoutService/results"
FITS_MAGIC = b"SIMPLE  ="


def log(msg: str) -> None:
    print(f"[casda-cutout] {msg}", flush=True)


def _dump(page, out_dir: Path, tag: str) -> None:
    """Best-effort diagnostic capture (screenshot + HTML) on failure; never raises."""
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(out_dir / f"FAIL_{tag}.png"), full_page=True)
        (out_dir / f"FAIL_{tag}.html").write_text(page.content())
        log(f"diagnostics written: {out_dir}/FAIL_{tag}.png|.html")
    except Exception as e:  # pragma: no cover - diagnostics must never mask the real error
        log(f"(could not write diagnostics: {e})")


def _dismiss_ack(page) -> None:
    """Dismiss the 'Acknowledgement of Country' modal that blocks clicks on load."""
    try:
        btn = page.query_selector('modal-container button:has-text("Continue")')
        if btn:
            btn.click(timeout=5000)
            page.wait_for_timeout(500)
    except Exception:
        pass


def login(page, username: str, password: str, out_dir: Path, step_timeout: int) -> None:
    """OPAL login through the DAP 'Sign In' modal. Raises RuntimeError if it can't confirm success."""
    page.goto(PORTAL, timeout=step_timeout, wait_until="domcontentloaded")
    page.wait_for_timeout(3000)
    _dismiss_ack(page)
    page.click("text=Sign In", timeout=step_timeout)
    page.wait_for_timeout(1500)
    page.click("text=Sign in with OPAL Account", timeout=step_timeout)
    page.wait_for_timeout(2500)
    # screenshot BEFORE typing the password (so the secret never lands in an artifact)
    try:
        page.screenshot(path=str(out_dir / "login_form.png"))
    except Exception:
        pass
    page.fill('input[placeholder="Username"]', username, timeout=step_timeout)
    page.fill('input[placeholder="Password"]', password, timeout=step_timeout)
    page.click('button:has-text("Sign in")', timeout=step_timeout)
    page.wait_for_timeout(6000)
    _dismiss_ack(page)
    signed_in = page.eval_on_selector_all(
        "a,button,span",
        "els => els.some(e => /sign out|log out|my account/i.test(e.innerText||''))",
    )
    err = page.eval_on_selector_all(
        "*",
        "els => els.some(e => /invalid|incorrect|denied|failed to/i.test((e.innerText||'').slice(0,80)))",
    )
    if not signed_in or err:
        _dump(page, out_dir, "login")
        raise RuntimeError("login not confirmed (no 'Sign out' marker or an error was shown)")
    log("OPAL login confirmed")


def request_results(
    page, ra: float, dec: float, surveys: list[str], size_arcmin: float, step_timeout: int
) -> None:
    """Navigate straight to the results URL (the cutout form just builds this URL)."""
    qs = "&".join(f"surveys={s}" for s in surveys)
    url = f"{RESULTS}?{qs}&ra={ra}&dec={dec}&size={size_arcmin}"
    log(f"requesting cutout: {url}")
    page.goto(url, timeout=step_timeout, wait_until="networkidle")
    page.wait_for_timeout(8000)
    _dismiss_ack(page)
    page.wait_for_timeout(2000)


def parse_results(page) -> tuple[int, str]:
    """Return (per-survey result count total, body text). Counts the 'RACS-* DR1 N' tallies."""
    body = page.inner_text("body")
    import re

    counts = [int(n) for n in re.findall(r"DR1\s+(\d+)", body)]
    return (sum(counts) if counts else 0), body


def download_matches(page, pol_regex: str, out_dir: Path, step_timeout: int) -> list[Path]:
    """Find result rows whose title matches the Stokes/pol filter and download their FITS.

    The results list renders one row per available image; each has a download affordance (a link to
    the CASDA data-access/SODA service or a 'Download' button). We match rows by visible text against
    ``pol_regex`` (default Stokes V), trigger the browser download, and validate the FITS magic bytes.
    """
    import re

    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    rows = page.query_selector_all("tr, li.result, div.result-row, [class*=result]")
    candidates = []
    for r in rows:
        try:
            txt = (r.inner_text() or "").strip()
        except Exception:
            continue
        if txt and re.search(pol_regex, txt, re.I):
            candidates.append((txt[:80], r))
    log(f"{len(candidates)} result row(s) match /{pol_regex}/i")
    for title, row in candidates:
        link = row.query_selector(
            "a[href*=data_access], a[href*=soda], a[href*='.fits'], "
            "a:has-text('Download'), button:has-text('Download')"
        )
        if link is None:
            continue
        try:
            with page.expect_download(timeout=step_timeout) as dl_info:
                link.click()
            dl = dl_info.value
            dest = out_dir / (dl.suggested_filename or "cutout.fits")
            dl.save_as(str(dest))
            head = dest.read_bytes()[:80]
            if FITS_MAGIC not in head or dest.stat().st_size < 2880:
                log(f"  REJECT {dest.name}: not a valid/complete FITS (head={head[:16]!r})")
                continue
            log(f"  saved {dest.name} ({dest.stat().st_size} bytes) for row '{title}'")
            saved.append(dest)
        except Exception as e:
            log(f"  download failed for row '{title}': {str(e).splitlines()[0][:80]}")
    return saved


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ra", type=float, required=True, help="J2000 RA in decimal degrees")
    ap.add_argument("--dec", type=float, required=True, help="J2000 Dec in decimal degrees")
    ap.add_argument(
        "--surveys",
        default="RACS-Mid,RACS-High",
        help="comma list of RACS surveys with Stokes V (default RACS-Mid,RACS-High; "
        "RACS-Low DR1 has no Stokes V)",
    )
    ap.add_argument("--size-arcmin", type=float, default=2.0, help="cutout radius in arcmin")
    ap.add_argument(
        "--pol-regex",
        default=r"stokes.?v|_v[._]| v ",
        help="regex matched against result-row text to pick the Stokes-V image",
    )
    ap.add_argument(
        "--out", default="casda_cutouts", help="output directory for FITS + diagnostics"
    )
    ap.add_argument(
        "--username",
        default=os.environ.get("CASDA_USERNAME", ""),
        help="OPAL username (or $CASDA_USERNAME)",
    )
    ap.add_argument(
        "--pw-file", default="~/.casda_pw", help="file holding the OPAL password (0600)"
    )
    ap.add_argument("--step-timeout", type=int, default=60000, help="per-step timeout (ms)")
    ap.add_argument("--budget-seconds", type=int, default=600, help="global wall-clock budget")
    ap.add_argument("--headed", action="store_true", help="show the browser (debugging)")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    if not args.username:
        log("ERROR: no --username and $CASDA_USERNAME unset")
        return 2
    pw_path = Path(args.pw_file).expanduser()
    if not pw_path.is_file() or pw_path.stat().st_size == 0:
        log(f"ERROR: password file {pw_path} missing or empty")
        return 2
    password = pw_path.read_text().strip()
    if not password:
        log(f"ERROR: password file {pw_path} is blank")
        return 2

    surveys = [s.strip() for s in args.surveys.split(",") if s.strip()]
    deadline = time.monotonic() + args.budget_seconds

    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        log(
            "ERROR: playwright not installed (uv pip install playwright && playwright install chromium)"
        )
        return 3

    with sync_playwright() as p:
        try:
            browser = p.chromium.launch(headless=not args.headed)
        except Exception as e:
            log(f"ERROR: could not launch chromium: {str(e).splitlines()[0]}")
            log("hint: run 'uv run playwright install chromium' (and install-deps if on a server)")
            return 3
        ctx = browser.new_context(accept_downloads=True)
        ctx.set_default_timeout(args.step_timeout)
        page = ctx.new_page()
        try:
            login(page, args.username, password, out_dir, args.step_timeout)
            if time.monotonic() > deadline:
                log("ERROR: budget exhausted after login")
                return 8
            request_results(page, args.ra, args.dec, surveys, args.size_arcmin, args.step_timeout)
            total, _ = parse_results(page)
            if total == 0:
                _dump(page, out_dir, "no_results")
                log(
                    "NO RESULTS for any requested survey. When CASDA's discovery backend is down the "
                    "cutout service returns 0 for EVERY position (including its own example) — verify "
                    "by re-running at a known-covered field; if that is also 0, it is the outage."
                )
                return 5
            log(f"{total} candidate image(s) reported across {surveys}")
            saved = download_matches(page, args.pol_regex, out_dir, args.step_timeout)
            if not saved:
                _dump(page, out_dir, "no_match")
                log(
                    "results existed but none matched the Stokes/pol filter (or had no download link)"
                )
                return 6
            log(f"DONE: {len(saved)} FITS downloaded to {out_dir}/")
            return 0
        except RuntimeError as e:
            log(f"LOGIN ERROR: {e}")
            return 4
        except Exception as e:
            log(f"UNEXPECTED ERROR: {str(e).splitlines()[0]}")
            _dump(page, out_dir, "unexpected")
            return 8
        finally:
            ctx.close()
            browser.close()


if __name__ == "__main__":
    raise SystemExit(main())
