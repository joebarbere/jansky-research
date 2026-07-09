/* ============================================================================
 * OVRO-LWA solar dynamic-spectrum (type II census, plan 50) — browser download helper
 * ============================================================================
 *
 * ⚠️  YOU PROBABLY DO NOT NEED THIS.  The daily dynamic-spectrum FITS are on
 *     AWS Open Data and are directly downloadable with NO login and NO CAPTCHA:
 *
 *        https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits/<YYYY>/<YYYYMMDD>.fits
 *
 *     The Cloudflare Turnstile on https://www.ovsa.njit.edu/lwadata-query only
 *     gates the interactive query UI, not the data.  So the simplest path is:
 *
 *        uv run python scripts/typeii_real.py --download 2024-05-14 2024-05-15
 *        # or:  curl -O https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits/2024/20240514.fits
 *
 * This script is a FALLBACK for when you want to browse/select days in the
 * portal UI by hand (e.g. only days flagged as having spectra) and let the
 * browser drive the downloads from your already-authenticated + CAPTCHA-passed
 * session.  It builds the same public S3 URLs and downloads them one by one.
 *
 * ----------------------------------------------------------------------------
 * HOW TO USE
 * ----------------------------------------------------------------------------
 * 1. Open  https://www.ovsa.njit.edu/lwadata-query  in your browser and, if the
 *    page shows a Cloudflare Turnstile checkbox, complete it. (Downloads below
 *    hit the public S3 bucket directly, so this step is only to satisfy you the
 *    session is live — it is not strictly required for the S3 GETs.)
 * 2. Open the browser devtools console (F12 → Console).
 * 3. Paste this whole file and press Enter to define the helpers.
 * 4. Pick your dates and call, e.g.:
 *
 *        // explicit list of days:
 *        await downloadOvroDays(["2024-05-14", "2024-05-15", "2024-05-16"]);
 *
 *        // OR every day in an inclusive range:
 *        await downloadOvroRange("2024-05-01", "2024-05-31");
 *
 *        // OR: only the days that actually HAVE spectra, scraped from the
 *        //     portal's own results after you run a query in the UI:
 *        await downloadOvroDays(scrapeDatesFromPortal());
 *
 * Files land in your browser's Downloads folder as <YYYYMMDD>.fits (each ~1.7 GB,
 * so pick a handful of days, not a whole year). A 3 s gap between downloads is
 * used to stay polite to the bucket; raise GAP_MS if you see throttling.
 * ----------------------------------------------------------------------------
 */

const OVRO_S3 = "https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits";
const GAP_MS = 3000; // pause between file downloads

/** "2024-05-14" (or "20240514") -> the public S3 URL for that day's dynamic spectrum. */
function ovroDspecUrl(date) {
  const ymd = String(date).replace(/-/g, "");
  return `${OVRO_S3}/${ymd.slice(0, 4)}/${ymd}.fits`;
}

/** Trigger a browser download of one URL as <name>, resolving when the click is dispatched. */
function triggerDownload(url, name) {
  const a = document.createElement("a");
  a.href = url;
  a.download = name;           // hint the filename; cross-origin S3 may ignore it (still saves)
  a.style.display = "none";
  document.body.appendChild(a);
  a.click();
  a.remove();
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/** Download an explicit list of days (["YYYY-MM-DD", ...]) sequentially. */
async function downloadOvroDays(dates) {
  if (!dates || !dates.length) {
    console.warn("[ovro] no dates given");
    return;
  }
  console.log(`[ovro] downloading ${dates.length} day(s), ~1.7 GB each ...`);
  for (let i = 0; i < dates.length; i++) {
    const ymd = String(dates[i]).replace(/-/g, "");
    const url = ovroDspecUrl(dates[i]);
    console.log(`[ovro] (${i + 1}/${dates.length}) ${url}`);
    triggerDownload(url, `${ymd}.fits`);
    if (i < dates.length - 1) await sleep(GAP_MS);
  }
  console.log("[ovro] all downloads dispatched — watch your browser's download manager.");
}

/** Download every day in an inclusive [start, end] range (both "YYYY-MM-DD"). */
async function downloadOvroRange(start, end) {
  const dates = [];
  for (let d = new Date(start + "T00:00:00Z"); d <= new Date(end + "T00:00:00Z");
       d.setUTCDate(d.getUTCDate() + 1)) {
    dates.push(d.toISOString().slice(0, 10));
  }
  await downloadOvroDays(dates);
}

/**
 * Best-effort scrape of dates from the portal's rendered results, so you only
 * grab days that have data. Run a query in the UI first, then call this. It
 * looks for YYYY-MM-DD / YYYYMMDD tokens in the results DOM; if the portal
 * markup changes, fall back to an explicit downloadOvroDays([...]) list.
 */
function scrapeDatesFromPortal() {
  const text = document.body.innerText || "";
  const iso = [...text.matchAll(/\b(20\d{2})-(\d{2})-(\d{2})\b/g)].map((m) => m[0]);
  const compact = [...text.matchAll(/\b(20\d{2})(\d{2})(\d{2})\b/g)].map(
    (m) => `${m[1]}-${m[2]}-${m[3]}`
  );
  const uniq = [...new Set([...iso, ...compact])].sort();
  console.log(`[ovro] scraped ${uniq.length} candidate date(s) from the page`, uniq);
  return uniq;
}

console.log(
  "[ovro] helpers ready: downloadOvroDays([...]), downloadOvroRange(a,b), " +
    "scrapeDatesFromPortal(). NOTE: the S3 bucket is public — you can also just " +
    "`curl -O` the URLs from ovroDspecUrl(date)."
);
