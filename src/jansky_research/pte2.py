"""Parkes Transient Events II: a per-source single-pulse heavy-tail census (plan 47, F10).

The PTE-II database (Yang+2025, ApJS; arXiv:2508.14403) holds 165,592 single pulses from 363 pulsars
reprocessed from Parkes Multibeam 1997--2001 archival data. The release paper already fits
log-normal/Gaussian *fluence* distributions to the 98 highest-count pulsars and reports a
population-level power-law, and HTRU-V (Burke-Spolaor+2012) is a prior 315-pulsar log-normal
single-pulse energy census --- so a generic "energy-distribution census" is NOT novel.

**The wedge this slice adds.** A uniform, per-source **giant-pulse test** across ALL 363 pulsars:
for each we ask whether its single-pulse S/N distribution has an extreme tail that its own
log-normal bulk cannot produce, via a Poisson excess test above a floor-robust threshold (a positive
giant-pulse detection, not merely a power-law fit; a Vuong power-law-vs-log-normal statistic and the
giant-tail index are reported as secondary diagnostics). We then rank the giant-pulse excess and
heavy-tailed fraction against spin-down luminosity (Edot) from the ATNF catalogue. HTRU-V tested
log-normal/Gaussian (never a heavy-tail model); PTE-II fit only 98 and its power-law is on the
*population* fluence-ratio distribution, not per-source tails --- so the per-source
heavy-tail-vs-Edot census is unclaimed.

**Honest scope.** (1) Per-pulse *fluence* is not stored in the DB; it is derived on the fly from raw
blobs by the paper's tool. We use per-pulse **S/N** (``snr_max``) as the energy proxy --- proportional
to fluence within a single source (fixed system temperature/gain), which is exactly the regime a
tail-shape test needs, but it forbids cross-source *absolute* energy comparison (only tail SHAPE and
the heavy/not-heavy classification are compared across sources). (2) The 1997--2001 sensitivity is
heterogeneous; the recover-a-known reports classification accuracy vs pulse count so the completeness
floor is explicit. (3) Detection imposes a low-S/N cut; the tail test operates above a
Clauset-selected lower bound, well clear of that cut.

Data (GATE-0 2026-07-10): open SQLite3 DB (GitHub LFS `Astroyx/Pulsar_collection`
`Pulsar_fits_database_v1.zip`, ~1.5 GB, no auth; mirror CSIRO DAP 10.25919/34am-zx04). Reuse:
`frbstats.fit_power_law`/`select_xmin`, `ppdot.spindown_luminosity`/`fetch_atnf_ppdot`, and the
`stacking` injection-recovery discipline.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .frbstats import fit_power_law

__all__ = [
    "fit_lognormal",
    "vuong_powerlaw_vs_lognormal",
    "fit_energy_tail",
    "synthetic_pulses",
    "inject_recover_tail",
    "census",
    "count_confound",
    "tail_vs_edot",
    "load_pulse_sn",
    "run",
]

# NB: the archetypal giant-pulse emitters (Crab J0534+2200, B1937+21 J1939+2134, Vela J0835-4510) do
# NOT survive the min-count cut in the 1997-2001 Parkes sample, so the recover-a-known is synthetic
# (inject_recover_tail) plus the one present tail-heavy source, B0950+08 (J0953+0755).
MIN_PULSES = 50  # per-pulsar minimum single-pulse count for a tail fit (pre-registered cut)
SIG = 0.05  # Vuong two-sided significance for preferring one model over the other

# PTE-II open data release (GATE-0 verified 2026-07-10).
PTE2_ZIP_URL = "https://github.com/Astroyx/Pulsar_collection/raw/main/Pulsar_fits_database_v1.zip"
PTE2_DB_NAME = "Pulsar_fits_database_v1.db"


def fit_lognormal(x: np.ndarray) -> tuple[float, float]:
    """MLE log-normal parameters ``(mu, sigma)`` of positive data ``x`` (mean/std of ``log x``)."""
    lx = np.log(np.asarray(x, float))
    return float(lx.mean()), float(lx.std(ddof=0))


def vuong_powerlaw_vs_lognormal(
    x_tail: np.ndarray, xmin: float, gamma: float, mu: float, sigma: float
) -> tuple[float, float]:
    """Normalized Vuong likelihood-ratio for power-law vs (truncated) log-normal on the tail.

    Both densities are normalized on ``x >= xmin``: the continuous power-law
    :math:`p(x)=(\\gamma-1)\\,x_{\\min}^{\\gamma-1}x^{-\\gamma}` and the log-normal truncated below
    ``xmin``. The per-point log-likelihood difference :math:`\\ell_i` gives the normalized statistic
    :math:`V=\\sqrt{n}\\,\\bar\\ell/\\mathrm{std}(\\ell)`, which is standard normal under the null of
    equally-good fits (Clauset et al. 2009). ``V>0`` favours the power-law (a heavy tail); the
    returned ``p`` is the two-sided normal tail probability.
    """
    from scipy.stats import norm

    x = np.asarray(x_tail, float)
    lx = np.log(x)
    lxmin = np.log(xmin)
    if gamma <= 1.0 or sigma <= 0 or x.size < 2:
        return 0.0, 1.0
    ll_pl = np.log(gamma - 1.0) - lxmin - gamma * (lx - lxmin)
    z = (lx - mu) / sigma
    trunc = 1.0 - norm.cdf((lxmin - mu) / sigma)
    trunc = max(trunc, 1e-300)
    ll_ln = -lx - np.log(sigma) - 0.5 * np.log(2 * np.pi) - 0.5 * z * z - np.log(trunc)
    ell = ll_pl - ll_ln
    s = ell.std(ddof=0)
    if s == 0:
        return 0.0, 1.0
    v = float(np.sqrt(ell.size) * ell.mean() / s)
    p = float(norm.sf(abs(v)) * 2.0)
    return v, p


def fit_energy_tail(
    sn: np.ndarray, *, min_pulses: int = MIN_PULSES, sig: float = SIG, k_sigma: float = 3.0
) -> dict:
    """Per-source giant-pulse test: does the extreme tail exceed the source's own log-normal bulk?

    A giant-pulse (heavy-tail) emitter has single pulses far above what a log-normal bulk allows. We
    count pulses above ``median(log S/N) + k_sigma * sigma_R`` (``n_giant``) and compare to the
    log-normal expectation (``n_exp``) with a Poisson excess test. Crucially, both the median and the
    right-side width ``sigma_R = P84.13(log S/N) - median`` are estimated from the UPPER half of the
    distribution only, which the single-pulse detection S/N floor does not truncate --- so a
    left-truncated bulk cannot fake an excess (the bias that inflates a naive whole-distribution fit
    or a bare Vuong test). A pulsar is ``heavy_tailed`` only on a significant, >=3-pulse excess: a
    positive detection of a giant-pulse population. The giant-tail power-law index and a Vuong
    power-law-vs-log-normal statistic are reported as secondary diagnostics. ``insufficient`` below
    ``min_pulses``.
    """
    from scipy.stats import norm, poisson

    x = np.asarray(sn, float)
    x = x[np.isfinite(x) & (x > 0)]
    n = int(x.size)
    if n < min_pulses:
        return {"n": n, "preferred": "insufficient", "heavy_tailed": False, "excess": float("nan")}
    lx = np.log(x)
    mu = float(
        np.median(lx)
    )  # median of log S/N (untruncated: the floor removes only the low tail)
    sigma = float(np.percentile(lx, 84.13) - mu)  # right-side width, floor-robust
    if sigma <= 0:
        sigma = float(lx.std(ddof=0)) or 1.0
    n_right = int((lx > mu).sum())  # the complete (untruncated) upper half
    x_hi = float(np.exp(mu + k_sigma * sigma))
    n_giant = int((x > x_hi).sum())
    frac_right = float(norm.sf(k_sigma) / 0.5)  # fraction of the upper half beyond k_sigma
    n_exp = float(n_right * frac_right)
    p_excess = float(poisson.sf(n_giant - 1, max(n_exp, 1e-9)))
    excess = float(np.log10((n_giant + 0.5) / (n_exp + 0.5)))
    heavy = bool(n_giant >= 3 and n_giant > n_exp and p_excess < sig)
    gamma = gamma_err = float("nan")
    vuong = vuong_p = float("nan")
    if n_giant >= 8:  # enough giants to characterise the tail's power-law index
        pl = fit_power_law(x, f_min=x_hi)
        gamma, gamma_err = float(pl.gamma), float(pl.gamma_err)
        vuong, vuong_p = vuong_powerlaw_vs_lognormal(x[x >= x_hi], x_hi, pl.gamma, mu, sigma)
    return {
        "n": n,
        "n_giant": n_giant,
        "n_exp": round(n_exp, 3),
        "excess": round(excess, 4),
        "p_excess": p_excess,
        "gamma": round(gamma, 4) if np.isfinite(gamma) else float("nan"),
        "gamma_err": round(gamma_err, 4) if np.isfinite(gamma_err) else float("nan"),
        "vuong": round(vuong, 4) if np.isfinite(vuong) else float("nan"),
        "vuong_p": vuong_p,
        "preferred": "power_law" if heavy else "lognormal",
        "heavy_tailed": heavy,
    }


def synthetic_pulses(
    *,
    kind: str = "lognormal",
    n: int = 400,
    mu: float = 2.5,  # median single-pulse S/N ~ e^2.5 ~ 12, comfortably above the detection floor
    sigma: float = 0.5,
    tail_index: float = 2.6,
    tail_frac: float = 0.04,
    tail_start: float = 4.0,
    sn_floor: float = 6.0,
    seed: int = 0,
) -> np.ndarray:
    """A synthetic single-pulse S/N set: a log-normal bulk, optionally with a power-law giant tail.

    ``kind='lognormal'`` is a pure log-normal bulk (no heavy tail); ``kind='heavy'`` replaces a
    fraction ``tail_frac`` of pulses with power-law (index ``tail_index``) giant pulses starting at
    ``tail_start`` x the bulk median --- the injected recover-a-known. A detection floor ``sn_floor``
    is applied (single-pulse searches only keep S/N above threshold).
    """
    rng = np.random.default_rng(seed)
    x = np.exp(rng.normal(mu, sigma, n))
    if kind == "heavy":
        k = int(round(tail_frac * n))
        if k > 0:
            x0 = tail_start * np.exp(mu)  # giant-pulse onset in S/N units
            u = rng.uniform(0.0, 1.0, k)
            giants = x0 * (1.0 - u) ** (-1.0 / (tail_index - 1.0))
            x[rng.choice(n, k, replace=False)] = giants
    return x[x >= sn_floor]


def inject_recover_tail(*, counts=(80, 150, 300, 600), n_each: int = 40, seed: int = 0) -> dict:
    """Recover-a-known: classify injected log-normal vs heavy-tailed sets, accuracy vs pulse count.

    For each pulse count we build ``n_each`` pure-log-normal and ``n_each`` heavy-tailed synthetic
    pulsars and record how often `fit_energy_tail` calls them correctly. Heavy-tail detection is
    honestly count-limited: a real giant-pulse tail in a low-count source is genuinely
    unrecoverable, so completeness rises with count (the reported floor). Also returns the
    false-positive rate (pure log-normals misclassified as heavy).
    """
    rng = np.random.default_rng(seed)
    curve = {}
    fp = 0
    fp_total = 0
    for c in counts:
        tp = 0
        for _ in range(n_each):
            s = int(rng.integers(1 << 30))
            heavy = synthetic_pulses(kind="heavy", n=c, seed=s)
            if fit_energy_tail(heavy)["heavy_tailed"]:
                tp += 1
        for _ in range(n_each):
            s = int(rng.integers(1 << 30))
            ln = synthetic_pulses(kind="lognormal", n=c, seed=s)
            fp_total += 1
            if fit_energy_tail(ln)["heavy_tailed"]:
                fp += 1
        curve[f"n{c}"] = round(tp / n_each, 3)
    return {
        "completeness_vs_count": curve,
        "false_positive_rate": round(fp / max(fp_total, 1), 4),
    }


def fit_energy_tail_truncated(sn: np.ndarray, *, k_sigma: float = 3.0, sig: float = 0.05) -> dict:
    """The floor-corrected giant-pulse excess test.

    The original estimator took mu = median and sigma = P84 - median of the OBSERVED log-S/N,
    then expected n_right * sf(3)/0.5 pulses above mu + 3 sigma. That is unbiased only if the
    detection floor removed nothing: a floor that removes a fraction f of the underlying
    pulses shifts the observed median up to the (1+f)/2 quantile and compresses the observed
    P84 - median, so the "3 sigma" threshold sits closer than 3 true sigma and every source
    shows a spurious excess (2.03x at f = 0.3; the committed census's aggregate obs/exp of
    1.85 implies f ~ 0.25-0.3 -- a referee derived this in closed form and the census
    confirmed it in every count bin).

    The fix: treat the floor as a KNOWN left-truncation at the observed minimum and fit
    (mu, sigma) by truncated-normal maximum likelihood. The expected count above
    mu + k sigma among the n observed pulses is then n * sf(k) / sf((a - mu)/sigma), which is
    correct for any f. Reports f_hat = Phi((a - mu)/sigma), the estimated truncated fraction.
    """
    from scipy import stats
    from scipy.optimize import minimize

    lx = np.log(np.asarray(sn, float))
    n = int(lx.size)
    a = float(lx.min()) - 1e-9  # the known truncation point: nothing below the floor survives

    def nll(theta: np.ndarray) -> float:
        mu, ls = theta
        sig_ = np.exp(ls)
        z = (lx - mu) / sig_
        za = (a - mu) / sig_
        return float(n * np.log(stats.norm.sf(za) + 1e-300) + 0.5 * np.sum(z**2) + n * ls)

    x0 = np.array(
        [float(np.median(lx)), np.log(max(float(np.percentile(lx, 84.13) - np.median(lx)), 1e-3))]
    )
    r = minimize(
        nll, x0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 2000}
    )
    mu, sig_ = float(r.x[0]), float(np.exp(r.x[1]))
    za = (a - mu) / sig_
    f_hat = float(stats.norm.cdf(za))
    thresh = mu + k_sigma * sig_
    n_giant = int((lx > thresh).sum())
    n_exp = float(n * stats.norm.sf(k_sigma) / max(stats.norm.sf(za), 1e-300))
    p_excess = float(stats.poisson.sf(n_giant - 1, n_exp)) if n_giant > 0 else 1.0
    heavy = bool(n_giant >= 3 and n_giant > n_exp and p_excess < sig)
    return {
        "mu": round(mu, 4),
        "sigma": round(sig_, 4),
        "f_hat": round(f_hat, 4),
        "n_giant": n_giant,
        "n_exp": round(n_exp, 3),
        "p_excess": p_excess,
        "heavy_tailed": heavy,
        "converged": bool(r.success),
    }


def census_statistics(rows: list[dict]) -> dict:
    """Population-level statistics over a completed census, from the committed rows alone.

    Added 2026-08-13 after a referee recomputed all of these from results/pte2_census.json and
    found the paper's headline did not survive them. Everything here derives from the census
    file, so it is committed evidence about committed evidence -- reproducible offline.

    - Multiple comparisons: 136 tests at p<0.05 expect ~6.8 chance positives, which the paper
      never stated. Benjamini-Hochberg and Bonferroni counts are what the headline can carry.
    - The aggregate observed/expected giant ratio: for a truly floor-robust estimator on
      log-normal sources this is 1; the census gives 1.85, the signature of the left-truncation
      bias in the median/P84 estimator (a floor removing fraction f of pulses shifts both).
    - The per-count-bin heavy fraction, which is NOT monotonic in count -- it peaks at 400-800
      pulses and falls to 0.14 in the highest-power bin, contradicting the "more power on more
      pulses" mechanism the paper asserted from a two-bin split.
    - The partial Spearman of excess vs log Edot CONTROLLING for log n: `excess` is count-biased
      (rho = -0.36 with log n) and Edot correlates with count, so the raw correlation carried a
      spurious negative term; controlling for count flips the sign to +0.10 (still null).
    - The Vuong table for every flagged source with >= 8 giants: a promised diagnostic that was
      computed and never reported. Its only significant value (J1243-6423, p = 3e-6) PREFERS
      the log-normal -- which supports the paper's null and should have been shown.
    """
    from scipy import stats

    if isinstance(rows, dict):
        rows = rows["all"]
    fitted = [r for r in rows if r.get("p_excess") is not None]
    n = len(fitted)
    pvals = np.array([r["p_excess"] for r in fitted])
    heavy = np.array([bool(r["heavy_tailed"]) for r in fitted])
    counts = np.array([r["n"] for r in fitted])
    n_giant = np.array([r["n_giant"] for r in fitted])
    n_exp = np.array([r["n_exp"] for r in fitted])

    # Benjamini-Hochberg at 0.05
    order = np.argsort(pvals)
    bh_thresh = 0.05 * (np.arange(1, n + 1) / n)
    passed = pvals[order] <= bh_thresh
    n_bh = int(np.max(np.nonzero(passed)[0]) + 1) if passed.any() else 0
    n_bonf = int((pvals < 0.05 / n).sum())

    bins = [(50, 100), (100, 200), (200, 400), (400, 800), (800, 1600), (1600, 10**9)]
    per_bin = []
    for lo, hi in bins:
        sel = (counts >= lo) & (counts < hi)
        per_bin.append(
            {
                "count_lo": lo,
                "count_hi": None if hi >= 10**9 else hi,
                "n": int(sel.sum()),
                "heavy_frac": round(float(heavy[sel].mean()), 3) if sel.any() else None,
                "obs_over_exp": (
                    round(float(n_giant[sel].sum() / max(n_exp[sel].sum(), 1e-9)), 2)
                    if sel.any()
                    else None
                ),
            }
        )

    # partial Spearman excess vs log Edot | log n
    edot_rows = [r for r in fitted if r.get("edot") is not None]
    partial = None
    if len(edot_rows) > 10:
        ex = np.array([r["excess"] for r in edot_rows])
        le = np.log10([r["edot"] for r in edot_rows])
        ln = np.log10([r["n"] for r in edot_rows])
        rex = stats.rankdata(ex)
        rle = stats.rankdata(le)
        rln = stats.rankdata(ln)

        def _resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
            k = np.polyfit(b, a, 1)
            return a - np.polyval(k, b)

        r_p, p_p = stats.pearsonr(_resid(rex, rln), _resid(rle, rln))
        partial = {
            "rho_partial": round(float(r_p), 3),
            "p_partial": round(float(p_p), 4),
            "rho_excess_logn": round(float(stats.spearmanr(ex, ln).statistic), 3),
            "rho_edot_logn": round(float(stats.spearmanr(le, ln).statistic), 3),
            "n": len(edot_rows),
        }

    vuong = [
        {
            "name": r["jname"],
            "n_giant": r["n_giant"],
            "gamma": r.get("gamma"),
            "gamma_err": r.get("gamma_err"),
            "vuong": r.get("vuong"),
            "vuong_p": r.get("vuong_p"),
        }
        for r in fitted
        if r["heavy_tailed"] and r["n_giant"] >= 8
    ]
    return {
        "n_fitted": n,
        "n_heavy_nominal": int(heavy.sum()),
        "n_heavy_bh05": n_bh,
        "n_heavy_bonferroni": n_bonf,
        "expected_chance_at_005": round(0.05 * n, 1),
        "aggregate_obs_over_exp": round(float(n_giant.sum() / n_exp.sum()), 3),
        "sum_n_giant": int(n_giant.sum()),
        "sum_n_exp": round(float(n_exp.sum()), 1),
        "heavy_frac_by_count_bin": per_bin,
        "edot_partial": partial,
        "vuong_flagged": vuong,
    }


def census(per_pulsar: dict, *, min_pulses: int = MIN_PULSES) -> list[dict]:
    """Run `fit_energy_tail` for every pulsar with >= ``min_pulses`` single pulses.

    ``per_pulsar`` maps J-name -> {"sn": array, ...metadata}. Returns one row per fitted pulsar
    (jname + tail fit + carried metadata), sorted by decreasing Vuong statistic (most heavy-tailed
    first). Sources below the count cut are dropped and reported via `run`'s counts.
    """
    rows = []
    for jname, d in per_pulsar.items():
        fit = fit_energy_tail(np.asarray(d["sn"], float), min_pulses=min_pulses)
        if fit["preferred"] == "insufficient":
            continue
        row = {"jname": jname, **fit}
        for k in ("p0", "s1400", "edot"):
            if k in d:
                row[k] = d[k]
        rows.append(row)
    rows.sort(
        key=lambda r: r["excess"] if np.isfinite(r.get("excess", np.nan)) else -1e9, reverse=True
    )
    return rows


def count_confound(rows: list[dict]) -> dict:
    """Is the heavy-tail classification driven by astrophysics or by pulse count (detection power)?

    The Poisson excess test has more power to flag a small departure in a source with more single
    pulses, so a real worry is that ``heavy_tailed`` just tracks pulse count. This tests it directly:
    the Mann-Whitney comparison of the pulse count of heavy vs non-heavy sources, the heavy fraction
    in the low- vs high-count halves, and the median giant-tail power-law index (a shallow ~2--3 index
    is a classic giant pulse; a steep ~10 index is a mild excess, not a giant-pulse population).
    ``count_limited`` is set when heavy sources have significantly more pulses --- the honest caveat
    that the flagged fraction is a detection-power floor, not a clean astrophysical incidence.
    """
    from scipy.stats import mannwhitneyu

    n = np.array([r["n"] for r in rows], float)
    heavy = np.array([r["heavy_tailed"] for r in rows], bool)
    gam = np.array(
        [r["gamma"] for r in rows if r["heavy_tailed"] and np.isfinite(r.get("gamma", np.nan))],
        float,
    )
    out: dict = {
        "n_fit": int(n.size),
        "n_heavy": int(heavy.sum()),
        "median_gamma_heavy": round(float(np.median(gam)), 2) if gam.size else float("nan"),
    }
    if heavy.any() and (~heavy).any():
        med = float(np.median(n))
        out["median_n_heavy"] = round(float(np.median(n[heavy])), 1)
        out["median_n_nonheavy"] = round(float(np.median(n[~heavy])), 1)
        out["heavy_frac_lowcount"] = round(float(heavy[n <= med].mean()), 3)
        out["heavy_frac_highcount"] = round(float(heavy[n > med].mean()), 3)
        if heavy.sum() >= 3 and (~heavy).sum() >= 3:
            _, pmw = mannwhitneyu(n[heavy], n[~heavy], alternative="greater")
            out["count_mw_p"] = float(pmw)
            out["count_limited"] = bool(pmw < 0.05)
    return out


def tail_vs_edot(rows: list[dict], edot_by_jname: dict) -> dict:
    """Rank-correlate the giant-pulse excess and heavy-tail incidence against spin-down luminosity.

    Attaches ``log10(Edot)`` to each fitted pulsar, then computes (a) the Spearman rank correlation
    of the giant-pulse ``excess`` statistic with ``log Edot`` and (b) the median ``log Edot`` of the
    heavy-tailed vs non-heavy-tailed sources with a Mann-Whitney test. Giant pulses are associated
    with high Edot (Crab/millisecond emitters), so a positive trend is the physical expectation ---
    reported honestly with its p-value, whatever the sign.
    """
    from scipy.stats import mannwhitneyu, spearmanr

    le_l: list[float] = []
    ex_l: list[float] = []
    hv_l: list[bool] = []
    for r in rows:
        ed = edot_by_jname.get(r["jname"])
        if ed is None or not np.isfinite(ed) or ed <= 0 or not np.isfinite(r.get("excess", np.nan)):
            continue
        le_l.append(float(np.log10(ed)))
        ex_l.append(r["excess"])
        hv_l.append(r["heavy_tailed"])
    le = np.asarray(le_l, float)
    ex = np.asarray(ex_l, float)
    heavy = np.asarray(hv_l, bool)
    out: dict = {"n_matched": int(le.size)}
    if le.size >= 5:
        rho, p = spearmanr(le, ex)
        out["spearman_excess_logedot"] = round(float(rho), 3)
        out["spearman_p"] = float(p)
    if heavy.sum() >= 3 and (~heavy).sum() >= 3:
        u, pu = mannwhitneyu(le[heavy], le[~heavy], alternative="two-sided")
        out["logedot_heavy_median"] = round(float(np.median(le[heavy])), 3)
        out["logedot_nonheavy_median"] = round(float(np.median(le[~heavy])), 3)
        out["mannwhitney_p"] = float(pu)
    return out


def load_pulse_sn(
    db_path: str | Path, *, sn_col: str = "snr_max", limit_pulsars: int | None = None
) -> dict:
    """Load per-pulsar single-pulse S/N arrays from the PTE-II SQLite database.

    Joins the segment table (one row per single pulse, carrying ``snr_max``) to ``file`` and
    ``pulsar`` to group by J-name, and pulls ``p0``/``s1400`` metadata. The segment table name and
    S/N column are auto-detected (the schema documents ``fileSegment``/``seg_file`` with
    ``snr_max``), so a minor release-to-release rename does not break the loader. Returns
    ``{jname: {"sn": np.ndarray, "p0": float, "s1400": float, "n": int}}``.
    """
    import sqlite3

    con = sqlite3.connect(str(db_path))
    try:
        con.row_factory = sqlite3.Row
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        seg = next((t for t in ("fileSegment", "seg_file", "segment") if t in tables), None)
        if seg is None:
            raise ValueError(f"no segment table found in {tables}")
        seg_cols = {r[1] for r in con.execute(f"PRAGMA table_info({seg})")}
        sn = sn_col if sn_col in seg_cols else next(c for c in seg_cols if "snr" in c.lower())
        q = (
            f"SELECT p.jname AS jname, p.p0 AS p0, p.s1400 AS s1400, s.{sn} AS sn "
            f"FROM {seg} s JOIN file f ON s.pfLinkID = f.pfLinkID "
            f"JOIN pulsar p ON f.pulsarID = p.pulsarID WHERE s.{sn} IS NOT NULL"
        )
        out: dict = {}
        for row in con.execute(q):
            d = out.setdefault(row["jname"], {"sn": [], "p0": row["p0"], "s1400": row["s1400"]})
            d["sn"].append(float(row["sn"]))
        for d in out.values():
            d["sn"] = np.asarray(d["sn"], float)
            d["n"] = int(d["sn"].size)
        if limit_pulsars is not None:
            out = dict(list(out.items())[:limit_pulsars])
        return out
    finally:
        con.close()


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known (classify injected heavy vs log-normal, accuracy vs count)."""
    import json

    if offline:
        metrics: dict = _synthetic_metrics()
    else:  # pragma: no cover - real leg reads the PTE-II SQLite DB + ATNF
        metrics = _real_census(out)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "pte2_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, op / "papers" / "pte2" / "figures")
    _write_macros(metrics, op / "papers" / "pte2" / "generated" / "macros.tex")
    return metrics


def _synthetic_metrics() -> dict:
    """Recover-a-known: the model-selection test classifies injected heavy vs log-normal sets."""
    rec = inject_recover_tail()
    # a single clean demonstration pair at high count
    heavy = fit_energy_tail(synthetic_pulses(kind="heavy", n=600, seed=1))
    lognormal = fit_energy_tail(synthetic_pulses(kind="lognormal", n=600, seed=2))
    return {
        "source": "synthetic single-pulse sets (log-normal bulk +/- injected power-law giant tail)",
        "is_real": False,
        "completeness_vs_count": rec["completeness_vs_count"],
        "false_positive_rate": rec["false_positive_rate"],
        "demo_heavy_preferred": heavy["preferred"],
        "demo_heavy_gamma": heavy["gamma"],
        "demo_lognormal_preferred": lognormal["preferred"],
        "recovered": bool(heavy["heavy_tailed"] and not lognormal["heavy_tailed"]),
    }


def _real_census(out: str) -> dict:  # pragma: no cover - network + multi-GB SQLite
    """Real leg: download the PTE-II DB, run the 363-pulsar census + Edot cross-match."""
    from .pte2_real import run_real_census

    return run_real_census(out)


def _figure(m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.8))
    # left: the two synthetic distributions (log-normal vs heavy-tailed) as complementary CDFs
    for kind, c in (("lognormal", "C0"), ("heavy", "C3")):
        x = np.sort(synthetic_pulses(kind=kind, n=4000, seed=7))
        ccdf = 1.0 - np.arange(x.size) / x.size
        ax1.loglog(x, ccdf, "-", color=c, lw=1.2, label=kind)
    ax1.set(xlabel="single-pulse S/N", ylabel="CCDF", title="Injected populations")
    ax1.legend(fontsize=8)
    # right: recover-a-known completeness vs pulse count
    cov = m.get("completeness_vs_count", {})
    if cov:
        counts = [int(k[1:]) for k in cov]
        ax2.plot(counts, [cov[k] for k in cov], "o-", color="C0")
    ax2.set(
        xlabel="pulses per source",
        ylabel="heavy-tail completeness",
        title=f"Recover-a-known (FP={m.get('false_positive_rate')})",
        ylim=(0, 1.05),
    )
    fig.tight_layout()
    fig.savefig(out / "pte2.pdf")
    plt.close(fig)


def _stats_macros() -> list[str]:
    """Macros from results/pte2_stats.json -- the bracketing statistics a referee demanded.

    Read from the committed file, never recomputed at write time, so what the paper quotes is
    what a reader can check.
    """
    import json

    path = Path("results/pte2_stats.json")
    if not path.is_file():
        return []
    st = json.loads(path.read_text())
    nv, tr = st["naive"], st["truncated"]
    return [
        "% Bracketing statistics (results/pte2_stats.json): the naive floor-ignorant null vs",
        "% the floor-marginalised truncated-normal null. Neither endpoint is the answer; the",
        "% bracket is.",
        rf"\newcommand{{\ptStatBH}}{{{nv['n_heavy_bh05']}}}",
        rf"\newcommand{{\ptStatBonf}}{{{nv['n_heavy_bonferroni']}}}",
        rf"\newcommand{{\ptStatChance}}{{{nv['expected_chance_at_005']:g}}}",
        rf"\newcommand{{\ptStatAggRatio}}{{{nv['aggregate_obs_over_exp']:.2f}}}",
        rf"\newcommand{{\ptStatSumGiant}}{{{nv['sum_n_giant']}}}",
        rf"\newcommand{{\ptStatSumExp}}{{{nv['sum_n_exp']:g}}}",
        rf"\newcommand{{\ptStatRhoPartial}}{{{nv['edot_partial']['rho_partial']:g}}}",
        rf"\newcommand{{\ptStatRhoPartialP}}{{{nv['edot_partial']['p_partial']:g}}}",
        rf"\newcommand{{\ptStatRhoExN}}{{{nv['edot_partial']['rho_excess_logn']:g}}}",
        rf"\newcommand{{\ptTruncNHeavy}}{{{tr['k_sweep']['3']['n_heavy']}}}",
        rf"\newcommand{{\ptTruncFMedian}}{{{tr['f_hat_summary']['median']:g}}}",
        rf"\newcommand{{\ptTruncFNine}}{{{tr['n_fhat_gt_09']}}}",
    ]


def _write_macros(m: dict, path: str | Path) -> None:
    """Emit both namespaces: ptSyn* (recover-a-known, always live) + ptReal* (the real census)."""
    from .report import _fmt_p

    def g(src: dict, key: str) -> str:
        v = src.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        return str(v)

    syn = _synthetic_metrics() if m.get("is_real") else m
    real = m if m.get("is_real") else {}
    cov = syn.get("completeness_vs_count", {})
    lines = [
        "% Auto-generated by jansky_research.pte2._write_macros -- do not edit.",
        "% ptSyn* = synthetic recover-a-known (always live); ptReal* = the real PTE-II census.",
        rf"\newcommand{{\ptSource}}{{{m['source']}}}",
        *_stats_macros(),
        rf"\newcommand{{\ptSynFP}}{{{g(syn, 'false_positive_rate')}}}",
        rf"\newcommand{{\ptSynRecovered}}{{{'yes' if syn.get('recovered') else 'no'}}}",
        rf"\newcommand{{\ptSynCompLow}}{{{cov.get('n80', '--')}}}",
        rf"\newcommand{{\ptSynCompHigh}}{{{cov.get('n600', '--')}}}",
    ]
    for macro, key in (
        ("NFit", "n_fit"),
        ("NHeavy", "n_heavy"),
        ("HeavyFrac", "heavy_fraction"),
        ("NMatched", "n_matched"),
        ("Spearman", "spearman_excess_logedot"),
        ("LogEdotHeavy", "logedot_heavy_median"),
        ("LogEdotNon", "logedot_nonheavy_median"),
        ("MedGamma", "median_gamma_heavy"),
        ("MedNHeavy", "median_n_heavy"),
        ("MedNNon", "median_n_nonheavy"),
        ("HeavyFracLow", "heavy_frac_lowcount"),
        ("HeavyFracHigh", "heavy_frac_highcount"),
    ):
        lines.append(rf"\newcommand{{\ptReal{macro}}}{{{g(real, key)}}}")
    lines.append(
        rf"\newcommand{{\ptRealCountLimited}}{{{'yes' if real.get('count_limited') else ('--' if not real else 'no')}}}"
    )
    for macro, key in (
        ("SpearmanP", "spearman_p"),
        ("MannWhitneyP", "mannwhitney_p"),
        ("CountMWP", "count_mw_p"),
    ):
        v = real.get(key)
        val = _fmt_p(v) if isinstance(v, (int, float)) and np.isfinite(v) else "--"
        lines.append(rf"\newcommand{{\ptReal{macro}}}{{{val}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    ap = argparse.ArgumentParser(description="PTE-II single-pulse heavy-tail census.")
    ap.add_argument("--out", default=".")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
