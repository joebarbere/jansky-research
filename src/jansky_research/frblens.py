"""Lensed-repeater test: recurring burst-to-burst delay patterns in CHIME Cat 2 (plan 42).

A strongly-lensed one-off FRB masquerades as a "repeater" whose bursts arrive in image
multiplets with a FIXED mutual delay (Dai & Lu 2017, 10.3847/1538-4357/aa8873; Li+2018, Nat.
Comm. 9:3833; review arXiv:2412.01536). No observational catalogue-level search exists (GATE-0
full-text pass 2026-07-06): CHIME's published lensing work is all intra-burst/us--ms regime ---
baseband echoes (Leung+2022/Kader+2022) and the Cat-2 microlensing search arXiv:2605.19653,
which autocorrelated single-burst dynamic spectra of NON-repeaters (ms delays, ~10^2-10^3 Msun
lenses). This slice is the untouched day-to-month galaxy-lens regime, searched across every
Cat-2 repeater. Expected yield ~ 0 (lensing optical depth ~1e-4): the deliverable is the first
empirical upper limit on the lensed fraction, framed that way from day one.

Statistic: for each repeater, scan candidate delays Delta (from the observed pairwise delays,
Delta > 1 d --- sub-transit pairs are indistinguishable from intrinsic clustering, stated) in
the BARYCENTRIC frame (mjd_400 is topocentric; a barycentre-fixed delay drifts by up to
~+-150 s/26 d in topocentric pair delays --- GATE-2 catch); for each Delta, greedily match
DISJOINT burst pairs (t, t+Delta) within 5 s and count matches M(Delta); the per-source
statistic is M_max. Consistency cut: matched pairs must agree in DM within 3*sqrt(2)*sigma of
the source's own fitted-DM scatter (the fitburst scatter is structure-driven, ~2-7 pc/cc for
active repeaters --- a naive 1 pc/cc "measurement floor" would reject most genuine image
pairs, a second GATE-2 catch). FAP: a PHASE scramble
--- keep each burst's sidereal day, permute the within-transit-window phases --- so all
day-scale structure (activity epochs, periodic windows, clustering) survives into the null and
only sub-window fixed-delay coherence, the actual lensing signature, can beat it. The first
real-catalogue contact proved this choice load-bearing: a day-scramble null (frbwait's, right
for periodicity) false-positived the three most clustered repeaters. A transit-survey geometry
fact the injection map makes explicit: BOTH images of a burst are only detectable when
Delta mod (sidereal day) falls inside the source's transit window, so sensitivity lives near
integer sidereal-day delays.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .frbwait import SIDEREAL_DAY, load_catalog2, repeater_trains

__all__ = [
    "all_pairs_delays",
    "match_count_at_delay",
    "recurring_delay_search",
    "scramble_fap",
    "inject_lensed_train",
    "sensitivity_map",
    "lensed_fraction_limit",
    "run",
]

MIN_DELAY_DAYS = 1.0  # below one transit spacing, image pairs and intrinsic clustering merge


def all_pairs_delays(mjds: np.ndarray, *, min_delay: float = MIN_DELAY_DAYS) -> np.ndarray:
    """All pairwise positive delays (days) above ``min_delay``, sorted."""
    t = np.sort(np.asarray(mjds, float))
    d = (t[None, :] - t[:, None])[np.triu_indices(t.size, k=1)]
    d = d[d >= min_delay]
    return np.sort(d)


def match_count_at_delay(
    mjds: np.ndarray,
    delay: float,
    *,
    tol_days: float,
    dm: np.ndarray | None = None,
    dm_tol: float = np.inf,
) -> int:
    """Greedy count of DISJOINT burst pairs separated by ``delay`` (+- ``tol_days``).

    A lensed one-off with image delay Delta turns each intrinsic burst into a (leading,
    trailing) pair: the pairs are disjoint by construction, so disjoint matching is the right
    counting. Earliest-first greedy matching is optimal for interval matching of this form.
    With ``dm`` given, a pair only matches if |DM_i - DM_j| <= ``dm_tol``.
    """
    t = np.asarray(mjds, float)
    order = np.argsort(t)
    t = t[order]
    dms = np.asarray(dm, float)[order] if dm is not None else None
    used = np.zeros(t.size, bool)
    count = 0
    for i in range(t.size):
        if used[i]:
            continue
        lo = np.searchsorted(t, t[i] + delay - tol_days, side="left")
        hi = np.searchsorted(t, t[i] + delay + tol_days, side="right")
        for j in range(lo, hi):
            if j == i or used[j]:
                continue
            if dms is not None and abs(dms[i] - dms[j]) > dm_tol:
                continue
            used[i] = used[j] = True
            count += 1
            break
    return count


TOL_DAYS = 5.0 / 86400.0  # >> TOA errors (median ~1.2 ms), DM-shift 26 ms/pc cc, bary residual


def barycentric_offset_fn(ra_deg: float, dec_deg: float):
    """Vectorised topocentric-MJD -> barycentric light-travel-time offset (days) for CHIME.

    Cat-2 `mjd_400` TOAs are TOPOCENTRIC (catalogue Appendix A); a lens delay is fixed at the
    BARYCENTER, and the annual Roemer term makes a fixed barycentric delay drift by up to
    ~+-150 s (at Delta ~ 26 d) in topocentric pair delays --- far beyond the 5-s tolerance.
    The search therefore runs on t_bary = t_topo + this offset. Uses astropy's built-in
    ephemeris (offline-safe); source-position uncertainty at arcminutes contributes < 0.2 s.
    """
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    loc = EarthLocation(lat=49.3208 * u.deg, lon=-119.6236 * u.deg, height=545.0 * u.m)
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)

    def fn(mjd: np.ndarray) -> np.ndarray:
        t = Time(np.atleast_1d(np.asarray(mjd, float)), format="mjd", scale="utc", location=loc)
        return np.asarray(t.light_travel_time(coord, kind="barycentric").to_value("day"))

    return fn


def phase_scramble(
    mjds: np.ndarray, rng: np.random.Generator, *, bary_offsets: np.ndarray | None = None
) -> np.ndarray:
    """Lensing null: keep each burst's sidereal DAY, permute the within-window phases.

    The lensing signature is a delay fixed to ~ms across every image pair; intrinsic bursts on
    the same days merely share the ~15-min transit window. Permuting the sidereal phases among
    bursts preserves the per-day burst counts exactly (activity epochs, periodic windows,
    clustering --- the things that made the day-scramble null flag 20220912A/20180916B/20201124A
    as false positives on first contact with the real catalogue) while destroying sub-window
    fixed-delay coherence --- precisely the signal.

    ``mjds`` are TOPOCENTRIC (the domain where the transit window lives). With
    ``bary_offsets`` (per-burst, day-keyed --- days are unchanged by the permutation, and the
    Roemer offset varies by < 0.01 s across a transit window), the returned null times are
    BARYCENTRIC, matching the search domain.
    """
    t = np.asarray(mjds, float)
    n = np.floor(t / SIDEREAL_DAY)
    phase = t - n * SIDEREAL_DAY
    scrambled = n * SIDEREAL_DAY + rng.permutation(phase)
    if bary_offsets is not None:
        scrambled = scrambled + np.asarray(bary_offsets, float)
    return np.sort(scrambled)


def recurring_delay_search(
    mjds: np.ndarray,
    *,
    dm: np.ndarray | None = None,
    tol_days: float = TOL_DAYS,
    dm_tol: float = 1.0,
    min_delay: float = MIN_DELAY_DAYS,
) -> dict:
    """Scan candidate delays (the observed pairwise delays) for the maximum recurring count.

    Returns M_max, the delay achieving it, and the number of candidates scanned. The input
    times must be BARYCENTRIC for real data (a lens delay is only fixed at the barycenter).
    tol default 5 s: catalogue TOA errors are ms-scale (repeater median ~1.2 ms, max 0.15 s),
    a Delta-DM of 1 pc/cc shifts the 400-MHz TOA by only 26 ms, residual barycentring error
    (arcminute positions, ephemeris) is < 0.2 s, and Dai & Lu-style secular delay drift is
    well under seconds across the span --- 5 s covers all of these with margin while staying
    far inside the ~15-min transit window, which is what separates fixed-delay coherence from
    intrinsic same-window clustering.
    """
    cands = all_pairs_delays(mjds, min_delay=min_delay)
    best_m, best_delay = 0, float("nan")
    if cands.size:
        # vectorised pre-filter: the sliding-window pair count over the sorted delay multiset
        # upper-bounds the disjoint match count, so the exact greedy matcher only needs to run
        # on the top window-count candidates (keeps the 476-burst repeater tractable)
        lo = np.searchsorted(cands, cands - tol_days, side="left")
        hi = np.searchsorted(cands, cands + tol_days, side="right")
        window = hi - lo
        top = np.argsort(window)[::-1][:32]
        top = top[window[top] >= 2]  # M_max >= 2 needs >= 2 co-delayed pairs
        for i in top:
            if window[i] <= best_m:
                break  # window count bounds the greedy count; sorted desc -> done
            m = match_count_at_delay(mjds, float(cands[i]), tol_days=tol_days, dm=dm, dm_tol=dm_tol)
            if m > best_m:
                best_m, best_delay = m, float(cands[i])
        if best_m == 0:
            # no recurring delay: report the (single-pair) maximum honestly
            j = int(np.argmax(window))
            best_m = match_count_at_delay(
                mjds, float(cands[j]), tol_days=tol_days, dm=dm, dm_tol=dm_tol
            )
            best_delay = float(cands[j])
    return {
        "m_max": int(best_m),
        "best_delay": best_delay,
        "n_candidates": int(cands.size),
        "tol_days": tol_days,
        "dm_tol": dm_tol,
    }


def scramble_fap(
    mjds: np.ndarray,
    m_obs: int,
    *,
    dm: np.ndarray | None = None,
    bary_offsets: np.ndarray | None = None,
    tol_days: float = TOL_DAYS,
    dm_tol: float = 1.0,
    n_scramble: int = 100,
    seed: int = 0,
) -> dict:
    """Phase-scramble FAP for M_max: p = (k+1)/(n+1) with the day-scale structure preserved.

    ``mjds`` are TOPOCENTRIC; the null (`phase_scramble`) keeps each burst's day and permutes
    the within-window phases, then re-applies the day-keyed ``bary_offsets`` so null and
    observed statistics live in the same (barycentric) domain. Intrinsic clustering, activity
    epochs, and periodic windows survive into every null realisation, so only genuine
    sub-window fixed-delay coherence can beat it. (The `frbwait`-style day scramble is the
    WRONG null here --- it flags intrinsic clustering; kept for the periodicity census only.)
    Under the permutation the DM-to-arrival association is scrambled too --- the correct H0.
    """
    rng = np.random.default_rng(seed)
    null_m = np.empty(n_scramble, dtype=int)
    for i in range(n_scramble):
        scr = phase_scramble(mjds, rng, bary_offsets=bary_offsets)
        null_m[i] = recurring_delay_search(scr, dm=dm, tol_days=tol_days, dm_tol=dm_tol)["m_max"]
    p = float((null_m >= m_obs).sum() + 1) / (n_scramble + 1)
    return {"p_value": p, "null_m": null_m, "m_obs": int(m_obs), "n_scramble": n_scramble}


def inject_lensed_train(
    mjds: np.ndarray,
    fluence: np.ndarray,
    *,
    delay: float,
    mag_ratio: float,
    fluence_limit: float | None = None,
    transit_window_min: float = 15.0,
    bary_fn=None,
    dm: np.ndarray | None = None,
    dm_sigma: float = 0.0,
    seed: int = 0,
) -> dict:
    """Turn a real burst train into a lensed one: add a trailing image per burst.

    ``delay`` is fixed in the BARYCENTRIC frame (the physical statement). With ``bary_fn``
    (from `barycentric_offset_fn`), the trailing image's topocentric arrival is
    t2 = (t + ltt(t) + delay) - ltt(t2) (one-step approximation; error << 1 s), so the annual
    Roemer drift the real data carries is modelled. The image is DETECTED only if (a) its
    fluence ``mag_ratio``-scaled clears ``fluence_limit`` (default: the faintest observed
    burst --- the source's empirical detection floor) and (b) it lands inside the source's
    TOPOCENTRIC transit window (width ``transit_window_min`` minutes) --- the transit-geometry
    selection that concentrates sensitivity near integer sidereal-day delays. With ``dm`` and
    ``dm_sigma``, each image gets an independently-"fitted" DM: the pair DM difference is
    drawn N(0, sqrt(2)*dm_sigma) --- the empirical per-source fit scatter --- so the search's
    DM-consistency cut costs the injections what it would cost real image pairs.
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(mjds, float)
    f = np.asarray(fluence, float)
    if fluence_limit is None:
        fluence_limit = float(np.nanmin(f))
    ltt = bary_fn(t) if bary_fn is not None else np.zeros_like(t)
    tb2 = t + ltt + delay
    t2 = tb2 - (bary_fn(tb2) if bary_fn is not None else 0.0)
    f2 = f * mag_ratio
    phase0 = np.median(t % SIDEREAL_DAY)
    w = transit_window_min / (24.0 * 60.0)
    dphase = np.abs(
        ((t2 % SIDEREAL_DAY) - phase0 + SIDEREAL_DAY / 2) % SIDEREAL_DAY - SIDEREAL_DAY / 2
    )
    det = (f2 >= fluence_limit) & (dphase <= w / 2)
    merged_topo = np.concatenate([t, t2[det]])
    order = np.argsort(merged_topo)
    out: dict = {
        "mjd": merged_topo[order],
        "n_images_detected": int(det.sum()),
        "detectable": bool(det.sum() >= 2),
    }
    if dm is not None:
        dm_arr = np.asarray(dm, float)
        dm_img = dm_arr[det] + rng.normal(0.0, np.sqrt(2.0) * dm_sigma, int(det.sum()))
        out["dm"] = np.concatenate([dm_arr, dm_img])[order]
    if bary_fn is not None:
        out["bary_offsets"] = bary_fn(out["mjd"])
    return out


def sensitivity_map(
    mjds: np.ndarray,
    fluence: np.ndarray,
    *,
    delays: np.ndarray,
    mag_ratios: np.ndarray,
    tol_days: float = TOL_DAYS,
    n_scramble: int = 50,
    detection_p: float | None = None,
    bary_fn=None,
    dm: np.ndarray | None = None,
    dm_sigma: float = 0.0,
    dm_tol: float = np.inf,
    seed: int = 0,
) -> dict:
    """Injection-recovery over a (delay, magnification-ratio) grid for one source.

    A cell is "sensitive" if the injected lensed train's M_max beats the phase-scramble null
    at p <= ``detection_p`` (default: the scramble resolution floor). Injections run through
    the SAME machinery as the real search: barycentric delay fixing (``bary_fn``), the
    DM-consistency cut (``dm``/``dm_tol``) with image-DM fit scatter (``dm_sigma``), and the
    transit-window selection. The map exposes the transit-geometry sensitivity structure
    honestly (cells away from integer sidereal-day delays are dark: the trailing images
    transit outside the beam).
    """
    if detection_p is None:
        detection_p = 2.0 / (n_scramble + 1)
    sens = np.zeros((delays.size, mag_ratios.size), bool)
    n_img = np.zeros_like(sens, dtype=int)
    for a, delta in enumerate(delays):
        for b, ratio in enumerate(mag_ratios):
            inj = inject_lensed_train(
                mjds,
                fluence,
                delay=float(delta),
                mag_ratio=float(ratio),
                bary_fn=bary_fn,
                dm=dm,
                dm_sigma=dm_sigma,
                seed=seed,
            )
            n_img[a, b] = inj["n_images_detected"]
            if not inj["detectable"]:
                continue
            offs = inj.get("bary_offsets")
            tb = inj["mjd"] + offs if offs is not None else inj["mjd"]
            found = recurring_delay_search(tb, dm=inj.get("dm"), tol_days=tol_days, dm_tol=dm_tol)
            if abs(found["best_delay"] - delta) > tol_days or found["m_max"] < 2:
                continue
            fap = scramble_fap(
                inj["mjd"],
                found["m_max"],
                dm=inj.get("dm"),
                bary_offsets=offs,
                tol_days=tol_days,
                dm_tol=dm_tol,
                n_scramble=n_scramble,
                seed=seed,
            )
            sens[a, b] = fap["p_value"] <= detection_p
    return {
        "delays": delays,
        "mag_ratios": mag_ratios,
        "sensitive": sens,
        "n_images": n_img,
        "sensitive_fraction": float(sens.mean()),
    }


def lensed_fraction_limit(n_searched: int, n_detected: int = 0, *, cl: float = 0.95) -> float:
    """Upper limit on the lensed-repeater fraction from ``n_detected`` among ``n_searched``.

    Poisson upper limit on the mean count (3.0 for 0 detected at 95%), divided by the number of
    repeaters searched. Quoted per searched repeater within the sensitivity region mapped by
    injections --- NOT an absolute optical-depth statement (the paper states the scope).
    """
    from scipy import stats

    if n_detected == 0:
        mu_up = -np.log(1.0 - cl)
    else:
        mu_up = 0.5 * stats.chi2.ppf(cl, 2 * (n_detected + 1))
    return float(mu_up / max(n_searched, 1))


def run(out: str = ".", *, offline: bool = True, n_scramble: int = 100) -> dict:
    """Offline: injected lensed train recovered + clean null control; real: all-repeater search."""
    import json

    detection_p = 2.0 / (n_scramble + 1)
    if offline:
        from .frbwait import synthetic_repeater_set

        base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
        t0 = base["mjd"]
        f0 = np.full(t0.size, 10.0)
        # injected: a lensed twin at 3 sidereal days, equal magnification
        delay = 3.0 * SIDEREAL_DAY
        inj = inject_lensed_train(t0, f0, delay=delay, mag_ratio=1.0)
        found = recurring_delay_search(inj["mjd"])
        fap = scramble_fap(inj["mjd"], found["m_max"], n_scramble=n_scramble)
        # control: the un-lensed train must be null
        found0 = recurring_delay_search(t0)
        fap0 = scramble_fap(t0, found0["m_max"], n_scramble=n_scramble, seed=1)
        rows = [
            {
                "name": "SYN-LENSED",
                "n_bursts": int(inj["mjd"].size),
                "m_max": found["m_max"],
                "best_delay": found["best_delay"],
                "p_value": fap["p_value"],
            },
            {
                "name": "SYN-CONTROL",
                "n_bursts": int(t0.size),
                "m_max": found0["m_max"],
                "best_delay": found0["best_delay"],
                "p_value": fap0["p_value"],
            },
        ]
        source = "synthetic transit-sampled trains"
        extra = {
            "true_delay": delay,
            "recovered_delay_err_s": round(abs(found["best_delay"] - delay) * 86400.0, 1),
        }
        n_searched, detections = 1, [r for r in rows[:1] if r["p_value"] <= detection_p]
        sens = sensitivity_map(
            t0,
            f0,
            delays=np.array([2.0 * SIDEREAL_DAY, 3.0 * SIDEREAL_DAY, 10.0, 45.7]),
            mag_ratios=np.array([0.2, 1.0]),
            n_scramble=max(20, n_scramble // 5),
        )
    else:  # pragma: no cover - needs the local Cat-2 mirror
        cat = load_catalog2()
        trains = repeater_trains(cat, min_bursts=5)
        rows = []
        for j, (name, tr) in enumerate(sorted(trains.items())):
            good = np.isfinite(tr["mjd"])
            t, dm_arr = tr["mjd"][good], tr["dm"][good]
            if t.size < 5 or (t.max() - t.min()) < 2 * MIN_DELAY_DAYS:
                continue
            # the search domain is BARYCENTRIC: a lens delay is only fixed at the barycenter
            # (mjd_400 is topocentric; the annual Roemer term drifts topocentric pair delays
            # by up to ~+-150 s -- far beyond the 5-s tolerance)
            bary_fn = barycentric_offset_fn(tr["ra"], tr["dec"])
            ltt = bary_fn(t)
            # DM tolerance from the source's own fitted-DM scatter (the fitburst scatter is
            # structure-driven, ~2-7 pc/cc for active repeaters -- a 1 pc/cc cut would reject
            # most GENUINE image pairs); pair ddm ~ N(0, sqrt(2) sigma) -> 3 sigma_pair keeps
            # ~99.7% of real pairs
            dm_sigma = float(1.4826 * np.nanmedian(np.abs(dm_arr - np.nanmedian(dm_arr))))
            dm_tol_src = max(1.0, 3.0 * np.sqrt(2.0) * dm_sigma)
            found = recurring_delay_search(t + ltt, dm=dm_arr, dm_tol=dm_tol_src)
            fap = scramble_fap(
                t,
                found["m_max"],
                dm=dm_arr,
                bary_offsets=ltt,
                dm_tol=dm_tol_src,
                n_scramble=n_scramble,
                seed=j,
            )
            rows.append(
                {
                    "name": name,
                    "n_bursts": int(t.size),
                    "m_max": found["m_max"],
                    "best_delay": round(found["best_delay"], 4),
                    "p_value": fap["p_value"],
                    "dm_tol": round(dm_tol_src, 2),
                }
            )
        source = "CHIME/FRB Catalog 2 repeaters (CANFAR DOI 10.11570/25.0066)"
        detections = [r for r in rows if r["p_value"] <= detection_p]
        n_searched = len(rows)
        # sensitivity map on the highest-count repeater (the census's deepest train) --
        # through the SAME machinery: barycentring, DM cut with image fit scatter, transit
        # selection. The other searched trains are shallower and necessarily less sensitive.
        big = max(trains.items(), key=lambda kv: kv[1]["n_bursts"])
        big_good = np.isfinite(big[1]["mjd"])
        big_dm = big[1]["dm"][big_good]
        big_sigma = float(1.4826 * np.nanmedian(np.abs(big_dm - np.nanmedian(big_dm))))
        sens = sensitivity_map(
            big[1]["mjd"][big_good],
            big[1]["fluence"][big_good],
            delays=np.concatenate([np.arange(2, 30, 4) * SIDEREAL_DAY, np.array([10.3, 33.7])]),
            mag_ratios=np.array([0.1, 0.3, 1.0]),
            n_scramble=max(20, n_scramble // 5),
            bary_fn=barycentric_offset_fn(big[1]["ra"], big[1]["dec"]),
            dm=big_dm,
            dm_sigma=big_sigma,
            dm_tol=max(1.0, 3.0 * np.sqrt(2.0) * big_sigma),
        )
        extra = {"sensitivity_source": big[0]}

    limit = lensed_fraction_limit(n_searched, len(detections))
    metrics = {
        "source": source,
        "is_real": not offline,
        "n_searched": n_searched,
        "n_detections": len(detections),
        "detection_names": [r["name"] for r in detections],
        "lensed_fraction_limit_95": round(limit, 4),
        "sensitive_fraction": round(sens["sensitive_fraction"], 3),
        "rows": rows,
        **extra,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .frbwait import _json_safe

    (op / "results" / "frblens_metrics.json").write_text(
        json.dumps(_json_safe(metrics), indent=2) + "\n"
    )
    _figure(rows, sens, op / "papers" / "frblens" / "figures")
    _write_macros(metrics, op / "papers" / "frblens" / "generated" / "macros.tex")
    return metrics


def _figure(rows: list[dict], sens: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ps = [r["p_value"] for r in rows]
    ns = [r["n_bursts"] for r in rows]
    ax1.scatter(ns, ps, s=16)
    ax1.set(
        xscale="log",
        yscale="log",
        xlabel="bursts",
        ylabel="M$_{max}$ p (sidereal scramble)",
        title="Recurring-delay search",
    )
    im = ax2.pcolormesh(
        sens["mag_ratios"],
        sens["delays"],
        sens["sensitive"].astype(float),
        cmap="Greys",
        vmin=0,
        vmax=1,
        shading="nearest",
    )
    fig.colorbar(im, ax=ax2, label="sensitive")
    ax2.set(
        xlabel="magnification ratio",
        ylabel="delay (days)",
        title="Injection sensitivity",
    )
    fig.tight_layout()
    fig.savefig(out / "frblens.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "flReal" if m.get("is_real") else "flSyn"
    lines = [
        "% Auto-generated by jansky_research.frblens._write_macros -- do not edit.",
        "% Synthetic (flSyn*) and real (flReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under flReal* (an offline rebuild resets flReal* to placeholders by design).",
        rf"\newcommand{{\flSource}}{{{m['source']}}}",
    ]
    keys = (
        ("NSearched", "n_searched"),
        ("NDet", "n_detections"),
        ("Limit", "lensed_fraction_limit_95"),
        ("SensFrac", "sensitive_fraction"),
    )
    for ns in ("flSyn", "flReal"):
        live = ns == pref
        for macro, key in keys:
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    from .frbwait import _json_safe

    p = argparse.ArgumentParser(description="Lensed-repeater delay-pattern search on Cat 2.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--n-scramble", type=int, default=100)
    args = p.parse_args(argv)
    m = run(args.out, offline=args.offline, n_scramble=args.n_scramble)
    m["rows"] = f"[{len(m['rows'])} rows in results/frblens_metrics.json]"
    print(json.dumps(_json_safe(m), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
