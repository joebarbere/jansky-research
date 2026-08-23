"""Independent reproduction of the Breakthrough Listen 3I/ATLAS GBT technosignature search.

On 2025-12-18 Breakthrough Listen observed the interstellar object 3I/ATLAS with the GBT at
1-12 GHz (four receivers, ABACAD on/off cadences), searched with turboSETI at +-4 Hz/s, and
reported a nondetection down to ~100 mW EIRP (Jacobson-Bell et al. 2025, RNAAS,
arXiv:2512.19763). The data are public (https://bldata.berkeley.edu/ATLAS/GB_ATLAS/), but no
independent reanalysis of the released files exists — the other 3I/ATLAS searches (ATA
arXiv:2512.18142, FAST arXiv:2603.19023) are teams analysing their own observations.

This module is an independent pipeline over those files: a physical-units (Hz/s) de-Doppler
search with per-channel bandpass normalisation and DC-spike excision — the two real-data
effects that defeated the teaching-grade detector in the merged ``driftsearch`` slice (plan 11
honestly reported that ``jansky.seti.drift_search`` recovers injected tones but *not* the real
Voyager-1 carrier). The ABACAD on/off filter then rejects anything not localised on the sky at
the target. Everything scientific is pure NumPy and offline-testable on a synthetic cadence;
file I/O and downloads are thin, network-marked wrappers.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np

__all__ = [
    "CADENCE_T0_SEC",
    "GB_ATLAS_BASE",
    "Hit",
    "L_BAND_NODES",
    "Scan",
    "cadence",
    "dedoppler",
    "eirp_limit_w",
    "find_hits",
    "normalize",
    "onoff_filter",
    "parse_index",
    "run",
    "synthetic_cadence",
]

GB_ATLAS_BASE = "https://bldata.berkeley.edu/ATLAS/GB_ATLAS/"
MJD = 61027  # 2025-12-18

# Cadence start (UT second of day) per receiver, pinned 2026-07-30 by remote header reads of
# the .0002 products (fsspec+h5py HTTP range requests; see plans/85).
CADENCE_T0_SEC = {"L": 17226, "S": 21817, "C": 26882, "X": 31308}

# Band -> node -> (f_lo, f_hi) in MHz, pinned from the same remote header reads. Each node
# spans 187.5 MHz. For C and X the recording used four banks with pairwise-duplicated coverage
# at bank boundaries (e.g. blc26 == blc30 in C); the duplicates are OMITTED here so a sweep
# searches each frequency chunk once — the skipped redundant nodes are listed per band in
# `DUPLICATE_NODES` (not a coverage gap: identical spectra recorded twice).
BAND_NODES: dict[str, dict[str, tuple[float, float]]] = {
    "L": {
        "blc21": (1876.5, 2064.0),
        "blc22": (1689.0, 1876.5),
        "blc23": (1501.5, 1689.0),
        "blc24": (1314.0, 1501.5),
        "blc25": (1126.5, 1314.0),
        "blc26": (939.0, 1126.5),
    },
    "S": {
        "blc22": (2589.0, 2776.5),
        "blc23": (2401.5, 2589.0),
        "blc24": (2214.0, 2401.5),
        "blc25": (2026.5, 2214.0),
        "blc26": (1839.0, 2026.5),
        "blc27": (1651.5, 1839.0),
    },
    "C": {
        "blc21": (8064.0, 8251.5),
        "blc22": (7876.5, 8064.0),
        "blc23": (7689.0, 7876.5),
        "blc24": (7501.5, 7689.0),
        "blc25": (7314.0, 7501.5),
        "blc26": (7126.5, 7314.0),
        "blc27": (6939.0, 7126.5),
        "blc32": (6751.5, 6939.0),
        "blc33": (6564.0, 6751.5),
        "blc34": (6376.5, 6564.0),
        "blc35": (6189.0, 6376.5),
        "blc36": (6001.5, 6189.0),
        "blc37": (5814.0, 6001.5),
        "blc62": (5626.5, 5814.0),
        "blc63": (5439.0, 5626.5),
        "blc64": (5251.5, 5439.0),
        "blc65": (5064.0, 5251.5),
        "blc66": (4876.5, 5064.0),
        "blc67": (4689.0, 4876.5),
        "blc72": (4501.5, 4689.0),
        "blc73": (4314.0, 4501.5),
        "blc74": (4126.5, 4314.0),
        "blc75": (3939.0, 4126.5),
    },
    "X": {
        "blc20": (12189.0, 12376.5),
        "blc21": (12001.5, 12189.0),
        "blc22": (11814.0, 12001.5),
        "blc23": (11626.5, 11814.0),
        "blc24": (11439.0, 11626.5),
        "blc25": (11251.5, 11439.0),
        "blc26": (11064.0, 11251.5),
        "blc27": (10876.5, 11064.0),
        "blc32": (10689.0, 10876.5),
        "blc33": (10501.5, 10689.0),
        "blc34": (10314.0, 10501.5),
        "blc35": (10126.5, 10314.0),
        "blc36": (9939.0, 10126.5),
        "blc37": (9751.5, 9939.0),
        "blc62": (9564.0, 9751.5),
        "blc63": (9376.5, 9564.0),
        "blc64": (9189.0, 9376.5),
        "blc65": (9001.5, 9189.0),
        "blc66": (8814.0, 9001.5),
        "blc67": (8626.5, 8814.0),
        "blc72": (8439.0, 8626.5),
        "blc73": (8251.5, 8439.0),
        "blc74": (8064.0, 8251.5),
        "blc75": (7876.5, 8064.0),
        "blc76": (7689.0, 7876.5),
    },
}
# Redundant recordings of chunks already covered above (identical band per the pinned headers).
DUPLICATE_NODES = {
    "C": {
        "blc30": "blc26",
        "blc31": "blc27",
        "blc60": "blc36",
        "blc61": "blc37",
        "blc70": "blc66",
        "blc71": "blc67",
    },
    "X": {
        "blc30": "blc26",
        "blc31": "blc27",
        "blc60": "blc36",
        "blc61": "blc37",
        "blc70": "blc66",
        "blc71": "blc67",
    },
}
L_BAND_NODES = BAND_NODES["L"]  # backward-compatible alias

AU_M = 1.495978707e11
# Geocentric distance of 3I/ATLAS at the L-cadence epoch, pinned 2026-07-30 from JPL Horizons
# (DES=C/2025 N1, center 500@399, 2025-12-18 05:00 UT: delta = 1.79801 au, rdot = -1.52 km/s).
DISTANCE_AU_DEFAULT = 1.798


@dataclass(frozen=True)
class Scan:
    """One rawspec product file in the GB_ATLAS archive listing."""

    node: str
    sec: int  # UT second of day (scan start)
    on: bool  # True = on-target, False = the _OFF reference position
    scan_id: str
    product: str  # "0000" (fine-frequency) / "0001" / "0002" (mid-resolution)
    filename: str


_FILE_RE = re.compile(
    r"(blc\d+)_guppi_(\d+)_(\d+)_DIAG_3I_ATLAS(_OFF)?_(\d+)\.rawspec\.(\d{4})\.h5"
)


def parse_index(text: str) -> list[Scan]:
    """Parse the GB_ATLAS HTTP index (or any text containing the filenames) into scans."""
    seen: dict[str, Scan] = {}
    for m in _FILE_RE.finditer(text):
        fn = m.group(0)
        seen[fn] = Scan(
            node=m.group(1),
            sec=int(m.group(3)),
            on=m.group(4) is None,
            scan_id=m.group(5),
            product=m.group(6),
            filename=fn,
        )
    return sorted(seen.values(), key=lambda s: (s.node, s.sec, s.product))


def cadence(scans: list[Scan], node: str, band: str = "L", product: str = "0000") -> list[Scan]:
    """Select one node's 6-scan ABACAD cadence for ``band``, validating the ON/OFF pattern.

    Scans are grouped by proximity to the pinned cadence start; a gap of more than 600 s to the
    block start ends the block (the next receiver's cadence starts ~1 h later).
    """
    t0 = CADENCE_T0_SEC[band]
    blk = sorted(
        (s for s in scans if s.node == node and s.product == product and 0 <= s.sec - t0 < 2400),
        key=lambda s: s.sec,
    )
    if len(blk) != 6:
        raise ValueError(f"expected 6 scans for {node} band {band}, found {len(blk)}")
    if [s.on for s in blk] != [True, False, True, False, True, False]:
        raise ValueError(f"scans for {node} band {band} are not an ABACAD on/off cadence")
    return blk


def normalize(wf: np.ndarray, *, clip: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Per-channel bandpass normalisation plus excision of grossly deviant channels.

    Each channel is divided by its median over time (flattening the bandpass), then channels
    whose median deviates from the band median by more than ``clip`` robust sigma (MAD) are
    masked to 1.0 — this removes the rawspec band-centre DC spike that plan 11 identified as
    the dominant artifact in real BL files. Returns ``(normalised, excised_mask)``.
    """
    wf = np.asarray(wf, float)
    med = np.median(wf, axis=0)
    band_med = np.median(med)
    mad = np.median(np.abs(med - band_med)) * 1.4826
    bad = np.abs(med - band_med) > clip * max(mad, 1e-300)
    safe = np.where(med > 0, med, 1.0)
    out = wf / safe
    out[:, bad] = 1.0
    return out, bad


def dedoppler(
    wf: np.ndarray, *, tsamp_s: float, foff_hz: float, drift_rates: np.ndarray
) -> np.ndarray:
    """Brute-force shift-and-sum de-Doppler in physical units.

    For each drift rate (Hz/s) every time row is shifted by ``round(rate * t / foff_hz)``
    channels (``foff_hz`` is signed — negative in BL files, where channel 0 is the highest
    frequency) and the rows are averaged. Returns shape ``(n_drift, n_freq)``; each row is the
    drift-aligned time-averaged spectrum, indexed by start channel.
    """
    wf = np.asarray(wf, float)
    n_time, n_freq = wf.shape
    t = np.arange(n_time) * tsamp_s
    out = np.zeros((len(drift_rates), n_freq))
    for i, rate in enumerate(np.asarray(drift_rates, float)):
        shifts = np.rint(rate * t / foff_hz).astype(int)
        acc = np.zeros(n_freq)
        for k in range(n_time):
            s = shifts[k]
            if s == 0:
                acc += wf[k]
            elif s > 0:
                acc[: n_freq - s] += wf[k, s:]
            else:
                acc[-s:] += wf[k, : n_freq + s]
        out[i] = acc / n_time
    return out


@dataclass(frozen=True)
class Hit:
    """One de-Doppler detection: start channel, drift rate, and integrated S/N."""

    chan: int
    drift_hz_s: float
    snr: float


def find_hits(
    wf: np.ndarray,
    *,
    tsamp_s: float,
    foff_hz: float,
    drift_rates: np.ndarray,
    threshold: float = 10.0,
    min_sep: int = 16,
) -> list[Hit]:
    """De-Doppler ``wf`` and return above-threshold peaks, strongest first.

    S/N is computed per drift row against that row's robust (median/MAD) noise. Peaks within
    ``min_sep`` channels of a stronger peak (at any drift) are suppressed, so one tone yields
    one hit.
    """
    dd = dedoppler(wf, tsamp_s=tsamp_s, foff_hz=foff_hz, drift_rates=drift_rates)
    med = np.median(dd, axis=1, keepdims=True)
    mad = np.median(np.abs(dd - med), axis=1, keepdims=True) * 1.4826
    snr = (dd - med) / np.maximum(mad, 1e-300)
    drift_rates = np.asarray(drift_rates, float)
    hits: list[Hit] = []
    for i, j in np.argwhere(snr > threshold):
        hits.append(Hit(chan=int(j), drift_hz_s=float(drift_rates[i]), snr=float(snr[i, j])))
    hits.sort(key=lambda h: -h.snr)
    kept: list[Hit] = []
    for h in hits:
        if all(abs(h.chan - k.chan) >= min_sep for k in kept):
            kept.append(h)
    return kept


def onoff_filter(
    hits_per_scan: list[list[Hit]],
    scan_secs: list[int],
    scan_on: list[bool],
    *,
    foff_hz: float,
    tol_chan: int = 32,
) -> list[Hit]:
    """ABACAD sky-localisation filter: keep hits present in every ON and absent in every OFF.

    A first-ON hit is extrapolated to each later scan along its own drift rate
    (``chan + rate * dt / foff_hz``); it survives only if a hit lies within ``tol_chan``
    channels of the prediction in every ON scan and in no OFF scan. This is what separates a
    sky signal at the target from local RFI (present in the OFFs too).
    """
    if not (len(hits_per_scan) == len(scan_secs) == len(scan_on)):
        raise ValueError("hits_per_scan, scan_secs and scan_on must have equal length")
    if not scan_on[0]:
        raise ValueError("the cadence must start with an ON scan")
    t0 = scan_secs[0]
    survivors = []
    for h in hits_per_scan[0]:
        ok = True
        for hits, sec, on in zip(hits_per_scan[1:], scan_secs[1:], scan_on[1:], strict=True):
            pred = h.chan + h.drift_hz_s * (sec - t0) / foff_hz
            matched = any(abs(g.chan - pred) <= tol_chan for g in hits)
            if on != matched:
                ok = False
                break
        if ok:
            survivors.append(h)
    return survivors


def synthetic_cadence(
    *,
    n_time: int = 16,
    n_freq: int = 4096,
    tsamp_s: float = 18.6,
    foff_hz: float = -2.79,
    tone_chan: int = 1024,
    tone_drift_hz_s: float = 0.4,
    tone_snr: float = 25.0,
    rfi_chan: int | None = 3000,
    dc_spike: bool = True,
    seed: int = 0,
) -> tuple[list[np.ndarray], list[int], list[bool]]:
    """Six synthetic scans mimicking one node's ABACAD cadence, with known contaminants.

    The drifting tone (the "signal") appears only in the ON scans and drifts coherently across
    the whole cadence; the RFI tone is a zero-drift carrier present in all six scans; the DC
    spike is a huge static spike at the band centre, as in real rawspec files. Returns
    ``(waterfalls, scan_secs, scan_on)`` ready for :func:`find_hits` + :func:`onoff_filter`.
    """
    rng = np.random.default_rng(seed)
    secs = [17226 + 318 * k for k in range(6)]
    on = [True, False, True, False, True, False]
    amp = tone_snr / math.sqrt(n_time)
    t_in_scan = np.arange(n_time) * tsamp_s
    wfs = []
    for sec, is_on in zip(secs, on, strict=True):
        wf = rng.normal(loc=10.0, scale=1.0, size=(n_time, n_freq))
        if is_on:
            dt = (sec - secs[0]) + t_in_scan
            chans = np.rint(tone_chan + tone_drift_hz_s * dt / foff_hz).astype(int)
            wf[np.arange(n_time), chans] += amp
        if rfi_chan is not None:
            wf[:, rfi_chan] += amp
        if dc_spike:
            wf[:, n_freq // 2] += 1e3
        wfs.append(wf)
    return wfs, secs, on


def eirp_limit_w(
    *,
    snr: float = 10.0,
    sefd_jy: float = 10.0,
    t_obs_s: float = 300.0,
    chan_hz: float = 2.79,
    distance_au: float = DISTANCE_AU_DEFAULT,
    n_pol: int = 2,
) -> float:
    """Minimum detectable EIRP (W) for a narrowband tone, EIRP = 4 pi d^2 * S_min * dnu.

    ``S_min = snr * SEFD / sqrt(n_pol * dnu * t)`` is the radiometer narrowband limit for one
    fine channel. With the paper's GBT L-band numbers (SEFD ~10 Jy, 5-min scans, ~3 Hz
    channels, d = 1.8 au) this reproduces their ~100 mW headline figure.
    """
    s_min_jy = snr * sefd_jy / math.sqrt(n_pol * chan_hz * t_obs_s)
    d_m = distance_au * AU_M
    return 4.0 * math.pi * d_m**2 * s_min_jy * 1e-26 * chan_hz


def run(out: str = ".", *, threshold: float = 10.0, seed: int = 0) -> dict:
    """Offline synthetic round-trip: recover the injected tone, reject RFI and the DC spike.

    Builds the synthetic cadence, runs normalisation + de-Doppler + the on/off filter, and
    writes metrics, a two-panel figure, and the paper macros. This is the tested, reproducible
    core; the real-data leg (``--node``, network) drives exactly the same functions.
    """
    import json
    from pathlib import Path

    tsamp, foff = 18.6, -2.79
    drifts = np.linspace(-1.0, 1.0, 41)
    wfs, secs, on = synthetic_cadence(tsamp_s=tsamp, foff_hz=foff, seed=seed)
    hits_per_scan = []
    dc_excised = []
    for wf in wfs:
        wfn, bad = normalize(wf)
        dc_excised.append(bool(bad[wf.shape[1] // 2]))
        hits_per_scan.append(
            find_hits(wfn, tsamp_s=tsamp, foff_hz=foff, drift_rates=drifts, threshold=threshold)
        )
    survivors = onoff_filter(hits_per_scan, secs, on, foff_hz=foff)
    recovered = any(abs(h.chan - 1024) < 32 for h in survivors)
    rfi_rejected = not any(abs(h.chan - 3000) < 32 for h in survivors)
    metrics = {
        "threshold": threshold,
        "n_hits_first_on": len(hits_per_scan[0]),
        "n_survivors": len(survivors),
        "recovered_injected": recovered,
        "rfi_rejected": rfi_rejected,
        "dc_spike_excised": all(dc_excised),
        "survivor_drift_hz_s": survivors[0].drift_hz_s if survivors else None,
        "eirp_limit_w_paper_params": eirp_limit_w(snr=16.0),
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "atlas3i_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    paper = op / "papers" / "atlas3i"
    _figure(wfs[0], tsamp, foff, drifts, paper / "figures")
    _write_macros(metrics, paper / "generated" / "macros_offline.tex")
    return metrics


def _figure(wf: np.ndarray, tsamp: float, foff: float, drifts: np.ndarray, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    wfn, _ = normalize(wf)
    dd = dedoppler(wfn, tsamp_s=tsamp, foff_hz=foff, drift_rates=drifts)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.2))
    ax1.imshow(wfn, origin="lower", aspect="auto", cmap="viridis")
    ax1.set(xlabel="channel", ylabel="time sample", title="Synthetic ON scan (normalised)")
    ax2.imshow(
        dd,
        origin="lower",
        aspect="auto",
        cmap="viridis",
        extent=[0, wf.shape[1], drifts.min(), drifts.max()],
    )
    ax2.set(xlabel="start channel", ylabel="drift rate (Hz/s)", title="De-Doppler plane")
    fig.tight_layout()
    fig.savefig(out / "atlas3i_dedoppler.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    """Emit LaTeX ``\\newcommand`` macros so the paper hard-codes no number."""
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.atlas3i._write_macros — do not edit by hand.",
        "% Synthetic-round-trip macros; the real-sweep macros come from sweep_macros.",
        rf"\newcommand{{\aiSynThreshold}}{{{m['threshold']:.0f}}}",
        rf"\newcommand{{\aiSynSurvivors}}{{{m['n_survivors']}}}",
        rf"\newcommand{{\aiSynEirpMw}}{{{m['eirp_limit_w_paper_params'] * 1e3:.0f}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a run may only ADD information, so an
    # offline rebuild can never blank a real value (report.preserve_live_macros).
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def vet_stamps(
    stamps: list[np.ndarray],
    secs: list[int],
    on: list[bool],
    hit: Hit,
    *,
    tsamp_s: float,
    foff_hz: float,
    stamp_center_chan: int,
    threshold: float = 8.0,
    tol_chan: int = 32,
) -> dict:
    """Vet one on/off survivor against small waterfall stamps from all six scans.

    The survivor's drift rate predicts where the tone must sit in every scan; each stamp is
    de-Dopplered at that rate and the S/N at the predicted channel is read off. A confirmed
    sky signal must track the prediction in every ON stamp (S/N above ``threshold`` within
    ``tol_chan``) and stay clean in every OFF — chance ON/OFF coincidences from dense RFI hit
    lists do not survive this coherence test. Zero-drift survivors are additionally flagged:
    a transmitter at the target must show nonzero barycentric drift, so an exactly-constant
    tone is terrestrial. Stamps are centred on ``stamp_center_chan`` (full-band channel).
    """
    if not (len(stamps) == len(secs) == len(on)):
        raise ValueError("stamps, secs and on must have equal length")
    t0 = secs[0]
    per_scan_snr = []
    for wf, sec in zip(stamps, secs, strict=True):
        wfn, _ = normalize(np.asarray(wf, float))
        dd = dedoppler(
            wfn, tsamp_s=tsamp_s, foff_hz=foff_hz, drift_rates=np.array([hit.drift_hz_s])
        )
        row = dd[0]
        med = np.median(row)
        mad = np.median(np.abs(row - med)) * 1.4826
        pred = (
            (hit.chan - stamp_center_chan)
            + wf.shape[1] // 2
            + hit.drift_hz_s * (sec - t0) / foff_hz
        )
        lo = max(0, int(pred) - tol_chan)
        hi = min(row.size, int(pred) + tol_chan + 1)
        peak = row[lo:hi].max() if hi > lo else med
        per_scan_snr.append(float((peak - med) / max(mad, 1e-300)))
    on_ok = all(s >= threshold for s, o in zip(per_scan_snr, on, strict=True) if o)
    off_clean = all(s < threshold / 2 for s, o in zip(per_scan_snr, on, strict=True) if not o)
    zero_drift = hit.drift_hz_s == 0.0
    return {
        "per_scan_snr": per_scan_snr,
        "tracks_in_ons": on_ok,
        "clean_in_offs": off_clean,
        "zero_drift": zero_drift,
        "confirmed": on_ok and off_clean and not zero_drift,
    }


def sweep_summary(results_dir: str = "results", band: str = "L") -> dict:
    """Aggregate one band's per-node sweep JSONs into the numbers the paper quotes."""
    import json
    from pathlib import Path

    rows = []
    for p in sorted(Path(results_dir).glob(f"atlas3i_blc*_{band}.json")):
        r = json.loads(p.read_text())
        rows.append(r)
    if not rows:
        raise FileNotFoundError(f"no atlas3i_blc*_{band}.json results in {results_dir}")
    return {
        "band": band,
        "nodes": [r["node"] for r in rows],
        "hits_per_scan": {r["node"]: r["n_hits_per_scan"] for r in rows},
        "survivors": {r["node"]: r["n_survivors"] for r in rows},
        "confirmed": {r["node"]: r["n_confirmed"] for r in rows},
        "total_hits": sum(sum(r["n_hits_per_scan"]) for r in rows),
        "total_survivors": sum(r["n_survivors"] for r in rows),
        "total_confirmed": sum(r["n_confirmed"] for r in rows),
        "sat_coherent": sum(
            1
            for r in rows
            for s in r["survivors"]
            if s.get("satellite_band")
            and s.get("tracks_in_ons")
            and s.get("clean_in_offs")
            and not s.get("zero_drift")
        ),
        "eirp_w": rows[0]["eirp_limit_w"],
        "distance_au": rows[0].get("distance_au", DISTANCE_AU_DEFAULT),
        "band_lo_mhz": min(lo for lo, _ in BAND_NODES[band].values()),
        "band_hi_mhz": max(hi for _, hi in BAND_NODES[band].values()),
    }


def sweep_macros(summary: dict, path) -> None:
    """Emit the real-sweep LaTeX macros so the RNAAS hard-codes no number."""
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.atlas3i.sweep_macros — do not edit by hand.",
        rf"\newcommand{{\aiNodes}}{{{len(summary['nodes'])}}}",
        rf"\newcommand{{\aiBandLo}}{{{summary['band_lo_mhz']:.0f}}}",
        rf"\newcommand{{\aiBandHi}}{{{summary['band_hi_mhz']:.0f}}}",
        rf"\newcommand{{\aiRawHits}}{{{summary['total_hits']:,}}}".replace(",", r"\,"),
        rf"\newcommand{{\aiSurvivors}}{{{summary['total_survivors']}}}",
        rf"\newcommand{{\aiConfirmed}}{{{summary['total_confirmed']}}}",
        rf"\newcommand{{\aiSatCoherent}}{{{summary['sat_coherent']}}}",
        rf"\newcommand{{\aiEirpMw}}{{{summary['eirp_w'] * 1e3:.1f}}}",
        rf"\newcommand{{\aiDistanceAu}}{{{summary['distance_au']:.3f}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a run may only ADD information, so an
    # offline rebuild can never blank a real value (report.preserve_live_macros).
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def sweep_figure(summary: dict, out_dir) -> None:
    """One-panel RNAAS figure: per-scan hit counts by node, with the vetting funnel."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 3.2))
    for i, node in enumerate(summary["nodes"]):
        hits = summary["hits_per_scan"][node]
        x = i + (np.arange(len(hits)) - (len(hits) - 1) / 2) * 0.11
        on = [k % 2 == 0 for k in range(len(hits))]
        ax.scatter(x, hits, c=["#1f77b4" if o else "#d62728" for o in on], s=22, zorder=3)
        lo, hi = BAND_NODES[summary.get("band", "L")][node]
        ax.annotate(
            f"{node}\n{lo:.0f}–{hi:.0f}",
            (i, 0.55),
            xycoords=("data", "axes fraction"),
            ha="center",
            fontsize=7,
        )
    ax.set_yscale("log")
    ax.set_xticks([])
    ax.set(
        ylabel="hits per 5-min scan (16$\\sigma$)",
        title=(
            f"GB_ATLAS {summary.get('band', 'L')} band: {summary['total_hits']:,} raw hits "
            f"$\\to$ {summary['total_survivors']} on/off survivors "
            f"$\\to$ {summary['total_confirmed']} confirmed"
        ),
    )
    ax.scatter([], [], c="#1f77b4", label="ON scans", s=22)
    ax.scatter([], [], c="#d62728", label="OFF scans", s=22)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "atlas3i_sweep.pdf")
    plt.close(fig)


def survey_summary(results_dir: str = "results") -> dict:
    """Aggregate every completed band's sweep into the full-survey numbers the paper quotes."""
    bands = {}
    for band in BAND_NODES:
        try:
            bands[band] = sweep_summary(results_dir, band=band)
        except FileNotFoundError:
            continue
    if not bands:
        raise FileNotFoundError(f"no per-band atlas3i results in {results_dir}")
    return {
        "bands": bands,
        "total_nodes": sum(len(b["nodes"]) for b in bands.values()),
        # Scans, not downloads. The paper previously hard-typed "366", which is the *transfer*
        # count from survey/atlas3i-findings.md (360 unique + the 6-scan blc25 L re-run kept as
        # atlas3i_blc25_L_pass1.json). Derived here so it cannot drift from the evidence again.
        "total_scans": sum(
            len(scans) for b in bands.values() for scans in b["hits_per_scan"].values()
        ),
        "total_hits": sum(b["total_hits"] for b in bands.values()),
        "total_survivors": sum(b["total_survivors"] for b in bands.values()),
        "total_confirmed": sum(b["total_confirmed"] for b in bands.values()),
        "sat_coherent": sum(b["sat_coherent"] for b in bands.values()),
        "eirp_w": bands[next(iter(bands))]["eirp_w"],
        "distance_au": bands[next(iter(bands))]["distance_au"],
    }


def survey_macros(survey: dict, path) -> None:
    """Emit full-survey + per-band LaTeX macros so the papers hard-code no number."""
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.atlas3i.survey_macros — do not edit by hand.",
        rf"\newcommand{{\aiTotNodes}}{{{survey['total_nodes']}}}",
        rf"\newcommand{{\aiTotScans}}{{{survey['total_scans']}}}",
        rf"\newcommand{{\aiTotBands}}{{{len(survey['bands'])}}}",
        rf"\newcommand{{\aiTotRawHits}}{{{survey['total_hits']:,}}}".replace(",", r"\,"),
        rf"\newcommand{{\aiTotSurvivors}}{{{survey['total_survivors']}}}",
        rf"\newcommand{{\aiTotConfirmed}}{{{survey['total_confirmed']}}}",
        rf"\newcommand{{\aiSatCoherent}}{{{survey['sat_coherent']}}}",
        rf"\newcommand{{\aiEirpMw}}{{{survey['eirp_w'] * 1e3:.1f}}}",
        rf"\newcommand{{\aiDistanceAu}}{{{survey['distance_au']:.3f}}}",
    ]
    for band, s in survey["bands"].items():
        lines += [
            rf"\newcommand{{\ai{band}Nodes}}{{{len(s['nodes'])}}}",
            rf"\newcommand{{\ai{band}BandLo}}{{{s['band_lo_mhz']:.0f}}}",
            rf"\newcommand{{\ai{band}BandHi}}{{{s['band_hi_mhz']:.0f}}}",
            rf"\newcommand{{\ai{band}RawHits}}{{{s['total_hits']:,}}}".replace(",", r"\,"),
            rf"\newcommand{{\ai{band}Survivors}}{{{s['total_survivors']}}}",
            rf"\newcommand{{\ai{band}Confirmed}}{{{s['total_confirmed']}}}",
        ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a run may only ADD information, so an
    # offline rebuild can never blank a real value (report.preserve_live_macros).
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def survey_figure(survey: dict, out_dir) -> None:
    """The full-survey figure: per-scan hit counts for all nodes, coloured by receiver band."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    colors = {"L": "#1f77b4", "S": "#2ca02c", "C": "#ff7f0e", "X": "#9467bd"}
    fig, ax = plt.subplots(figsize=(8.5, 3.4))
    i = 0
    for band, s in survey["bands"].items():
        first = True
        for node in s["nodes"]:
            hits = np.maximum(s["hits_per_scan"][node], 0.7)  # log axis: show zeros at floor
            x = i + (np.arange(len(hits)) - (len(hits) - 1) / 2) * 0.1
            ax.scatter(
                x,
                hits,
                c=colors.get(band, "0.4"),
                s=10,
                zorder=3,
                label=(
                    f"{band} ({s['band_lo_mhz']:.0f}–{s['band_hi_mhz']:.0f} MHz)" if first else None
                ),
            )
            first = False
            i += 1
    ax.set_yscale("log")
    ax.set_xticks([])
    ax.set_xlabel("compute node (grouped by receiver band, increasing frequency)")
    ax.set(
        ylabel="hits per 5-min scan (16$\\sigma$)",
        title=(
            f"GB_ATLAS 1–12 GHz: {survey['total_hits']:,} raw hits "
            f"$\\to$ {survey['total_survivors']} on/off survivors "
            f"$\\to$ {survey['total_confirmed']} confirmed"
        ),
    )
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "atlas3i_survey.pdf")
    plt.close(fig)


# L-band satellite downlink allocations (MHz): transmitters genuinely in the sky, so they
# defeat a two-position ON/OFF filter by design and must be flagged by frequency. Ranges per
# ITU allocations; the same lists underlie the exclusion masks in BL L-band searches.
SATELLITE_BANDS_MHZ = {
    "Inmarsat/MSS downlink": (1525.0, 1559.0),
    "GPS L1 / Galileo E1": (1559.0, 1591.0),
    "GLONASS L1": (1592.0, 1610.0),
    "Iridium downlink": (1616.0, 1626.5),
    "GPS L2 (+GLONASS L2)": (1215.0, 1254.0),
    "GPS L5 / Galileo E5": (1164.0, 1215.0),
    "S-band MSS downlink": (2170.0, 2200.0),
    "S-DARS broadcast": (2320.0, 2345.0),
    "Ku FSS/DBS downlink": (11700.0, 12750.0),
}


def classify_band(freq_mhz: float) -> str | None:
    """Name the known satellite downlink allocation containing ``freq_mhz``, else ``None``."""
    for name, (lo, hi) in SATELLITE_BANDS_MHZ.items():
        if lo <= freq_mhz < hi:
            return name
    return None


# ----------------------------------------------------------------------------- real-data leg


class _HttpRangeFile:  # pragma: no cover - network
    """Minimal read-only file-like over HTTP Range requests, for remote h5py access.

    Reads are served from 256 KB-aligned cached blocks so h5py's many small metadata reads
    do not each cost a round trip. Dependency-free (urllib only).
    """

    BLOCK = 256 * 1024

    def __init__(self, url: str):
        from urllib.request import Request, urlopen

        self.url = url
        self.pos = 0
        self._cache: dict[int, bytes] = {}
        req = Request(url, headers={"Range": "bytes=0-0"})
        with urlopen(req) as r:
            rng = r.headers.get("Content-Range", "")
        self.size = int(rng.rsplit("/", 1)[-1])

    def _block(self, idx: int) -> bytes:
        if idx not in self._cache:
            from urllib.request import Request, urlopen

            start = idx * self.BLOCK
            end = min(start + self.BLOCK, self.size) - 1
            req = Request(self.url, headers={"Range": f"bytes={start}-{end}"})
            with urlopen(req) as r:
                self._cache[idx] = r.read()
        return self._cache[idx]

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self.pos
        out = bytearray()
        while n > 0 and self.pos < self.size:
            idx, off = divmod(self.pos, self.BLOCK)
            chunk = self._block(idx)[off : off + n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def seek(self, pos: int, whence: int = 0) -> int:
        self.pos = {0: pos, 1: self.pos + pos, 2: self.size + pos}[whence]
        return self.pos

    def tell(self) -> int:
        return self.pos

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def writable(self) -> bool:
        return False


def fetch_stamp(  # pragma: no cover - network
    scan: Scan, center_chan: int, half_width: int = 2048
) -> tuple[np.ndarray, dict]:
    """Remote-read a small frequency window around ``center_chan`` from one archive file."""
    import h5py
    import hdf5plugin  # noqa: F401 - registers the bitshuffle filter

    with h5py.File(_HttpRangeFile(GB_ATLAS_BASE + scan.filename), "r") as f:
        d = f["data"]
        hdr = {
            k: (v.item() if hasattr(v, "item") and getattr(v, "size", 1) == 1 else v)
            for k, v in d.attrs.items()
        }
        n_freq = d.shape[-1]
        lo = max(0, center_chan - half_width)
        hi = min(n_freq, center_chan + half_width)
        block = d[:, 0, lo:hi] if d.ndim == 3 else d[:, lo:hi]
        return np.asarray(block, float), hdr


def vet_node(  # pragma: no cover - network
    node: str, survivors: list[dict], *, band: str = "L", half_width: int = 2048
) -> list[dict]:
    """Remote-vet a node's on/off survivors: stamp reads + drift-coherence check per hit."""
    scans = cadence(fetch_index(), node, band=band)
    out = []
    for sv in survivors:
        hit = Hit(chan=int(sv["chan"]), drift_hz_s=float(sv["drift_hz_s"]), snr=float(sv["snr"]))
        stamps: list[np.ndarray] = []
        hdr: dict = {}
        for s in scans:
            wf, hdr = fetch_stamp(s, hit.chan, half_width)
            stamps.append(wf)
        verdict = vet_stamps(
            stamps,
            [s.sec for s in scans],
            [s.on for s in scans],
            hit,
            tsamp_s=float(hdr["tsamp"]),
            foff_hz=float(hdr["foff"]) * 1e6,
            stamp_center_chan=max(0, hit.chan - half_width) + half_width,
        )
        freq_mhz = float(hdr["fch1"]) + float(hdr["foff"]) * hit.chan
        out.append({**sv, "freq_mhz": freq_mhz, **verdict})
        print(f"[atlas3i] vet chan={hit.chan} f={freq_mhz:.3f} MHz: {verdict}", flush=True)
    return out


def fetch_index() -> list[Scan]:  # pragma: no cover - network
    """Download and parse the live GB_ATLAS archive index."""
    from urllib.request import urlopen

    with urlopen(GB_ATLAS_BASE) as r:
        return parse_index(r.read().decode("utf-8", "replace"))


def fetch_scan(scan: Scan, dest_dir: str) -> str:  # pragma: no cover - network, large
    """Download one scan file (resumable via a .part file); returns the local path."""
    from pathlib import Path
    from urllib.request import Request, urlopen

    dest = Path(dest_dir) / scan.filename
    if dest.exists():
        return str(dest)
    part = dest.with_suffix(dest.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    req = Request(GB_ATLAS_BASE + scan.filename)
    if offset:
        req.add_header("Range", f"bytes={offset}-")
    with urlopen(req) as r, open(part, "ab" if offset else "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    part.rename(dest)
    return str(dest)


def read_scan(path: str):  # pragma: no cover - optional voyager extra (h5py)
    """Read a rawspec .h5 product -> ``(waterfall, header)`` (requires the ``voyager`` extra)."""
    import h5py
    import hdf5plugin  # noqa: F401 - registers the bitshuffle filter

    with h5py.File(path, "r") as f:
        d = f["data"]
        hdr = {
            k: (v.item() if hasattr(v, "item") and getattr(v, "size", 1) == 1 else v)
            for k, v in d.attrs.items()
        }
        return np.asarray(d[:]).squeeze().astype(float), hdr


def search_node(  # pragma: no cover - network + large data
    node: str,
    *,
    band: str = "L",
    dest_dir: str = "data/atlas3i",
    threshold: float = 16.0,
    threshold_off: float = 10.0,
    max_drift_hz_s: float = 4.0,
    n_drift: int = 129,
    chunk: int = 1 << 20,
    delete_after: bool = True,
) -> dict:
    """Full real-data pass for one node: fetch its cadence, search, on/off filter, report.

    Fine-frequency files are searched in frequency chunks (~1M channels, overlapping by the
    maximum drift excursion) to bound memory. With ``delete_after`` the 10 GB scan files are
    removed once searched, so peak disk stays at one cadence (~60 GB).
    """
    import os

    import h5py
    import hdf5plugin  # noqa: F401 - registers the bitshuffle filter

    drifts = np.linspace(-max_drift_hz_s, max_drift_hz_s, n_drift)
    scans = cadence(fetch_index(), node, band=band)
    hits_per_scan: list[list[Hit]] = []
    hdr: dict = {}
    paths = []
    for s in scans:
        path = fetch_scan(s, dest_dir)
        paths.append(path)
        # Chunked lazy reads: a full fine-res waterfall is ~10 GB on disk and would double
        # twice in RAM as float64; slicing the dataset per frequency chunk bounds memory.
        hits: list[Hit] = []
        with h5py.File(path, "r") as f:
            d = f["data"]
            hdr = {
                k: (v.item() if hasattr(v, "item") and getattr(v, "size", 1) == 1 else v)
                for k, v in d.attrs.items()
            }
            tsamp = float(hdr["tsamp"])
            foff = float(hdr["foff"]) * 1e6  # MHz -> Hz
            n_freq = d.shape[-1]
            n_time = d.shape[0]
            pad = int(max_drift_hz_s * tsamp * n_time / abs(foff)) + 1
            for lo in range(0, n_freq, chunk):
                hi = min(n_freq, lo + chunk + pad)
                block = d[:, 0, lo:hi] if d.ndim == 3 else d[:, lo:hi]
                wfn, _ = normalize(np.asarray(block, float))
                # Paper convention (arXiv:2512.19763 §III): a lower OFF-scan threshold so a
                # candidate slightly weaker in the OFFs still registers there and is vetoed.
                thr = threshold if s.on else threshold_off
                hits += [
                    Hit(chan=h.chan + lo, drift_hz_s=h.drift_hz_s, snr=h.snr)
                    for h in find_hits(
                        wfn, tsamp_s=tsamp, foff_hz=foff, drift_rates=drifts, threshold=thr
                    )
                    if h.chan < (hi - lo) - (pad if hi < n_freq else 0)
                ]
        print(f"[atlas3i] {s.filename}: {len(hits)} hits", flush=True)
        hits_per_scan.append(hits)
    foff_hz = float(hdr["foff"]) * 1e6
    tsamp = float(hdr["tsamp"])
    survivors = onoff_filter(
        hits_per_scan,
        [s.sec for s in scans],
        [s.on for s in scans],
        foff_hz=foff_hz,
    )
    # Vet inline while the scan files are still on disk: drift-coherence stamps are cheap
    # locally, and the remote-stamp path proved slow/fragile (server throttling).
    half_width = 2048
    vetted = []
    for h in survivors:
        stamps = []
        for path in paths:
            with h5py.File(path, "r") as f:
                d = f["data"]
                lo = max(0, h.chan - half_width)
                hi = min(d.shape[-1], h.chan + half_width)
                block = d[:, 0, lo:hi] if d.ndim == 3 else d[:, lo:hi]
                stamps.append(np.asarray(block, float))
        verdict = vet_stamps(
            stamps,
            [s.sec for s in scans],
            [s.on for s in scans],
            h,
            tsamp_s=tsamp,
            foff_hz=foff_hz,
            stamp_center_chan=max(0, h.chan - half_width) + half_width,
        )
        freq_mhz = float(hdr["fch1"]) + float(hdr["foff"]) * h.chan
        sat = classify_band(freq_mhz)
        verdict["confirmed"] = verdict["confirmed"] and sat is None
        vetted.append({**vars(h), "freq_mhz": freq_mhz, "satellite_band": sat, **verdict})
        print(f"[atlas3i] vet chan={h.chan} f={freq_mhz:.3f} MHz sat={sat}: {verdict}", flush=True)
    if delete_after:
        for p in paths:
            os.remove(p)
    return {
        "node": node,
        "band": band,
        "n_hits_per_scan": [len(h) for h in hits_per_scan],
        "n_survivors": len(survivors),
        "n_confirmed": sum(bool(v["confirmed"]) for v in vetted),
        "survivors": vetted,
        "eirp_limit_w": eirp_limit_w(snr=threshold),
    }


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="BL 3I/ATLAS GBT reproduction (plan 85).")
    p.add_argument("--out", default=".")
    p.add_argument("--node", help="run the real-data leg on this node (e.g. blc25)")
    p.add_argument("--band", default="L")
    p.add_argument("--threshold", type=float, default=16.0, help="ON-scan S/N threshold")
    p.add_argument(
        "--threshold-off", type=float, default=10.0, help="OFF-scan S/N threshold (paper: 10)"
    )
    p.add_argument("--keep", action="store_true", help="keep downloaded scan files")
    p.add_argument("--vet", help="vet the survivors in this search-result JSON (remote stamps)")
    p.add_argument("--sweep", action="store_true", help="search+vet every node of --band serially")
    p.add_argument(
        "--paper", action="store_true", help="write the RNAAS figure + macros from results/"
    )
    args = p.parse_args(argv)
    if args.paper:
        from pathlib import Path

        paper = Path(args.out) / "papers" / "atlas3i"
        sv = survey_summary(str(Path(args.out) / "results"))
        survey_figure(sv, paper / "figures")
        survey_macros(sv, paper / "generated" / "macros.tex")
        # keep the single-band L figure for the RNAAS history / comparisons
        sweep_figure(sv["bands"]["L"], paper / "figures") if "L" in sv["bands"] else None
        print(json.dumps({k: v for k, v in sv.items() if k != "bands"}, indent=2))
        return 0
    if args.sweep:
        from pathlib import Path

        rp = Path(args.out) / "results"
        rp.mkdir(parents=True, exist_ok=True)
        for node in BAND_NODES[args.band]:
            dest = rp / f"atlas3i_{node}_{args.band}.json"
            if dest.exists():
                print(f"[atlas3i] {node}: already done, skipping", flush=True)
                continue
            r = search_node(
                node,
                band=args.band,
                threshold=args.threshold,
                threshold_off=args.threshold_off,
                delete_after=not args.keep,
            )
            dest.write_text(json.dumps(r, indent=2) + "\n")
            print(
                f"[atlas3i] {node}: {r['n_survivors']} survivors, {r['n_confirmed']} confirmed",
                flush=True,
            )
        return 0
    if args.vet:
        from pathlib import Path

        prev = json.loads(Path(args.vet).read_text())
        res = {**prev, "survivors": vet_node(prev["node"], prev["survivors"], band=prev["band"])}
        res["n_confirmed"] = sum(bool(s["confirmed"]) for s in res["survivors"])
        Path(args.vet.replace(".json", "_vetted.json")).write_text(json.dumps(res, indent=2) + "\n")
    elif args.node:
        from pathlib import Path

        res = search_node(
            args.node,
            band=args.band,
            threshold=args.threshold,
            threshold_off=args.threshold_off,
            delete_after=not args.keep,
        )
        rp = Path(args.out) / "results"
        rp.mkdir(parents=True, exist_ok=True)
        (rp / f"atlas3i_{args.node}_{args.band}.json").write_text(json.dumps(res, indent=2) + "\n")
    else:
        res = run(args.out)
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
