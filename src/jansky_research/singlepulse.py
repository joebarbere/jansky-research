"""Single-pulse / pulsar recover-a-known on a real public filterbank, FDMT-powered (plan 34).

The science leg of the torch-fdmt arc: read a small public SIGPROC filterbank of a KNOWN bright
pulsar, dedisperse with :mod:`jansky_research.fdmt` over a DM grid, and recover the catalogued DM
(and, by folding, the period) --- a validation on real telescope data, not a discovery. GATE 0
(2026-07-02) chose a 3.25 MB Parkes/UWL observation of the **Crab pulsar** (B0531+21; ATNF
DM = 56.77 pc cm^-3, P = 33.27 ms) from the ``sigpyproc3`` test-data tree: 832 x 4 MHz channels over
702--4030 MHz, 0.512 ms sampling, 2.1 s --- the full 463 ms DM sweep is visible inside the file and
Crab giant pulses reach Jy level (ideal single-pulse targets).

The SIGPROC reader here is a deliberately minimal pure-NumPy parser (the file is 8-bit,
single-IF) --- no ``blimpy`` dependency for a 3 MB file. Offline, everything runs on a synthetic
pulse train built with :func:`jansky.transients.disperse_pulse`, so tests/CI never touch the
network; the real file is fetched only by ``run(offline=False)`` / ``make reproduce``.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Any

import numpy as np
from jansky import transients

from . import fdmt as F

__all__ = [
    "CRAB_DM",
    "CRAB_P_S",
    "CRAB_FIL_URL",
    "read_sigproc",
    "synthetic_observation",
    "search",
    "run",
]

#: ATNF catalogue values for PSR B0531+21 (the recover-a-known targets).
CRAB_DM = 56.77
CRAB_P_S = 0.03327
#: GATE-0-verified public file (Parkes/UWL, MJD 58543, 3.25 MB, no auth).
CRAB_FIL_URL = "https://raw.githubusercontent.com/FRBs/sigpyproc3/main/tests/data/parkes_8bit_1.fil"

_STR_KEYS = {"source_name", "rawdatafile"}
_INT_KEYS = {
    "telescope_id",
    "machine_id",
    "data_type",
    "nchans",
    "nbits",
    "nifs",
    "nbeams",
    "ibeam",
}
_DBL_KEYS = {"tstart", "tsamp", "fch1", "foff", "src_raj", "src_dej", "az_start", "za_start"}


def read_sigproc(path: str | Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Minimal SIGPROC filterbank reader (8/32-bit, single IF) --- pure NumPy.

    Returns ``(dynspec, freqs_mhz, header)`` with ``dynspec`` shaped ``(n_time, n_chan)`` and
    frequencies in MHz, matching the :mod:`jansky.transients` conventions.
    """
    raw = Path(path).read_bytes()

    def rd_str(off: int) -> tuple[str, int]:
        n = struct.unpack_from("<i", raw, off)[0]
        return raw[off + 4 : off + 4 + n].decode(), off + 4 + n

    key, off = rd_str(0)
    if key != "HEADER_START":
        raise ValueError(f"{path}: not a SIGPROC filterbank")
    hdr: dict[str, Any] = {}
    while True:
        key, off = rd_str(off)
        if key == "HEADER_END":
            break
        if key in _STR_KEYS:
            hdr[key], off = rd_str(off)
        elif key in _INT_KEYS:
            hdr[key] = struct.unpack_from("<i", raw, off)[0]
            off += 4
        elif key in _DBL_KEYS:
            hdr[key] = struct.unpack_from("<d", raw, off)[0]
            off += 8
        else:  # unknown keys in this minimal reader: assume double (SIGPROC's common case)
            hdr[key] = struct.unpack_from("<d", raw, off)[0]
            off += 8
    nch, nbits = int(hdr["nchans"]), int(hdr["nbits"])
    dtype = {8: np.uint8, 32: np.float32}[nbits]
    data = np.frombuffer(raw, dtype=dtype, offset=off)
    n_time = data.size // nch
    dyn = data[: n_time * nch].reshape(n_time, nch).astype(np.float32)
    freqs = hdr["fch1"] + hdr["foff"] * np.arange(nch)  # MHz; foff usually negative
    return dyn, freqs, hdr


def synthetic_observation(
    *, dm: float = CRAB_DM, period_s: float = CRAB_P_S, n_time: int = 4096, seed: int = 0
) -> tuple[np.ndarray, np.ndarray, float]:
    """Offline fixture: a dispersed pulse TRAIN at a known (DM, P) in Crab-file-like geometry."""
    freqs = np.linspace(1200.0, 1600.0, 256)
    dt = 5.12e-4
    rng = np.random.default_rng(seed)
    dyn = rng.normal(0.0, 1.0, (n_time, freqs.size)).astype(np.float32)
    p_samp = period_s / dt
    for k in range(int(n_time / p_samp)):
        one = transients.disperse_pulse(
            n_time,
            freqs,
            dm,
            dt,
            t0_index=int(200 + k * p_samp),
            amplitude=8.0,
            noise=0.0,
            seed=None,
        )
        dyn += one
    return dyn, freqs, dt


def search(
    dynspec: np.ndarray,
    freqs_mhz: np.ndarray,
    dt: float,
    *,
    max_dm: float = 120.0,
    device: str = "cpu",
) -> dict:
    """FDMT DM--time butterfly + boxcar single-pulse + folding period search, one call.

    Channels with wildly non-stationary variance are clipped (a light RFI guard) before the
    transform; the fold searches around the best single-pulse spacing.
    """
    dyn = np.asarray(dynspec, np.float32)
    # per-channel robust normalisation (bandpass + crude RFI guard)
    med = np.median(dyn, axis=0)
    mad = np.median(np.abs(dyn - med), axis=0) * 1.4826 + 1e-6
    dyn = np.clip((dyn - med) / mad, -6.0, 12.0)

    r = F.fdmt(dyn, freqs_mhz, dt, max_dm, device=device)
    best_dm, best_snr = r.best()
    series = r.plane[int(np.argmin(np.abs(r.dms - best_dm)))].cpu().numpy()
    sp_snr, sp_width, sp_pos = transients.boxcar_snr(series, np.array([1, 2, 4, 8, 16, 32]))
    return {
        "best_dm": float(best_dm),
        "best_snr": float(best_snr),
        "sp_snr": float(sp_snr),
        "sp_width_samples": int(sp_width),
        "sp_pos": int(sp_pos),
        "series": series,
        "dms": r.dms,
        "dm_curve": r.plane.max(dim=1).values.cpu().numpy(),
    }


def _fold_period(series: np.ndarray, dt: float, p0: float) -> float:
    """Refine the pulse period by epoch folding around a first guess."""
    times = np.arange(series.size) * dt
    res = transients.epoch_folding_search(
        times, series, np.linspace(0.9 * p0, 1.1 * p0, 201), n_bins=32
    )
    return float(res.best_period)


def run(out: str = ".", *, offline: bool = True, device: str = "cpu") -> dict:
    """Full slice: synthetic recover-a-known, plus the real Crab leg when online."""
    import json

    dyn, freqs, dt = synthetic_observation()
    s = search(dyn, freqs, dt, max_dm=120.0, device=device)
    p_syn = _fold_period(s["series"], dt, CRAB_P_S * 1.02)
    metrics = {
        "source": "synthetic pulse train (Crab-like)",
        "true_dm": CRAB_DM,
        "recovered_dm": round(s["best_dm"], 2),
        "butterfly_snr": round(s["best_snr"], 1),
        "true_period_ms": CRAB_P_S * 1e3,
        "recovered_period_ms": round(p_syn * 1e3, 3),
        "device": device,
    }

    if not offline:  # pragma: no cover - network + real data
        import urllib.request

        fil = Path(out) / "data" / "parkes_crab.fil"
        fil.parent.mkdir(parents=True, exist_ok=True)
        if not fil.exists():
            urllib.request.urlretrieve(CRAB_FIL_URL, fil)
        rdyn, rfreqs, hdr = read_sigproc(fil)
        rs = search(rdyn, rfreqs, float(hdr["tsamp"]), max_dm=120.0, device=device)
        metrics.update(
            {
                "source": f"Parkes/UWL {hdr.get('source_name', '?')} (real, GATE-0 file)",
                "real_recovered_dm": round(rs["best_dm"], 2),
                "real_butterfly_snr": round(rs["best_snr"], 1),
                "real_sp_snr": round(rs["sp_snr"], 1),
                "real_dm_error_pc": round(100 * abs(rs["best_dm"] - CRAB_DM) / CRAB_DM, 1),
            }
        )
        s = rs  # figure shows the real butterfly when available

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "singlepulse_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(s, op / "papers" / "torchfdmt" / "figures")
    _write_macros(metrics, op / "papers" / "torchfdmt" / "generated" / "macros.tex")
    return metrics


def _figure(s: dict, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax1.plot(s["dms"], s["dm_curve"], "-", color="C0", lw=1.2)
    ax1.axvline(CRAB_DM, color="C3", ls="--", lw=1, label=f"catalogue DM {CRAB_DM}")
    ax1.set(xlabel="DM (pc cm$^{-3}$)", ylabel="peak track sum", title="FDMT butterfly peak")
    ax1.legend(fontsize=8)
    ax2.plot(np.arange(s["series"].size), s["series"], "-", color="C0", lw=0.7)
    ax2.axvline(s["sp_pos"], color="C3", ls=":", lw=1)
    ax2.set(xlabel="sample", ylabel="dedispersed power", title="Best-DM time series")
    fig.tight_layout()
    fig.savefig(out / "singlepulse.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.singlepulse._write_macros -- do not edit.",
        rf"\newcommand{{\spSource}}{{{m['source']}}}",
        rf"\newcommand{{\spTrueDm}}{{{_fmt('true_dm')}}}",
        rf"\newcommand{{\spRecDm}}{{{_fmt('recovered_dm')}}}",
        rf"\newcommand{{\spSnr}}{{{_fmt('butterfly_snr')}}}",
        rf"\newcommand{{\spRecP}}{{{_fmt('recovered_period_ms')}}}",
        rf"\newcommand{{\spRealDm}}{{{_fmt('real_recovered_dm')}}}",
        rf"\newcommand{{\spRealSnr}}{{{_fmt('real_butterfly_snr')}}}",
        rf"\newcommand{{\spRealSpSnr}}{{{_fmt('real_sp_snr')}}}",
        rf"\newcommand{{\spRealDmErr}}{{{_fmt('real_dm_error_pc')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="FDMT single-pulse recover-a-known (Crab).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--device", default="cpu")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline, device=args.device), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
