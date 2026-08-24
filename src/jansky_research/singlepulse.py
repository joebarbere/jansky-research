"""Single-pulse / pulsar recover-a-known on a real public filterbank, FDMT-powered (plan 34).

The science leg of the torch-fdmt arc: read a small public SIGPROC filterbank of a KNOWN bright
pulsar, dedisperse with :mod:`jansky_research.fdmt` over a DM grid, and recover the catalogued DM
(and, by folding, the period) --- a validation on real telescope data, not a discovery. GATE 0
(2026-07-02) chose a 3.25 MB Parkes/UWL observation of the **Crab pulsar** (B0531+21; ATNF
DM = 56.77 pc cm^-3, P ~ 33.39 ms (epoch-dependent)) from the ``sigpyproc3`` test-data tree: 832 x 4 MHz channels over
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
    "shift_null",
    "run",
]

#: ATNF catalogue values for PSR B0531+21 (the recover-a-known targets).
CRAB_DM = 56.77
CRAB_P_S = 0.033392  # ATNF P0 (epoch-dependent: Crab Pdot~4.2e-13; value near the file MJD)
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
    dyn = _normalise(dynspec)
    r = F.fdmt(dyn, freqs_mhz, dt, max_dm, device=device)
    best_dm, best_snr = r.best()
    series = r.plane[int(np.argmin(np.abs(r.dms - best_dm)))].cpu().numpy()
    sp_snr, sp_width, sp_pos = transients.boxcar_snr(series, np.array([1, 2, 4, 8, 16, 32]))
    # The per-row z-scored peak — the SAME statistic best() maximises. An earlier figure
    # plotted the raw track sum, which grows with delay row (more samples integrated), so its
    # maximum sat ~60 pc/cm^3 from the value the caption quoted.
    plane = r.plane
    med = plane.median(dim=1, keepdim=True).values
    mad = (plane - med).abs().median(dim=1, keepdim=True).values * 1.4826 + 1e-12
    snr_curve = ((plane - med) / mad).max(dim=1).values.cpu().numpy()
    return {
        "best_dm": float(best_dm),
        "best_snr": float(best_snr),
        "sp_snr": float(sp_snr),
        "sp_width_samples": int(sp_width),
        "sp_pos": int(sp_pos),
        "n_time": int(dyn.shape[0]),
        "series": series,
        "dms": r.dms,
        "dm_curve": snr_curve,
        "dm_curve_raw": plane.max(dim=1).values.cpu().numpy(),
    }


def _normalise(dynspec: np.ndarray) -> np.ndarray:
    """Per-channel robust normalisation (bandpass + crude RFI guard), shared by search and null."""
    dyn = np.asarray(dynspec, np.float32)
    med = np.median(dyn, axis=0)
    mad = np.median(np.abs(dyn - med), axis=0) * 1.4826 + 1e-6
    return np.clip((dyn - med) / mad, -6.0, 12.0)


def shift_null(
    dynspec: np.ndarray,
    freqs_mhz: np.ndarray,
    dt: float,
    *,
    max_dm: float = 120.0,
    n_reps: int = 200,
    device: str = "cpu",
    seed: int = 0,
) -> dict:
    """Null distribution of the butterfly peak under per-channel circular time shifts.

    Each repetition rolls every channel by an independent uniform offset, destroying any
    dispersed track while preserving each channel's own statistics and RFI, then records the
    maximum of the same per-row z-scored statistic :func:`search` reports. This is the trials
    accounting a bare "S/N = 6.0" lacks: the butterfly plane holds millions of cells, so its
    noise maximum is itself several sigma, and only this distribution says where a measured
    peak height sits against it.
    """
    rng = np.random.default_rng(seed)
    dyn = _normalise(dynspec)
    n_t, n_c = dyn.shape
    vals = []
    for _ in range(n_reps):
        shifted = np.empty_like(dyn)
        offs = [int(x) for x in rng.integers(0, n_t, size=(n_c,))]
        for c in range(n_c):
            shifted[:, c] = np.roll(dyn[:, c], offs[c])
        r = F.fdmt(shifted, freqs_mhz, dt, max_dm, device=device)
        vals.append(r.best()[1])
    arr = np.asarray(vals)
    return {
        "n_reps": int(n_reps),
        "mean": float(arr.mean()),
        "p50": float(np.percentile(arr, 50)),
        "p99": float(np.percentile(arr, 99)),
        "max": float(arr.max()),
        "best_snrs": arr,
    }


def _dm_step(freqs_mhz: np.ndarray, dt: float) -> float:
    """DM change corresponding to one sample of dispersive delay across the band.

    The FDMT butterfly indexes rows by integer delay samples, so the recovered DM is
    quantised at this step and an offset should be read in trials, not in percent.
    """
    f = np.asarray(freqs_mhz, dtype=float)
    lo, hi = float(f.min()), float(f.max())
    return dt / (4.148808e3 * (lo**-2 - hi**-2))


def _fold_period(series: np.ndarray, dt: float, p0: float) -> float:
    """Refine the pulse period by epoch folding around a first guess."""
    times = np.arange(series.size) * dt
    res = transients.epoch_folding_search(
        times, series, np.linspace(0.9 * p0, 1.1 * p0, 201), n_bins=32
    )
    return float(res.best_period)


def run(
    out: str = ".",
    *,
    offline: bool = True,
    device: str = "cpu",
    bench: bool = False,
    bench_devices: tuple[str, ...] | None = None,
    null_reps: int = 200,
) -> dict:
    """Full slice: synthetic recover-a-known, plus the real Crab leg when online.

    ``bench_devices`` decouples the benchmark's device set from the science leg's ``device``:
    previously ``devs`` was derived from ``device``, so the combination the paper reports
    (CPU science + both benchmark columns) was unreachable by any single invocation — which is
    exactly how the committed row came to be assembled by hand from two runs.
    """

    dyn, freqs, dt = synthetic_observation()
    s = search(dyn, freqs, dt, max_dm=120.0, device=device)
    p_syn = _fold_period(s["series"], dt, CRAB_P_S * 1.02)
    metrics = {
        "source": "synthetic pulse train (Crab-like)",
        "true_dm": CRAB_DM,
        "catalogue_dm": CRAB_DM,
        "recovered_dm": round(s["best_dm"], 2),
        "butterfly_snr": round(s["best_snr"], 1),
        "true_period_ms": round(CRAB_P_S * 1e3, 3),
        "recovered_period_ms": round(p_syn * 1e3, 3),
        "device": device,
    }

    if bench:  # pragma: no cover - timing-dependent; reproduced via --benchmark
        from . import fdmt as _F

        devs = bench_devices or (("cpu", "cuda") if device != "cpu" else ("cpu",))
        b = _F.benchmark(n_time=8192, n_chan=1024, max_dm=800.0, devices=devs, repeats=3)
        hw = "unknown"
        try:
            import platform

            import torch

            if "cuda" in devs and torch.cuda.is_available():
                hw = f"{torch.cuda.get_device_name(0)}, torch {torch.__version__}"
            else:
                hw = f"{platform.processor() or platform.machine()}, torch {torch.__version__}"
        except Exception:
            pass
        metrics.update(
            {
                "bench_brute_cpu_s": round(b["brute_cpu_s"], 2),
                "bench_fdmt_cpu_s": round(b["fdmt_cpu_s"], 2),
                "bench_brute_gpu_s": round(b.get("brute_cuda_s", float("nan")), 2)
                if "brute_cuda_s" in b
                else None,
                "bench_fdmt_gpu_s": round(b.get("fdmt_cuda_s", float("nan")), 2)
                if "fdmt_cuda_s" in b
                else None,
                "bench_numpy_oracle_s": round(b["numpy_oracle_reduced_s"], 2),
                "bench_n_dm": int(b["n_dm_trials"]),
                # Which invocation produced this row. A CPU-only run leaves the GPU
                # columns null; patching GPU numbers into that JSON later yields a row no
                # single run could produce, which is how the committed row came to pair
                # device="cpu" with both columns filled.
                "bench_devices": list(devs),
                # Emitted by torch introspection, never typed: a hand-entered hardware string
                # would be inherited by future runs through the results merge and misattribute
                # them (the referee-caught failure).
                "benchmark_device": "cuda" if "cuda" in devs else "cpu",
                "benchmark_hardware": hw,
            }
        )

    if not offline:  # pragma: no cover - network + real data
        import urllib.request

        fil = Path(out) / "data" / "parkes_crab.fil"
        fil.parent.mkdir(parents=True, exist_ok=True)
        if not fil.exists():
            urllib.request.urlretrieve(CRAB_FIL_URL, fil)
        rdyn, rfreqs, hdr = read_sigproc(fil)
        rdt = float(hdr["tsamp"])
        rs = search(rdyn, rfreqs, rdt, max_dm=120.0, device=device)
        dm_step = _dm_step(rfreqs, rdt)
        n_trials = int(rs["dms"].size)
        offset_trials = abs(rs["best_dm"] - CRAB_DM) / dm_step
        # The trials accounting for the butterfly height: per-channel circular-shift null.
        null = shift_null(rdyn, rfreqs, rdt, max_dm=120.0, n_reps=null_reps, device=device)
        n_ge = int(np.sum(null["best_snrs"] >= rs["best_snr"]))
        metrics.update(
            {
                # Name BOTH legs. The unprefixed keys (recovered_dm, butterfly_snr,
                # recovered_period_ms) are the SYNTHETIC injection; the real ones all carry
                # a real_ prefix. A `source` naming only the Parkes file sat directly above
                # the synthetic recovered_dm, so an auditor reading the evidence file top-
                # down would take 56.63 for the Crab recovery when the Crab value is 56.59.
                "source": (
                    "synthetic injection (unprefixed keys) + Parkes/UWL "
                    f"{hdr.get('source_name', '?')} (real, GATE-0 file; real_* keys)"
                ),
                "real_recovered_dm": round(rs["best_dm"], 2),
                "real_butterfly_snr": round(rs["best_snr"], 1),
                "real_sp_snr": round(rs["sp_snr"], 1),
                "real_sp_pos": int(rs["sp_pos"]),
                "real_sp_width_samples": int(rs["sp_width_samples"]),
                "real_n_time": int(rs["n_time"]),
                "real_dm_error_pc": round(100 * abs(rs["best_dm"] - CRAB_DM) / CRAB_DM, 1),
                # One FDMT row is one delay sample, so this is the DM quantisation of the
                # recovery: without it a "0.3%" offset cannot be read as few-trial agreement.
                "real_dm_step_pc": round(dm_step, 4),
                "real_n_dm_trials": n_trials,
                "real_dm_offset_trials": round(float(offset_trials), 1),
                # p for a uniformly placed noise maximum to land this close to the catalogue
                # DM: the honest positional statement behind "recovers the Crab DM".
                "real_p_positional": round((2.0 * offset_trials + 1.0) / n_trials, 5),
                "real_null_reps": null["n_reps"],
                "real_null_p50": round(null["p50"], 2),
                "real_null_p99": round(null["p99"], 2),
                "real_null_max": round(null["max"], 2),
                # (k+1)/(n+1): the butterfly peak height against the shift null.
                "real_p_null": round((n_ge + 1) / (null["n_reps"] + 1), 4),
            }
        )
        s = rs  # figure shows the real butterfly when available

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "singlepulse_metrics.json")
    _figure(s, op / "papers" / "torchfdmt" / "figures")
    _write_macros(metrics, op / "papers" / "torchfdmt" / "generated" / "macros.tex")
    return metrics


def _figure(s: dict, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    # Plot the statistic best() maximises (per-row z-scored peak). The raw track sum rises
    # with delay row and put the displayed maximum ~60 pc/cm^3 from the quoted DM.
    ax1.plot(s["dms"], s["dm_curve"], "-", color="C0", lw=1.2, label="row-normalised peak")
    ax1.axvline(CRAB_DM, color="C3", ls="--", lw=1, label=f"catalogue DM {CRAB_DM}")
    ax1.set(
        xlabel="DM (pc cm$^{-3}$)",
        ylabel="peak of z-scored series",
        title="FDMT butterfly peak",
    )
    ax1.legend(fontsize=8)
    ax2.plot(np.arange(s["series"].size), s["series"], "-", color="C0", lw=0.7)
    ax2.axvline(s["sp_pos"], color="C3", ls=":", lw=1, label="best boxcar detection")
    ax2.legend(fontsize=8)
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
        rf"\newcommand{{\spTrueP}}{{{_fmt('true_period_ms')}}}",
        rf"\newcommand{{\spRecDm}}{{{_fmt('recovered_dm')}}}",
        rf"\newcommand{{\spSnr}}{{{_fmt('butterfly_snr')}}}",
        rf"\newcommand{{\spRecP}}{{{_fmt('recovered_period_ms')}}}",
        rf"\newcommand{{\spRealDm}}{{{_fmt('real_recovered_dm')}}}",
        rf"\newcommand{{\spRealSnr}}{{{_fmt('real_butterfly_snr')}}}",
        rf"\newcommand{{\spRealSpSnr}}{{{_fmt('real_sp_snr')}}}",
        rf"\newcommand{{\spRealDmErr}}{{{_fmt('real_dm_error_pc')}}}",
        rf"\newcommand{{\spCatDm}}{{{_fmt('catalogue_dm')}}}",
        rf"\newcommand{{\spRealDmStep}}{{{_fmt('real_dm_step_pc')}}}",
        rf"\newcommand{{\spRealTrials}}{{{_fmt('real_n_dm_trials')}}}",
        rf"\newcommand{{\spRealOffTrials}}{{{_fmt('real_dm_offset_trials')}}}",
        rf"\newcommand{{\spRealPpos}}{{{_fmt('real_p_positional')}}}",
        rf"\newcommand{{\spRealNullReps}}{{{_fmt('real_null_reps')}}}",
        rf"\newcommand{{\spRealNullMedian}}{{{_fmt('real_null_p50')}}}",
        rf"\newcommand{{\spRealNullNinetyNine}}{{{_fmt('real_null_p99')}}}",
        rf"\newcommand{{\spRealNullMax}}{{{_fmt('real_null_max')}}}",
        rf"\newcommand{{\spRealPnull}}{{{_fmt('real_p_null')}}}",
        rf"\newcommand{{\spBruteCpu}}{{{_fmt('bench_brute_cpu_s')}}}",
        rf"\newcommand{{\spBruteGpu}}{{{_fmt('bench_brute_gpu_s')}}}",
        rf"\newcommand{{\spFdmtCpu}}{{{_fmt('bench_fdmt_cpu_s')}}}",
        rf"\newcommand{{\spFdmtGpu}}{{{_fmt('bench_fdmt_gpu_s')}}}",
        rf"\newcommand{{\spOracle}}{{{_fmt('bench_numpy_oracle_s')}}}",
        rf"\newcommand{{\spBenchNdm}}{{{_fmt('bench_n_dm')}}}",
        rf"\newcommand{{\spBenchHw}}{{{_fmt('benchmark_hardware')}}}",
    ]
    # Ratios are DERIVED, never typed. A hand-written "~24x" survived in the paper next to
    # the two macros it is computed from, long after a re-benchmark made it 29x.
    for name, num, den in (
        ("spBruteSpeedup", "bench_brute_cpu_s", "bench_brute_gpu_s"),
        ("spFdmtSpeedup", "bench_fdmt_cpu_s", "bench_fdmt_gpu_s"),
    ):
        a_, b_ = m.get(num), m.get(den)
        ratio = "--"
        if isinstance(a_, int | float) and isinstance(b_, int | float) and b_:
            ratio = f"{a_ / b_:.0f}"
        lines.append(rf"\newcommand{{\{name}}}{{{ratio}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a CPU run has no GPU benchmark keys and would blank
    # the GPU macros a previous run wrote. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="FDMT single-pulse recover-a-known (Crab).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--device", default="cpu")
    p.add_argument("--benchmark", action="store_true")
    p.add_argument(
        "--bench-devices",
        default=None,
        help="comma-separated benchmark device set, decoupled from --device (e.g. cpu,cuda)",
    )
    p.add_argument("--null-reps", type=int, default=200)
    args = p.parse_args(argv)
    bdevs = tuple(args.bench_devices.split(",")) if args.bench_devices else None
    print(
        json.dumps(
            run(
                args.out,
                offline=args.offline,
                device=args.device,
                bench=args.benchmark,
                bench_devices=bdevs,
                null_reps=args.null_reps,
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
