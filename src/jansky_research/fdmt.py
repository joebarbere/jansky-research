"""A pure-PyTorch Fast Dispersion Measure Transform (FDMT) --- GPU dedispersion without CUDA.

The FDMT (Zackay & Ofek 2017, ApJ 835, 11) computes the full DM--time plane of a dynamic spectrum in
:math:`O(N_f N_t \\log N_f)` operations instead of the brute-force :math:`O(N_f N_t N_\\mathrm{DM})`,
by recursively merging sub-band transforms (the dispersion curve :math:`t \\propto \\nu^{-2}` is
additive across sub-bands). Every production GPU dedisperser (Heimdall/``dedisp``, astro-accelerate,
FREDDA) ships CUDA kernels; this module instead expresses both the FDMT recursion and the brute-force
baseline as **pure tensor ops** (gather/roll/add), so one ``device=`` argument runs them on CPU,
CUDA, *or* ROCm (verified on an RX 7600 XT / gfx1102) --- no custom kernels.

Correctness is anchored three ways (see ``tests/test_fdmt.py``): the delay-0 row must equal the plain
channel sum exactly; the recovered DM must agree with both :func:`brute_dedisperse` and the tested
:func:`jansky.transients.dm_search` oracle on a synthetic dispersed pulse. One semantic difference is
intentional: FDMT sums the **full track** including the intra-channel dispersion smear (its
per-channel partial sums span the in-channel delay), while brute roll-and-sum picks one sample per
channel --- so for a smeared pulse the FDMT peak amplitude and S/N are *higher* (it integrates all
the pulse energy). Peak location is the like-for-like check. Conventions follow
``jansky.transients``: dynamic spectra are ``(n_time, n_chan)``, frequencies in MHz, delays measured
relative to the **highest** frequency.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from jansky.constants import DM_CONST

if TYPE_CHECKING:  # pragma: no cover - typing only
    pass

__all__ = [
    "FDMTResult",
    "benchmark",
    "brute_dedisperse",
    "delay_samples",
    "dm_from_delay",
    "fdmt",
]


def _require_torch() -> Any:
    """Import torch lazily so the core package works without the ``fdmt`` extra."""
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "jansky_research.fdmt needs PyTorch: install with `uv sync --extra fdmt` "
            "(CPU wheel; for GPU use a ROCm/CUDA venv)"
        ) from exc
    return torch


def delay_samples(dm: float, f_lo_mhz: float, f_hi_mhz: float, dt: float) -> int:
    r"""Total dispersion delay across the band, in whole samples.

    :math:`\Delta t = k_\mathrm{DM}\,\mathrm{DM}\,(f_\mathrm{lo}^{-2} - f_\mathrm{hi}^{-2})`,
    rounded to samples of length ``dt`` --- the FDMT's integer DM-axis unit.
    """
    return int(round(DM_CONST * dm * (f_lo_mhz**-2.0 - f_hi_mhz**-2.0) / dt))


def dm_from_delay(delay: np.ndarray | int, f_lo_mhz: float, f_hi_mhz: float, dt: float) -> Any:
    """Invert :func:`delay_samples`: the DM (pc cm^-3) whose band delay is ``delay`` samples."""
    return np.asarray(delay, float) * dt / (DM_CONST * (f_lo_mhz**-2.0 - f_hi_mhz**-2.0))


@dataclass
class FDMTResult:
    """The DM--time plane plus its axes: ``plane[d, t]`` sums the track with band delay ``d``."""

    plane: Any  # torch.Tensor (n_delay, n_time)
    delays: np.ndarray  # integer band delays (samples)
    dms: np.ndarray  # the same axis in pc cm^-3

    def best(self) -> tuple[float, float]:
        """(best_dm, peak value of a z-scored time series at that delay) --- the butterfly peak."""
        plane = self.plane
        med = plane.median(dim=1, keepdim=True).values
        mad = (plane - med).abs().median(dim=1, keepdim=True).values * 1.4826 + 1e-12
        snr = ((plane - med) / mad).max(dim=1).values
        i = int(snr.argmax())
        return float(self.dms[i]), float(snr[i])


def _subband_delay(d_total: int, g_mid: float, g_hi: float, g_lo: float) -> int:
    """The high-sub-band share of a track's total delay: d * (g_mid - g_hi) / (g_lo - g_hi)."""
    if g_lo == g_hi:
        return 0
    return int(round(d_total * (g_mid - g_hi) / (g_lo - g_hi)))


def fdmt(
    dynspec: np.ndarray,
    freqs_mhz: np.ndarray,
    dt: float,
    max_dm: float,
    *,
    device: str = "cpu",
) -> FDMTResult:
    """The Fast DM Transform of a dynamic spectrum, on any torch device.

    Parameters follow :func:`jansky.transients.dedisperse`: ``dynspec`` is ``(n_time, n_chan)``,
    ``freqs_mhz`` the channel centres, ``dt`` the sample time (s). Returns every integer band delay
    ``0..delay_samples(max_dm)`` --- the full DM--time "butterfly" plane.

    Implementation: the classic recursion on channel-count halvings. State is a tensor
    ``A[band, delay, time]`` of partial track sums; each iteration merges adjacent band pairs,
    splitting each output delay at the shared frequency edge with the :math:`\\nu^{-2}` map and
    adding the high band shifted by the low band's delay. Channels are padded to a power of two
    with zero rows (standard practice; zero rows add nothing to any track).
    """
    torch = _require_torch()
    dynspec = np.asarray(dynspec, dtype=np.float32)
    freqs = np.asarray(freqs_mhz, dtype=float)
    n_time = dynspec.shape[0]

    # Channels ordered by DESCENDING frequency: row 0 arrives first, delays grow downward.
    order = np.argsort(freqs)[::-1]
    freqs = freqs[order]
    data = torch.as_tensor(np.ascontiguousarray(dynspec[:, order].T), device=device)

    # Channel edges from centres (uniform grid assumed --- filterbank convention).
    df = abs(freqs[0] - freqs[1]) if freqs.size > 1 else abs(freqs[0]) * 1e-3
    edges_hi = freqs + df / 2.0
    edges_lo = freqs - df / 2.0
    f_top, f_bot = float(edges_hi[0]), float(edges_lo[-1])
    d_max = delay_samples(max_dm, f_bot, f_top, dt)

    # Pad the channel axis to a power of two with zero rows below the band.
    n_chan = len(freqs)
    n_pow2 = 1 << (n_chan - 1).bit_length()
    if n_pow2 != n_chan:
        pad = torch.zeros(n_pow2 - n_chan, n_time, device=device, dtype=data.dtype)
        data = torch.cat([data, pad], dim=0)
        # pad rows carry ZERO dispersion span (both edges at the real band bottom), so the
        # merge tree assigns them no share of any track's delay budget
        extra = np.full(n_pow2 - n_chan, f_bot)
        edges_hi = np.concatenate([edges_hi, extra])
        edges_lo = np.concatenate([edges_lo, extra])

    g = lambda f: float(f) ** -2.0  # noqa: E731 - the dispersion kernel
    g_top = g(f_top)
    span = g(f_bot) - g_top

    def band_dmax(lo: float, hi: float) -> int:
        """Max in-band delay a d_max track accumulates across sub-band [lo, hi]."""
        return int(np.ceil(d_max * (g(lo) - g(hi)) / span)) if span > 0 else 0

    # --- initialisation: per-channel partial sums over 0..band_dmax(channel) samples ---
    d_init = max(band_dmax(float(edges_lo[i]), float(edges_hi[i])) for i in range(n_pow2))
    state = torch.zeros(n_pow2, d_init + 1, n_time, device=device, dtype=data.dtype)
    state[:, 0, :] = data
    for d in range(1, d_init + 1):
        # A[c, d, t] = A[c, d-1, t] + I[c, t+d]: the track's sum extends one later sample.
        state[:, d, :-d] = state[:, d - 1, :-d] + data[:, d:]
        state[:, d, -d:] = state[:, d - 1, -d:]
    lo_edges = edges_lo.copy()
    hi_edges = edges_hi.copy()

    # --- recursion: merge adjacent band pairs until one band spans everything ---
    while state.shape[0] > 1:
        n_band = state.shape[0]
        n_out = n_band // 2
        new_lo = lo_edges[1::2]  # lower edge of each merged pair (bands are freq-descending)
        new_hi = hi_edges[0::2]
        out_dmax = max(band_dmax(float(new_lo[b]), float(new_hi[b])) for b in range(n_out))
        out = torch.zeros(n_out, out_dmax + 1, n_time, device=device, dtype=data.dtype)
        for b in range(n_out):
            hi_b, mid_b, lo_b = float(hi_edges[2 * b]), float(lo_edges[2 * b]), float(new_lo[b])
            b_dmax = band_dmax(lo_b, hi_b)
            for d in range(b_dmax + 1):
                # delay accumulated in the HIGH sub-band [mid_b, hi_b]; the nu^-2 curve is
                # additive, so its share of d is (g(mid)-g(hi)) / (g(lo)-g(hi)).
                d_hi = _subband_delay(d, g(mid_b), g(hi_b), g(lo_b))
                d_lo = d - d_hi
                d_hi_c = min(d_hi, state.shape[1] - 1)
                d_lo_c = min(d_lo, state.shape[1] - 1)
                a = state[2 * b, d_hi_c, :]
                bl = state[2 * b + 1, d_lo_c, :]
                if d_hi > 0:
                    # the low band starts d_hi samples after the track enters the high band
                    shifted = torch.zeros_like(bl)
                    shifted[:-d_hi] = bl[d_hi:]
                    out[b, d, :] = a + shifted
                else:
                    out[b, d, :] = a + bl
        state = out
        lo_edges, hi_edges = new_lo, new_hi

    plane = state[0, : d_max + 1, :]
    delays = np.arange(plane.shape[0])
    return FDMTResult(plane=plane, delays=delays, dms=dm_from_delay(delays, f_bot, f_top, dt))


def brute_dedisperse(
    dynspec: np.ndarray,
    freqs_mhz: np.ndarray,
    dt: float,
    dms: np.ndarray,
    *,
    device: str = "cpu",
) -> Any:
    """Brute-force roll-and-sum dedispersion over a DM grid, as one batched tensor gather.

    The torch twin of :func:`jansky.transients.dedisperse` (same per-channel integer shifts
    relative to the highest frequency), vectorised over all trial DMs --- the baseline the FDMT
    is benchmarked against, and the always-correct fallback engine.
    """
    torch = _require_torch()
    dynspec = np.asarray(dynspec, dtype=np.float32)
    freqs = np.asarray(freqs_mhz, dtype=float)
    dms = np.asarray(dms, dtype=float)
    n_time = dynspec.shape[0]
    f_ref = freqs.max()

    data = torch.as_tensor(np.ascontiguousarray(dynspec.T), device=device)  # (n_chan, n_time)
    # integer shift per (dm, chan), matching the oracle's round(delay/dt)
    shifts = np.rint(DM_CONST * dms[:, None] * (freqs[None, :] ** -2.0 - f_ref**-2.0) / dt)
    sh = torch.as_tensor(shifts.astype(np.int64), device=device)
    t_idx = torch.arange(n_time, device=device)
    out = torch.empty(len(dms), n_time, device=device, dtype=data.dtype)
    batch = max(1, int(2e8 // (data.numel() or 1)))  # keep gather buffers ~sub-GB
    for i in range(0, len(dms), batch):
        idx = (t_idx[None, None, :] + sh[i : i + batch, :, None]) % n_time  # roll(-shift)
        out[i : i + batch] = data[None].expand(idx.shape[0], -1, -1).gather(2, idx).sum(1)
    return out


def benchmark(
    *,
    n_time: int = 4096,
    n_chan: int = 512,
    max_dm: float = 500.0,
    dt: float = 1e-3,
    f_lo: float = 1200.0,
    f_hi: float = 1600.0,
    devices: tuple[str, ...] = ("cpu",),
    repeats: int = 3,
) -> dict:
    """Wall-time table: numpy oracle vs torch brute vs torch FDMT, per device.

    The paper's headline table. Honest scope: torch-vs-CPU-baselines on this hardware --- not a
    comparison against tuned CUDA dedispersers on datacentre GPUs.
    """
    import time

    from jansky import transients

    torch = _require_torch()
    freqs = np.linspace(f_lo, f_hi, n_chan)
    dyn = transients.disperse_pulse(n_time, freqs, max_dm / 2.0, dt, t0_index=n_time // 8, seed=0)
    n_dm = delay_samples(max_dm, float(freqs.min()), float(freqs.max()), dt) + 1
    dms = np.asarray(dm_from_delay(np.arange(n_dm), float(freqs.min()), float(freqs.max()), dt))

    rows: dict[str, float] = {}
    t0 = time.perf_counter()
    transients.dm_search(dyn, freqs, dt, dms[:: max(1, n_dm // 256)])  # oracle at reduced grid
    rows["numpy_oracle_reduced_s"] = time.perf_counter() - t0

    for dev in devices:
        for name, fn in (
            ("brute", lambda d=dev: brute_dedisperse(dyn, freqs, dt, dms, device=d)),
            ("fdmt", lambda d=dev: fdmt(dyn, freqs, dt, max_dm, device=d)),
        ):
            fn()  # warm-up (JIT/alloc)
            if dev != "cpu":  # pragma: no cover - GPU only
                torch.cuda.synchronize()
            t0 = time.perf_counter()
            for _ in range(repeats):
                fn()
            if dev != "cpu":  # pragma: no cover - GPU only
                torch.cuda.synchronize()
            rows[f"{name}_{dev}_s"] = (time.perf_counter() - t0) / repeats
    rows["n_dm_trials"] = float(n_dm)
    return rows
