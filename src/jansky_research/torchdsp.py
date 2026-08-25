"""torch-dsp: coherent dedispersion, RFI excision, and an FFA in pure PyTorch (plan 43).

Extends the merged `torch-fdmt` arc into a coherent-DSP suite. GATE-0 (2026-07-06, per-kernel
repo/full-text sweep): **no pure-PyTorch or JAX implementation exists of any of the three
kernels**. Coherent dedispersion lives in dspsr (C++/CUDA), CDMT (98.7% CUDA), a Julia package
(CoherentDedispersion.jl) and a SYCL prototype (fxzjshm) --- so the claim is "pure-PyTorch /
pip-installable / autograd-native", NOT "first portable". SumThreshold/spectral-kurtosis GPU
kernels exist only CuPy-locked (`jess`, Kania et al. 2026, AJ 171, 73); IQRM and AOFlagger are
CPU. No pure-torch FFA exists (riptide is C++/CPU, its GPU issue closed unimplemented
2024-03-27; `gaffa`, a CUDA FFA scaffold, appeared 2026-06-12 --- cite as concurrent CUDA work,
claim only the torch niche). LPT fence from fable-ideas: the FFA targets classical pulsar
reprocessing, not LPT discovery.

Anchors (plan 43): a synthetic dispersed impulse must re-collapse to its impulse (round-trip);
torch SK/SumThreshold must match the tested CPU oracles in `jansky.rfi`; the FFA must recover
injected periods against a brute-force folding oracle. Real legs: a CHIME/FRB baseband event
(DOI 10.11570/23.0029, tied-beam voltages at 2.56 us) re-dedispersed to its catalogue DM, and
the vendored Crab filterbank for real-RFI masks + the FFA period re-find. Pure tensor ops
throughout: one ``device=`` argument covers CPU, CUDA, and ROCm (RX 7600 XT / gfx1102).

SumThreshold evaluation-order note (stated, not hidden): the `jansky.rfi` oracle flags windows
SEQUENTIALLY within a pass (flags set by earlier windows change later window means); the torch
performance path evaluates all windows of a pass in PARALLEL from the pass-start mask. Tests
pin the kernel maths by running the torch kernel in sequential mode (exact mask equality with
the oracle) and quantify the parallel/sequential agreement separately.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from jansky.constants import DM_CONST  # s MHz^2 (pc cm^-3)^-1

__all__ = [
    "chirp_group_delay_residual",
    "coherent_dedisperse",
    "dedisperse_channelized",
    "synthetic_dispersed_voltage",
    "spectral_kurtosis",
    "sumthreshold",
    "sumthreshold2d",
    "ffa_search",
    "fold_snr",
    "hardware_string",
    "load_chime_baseband",
    "run",
    "run_benchmark_only",
]

CHIME_BASEBAND_LOCAL = Path("data/FRB20181231C_24366209_beamformed.h5")
CHIME_BASEBAND_DOI = "10.11570/23.0029"
CHIME_CHAN_BW_MHZ = 400.0 / 1024.0  # 2.56 us complex sampling per channel


def _require_torch() -> Any:
    """Import torch lazily so the core package works without the ``fdmt`` extra."""
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "jansky_research.torchdsp needs PyTorch: install with `uv sync --extra fdmt` "
            "(CPU wheel; for GPU use a ROCm/CUDA venv)"
        ) from exc
    return torch


# --------------------------------------------------------------------------------------------
# Kernel A: coherent dedispersion
# --------------------------------------------------------------------------------------------


def _chirp_phase(dm: float, f_base_mhz: np.ndarray, f0_mhz: float) -> np.ndarray:
    r"""Dedispersion chirp phase (radians) at baseband offsets ``f_base_mhz`` from ``f0_mhz``.

    :math:`\phi(f) = 2\pi\,k_\mathrm{DM}\,\mathrm{DM}\,f^2 / (f_0^2 (f_0 + f))` (Hankins &
    Rickett 1975), with :math:`k_\mathrm{DM}` = ``jansky.constants.DM_CONST`` in
    s MHz\ :sup:`2`; the 1e6 converts s MHz to radians-per-:math:`2\pi`. SIGN CONVENTION:
    which sign of the kernel dedisperses depends on the backend's sideband/conjugation
    convention (dspsr parameterises exactly this); the sign used here is anchored
    EMPIRICALLY by the CHIME real leg (burst S/N maximised at the catalogue DM) --- other
    backends may need the conjugate. Computed in float64: in a narrow channelized-voltage
    channel the phase reaches ~1e4 rad (DM ~ 556 at 400 MHz, 0.39-MHz channel); for WIDEBAND
    (unchannelized) use it exceeds 1e6 rad and float32 accumulation would genuinely smear.
    """
    f = np.asarray(f_base_mhz, dtype=np.float64)
    return 2.0 * np.pi * 1.0e6 * DM_CONST * dm * f**2 / (f0_mhz**2 * (f0_mhz + f))


def coherent_dedisperse(
    voltages: Any, dm: float, f0_mhz: float, *, chan_bw_mhz: float, device: str = "cpu"
) -> Any:
    """Coherently dedisperse one channel of complex baseband voltages (pure torch).

    Removes the intra-channel dispersion smear exactly: FFT the voltage series, multiply by
    the conjugate dispersion chirp for centre frequency ``f0_mhz``, inverse FFT. The chirp is
    built in float64 and applied as complex128, then cast back --- phase error stays << 1 rad
    across the band even at CHIME DMs.
    """
    torch = _require_torch()
    v = torch.as_tensor(np.asarray(voltages), device=device)
    n = v.shape[-1]
    f_base = np.fft.fftfreq(n, d=1.0 / chan_bw_mhz)  # MHz offsets within the channel
    phase = _chirp_phase(dm, f_base, f0_mhz)
    kernel = torch.as_tensor(np.exp(1j * phase), device=device)  # complex128
    spec = torch.fft.fft(v.to(torch.complex128), dim=-1)
    return torch.fft.ifft(spec * kernel, dim=-1).to(torch.complex64)


def chirp_group_delay_residual(
    dm: float = 556.0, f0_mhz: float = 400.0, *, chan_bw_mhz: float = CHIME_CHAN_BW_MHZ
) -> dict:
    r"""Check the chirp against the cold-plasma dispersion law it claims to implement.

    The round-trip validation composes :func:`_chirp_phase` with its own negation, so it can
    never detect a wrong dispersion constant, frequency exponent, or band convention. This
    can: the chirp's numerical group delay :math:`\tau(f)=(1/2\pi)\,d\phi/df` must equal the
    relative cold-plasma delay :math:`k_\mathrm{DM}\,\mathrm{DM}\,(\nu^{-2}-f_0^{-2})` at
    every offset across the channel. Returns the max absolute residual in microseconds and
    in complex samples at ``chan_bw_mhz``.
    """
    f = np.linspace(-chan_bw_mhz / 2, chan_bw_mhz / 2, 4097)
    df = f[1] - f[0]
    phi = _chirp_phase(dm, f, f0_mhz)
    tau_num_us = np.gradient(phi, df) / (2.0 * np.pi)  # rad/MHz -> us
    tau_law_us = -1.0e6 * DM_CONST * dm * ((f0_mhz + f) ** -2 - f0_mhz**-2.0)
    resid_us = float(np.max(np.abs(tau_num_us - tau_law_us)[1:-1]))  # gradient edges excluded
    return {
        "max_resid_us": resid_us,
        "max_resid_samples": resid_us * chan_bw_mhz,  # 1/chan_bw us per complex sample
        "dm": dm,
        "f0_mhz": f0_mhz,
    }


def dedisperse_channelized(
    voltages: Any,
    dm: float,
    freqs_mhz: np.ndarray,
    *,
    chan_bw_mhz: float,
    ref_freq_mhz: float | None = None,
    conjugate: bool = False,
    device: str = "cpu",
) -> Any:
    """Coherent (intra-channel) + integer-sample (inter-channel) dedispersion of (nchan, ntime).

    Each channel is chirp-dedispersed at its own centre frequency, then rolled by the whole-
    sample inter-channel delay relative to ``ref_freq_mhz`` (default: the highest channel) ---
    the standard coherent-dedispersion-of-channelized-voltages ladder. BOUNDARY SEMANTICS
    (stated, not hidden): both the chirp multiply (circular convolution, no overlap-save
    discard region) and the inter-channel gather are CIRCULAR --- energy within a dispersion
    sweep of the buffer edge wraps around, unlike dspsr's discard-region handling. Fine for
    a burst captured mid-buffer (the shipped legs); a wideband streaming pipeline would need
    overlap-save on top. ``conjugate=True`` applies the opposite chirp sign (the other
    backend sideband convention) -- used by the real leg to measure, not assert, which sign
    dedisperses CHIME baseband.
    """
    torch = _require_torch()
    freqs = np.asarray(freqs_mhz, dtype=np.float64)
    if ref_freq_mhz is None:
        ref_freq_mhz = float(freqs.max())
    dt_s = 1.0e-6 / chan_bw_mhz  # complex sample interval
    v = torch.as_tensor(np.asarray(voltages), device=device)
    n = v.shape[-1]
    f_base = np.fft.fftfreq(n, d=1.0 / chan_bw_mhz)
    shifts = np.array(
        [int(round(DM_CONST * dm * (f0**-2 - ref_freq_mhz**-2) / dt_s)) for f0 in freqs]
    )
    idx = torch.arange(n, device=device)
    out = torch.empty(v.shape, dtype=torch.complex64, device=device)
    # batched per-chunk: ONE FFT call per chunk of channels (the per-channel host loop was the
    # GPU bottleneck -- the torch-fdmt lesson); chunking bounds the working set
    chunk = max(1, min(len(freqs), 16))
    fb = torch.as_tensor(f_base, device=device)  # float64, on device
    f0s = torch.as_tensor(freqs, device=device)
    for a in range(0, len(freqs), chunk):
        b = min(a + chunk, len(freqs))
        f0 = f0s[a:b, None]
        # chirp built ON DEVICE (host-side numpy exp of nchan x n complex128 was the bottleneck);
        # the phase is accumulated in float64 (~1e4 rad per channel at CHIME DMs; >1e6 rad in
        # wideband use), WRAPPED mod 2pi, and only then cast to float32 -- so the FFT/multiply
        # run in complex64, where consumer GPUs (1/32-rate f64) actually shine, at < 1e-6 rad
        # phase error
        phase = 2.0 * np.pi * 1.0e6 * DM_CONST * dm * fb[None, :] ** 2 / (f0**2 * (f0 + fb))
        if conjugate:
            phase = -phase
        phase32 = torch.remainder(phase, 2.0 * np.pi).to(torch.float32)
        kernel = torch.polar(torch.ones_like(phase32), phase32)  # complex64
        spec = torch.fft.fft(v[a:b].to(torch.complex64), dim=-1)
        ded = torch.fft.ifft(spec * kernel, dim=-1)
        sh = torch.as_tensor(shifts[a:b], device=device)[:, None]
        gather_idx = (idx[None, :] + sh) % n  # roll each channel by its own delay
        out[a:b] = torch.gather(ded, 1, gather_idx)
    return out


def synthetic_dispersed_voltage(
    *,
    n_time: int = 4096,
    dm: float = 100.0,
    f0_mhz: float = 600.0,
    chan_bw_mhz: float = CHIME_CHAN_BW_MHZ,
    noise: float = 0.05,
    seed: int = 0,
) -> dict:
    """One channel of complex voltage carrying an impulse DISPERSED with the exact inverse chirp.

    The recover-a-known: :func:`coherent_dedisperse` at the same DM must re-collapse the smear
    to (near) the original impulse --- peak power back within a few samples' width and >~90% of
    the impulse energy re-concentrated.
    """
    rng = np.random.default_rng(seed)
    v = np.zeros(n_time, dtype=np.complex128)
    v[n_time // 2] = 1.0 * n_time**0.5  # unit-energy impulse (spread across the band)
    spec = np.fft.fft(v)
    f_base = np.fft.fftfreq(n_time, d=1.0 / chan_bw_mhz)
    spec *= np.exp(-1j * _chirp_phase(dm, f_base, f0_mhz))  # inverse of the dedispersion kernel
    v_disp = np.fft.ifft(spec)
    v_disp += noise * (rng.standard_normal(n_time) + 1j * rng.standard_normal(n_time))
    return {
        "voltage": v_disp.astype(np.complex64),
        "dm": dm,
        "f0_mhz": f0_mhz,
        "chan_bw_mhz": chan_bw_mhz,
        "impulse_index": n_time // 2,
    }


# --------------------------------------------------------------------------------------------
# Kernel B: RFI excision (spectral kurtosis + SumThreshold)
# --------------------------------------------------------------------------------------------


def spectral_kurtosis(power: Any, *, axis: int = 0, device: str = "cpu") -> Any:
    """Generalised SK estimator, pure torch; same maths as ``jansky.rfi.spectral_kurtosis``.

    SK = (M+1)/(M-1) (M S2/S1^2 - 1) over the ``M`` samples along ``axis``; ~1 for clean
    Gaussian noise, <1 for steady (CW) interference, >1 for spiky interference (Nita & Gary
    2010). Tested for exact agreement with the CPU oracle.
    """
    torch = _require_torch()
    p = torch.as_tensor(np.asarray(power, dtype=np.float64), device=device)
    m = p.shape[axis]
    s1 = p.sum(dim=axis)
    s2 = (p**2).sum(dim=axis)
    return (m + 1) / (m - 1) * (m * s2 / s1**2 - 1)


def _mad_sigma_np(x: np.ndarray) -> float:
    # numpy medians (even-length = midpoint average), matching jansky.rfi.mad_sigma exactly;
    # torch.median returns the LOWER middle element and would break oracle parity
    med = np.median(x)
    return float(1.4826 * np.median(np.abs(x - med)))


def sumthreshold(
    values: Any,
    *,
    max_window: int = 8,
    threshold: float = 3.5,
    rho: float = 1.5,
    sigma: float | None = None,
    mask: Any = None,
    sequential: bool = False,
    device: str = "cpu",
) -> np.ndarray:
    """SumThreshold RFI flagger (Offringa et al. 2010), pure torch.

    Same parameters and thresholds as the tested CPU oracle ``jansky.rfi.sumthreshold``
    (per-window threshold chi_M = threshold*sigma / rho^log2(M); median-removed residual;
    robust sigma from the unflagged samples). Two evaluation modes:

    - ``sequential=True``: windows are processed in order within each pass and flags take
      effect immediately --- byte-identical to the CPU oracle (used by the tests to pin the
      kernel maths; O(n) python loop, not the performance path).
    - ``sequential=False`` (default): all windows of a pass are evaluated in parallel from
      the pass-start mask via cumulative sums, and their flags are unioned at pass end ---
      the tensor-friendly variant (the same design choice as GPU implementations of the
      algorithm). Agreement with the sequential oracle is quantified, not assumed.
    """
    torch = _require_torch()
    resid_np = np.asarray(values, dtype=np.float64).ravel()
    resid_np = resid_np - np.median(resid_np)
    out_np = (
        np.zeros(resid_np.size, dtype=bool)
        if mask is None
        else np.asarray(mask, dtype=bool).ravel().copy()
    )
    if sigma is not None:
        s = float(sigma)
    else:
        unflagged = resid_np[~out_np]
        s = _mad_sigma_np(unflagged) if unflagged.size else 0.0
    resid = torch.as_tensor(resid_np, device=device)
    out = torch.as_tensor(out_np, device=device)
    n = resid.shape[0]
    if s <= 0:
        return out.cpu().numpy()
    m = 1
    while m <= max_window:
        chi = threshold * s / (rho ** np.log2(m))
        if m == 1:
            out |= resid.abs() > chi
        elif sequential:
            for j in range(n - m + 1):
                wm = out[j : j + m]
                if bool(wm.all()):
                    continue
                if float(resid[j : j + m][~wm].mean().abs()) > chi:
                    out[j : j + m] = True
        else:
            # parallel pass: window sums of unflagged residuals / unflagged counts via cumsum
            r0 = torch.where(out, torch.zeros((), dtype=resid.dtype, device=device), resid)
            c0 = (~out).to(resid.dtype)
            rc = torch.cat([torch.zeros(1, dtype=resid.dtype, device=device), r0.cumsum(0)])
            cc = torch.cat([torch.zeros(1, dtype=resid.dtype, device=device), c0.cumsum(0)])
            wsum = rc[m:] - rc[:-m]
            wcnt = cc[m:] - cc[:-m]
            hit = (wcnt > 0) & (wsum.abs() / wcnt.clamp(min=1.0) > chi)
            if bool(hit.any()):
                # union of all hit windows: mark [j, j+m) for each hit start j
                starts = torch.nonzero(hit).ravel()
                delta = torch.zeros(n + 1, dtype=torch.int32, device=device)
                delta.index_add_(
                    0, starts, torch.ones(starts.numel(), dtype=torch.int32, device=device)
                )
                ends = (starts + m).clamp(max=n)
                delta.index_add_(
                    0, ends, -torch.ones(starts.numel(), dtype=torch.int32, device=device)
                )
                out |= delta[:-1].cumsum(0) > 0
        m *= 2
    return out.cpu().numpy()


def sumthreshold2d(
    dynspec: Any,
    *,
    max_window: int = 8,
    threshold: float = 3.5,
    rho: float = 1.5,
    n_iter: int = 1,
    sequential: bool = False,
    device: str = "cpu",
) -> np.ndarray:
    """2-D SumThreshold over ``(n_time, n_chan)``; mirrors ``jansky.rfi.sumthreshold2d``.

    Time pass per channel, frequency pass per time sample, the accumulating mask threaded
    through, ``n_iter`` outer iterations --- the same structure as the CPU oracle so that
    ``sequential=True`` reproduces it exactly.
    """
    data = np.asarray(dynspec, dtype=np.float64)
    mask = np.zeros(data.shape, dtype=bool)
    for _ in range(max(1, n_iter)):
        seed_mask = mask
        acc = seed_mask.copy()
        for c in range(data.shape[1]):
            acc[:, c] |= sumthreshold(
                data[:, c],
                max_window=max_window,
                threshold=threshold,
                rho=rho,
                mask=seed_mask[:, c],
                sequential=sequential,
                device=device,
            )
        for t in range(data.shape[0]):
            acc[t, :] |= sumthreshold(
                data[t, :],
                max_window=max_window,
                threshold=threshold,
                rho=rho,
                mask=seed_mask[t, :],
                sequential=sequential,
                device=device,
            )
        mask = acc
    return mask


# --------------------------------------------------------------------------------------------
# Kernel C: FFA (stretch goal)
# --------------------------------------------------------------------------------------------


def _ffa_transform(rows: Any) -> Any:
    """Radix-2 FFA merge of ``(m, p)`` folded rows -> ``(m, p)`` drift profiles (pure torch).

    Output row ``j`` is the fold with a total phase drift of ``j`` samples across the ``m``
    rows (candidate period ``p + j/(m-1)`` samples). Classic Staelin (1969) recursion,
    vectorised per stage with gather (the same pattern as `torch-fdmt`'s merges).
    """
    torch = _require_torch()
    x = rows
    m, p = x.shape
    idx = torch.arange(p, device=x.device)
    g = 1
    while g < m:
        x = x.reshape(m // (2 * g), 2, g, p)
        top, bot = x[:, 0], x[:, 1]  # (n_groups, g, p)
        j = torch.arange(2 * g, device=x.device)
        h = j // 2  # both halves use drift-h profiles (Staelin/Morello recursion)
        shift = h + (j & 1)  # total output drift = h (top) + h (bottom) + (j&1) = j
        gather_idx = (idx[None, :] + shift[:, None]) % p  # (2g, p)
        merged = top[:, h, :] + torch.gather(
            bot[:, h, :], 2, gather_idx[None, :, :].expand(x.shape[0], -1, -1)
        )
        x = merged
        g *= 2
        x = x.reshape(m // g, g, p) if g < m else x
    return x.reshape(m, p)


def _profile_snr(profiles: Any) -> Any:
    """(max - median) / MAD-sigma per profile row --- the same S/N used by the brute oracle."""
    torch = _require_torch()
    med = profiles.median(dim=-1, keepdim=True).values
    mad = (profiles - med).abs().median(dim=-1).values * 1.4826
    return (profiles.max(dim=-1).values - med.squeeze(-1)) / torch.clamp(mad, min=1e-12)


def fold_snr(series: Any, period_samples: float, *, device: str = "cpu") -> float:
    """Brute-force folding oracle: fold at one (possibly fractional) trial period, return S/N."""
    torch = _require_torch()
    x = torch.as_tensor(np.asarray(series, dtype=np.float32), device=device)
    n = x.shape[0]
    nbin = int(np.floor(period_samples))
    phase = torch.arange(n, device=device, dtype=torch.float64) % period_samples / period_samples
    bins = (phase * nbin).long().clamp(max=nbin - 1)
    prof = torch.zeros(nbin, dtype=torch.float32, device=device)
    cnt = torch.zeros(nbin, dtype=torch.float32, device=device)
    prof.index_add_(0, bins, x)
    cnt.index_add_(0, bins, torch.ones_like(x))
    prof = prof / cnt.clamp(min=1.0)
    return float(_profile_snr(prof[None, :])[0])


def ffa_search(
    series: Any,
    *,
    pmin_samples: int,
    pmax_samples: int,
    device: str = "cpu",
) -> dict:
    """Pure-torch FFA period search over integer base periods [pmin, pmax] samples.

    For each base period p the series is folded into the largest power-of-two number of rows,
    FFA-merged, and scored with :func:`_profile_snr`; candidate periods are p + j/(m-1)
    samples. Returns the best period (samples), its S/N, and the per-base-period best S/N
    curve. Honest scope (the `torch-fdmt` lesson): the python loop is over base periods only;
    each transform is fully vectorised tensor ops.
    """
    torch = _require_torch()
    x = torch.as_tensor(np.asarray(series, dtype=np.float32), device=device)
    x = x - x.median()
    n = x.shape[0]
    best: dict[str, Any] = {"period_samples": float("nan"), "snr": -np.inf}
    curve = []
    for p in range(int(pmin_samples), int(pmax_samples) + 1):
        m = n // p
        if m < 2:
            break
        m2 = 1 << (m.bit_length() - 1)  # largest power of two <= m
        rows = x[: m2 * p].reshape(m2, p)
        profiles = _ffa_transform(rows)
        snr = _profile_snr(profiles)
        k = int(torch.argmax(snr))
        s = float(snr[k])
        period = p + (k / (m2 - 1) if m2 > 1 else 0.0)
        curve.append((p, s))
        if s > best["snr"]:
            best = {"period_samples": float(period), "snr": s, "n_rows": m2, "profile_j": k}
    best["curve"] = np.asarray(curve)
    return best


def benchmark(device: str = "cpu", *, n_chan: int = 64, n_time: int = 1 << 20) -> dict:
    """Wall-clock the three kernels at survey-ish sizes on ``device`` (honest, single-run).

    Portability first, speed second (the `torch-fdmt` framing): the point is that the SAME
    code runs on CPU/CUDA/ROCm, not a shoot-out against tuned CUDA kernels on datacentre
    hardware. Sizes: ``n_chan`` channels x ``n_time`` complex samples for the chirp;
    (8192 x 256) for SumThreshold; a 2^22-sample series over 64 base periods for the FFA.
    """
    import time

    torch = _require_torch()
    rng = np.random.default_rng(0)

    v = (rng.standard_normal((n_chan, n_time)) + 1j * rng.standard_normal((n_chan, n_time))).astype(
        np.complex64
    )
    freqs = np.linspace(500.0, 700.0, n_chan)
    t0 = time.perf_counter()
    out = dedisperse_channelized(v, 100.0, freqs, chan_bw_mhz=CHIME_CHAN_BW_MHZ, device=device)
    if device != "cpu":  # pragma: no cover - GPU only
        torch.cuda.synchronize()
    t_chirp = time.perf_counter() - t0
    del out

    dyn = rng.standard_normal((8192, 256))
    dyn[:, 100] += 6.0
    t0 = time.perf_counter()
    sumthreshold2d(dyn, device=device)
    t_st = time.perf_counter() - t0

    ts = rng.standard_normal(1 << 22).astype(np.float32)
    t0 = time.perf_counter()
    ffa_search(ts, pmin_samples=200, pmax_samples=263, device=device)
    if device != "cpu":  # pragma: no cover - GPU only
        torch.cuda.synchronize()
    t_ffa = time.perf_counter() - t0

    return {
        "device": device,
        "chirp_s": round(t_chirp, 2),
        "chirp_msamples": round(n_chan * n_time / 1e6, 1),
        "chirp_n_chan": n_chan,
        "chirp_n_time": n_time,
        "sumthreshold_s": round(t_st, 2),
        "st_shape": [8192, 256],
        "ffa_s": round(t_ffa, 2),
        "ffa_n_samples": 1 << 22,
        "ffa_n_base_periods": 64,
    }


def hardware_string(device: str = "cpu") -> str:
    """Hardware + torch build for one device, from torch introspection -- never typed by hand.

    A hand-entered hardware string is inherited by future runs through the results merge and
    misattributes them (the referee-caught failure on both this slice and ``torchfdmt``).
    """
    try:
        import platform

        torch = _require_torch()

        if device != "cpu" and torch.cuda.is_available():
            return f"{torch.cuda.get_device_name(0)}, torch {torch.__version__}"
        return f"{platform.processor() or platform.machine()}, torch {torch.__version__}"
    except Exception:  # pragma: no cover - introspection is best-effort
        return "unknown"


def run_benchmark_only(out: str = ".", *, device: str = "cpu") -> dict:  # pragma: no cover
    """One single-session benchmark producing BOTH timing columns, merged into the real file.

    This is the command that makes the committed benchmark table producible: the GPU and CPU
    legs run back-to-back in the same interpreter, the hardware strings for both come from
    introspection, and `preserve_live_results` merges the block into the science metrics with
    a ``_merge`` record instead of a hand-edit. On a GPU device it also re-runs the offline
    kernel checks on that device and commits the cross-device comparison -- the portability
    claim measured. Touches no science key.
    """
    torch = _require_torch()

    if device != "cpu":
        benchmark(device, n_chan=4, n_time=1 << 16)  # warm-up: first call pays GPU init
    bench_dev = benchmark(device)
    bench_cpu = benchmark("cpu") if device != "cpu" else bench_dev
    payload: dict = {
        "is_real": True,  # merges as real without touching the science `source` string
        "benchmark": bench_dev,
        "benchmark_cpu": bench_cpu,
        "benchmark_device": device,
        "benchmark_hardware": hardware_string(device),
        "benchmark_hardware_cpu": hardware_string("cpu"),
        "benchmark_torch": str(torch.__version__),
        "benchmark_single_session": True,
    }
    if device != "cpu":
        gpu_checks = {
            k: v for k, v in _kernel_checks(device).items() if k not in ("ffa_curve", "device")
        }
        payload["cross_device_check"] = {"device": device, **gpu_checks}

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    merged = write_results(payload, op / "results" / "torchdsp_metrics.json")
    _write_macros(merged, op / "papers" / "torchdsp" / "generated" / "macros.tex")
    return merged


# --------------------------------------------------------------------------------------------
# Real-data loaders
# --------------------------------------------------------------------------------------------


def load_chime_baseband(
    path: str | Path = CHIME_BASEBAND_LOCAL,
) -> dict:  # pragma: no cover - needs the 150 MB CANFAR file
    """Load a CHIME/FRB tied-beam baseband release file (DOI 10.11570/23.0029).

    Returns the good-channel complex voltages (nchan, npol, ntime), channel centre
    frequencies (MHz), and the per-channel start-time alignment in 2.56-us frames
    (channels are rolled to a common start).
    """
    import h5py

    with h5py.File(path, "r") as f:
        bb = f["tiedbeam_baseband"][...]
        freq_map = f["index_map/freq"][...]
        good_ids = f["index_map/good_channels"][...]
        fpga0 = f["time0"]["fpga_count"][...].astype(np.int64)
    freqs = freq_map["centre"]
    # good_channels lists GLOBAL channel ids (0-1023); map them onto the saved channel axis
    good = np.isin(freq_map["id"], good_ids)
    bb = bb[good]
    freqs = freqs[good]
    fpga0 = fpga0[good]
    off = fpga0 - fpga0.min()
    aligned = np.zeros_like(bb)
    for k in range(bb.shape[0]):  # roll each channel to the common start
        aligned[k] = np.roll(bb[k], int(off[k]), axis=-1)
    aligned = np.nan_to_num(aligned)
    return {"voltages": aligned, "freqs_mhz": freqs, "offsets_frames": off}


# --------------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------------


def _dedisp_recovery_metrics(device: str) -> dict:
    syn = synthetic_dispersed_voltage()
    n = syn["voltage"].size
    power_before = np.abs(syn["voltage"]) ** 2

    def _peak_and_frac(power_after: np.ndarray) -> tuple[int, float]:
        peak = int(np.argmax(power_after))
        w = 4  # energy within +-4 samples of the impulse, as a fraction of total
        frac = float(
            power_after[max(0, peak - w) : peak + w + 1].sum() / max(power_after.sum(), 1e-12)
        )
        return peak, frac

    v_fixed = coherent_dedisperse(
        syn["voltage"], syn["dm"], syn["f0_mhz"], chan_bw_mhz=syn["chan_bw_mhz"], device=device
    )
    peak, frac = _peak_and_frac((v_fixed.abs() ** 2).cpu().numpy())
    # The same round trip through the SHIPPED path: every science and benchmark result uses
    # dedisperse_channelized (phase wrapped mod 2pi, complex64), not coherent_dedisperse
    # (complex128), so the validated path must be the shipped one.
    v_chan = dedisperse_channelized(
        syn["voltage"][None, :],
        syn["dm"],
        np.array([syn["f0_mhz"]]),
        chan_bw_mhz=syn["chan_bw_mhz"],
        device=device,
    )
    peak_c, frac_c = _peak_and_frac((v_chan[0].abs() ** 2).cpu().numpy())
    gd = chirp_group_delay_residual(dm=syn["dm"], f0_mhz=syn["f0_mhz"])
    return {
        "peak_offset_samples": abs(peak - syn["impulse_index"]),
        "reconcentrated_energy_frac": round(frac, 3),
        "peak_offset_samples_channelized": abs(peak_c - syn["impulse_index"]),
        "reconcentrated_energy_frac_channelized": round(frac_c, 3),
        # The chirp vs the cold-plasma law itself -- the check the self-composing round trip
        # cannot perform (it would pass with a wrong dispersion constant).
        "group_delay_max_resid_samples": float(f"{gd['max_resid_samples']:.2e}"),
        "dispersed_peak_frac": round(
            float(np.sort(power_before)[-9:].sum() / max(power_before.sum(), 1e-12)),
            3,
        ),
        "n_time": n,
    }


def _mask_agreement(mask_par: np.ndarray, mask_cpu: np.ndarray) -> dict:
    """Jaccard plus its DIRECTION: Jaccard is symmetric, but over-flagging (harmless data
    loss) and under-flagging (missed RFI) are operationally different for a user."""
    inter = int(np.logical_and(mask_par, mask_cpu).sum())
    union = int(np.logical_or(mask_par, mask_cpu).sum())
    return {
        "jaccard": round(float(inter / union) if union else 1.0, 4),
        "n_par_only": int(np.logical_and(mask_par, ~mask_cpu).sum()),  # over-flagged
        "n_cpu_only": int(np.logical_and(~mask_par, mask_cpu).sum()),  # under-flagged
    }


def _st_agreement_sweep(dyn: np.ndarray, device: str) -> dict:
    """Parallel-vs-sequential SumThreshold agreement across the settings that grow divergence.

    The divergence is generated by flags taking effect mid-pass, so it must grow with flag
    occupancy and with ``n_iter``; the shipped default (threshold 3.5, n_iter 1) is the
    configuration in which it is smallest. This sweep bounds it elsewhere.
    """
    from jansky import rfi as rfi_cpu

    out = {}
    for thr in (3.0, 3.5, 4.0):
        for n_iter in (1, 2):
            mask_cpu = rfi_cpu.sumthreshold2d(dyn, threshold=thr, n_iter=n_iter)
            mask_par = sumthreshold2d(dyn, threshold=thr, n_iter=n_iter, device=device)
            out[f"thr{thr}_iter{n_iter}"] = _mask_agreement(mask_par, mask_cpu)
    return out


def _kernel_checks(device: str) -> dict:
    """The three per-kernel recover-a-knowns + oracle matches (device-parameterised).

    Factored out of :func:`run` so a GPU session can run the identical checks and commit a
    cross-device comparison -- the paper's portability claim measured, not asserted.
    """
    from jansky import rfi as rfi_cpu

    rng = np.random.default_rng(0)

    # Kernel A: dispersed-impulse round trip (both code paths) + group-delay law check
    dedisp = _dedisp_recovery_metrics(device)

    # Kernel B: torch SK + SumThreshold vs the jansky.rfi CPU oracles on synthetic RFI
    n_t, n_c = 512, 64
    dyn = rng.standard_normal((n_t, n_c))
    # a full-length constant line vanishes in the (median-removed) time pass by construction;
    # it is the FREQUENCY pass that must catch it, per time sample -- inject clearly above chi_1
    dyn[:, 17] += 6.0  # narrowband CW line
    dyn[300:308, :] += 3.0  # broadband burst
    sk_t = spectral_kurtosis(dyn**2, axis=0, device=device)
    sk_cpu = rfi_cpu.spectral_kurtosis(dyn**2, axis=0)
    sk_max_diff = float(np.max(np.abs(np.asarray(sk_t.cpu()) - sk_cpu)))
    mask_seq = sumthreshold2d(dyn, sequential=True, device=device)
    mask_cpu = rfi_cpu.sumthreshold2d(dyn)
    mask_par = sumthreshold2d(dyn, device=device)
    seq_equal = bool(np.array_equal(mask_seq, mask_cpu))
    agree = _mask_agreement(mask_par, mask_cpu)
    line_caught = bool(mask_par[:, 17].mean() > 0.9)
    burst_caught = bool(mask_par[300:308, :].mean() > 0.9)
    st_sweep = _st_agreement_sweep(dyn, device)

    # Kernel C: FFA recovers an injected period; brute-fold oracle agrees
    p_true = 233.7  # samples
    n = 1 << 16
    ts = rng.standard_normal(n).astype(np.float32)
    pulse_idx = (np.arange(0, n, p_true)).astype(int)
    ts[pulse_idx[pulse_idx < n]] += 6.0
    ffa = ffa_search(ts, pmin_samples=200, pmax_samples=270, device=device)
    oracle_snr = fold_snr(ts, p_true, device=device)
    ffa_err = abs(ffa["period_samples"] - p_true)
    m2 = int(ffa.get("n_rows", 2))
    m_full = n // int(p_true)
    return {
        "device": device,
        "dedisp": dedisp,
        "sk_max_diff": sk_max_diff,
        "sumthreshold_sequential_equals_oracle": seq_equal,
        "sumthreshold_parallel_jaccard": agree["jaccard"],
        "sumthreshold_par_only_flags": agree["n_par_only"],
        "sumthreshold_cpu_only_flags": agree["n_cpu_only"],
        "sumthreshold_agreement_sweep": st_sweep,
        "rfi_line_caught": line_caught,
        "rfi_burst_caught": burst_caught,
        "ffa_period_true": p_true,
        "ffa_period_found": round(ffa["period_samples"], 2),
        "ffa_period_err_samples": round(ffa_err, 2),
        # the injected 233.7 is NOT on the drift grid: candidates are p + k/(m2-1), so the
        # honest statement is "to the nearest grid point", never "exactly"
        "ffa_drift_resolution_samples": round(1.0 / max(m2 - 1, 1), 4),
        "ffa_snr": round(ffa["snr"], 1),
        "fold_oracle_snr": round(oracle_snr, 1),
        # the FFA/oracle S/N gap, predicted (power-of-two row truncation sqrt(m2/m) times
        # ~2-bin smearing 1/sqrt(2)) vs measured -- an assertion in the paper becomes a check
        "ffa_snr_ratio_predicted": round(float(np.sqrt(m2 / m_full) / np.sqrt(2.0)), 2),
        "ffa_snr_ratio_measured": round(float(ffa["snr"] / max(oracle_snr, 1e-9)), 2),
        "ffa_curve": ffa["curve"],  # array; stays out of the JSON (json-safe drops it)
    }


def run(out: str = ".", *, offline: bool = True, device: str = "cpu", bench: bool = False) -> dict:
    """Offline: per-kernel recover-a-knowns + oracle matches; real: CHIME + Crab legs."""

    metrics: dict = {
        "source": "synthetic per-kernel recover-a-knowns",
        "is_real": not offline,
    }
    metrics.update(_kernel_checks(device))

    if not offline:  # pragma: no cover - needs the local CANFAR file + network Crab fetch
        metrics.update(_real_legs(device))
        metrics["source"] = f"CHIME baseband (DOI {CHIME_BASEBAND_DOI}) + Parkes Crab filterbank"
    if bench:  # pragma: no cover - wall-clock, not run in CI
        if device != "cpu":
            benchmark(device, n_chan=4, n_time=1 << 16)  # warm-up: first call pays GPU init
        metrics["benchmark"] = benchmark(device)
        if device != "cpu":
            metrics["benchmark_cpu"] = benchmark("cpu")  # same venv, same code, CPU reference
            metrics["benchmark_hardware_cpu"] = hardware_string("cpu")
        metrics["benchmark_device"] = device
        metrics["benchmark_hardware"] = hardware_string(device)
        metrics["benchmark_single_session"] = True

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    json_metrics = {k: v for k, v in metrics.items() if k != "ffa_curve"}
    from .report import write_results

    write_results(json_metrics, op / "results" / "torchdsp_metrics.json")
    _figure(metrics, op / "papers" / "torchdsp" / "figures")
    _write_macros(metrics, op / "papers" / "torchdsp" / "generated" / "macros.tex")
    return metrics


def _real_legs(device: str) -> dict:  # pragma: no cover - real data + network
    """CHIME baseband re-dedispersion + Crab filterbank RFI/FFA legs."""
    import csv

    torch = _require_torch()

    out: dict = {}
    # --- CHIME: coherent dedispersion sharpens the burst; S/N peaks at the catalogue DM
    bb = load_chime_baseband()
    with open("data/chimefrbcat2.csv") as fh:
        row = next(r for r in csv.DictReader(fh) if r["tns_name"].startswith("FRB20181231C"))
    dm_cat = float(row["dm_fitb"])
    v0 = bb["voltages"][:, 0, :]
    v1 = bb["voltages"][:, 1, :]
    freqs = bb["freqs_mhz"]

    def burst_snr(dm: float, *, conjugate: bool = False) -> tuple[float, int]:
        # dedisperse each polarization coherently, sum the powers (total intensity)
        d0 = dedisperse_channelized(
            v0, dm, freqs, chan_bw_mhz=CHIME_CHAN_BW_MHZ, conjugate=conjugate, device=device
        )
        d1 = dedisperse_channelized(
            v1, dm, freqs, chan_bw_mhz=CHIME_CHAN_BW_MHZ, conjugate=conjugate, device=device
        )
        power = (d0.abs() ** 2 + d1.abs() ** 2).sum(dim=0)
        # 100-sample (0.256 ms) boxcar S/N; the argmax is recorded so the committed evidence
        # can distinguish "the burst sharpens at the catalogue DM" from "the largest noise
        # excursion happened to be biggest at that trial" (the boxcar max is a search)
        box = torch.nn.functional.avg_pool1d(power[None, None, :], 100, stride=1)[0, 0]
        med = box.median()
        mad = (box - med).abs().median() * 1.4826
        snr = float((box.max() - med) / torch.clamp(mad, min=1e-12))
        return snr, int(torch.argmax(box))

    dms = [0.0, dm_cat - 20.0, dm_cat - 5.0, dm_cat, dm_cat + 5.0, dm_cat + 20.0]
    snrs: dict[str, float] = {}
    peaks: dict[str, int] = {}
    for d in dms:
        s_val, arg = burst_snr(d)
        snrs[f"{d:.1f}"] = round(s_val, 1)
        peaks[f"{d:.1f}"] = arg
    out["chime_dm_catalogue"] = round(dm_cat, 2)
    out["chime_snr_vs_dm"] = snrs
    out["chime_peak_sample_vs_dm"] = peaks
    out["chime_peaks_at_catalogue_dm"] = bool(
        snrs[f"{dm_cat:.1f}"] >= max(v2 for k2, v2 in snrs.items())
    )
    # The chirp-sign anchor, measured rather than asserted: the conjugate kernel at the same
    # DM must NOT sharpen the burst, or the "anchored empirically" sentence has no arm.
    s_conj, _ = burst_snr(dm_cat, conjugate=True)
    out["chime_snr_conjugate_sign"] = round(s_conj, 1)
    # descriptors, emitted rather than hand-typed in prose
    out["chime_n_chan_good"] = int(len(freqs))
    out["chime_bandwidth_mhz"] = round(float(len(freqs) * CHIME_CHAN_BW_MHZ), 1)
    out["chime_n_time"] = int(v0.shape[-1])

    # --- Crab filterbank: real-RFI masks + FFA period re-find (dedispersed at the Crab DM)
    import urllib.request

    from jansky import rfi as rfi_cpu

    from .fdmt import delay_samples
    from .singlepulse import CRAB_DM, CRAB_FIL_URL, read_sigproc

    fil = Path("data/parkes_crab.fil")
    if not fil.exists():
        urllib.request.urlretrieve(CRAB_FIL_URL, fil)
    dyn, rfreqs, hdr = read_sigproc(fil)  # (n_time, n_chan), MHz, header
    dt = float(hdr["tsamp"])
    # descriptors: the paper's "2.1 s" was hand-typed; emit what the file actually holds
    out["crab_tsamp_us"] = round(dt * 1e6, 3)
    out["crab_n_samples"] = int(dyn.shape[0])
    out["crab_n_chan"] = int(dyn.shape[1])
    out["crab_duration_s"] = round(dyn.shape[0] * dt, 3)
    mask_cpu = rfi_cpu.sumthreshold2d(dyn[:2048])
    mask_par = sumthreshold2d(dyn[:2048], device=device)
    agree = _mask_agreement(mask_par, mask_cpu)
    out["crab_mask_jaccard"] = agree["jaccard"]
    out["crab_mask_par_only_flags"] = agree["n_par_only"]
    out["crab_mask_cpu_only_flags"] = agree["n_cpu_only"]
    f_top = float(np.max(rfreqs))
    ded = np.zeros(dyn.shape[0], dtype=np.float32)
    for c, fc in enumerate(rfreqs):  # incoherent dedispersion at the Crab DM, then sum
        ded += np.roll(dyn[:, c], -delay_samples(CRAB_DM, float(fc), f_top, dt))
    p_crab = 0.0337  # s, approximate topocentric Crab period
    pmin = int(0.9 * p_crab / dt)
    pmax = int(1.1 * p_crab / dt)
    ffa = ffa_search(ded, pmin_samples=pmin, pmax_samples=pmax, device=device)
    out["crab_period_found_ms"] = round(ffa["period_samples"] * dt * 1e3, 3)
    out["crab_ffa_snr"] = round(ffa["snr"], 1)
    out["crab_period_published_ms"] = 33.7
    # honesty diagnostics: the FFA's period resolution at this data length, and the
    # brute-fold oracle at the published period -- if that too is marginal, the file supports
    # only a marginal periodicity detection and the synthetics carry the algorithm validation
    m2 = ffa.get("n_rows", 2)
    out["crab_ffa_resolution_ms"] = round((p_crab / dt) / max(m2 - 1, 1) * dt * 1e3, 2)
    out["crab_fold_snr_at_published"] = round(fold_snr(ded, p_crab / dt, device=device), 1)
    # The null's threshold, measured: the same search volume on pure noise. Without this,
    # "reported as not significant" rests on a threshold that exists nowhere in the evidence.
    rng = np.random.default_rng(1)
    noise = rng.standard_normal(ded.shape[0]).astype(np.float32)
    ffa_noise = ffa_search(noise, pmin_samples=pmin, pmax_samples=pmax, device=device)
    out["crab_noise_ffa_snr"] = round(ffa_noise["snr"], 1)
    # The null's selection function, measured: inject a pulse train at the published period
    # into the REAL dedispersed series at a ladder of amplitudes (units of the series' robust
    # sigma) and record the weakest the FFA recovers -- that is what the 2.1-s file could
    # and could not have seen.
    sigma = float(1.4826 * np.median(np.abs(ded - np.median(ded))))
    p_samples = p_crab / dt
    res_ms = out["crab_ffa_resolution_ms"]
    ladder: dict[str, dict] = {}
    min_amp = None
    for amp in (0.1, 0.2, 0.5, 1.0, 2.0, 5.0):
        inj = ded.copy()
        idx = np.arange(0, inj.shape[0], p_samples).astype(int)
        inj[idx[idx < inj.shape[0]]] += amp * sigma
        f_inj = ffa_search(inj, pmin_samples=pmin, pmax_samples=pmax, device=device)
        period_ms = f_inj["period_samples"] * dt * 1e3
        recovered = bool(
            abs(period_ms - out["crab_period_published_ms"]) < 2.0 * res_ms
            and f_inj["snr"] > out["crab_noise_ffa_snr"] + 2.0
        )
        ladder[str(amp)] = {
            "period_ms": round(period_ms, 3),
            "snr": round(f_inj["snr"], 1),
            "recovered": recovered,
        }
        if recovered and min_amp is None:
            min_amp = amp
    out["crab_injection_ladder_sigma"] = ladder
    out["crab_injection_min_amp_sigma"] = min_amp
    return out


def _figure(m: dict, out_dir: str | Path) -> None:
    try:
        from .report import _agg
    except ImportError:  # pragma: no cover - the ROCm bench venv has no matplotlib
        return

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    syn = synthetic_dispersed_voltage()
    v_fixed = coherent_dedisperse(
        syn["voltage"], syn["dm"], syn["f0_mhz"], chan_bw_mhz=syn["chan_bw_mhz"]
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax1.plot(np.abs(syn["voltage"]) ** 2, lw=0.5, label="dispersed")
    ax1.plot(np.asarray((v_fixed.abs() ** 2).cpu()), lw=0.8, label="coherently dedispersed")
    ax1.set(xlabel="sample", ylabel="power", title="Chirp round-trip")
    ax1.legend(fontsize=8)
    curve = m.get("ffa_curve")  # the run's own periodogram, not a re-simulation
    if curve is not None:
        ax2.plot(curve[:, 0], curve[:, 1], lw=0.8)
        ax2.axvline(m["ffa_period_true"], color="C3", ls="--", label="injected period")
        ax2.set(xlabel="base period (samples)", ylabel="best S/N", title="FFA periodogram")
        ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "torchdsp.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    def g_sci(key: str, src: dict | None = None) -> str:
        """A small float as a LaTeX math body, e.g. 1.5e-14 -> ``1.5\\times10^{-14}``."""
        v = (src if src is not None else m).get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        mant, exp = f"{float(v):.1e}".split("e")
        return rf"{mant}\times10^{{{int(exp)}}}"

    pref = "tdReal" if m.get("is_real") else "tdSyn"
    dedisp = m.get("dedisp") or {}
    sweep = m.get("sumthreshold_agreement_sweep") or {}
    sweep_j = [v.get("jaccard") for v in sweep.values() if v.get("jaccard") is not None]
    m = dict(m)
    if sweep_j:
        m["st_sweep_jaccard_lo"] = min(sweep_j)
        m["st_sweep_jaccard_hi"] = max(sweep_j)
    snrs = m.get("chime_snr_vs_dm") or {}
    cat = m.get("chime_dm_catalogue")
    off = [v for k, v in snrs.items() if k != "0.0" and (cat is None or k != f"{cat:.1f}")]
    if off:
        m["chime_off_dm_max_snr"] = max(off)
    lines = [
        "% Auto-generated by jansky_research.torchdsp._write_macros -- do not edit.",
        "% Synthetic (tdSyn*) and real (tdReal*) namespaces are both always emitted. The",
        "% prefix labels the RUN MODE that produced the value, not the quantity's own",
        "% provenance: kernel-oracle checks are synthetic measurements in either mode, and",
        "% the prose must say so where it cites them. The inactive namespace is emitted as",
        "% placeholders; preserve_live_macros then keeps any real value already on disk and",
        "% refuses synthetic-over-real downgrades (it does NOT reset tdReal* on a rebuild).",
        rf"\newcommand{{\tdSource}}{{{m.get('source', '--')}}}",
        rf"\newcommand{{\tdDevice}}{{{m.get('device', '--')}}}",
        # benchmark provenance: introspected, never typed (underscores TeX-escaped)
        rf"\newcommand{{\tdBenchHwGpu}}{{{g('benchmark_hardware').replace('_', chr(92) + '_')}}}",
        rf"\newcommand{{\tdBenchHwCpu}}{{{g('benchmark_hardware_cpu').replace('_', chr(92) + '_')}}}",
    ]
    keys = (
        ("SkMaxDiff", "sk_max_diff"),
        ("StJaccard", "sumthreshold_parallel_jaccard"),
        ("StJaccardLo", "st_sweep_jaccard_lo"),
        ("StJaccardHi", "st_sweep_jaccard_hi"),
        ("StParOnly", "sumthreshold_par_only_flags"),
        ("StCpuOnly", "sumthreshold_cpu_only_flags"),
        ("FfaFound", "ffa_period_found"),
        ("FfaErr", "ffa_period_err_samples"),
        ("FfaDriftRes", "ffa_drift_resolution_samples"),
        ("FfaSnr", "ffa_snr"),
        ("OracleSnr", "fold_oracle_snr"),
        ("FfaRatioPred", "ffa_snr_ratio_predicted"),
        ("FfaRatioMeas", "ffa_snr_ratio_measured"),
        ("CrabPeriodMs", "crab_period_found_ms"),
        ("CrabFoldPub", "crab_fold_snr_at_published"),
        ("CrabJaccard", "crab_mask_jaccard"),
        ("CrabParOnly", "crab_mask_par_only_flags"),
        ("CrabCpuOnly", "crab_mask_cpu_only_flags"),
        ("CrabNoiseSnr", "crab_noise_ffa_snr"),
        ("CrabInjMin", "crab_injection_min_amp_sigma"),
        ("CrabDur", "crab_duration_s"),
        ("ChimeDm", "chime_dm_catalogue"),
        ("ChimeOffMax", "chime_off_dm_max_snr"),
        ("ChimeConjSnr", "chime_snr_conjugate_sign"),
        ("ChimeNchan", "chime_n_chan_good"),
        ("ChimeBw", "chime_bandwidth_mhz"),
    )
    bench_keys = (("Chirp", "chirp_s"), ("St", "sumthreshold_s"), ("Ffa", "ffa_s"))
    for ns in ("tdSyn", "tdReal"):
        live = ns == pref
        for macro, key in keys:
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
        lines.append(
            rf"\newcommand{{\{ns}SkMaxDiffSci}}{{{g_sci('sk_max_diff') if live else '--'}}}"
        )
        lines.append(
            rf"\newcommand{{\{ns}GroupDelayResid}}"
            rf"{{{g_sci('group_delay_max_resid_samples', dedisp) if live else '--'}}}"
        )
        # `device` labels the SCIENCE leg; `benchmark_device` labels the timing run.
        # They differ legitimately: the real CHIME leg ran on CPU while the benchmark was
        # measured on the RX 7600 XT. Overloading one field for both would either mislabel
        # CPU science as GPU or discard a real GPU measurement.
        bench_dev = m.get("benchmark_device") or m.get("device")
        gpu_src = m.get("benchmark") if bench_dev != "cpu" else None
        # the CPU column comes from benchmark_cpu, or from `benchmark` ONLY when the run
        # itself was on CPU -- never silently label GPU timings as CPU
        cpu_src = m.get("benchmark_cpu") or (
            m.get("benchmark") if m.get("device") == "cpu" else None
        )
        for macro, key in bench_keys:
            gpu = (gpu_src or {}).get(key) if live else None
            cpu = (cpu_src or {}).get(key) if live else None
            lines.append(
                rf"\newcommand{{\{ns}Bench{macro}Gpu}}{{{gpu if gpu is not None else '--'}}}"
            )
            lines.append(
                rf"\newcommand{{\{ns}Bench{macro}Cpu}}{{{cpu if cpu is not None else '--'}}}"
            )
            # the speedup DERIVED from the two committed timings, so the abstract can never
            # again carry a hand-typed ratio that contradicts its own parenthetical
            speedup = None
            if gpu is not None and cpu is not None and float(gpu) > 0:
                speedup = round(float(cpu) / float(gpu), 1)
            lines.append(
                rf"\newcommand{{\{ns}Speedup{macro}}}{{{speedup if speedup is not None else '--'}}}"
            )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: a CPU run has no GPU benchmark keys and would blank the
    # GPU macros a previous ROCm run wrote. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="torch-dsp: coherent dedispersion/RFI/FFA suite.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--device", default="cpu")
    p.add_argument("--benchmark", action="store_true")
    p.add_argument(
        "--benchmark-only",
        action="store_true",
        help="single-session GPU+CPU benchmark (and, on GPU, the cross-device kernel checks) "
        "merged into the existing results file; touches no science leg",
    )
    args = p.parse_args(argv)
    if args.benchmark_only:
        m = run_benchmark_only(args.out, device=args.device)
    else:
        m = run(args.out, offline=args.offline, device=args.device, bench=args.benchmark)
    m.pop("ffa_curve", None)  # array; the JSON artifact already excludes it
    print(json.dumps(m, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
