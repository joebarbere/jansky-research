"""hiblend R7 prediction (findings step 10): beta(u) from the injection pool alone -- separations and flux ratios, no beta measured. Prints only; writes nothing."""
import numpy as np
from jansky_research import fashienv as fe, hiblend as h
raw = fe.fetch_fashi_dr2()
ok = (np.isfinite(raw["flux"]) & (raw["flux"] > 0) & np.isfinite(raw["flux_err"]) & (raw["flux_err"] > 0)
      & np.isfinite(raw["cz"]) & np.isfinite(raw["w50"]) & np.isfinite(raw["ra"]) & np.isfinite(raw["dec"]))
fashi = {k: np.asarray(v)[ok] for k, v in raw.items()}
f = h.build_field(fashi, h.fetch_alfalfa())
m = h.sample_masks(len(f["flux_f"]), f["nb"]); nb = f["nb"]
ring = m["primary"][nb.target] & nb.overlap & (nb.sep_arcmin >= h.NEIGHBOUR_RMIN_ARCMIN)
sep, rho = nb.sep_arcmin[ring], nb.flux[ring] / f["flux_f"][nb.target[ring]]
bA, bF = h.beam_response(sep, h.ALFA_FWHM_ARCMIN), h.beam_response(sep, h.FAST_FWHM_ARCMIN)
# first order with finite rho: resid = log10((1+u rho bA)/(1+rho bF)), R = log10((1+rho bA)/(1+rho bF))
R = np.log10((1 + rho * bA) / (1 + rho * bF))
for u in (0.6, 0.86, 1.0):
    res = np.log10((1 + u * rho * bA) / (1 + rho * bF))
    lin = np.sum(R * res) / np.sum(R * R)  # OLS through the origin on the pool
    print(u, "pool-OLS beta", round(float(lin), 3))
d = bA - bF; w = rho**2
print("small-rho A,B", round(float(np.sum(w*d*bA)/np.sum(w*d*d)),3), round(float(np.sum(w*d*bF)/np.sum(w*d*d)),3))
print("pool n", sep.size, "median sep", round(float(np.median(sep)),2))
