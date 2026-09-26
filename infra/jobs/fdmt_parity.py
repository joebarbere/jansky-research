"""CPU vs CUDA parity for the torch-fdmt kernels (plan 96, phase 2).

Runs FDMT and brute-force dedispersion on the repo's offline Crab-like fixture on both devices
and reports the largest difference, plus whether both devices pick the same best DM. The CPU
path is itself checked against the NumPy oracle by the test suite, so CPU == CUDA here closes
the chain oracle -> CPU -> CUDA.
"""

import json

import numpy as np
import torch

from jansky_research.fdmt import brute_dedisperse, fdmt
from jansky_research.singlepulse import CRAB_DM, synthetic_observation

dyn, freqs, dt = synthetic_observation()
out: dict = {"torch": torch.__version__, "cuda": torch.version.cuda}
out["gpu"] = torch.cuda.get_device_name(0)

a = fdmt(dyn, freqs, dt, 120.0, device="cpu")
b = fdmt(dyn, freqs, dt, 120.0, device="cuda")
pa, pb = a.plane.cpu().numpy(), b.plane.cpu().numpy()
scale = float(np.abs(pa).max())
out["fdmt_shape"] = list(pa.shape)
out["fdmt_max_abs_diff"] = float(np.abs(pa - pb).max())
out["fdmt_max_rel_diff"] = float(np.abs(pa - pb).max() / scale)
out["fdmt_best_cpu"] = a.best()
out["fdmt_best_cuda"] = b.best()
out["fdmt_same_best_dm"] = a.best()[0] == b.best()[0]

dms = np.linspace(0.0, 120.0, 241)
bc = brute_dedisperse(dyn, freqs, dt, dms, device="cpu").cpu().numpy()
bg = brute_dedisperse(dyn, freqs, dt, dms, device="cuda").cpu().numpy()
out["brute_max_abs_diff"] = float(np.abs(bc - bg).max())
out["brute_max_rel_diff"] = float(np.abs(bc - bg).max() / float(np.abs(bc).max()))
out["brute_best_dm_cpu"] = float(dms[bc.max(axis=1).argmax()])
out["brute_best_dm_cuda"] = float(dms[bg.max(axis=1).argmax()])
out["catalogue_dm"] = CRAB_DM
print(json.dumps(out, indent=2))
