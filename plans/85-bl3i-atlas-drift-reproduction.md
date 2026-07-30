# 85 — Independent reproduction of the Breakthrough Listen 3I/ATLAS GBT nondetection

Status: 🚧 in progress — GATE 0 done 2026-07-30: data verified public + novelty PASS (details
below). Step 1 done 2026-07-30 — node→band mapping pinned by remote header reads (fsspec+h5py
HTTP range requests on the `.0002` products, ~KB per file): MJD 61027, four cadences by start
second — **L** t0=17226: blc21–26 = 939.0–2064.0 MHz; **S** t0=21817: blc22–27 =
1651.5–2776.5 MHz; **C** t0=26882: blc21–27/30–37/60–67/70–75 = 3939.0–8251.5 MHz; **X**
t0=31308: blc20–27/30–37/60–67/70–76 = 7689.0–12376.5 MHz. Each node = 187.5 MHz; `.0002` is
the 65536-channel (~2.86 kHz) mid-res product. Minimum-unit node: **blc25** (1126.5–1314 MHz,
inside FAST's 1.05–1.45 GHz comparison band). Disk checked: 348 GB free vs ~60 GB node-serial
peak.

## Context

On 2025-12-18 Breakthrough Listen observed the interstellar object 3I/ATLAS with the GBT across
1–12 GHz (four receivers L/S/C/X, ABACAD 6×5-min on/off cadences), searched with turboSETI at
±4 Hz s⁻¹ / ~16σ, and reported a nondetection down to ~100 mW EIRP (RNAAS, arXiv:2512.19763).
The data are public, but — as with most technosignature nulls — **no one has independently
reproduced the search from the released files**: the other 3I/ATLAS searches (ATA
arXiv:2512.18142; FAST arXiv:2603.19023, L-band 1.05–1.45 GHz; MeerKAT/Parkes on the same
portal) are each teams analysing their *own* observations. That leaves the reproducibility cell
empty, and it is exactly our pattern: the merged `driftsearch` slice (plan 11) already gives us
a tested pure-NumPy dedoppler search, an injection-recovery-calibrated threshold, and the
h5py/hdf5plugin plumbing validated on the Voyager-1 BL file — and plan 11's DC-spike lesson
("the apparent Voyager detection was an artifact") is precisely the kind of check an
independent pipeline contributes to someone else's candidate list. Deliverable framing: an
independent-pipeline confirmation (or not) of the null + our own EIRP limit on the same haystack
axes as the GBT/ATA/FAST trio, honest about the sub-band we actually processed.

**GATE 0 (done 2026-07-30).** *Data:* `https://bldata.berkeley.edu/ATLAS/GB_ATLAS/` (direct
HTTP index, egress verified) — files
`blc{20..35}_guppi_61027_*_DIAG_3I_ATLAS[_OFF]_*.rawspec.{0000,0001,0002}.h5`; per scan the
three rawspec products are ~10 GB (fine-frequency, ~3 Hz — the SETI product), ~1 GB, and
~50 MB; timestamps 09:02–12:06 UT spanning the four receiver cadences. Full fine-res volume is
multi-TB → scope to **L band** (matches the FAST comparison band; most sensitive receiver),
processed node-serially (download → search → delete, resumable), minimum publishable unit = one
full ABACAD cadence for one compute node (~187.5 MHz, 6×10 GB on disk at peak), stretch = all
L-band nodes. *Novelty:* web pass 2026-07-30 finds no independent reanalysis of the released
GB_ATLAS data; the RNAAS is 7 months old and the follow-ups above are all original
observations. Standing fable-ideas caveat satisfied (full-text + data-URL check done here).

## Deliverables

- `src/jansky_research/atlas3i.py`: `fetch_cadence` (portal HTTP download of one node's 6-scan
  cadence, resumable, delete-after option, `# pragma`), `read_scan` (BL rawspec `.h5` →
  frequency/time grid + header; shares the plan-11 `voyager` extra: h5py + hdf5plugin),
  `dedoppler_cadence` (drift search ±4 Hz s⁻¹ over the fine-res product, reusing
  `jansky.seti.drift_search` / plan-11 threshold calibration; optional torch backend via the
  merged `torchdsp` suite if NumPy is too slow at 67M channels), `onoff_filter` (ABACAD
  rejection: candidate must appear in all ONs, no OFFs), `eirp_limit` (EIRP = 4πd²·S·Δν from
  the measured noise + the 2025-12-18 geocentric distance; compare against the paper's
  ~100 mW and the ATA/FAST limits), `synthetic_cadence` (injected drifting tone + an
  on/off-rejectable RFI tone in a small synthetic cadence), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor (offline, synthetic cadence); `papers/atlas3i/`;
  `survey/atlas3i-findings.md`; wiring (`data.py`, Snakefile, README table).

## Approach

1. Pin the node→band mapping by reading header attrs (`fch1`/`foff`/`tstart`) from the ~50 MB
   `.0002` products of every node — cheap, and it fixes which blc nodes are the L-band scope.
2. Tooling + synthetic round-trip: injected tone recovered through
   `dedoppler_cadence`+`onoff_filter`; injected always-on RFI tone rejected; threshold false
   positive rate calibrated as in plan 11.
3. Real leg (network, `# pragma`): one L-band node full cadence → candidate list → on/off
   filter → compare with the paper's (empty) candidate table; extend node-serially across
   L band as disk/time allow. Record every sub-band processed vs skipped — no silent caps.
4. EIRP limit from our own noise measurement + Horizons distance at the epoch; place on the
   haystack axes vs GBT/ATA/FAST.
5. GATE-2 science review: threshold/unit honesty (our σ vs their 16σ convention — Choza+2024),
   the "reproduction ≠ same-data-independence" framing (same data, independent pipeline),
   coverage honesty (which GHz actually searched).
6. Paper (`papers/atlas3i/`, RNAAS-length target), `\software{}`-cites `jansky-research` +
   `jansky`.

## Verification

Synthetic cadence round-trip (recover injected, reject RFI, FPR-calibrated threshold); real-leg
sanity: our noise level consistent with the paper's implied sensitivity before any limit claim;
`make cov` ≥85%, ruff+mypy clean; GATE-2 sign-off.

## Risks & mitigations

- **Data volume** (10 GB × 6 per node): node-serial download→search→delete with a resumable
  manifest; minimum publishable unit is one node, stated plainly in the paper.
- **Portal file layout surprises** (segment meaning, node→band): step 1 pins everything from
  real headers before any bulk download.
- **NumPy too slow at 67M channels × drift trials**: fall back to the merged `torchdsp` GPU
  path (ROCm leg already validated), or coarser drift grid with the sensitivity cost stated.
- **A "candidate" that survives our filter but not theirs (or vice versa)**: that is a result,
  not a failure — report the disagreement with the plan-11 artifact checklist applied.
