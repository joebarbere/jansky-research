# 78 — Open H-line pipeline: calibrated rooftop spectra recovering LAB/HI4PI per pointing

Status: 📋 planned (hardware-gated) — needs the H-line receiver chain assembled and producing
first-light captures (in acquisition; first light targeted late summer 2026); the software +
synthetic leg and the LAB/HI4PI comparison harness can be built and CI'd now, offline

## Context

The campus-telescope/CHART genre (arXiv:2307.11173, 2404.17893) publishes small-dish H-line
detections but ships no tested, CI'd pipeline, and validation against the public surveys is
qualitative ("the bump is there"). The station sweep's common finding: the amateur literature
stops at qualitative detections, so rigor — per-pointing calibrated spectra reproducing LAB/HI4PI
profiles to *stated accuracy* — is the novelty axis. This slice builds that pipeline for the
Philadelphia rooftop receiver (`station/hydrogen-line-receiver.md`), extending `hi.py`'s LAB
reader (`read_lab_slice`/`fetch_lab_longitude`) into a survey-comparison harness — and it is the
calibration substrate for plans 79 (annual Doppler) and 80 (drift-scan), which sequence after it.

## Deliverables

- `src/jansky_research/hline.py`: `read_capture` (averaged-spectra format from the station's
  capture service), `bandpass_calibrate` (reference-load/frequency-switched gain removal),
  `freq_to_vlsr` (barycentric + LSR correction), `survey_profile` (LAB/HI4PI profile at a
  pointing, extending `hi.py`), `compare_pointing` (residual + accuracy statistic),
  `synthetic_capture` (injected profile + gain shape + noise), `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic/offline fixtures — no sky data needed for tests);
  `papers/hline/`; `survey/hline-findings.md`; wiring.

## Approach

0. GATE 0: receiver chain per `station/hydrogen-line-receiver.md` mounted and producing raw
   averaged spectra (commissioning waterfall done, gain set below overload).
1. Tooling + synthetic recover-a-known, ahead of hardware: inject a known HI profile through a
   simulated bandpass + noise; the pipeline must recover the profile's amplitude, centroid, and
   width to stated tolerance. Build the HI4PI/LAB per-pointing comparison harness offline too.
2. First-light leg (`# pragma: no cover`): pointed 60–300 s integrations toward the Galactic
   plane in Cygnus + off-plane reference; calibrate; compare each pointing to LAB/HI4PI.
3. Extend to a handful of longitudes (the rotation-curve pointings); report per-pointing
   accuracy (K-scale or relative), centroid error in km/s, and the urban RFI environment.
4. GATE-2 science review → paper: an open, tested, CI'd small-dish H-line pipeline with
   quantified survey agreement — the artifact the genre lacks.

## Verification

Synthetic round-trip recovers injected profile parameters; real leg reproduces LAB/HI4PI
profiles per pointing to a stated, honest accuracy (that number *is* the result); checks green;
GATE-2 sign-off.

## Risks & mitigations

- **Hardware/first-light delay** → the entire software + synthetic + survey-harness stack lands
  first; only step 2 waits on the roof.
- **Urban RFI in the protected band's shoulders** (Philadelphia rooftop) → RFI flagging in
  `bandpass_calibrate`; off-plane reference pointings; report occupancy honestly.
- **Absolute calibration is hard with hobby parts** → lead with relative/shape accuracy and a
  clearly stated temperature-scale uncertainty; plan 84's Cas A/Cyg A work sharpens it later.
