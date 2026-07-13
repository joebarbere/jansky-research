---
name: pull-station-data
description: Pull codified observation bundles (averaged HI spectra + full provenance) from the jansky-observe rooftop station into data/station/, over the station's MCP tools or its HTTP+zip API. Use when you need self-collected station spectra for plan 78's hline pipeline (read_capture) or plans 79/80. Read-only; the station carries no write/delete verbs.
---

# Pull station data

Fetch **codified observation bundles** from the [`jansky-observe`](https://github.com/joebarbere/jansky-observe)
rooftop station into this repo. Each bundle is one documented, machine-recoverable
export per observation — schema **`jansky-observe.observation-bundle/1`**:

- **`bundle.json`** — the provenance manifest: the station **UUID**, the observation
  with its pointing, and a per-capture block carrying SDR settings (gain, center
  freq, sample rate), **LST at the capture's start**, timestamps, cal-epoch
  reference, and classifier verdicts.
- **`capture-<id>.npz`** per npz capture — the averaged spectrum (`frequency_hz`,
  `power_db`) plus scalar provenance (`station_uuid`, `az_deg`, `el_deg`,
  `lst_hours`, `gain`, `center_freq_hz`, `sample_rate_hz`) so each file is usable
  standalone.

This is the **"averaged-spectra format from the station's capture service"** that
[plan 78](../../plans/78-station-hline-pipeline.md)'s `hline.read_capture` consumes,
and the calibration substrate plans 79/80 build on. The station is the single owner
of its data; this skill only ever **reads**.

## When the station has real data

⚠️ **Hardware-gated.** Until the H-line receiver chain is assembled and producing
first light (`station/hydrogen-line-receiver.md`; targeted late summer 2026), the
station runs its synthetic source — bundles pulled now carry synthetic/fake-HI
spectra, useful for wiring `read_capture` and the comparison harness offline but
not real sky. Pull real spectra once the roof is live.

## Prerequisites

- The station reachable on the LAN. Default base URL `http://raspberrypi.local:8000`
  (override with `--station-url` or `$JANSKY_STATION_URL`; the Pi is also at
  `http://10.3.1.106:8000`). Check: `curl -fsS $URL/healthz`.
- No Python dependency — the bulk script is stdlib only.

## Two ways to pull

### 1. Interactive, via the station's MCP (a few observations)

When the station's MCP is connected
(`claude mcp add --transport http jansky-observe http://raspberrypi.local:8000/mcp`),
use its read tools directly:

- `get_station_identity` → the station UUID (the per-station key everything is filed under).
- `list_observations` (optionally `status="done"`) → pick observations.
- `get_observation_bundle(observation_id)` → the **manifest** (JSON provenance) for one
  observation, in one call.
- `get_spectrum(capture_id, axis="mhz"|"vlsr")` / `get_capture_meta(capture_id)` /
  `export_capture(...)` → per-capture spectra + metadata when you want a single trace.

The manifest is enough to decide *what* to fetch; the averaged-spectrum arrays themselves
come down as the `.npz` files in the zip (below) or via `get_spectrum`.

### 2. Bulk, via the HTTP+zip API (many observations)

`pull_bundles.py` hits the same data over the station's plain JSON+zip endpoints and
unpacks each bundle into `data/station/<station-uuid>/observation-<id>/`:

```bash
# List what would be pulled (downloads nothing):
uv run python .claude/skills/pull-station-data/pull_bundles.py --list

# Pull every 'done' observation into data/station/:
uv run python .claude/skills/pull-station-data/pull_bundles.py

# One observation, from the Pi's IP, into a scratch dir:
uv run python .claude/skills/pull-station-data/pull_bundles.py \
    --station-url http://10.3.1.106:8000 \
    --observation-id 12 --out data/station
```

Flags: `--station-url` (env `JANSKY_STATION_URL`), `--out` (default `data/station`),
`--observation-id` (repeatable; default pulls the list), `--status` (default `done`,
`any` disables the filter), `--limit`, `--list`, `--timeout`.

## Layout it writes

```
data/station/<station-uuid>/
  observation-<id>/
    bundle.json          # provenance manifest (schema jansky-observe.observation-bundle/1)
    capture-<id>.npz     # one per npz capture: frequency_hz, power_db + scalar provenance
```

Filing under the station UUID keeps a future second station's data separate. `data/` is
ruff/CI-excluded, so pulled bundles never enter lint/coverage.

## Exit codes

| code | meaning |
|---|---|
| 0 | success — every requested bundle fetched and unpacked (or nothing matched) |
| 2 | bad usage |
| 3 | station unreachable / identity fetch failed |
| 4 | the observation list could not be fetched |
| 5 | one or more requested observations failed |
| 6 | a bundle downloaded but was not a valid zip |

## The contract

jansky-observe owns the bundle format + the MCP tools this skill calls
(`get_station_identity`, `list_observations`, `get_observation_bundle`, `get_spectrum`,
`get_capture_meta`, `export_capture`) — shipped in its **v0.9.0 / M8** milestone. This
skill is jansky-research's consumer side; `hline.read_capture` (plan 78) reads the
`.npz` files it lands.
