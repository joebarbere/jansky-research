# Long-Duration Operations

The station is designed to run continuously for a year or more, because for small apertures **integration time substitutes for collecting area**, and several of the target measurements only exist at year scale:

- **Earth's orbital velocity from hydrogen.** Daily spectra of a fixed galactic region show the ±30 km/s orbital velocity as a ±142 kHz annual Doppler sinusoid — a measurement of Earth's orbit made with a rooftop dish, and one that requires a full year by definition.
- **Drift-scan sky mapping.** A fixed dish sweeps the same declination strip every sidereal day (23h56m); months of stacked passes yield a high-SNR hydrogen strip, and periodic elevation steps build a crude all-sky map over one to two years.
- **Solar radio monitoring** as the Sun transits the beam daily.
- **Sidereal/solar separation as data quality control:** astronomical signals arrive ~4 minutes earlier each day, while interference follows the 24-hour human clock — after months of data the two populations separate cleanly.

## Reliability target

The tracked service-level objective is **data continuity** — the percentage of scheduled spectra successfully captured — rather than raw process uptime. Data continuity is the honest metric: it survives planned maintenance windows, and it measures what the science actually needs.

## Architecture

Outdoor, at the dish: an active hydrogen-line feed with integrated LNAs and SAW filters, powered up the coax by a USB bias-tee injector (~5 V, ~120 mA); the SDR's own bias tee stays disabled and the feed's internal DC block keeps DC off the SDR. Drip loops and self-amalgamating tape on all connections; short coax to the SDR enclosure on the same mast; a grounding block at building entry; non-penetrating ballasted roof mount.

At the dish, in a weatherproof enclosure: the SDR and a Raspberry Pi 5, powered and networked over Power-over-Ethernet from an indoor switch on a small UPS. The Pi runs:

- a capture service producing **averaged spectra only** (megabytes per day; raw IQ is never logged),
- systemd units with restart policies plus a hardware watchdog,
- nightly rsync of data off-box, with monthly cold backups,
- Tailscale for remote administration,
- Prometheus + Grafana for both station health (temperatures, disk, capture count) and live science dashboards (SNR trends, latest spectra).

## Calibration schedule

Weekly reference measurements (50 Ω load or cold-sky pointing) anchor the year-long gain history; comparing January's spectra to July's is fundamentally a calibration-stability problem, which is also why a Dicke-switched second dish (see [interferometry](interferometry.md)) is under consideration.

## Year-1 plan

Months 0–1: commissioning and RFI baseline. Months 1–12: daily fixed-pointing spectra for the annual Doppler curve, with continuous drift-scan mapping in parallel and monthly elevation steps. Results published quarterly.
