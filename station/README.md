# The rooftop station

SDR-based radio-astronomy instruments built and operated from a rooftop in Philadelphia,
documented end to end — hardware, software, and operations. These are build guides for the
physical station that is meant to feed **self-collected** data into this repo's research slices,
alongside the public archives the existing slices draw on.

Named for Karl Jansky, who detected the first cosmic radio emission with homemade equipment in New
Jersey, about 60 miles from here.

## Instruments

| Instrument | Target | Status |
|---|---|---|
| [Hydrogen line receiver](hydrogen-line-receiver.md) | 21 cm neutral hydrogen — Milky Way structure and rotation | Acquiring hardware |
| [Meteor scatter station](meteor-scatter-station.md) | FM forward-scatter meteor detection, 24/7 | Planned |
| [Two-dish interferometer](interferometry.md) | Solar fringes; additive interferometry | Future |

Supporting references: [test equipment](test-equipment.md) · [long-duration operations](operations.md)

## Design principles

- **Integration time substitutes for aperture.** Small dishes, long uptimes.
- **The LNA lives at the antenna.** Signal-chain losses before amplification are unrecoverable.
- **Buy once, reuse everywhere.** Wideband test gear and swappable antennas over single-purpose parts.
- **Treat the station as a service.** Watchdogs, dashboards, and a data-continuity SLO rather than raw uptime.

## Status

Station started July 2026. First light targeted for late summer 2026, with the meteor station's
first shower target being the Perseids (Aug 12–13). Once the receiver is producing calibrated
spectra, self-collected data joins the public-archive slices under `../papers/`.

> These are the public-facing build guides. The owner's working notes — purchase log with prices,
> per-part rationale, and personal reminders — live in an Obsidian vault, not this repo.
