# Two-Dish Interferometry

A planned extension: adding a second dish converts the station from a single-aperture radiometer into an interferometer — with an important physical caveat that shapes the whole plan.

## What a small interferometer can and cannot see

An interferometer responds to structure smaller than its fringe spacing and resolves *out* smooth extended emission. Galactic 21 cm hydrogen fills degrees of sky, so on any baseline longer than a couple of meters most of it disappears from the correlated signal. Two-dish amateur interferometry therefore targets **compact continuum sources**: the Sun first (strong and easy), then Cassiopeia A and Cygnus A as stretch goals. The hydrogen line hardware is reused, but the interferometric science is continuum.

## Tier 1 — Additive (phase-switched) interferometer

The achievable configuration, and the plan of record:

```
dish A → LNA + 1420 MHz filter ─┐
                                 ├─ matched-length coax → RF power combiner → single SDR
dish B → LNA + 1420 MHz filter ─┘
```

With one receiver there is no inter-channel coherence problem. As the Sun drifts through the overlapped beams, total power oscillates — fringes — and the fringe spacing yields the Sun's angular size. This reproduces 1950s-era radio astronomy techniques faithfully, for roughly $200 of additional hardware (second dish, second LNA, combiner, matched cables). The two bandpass filters were purchased from the same production batch for this purpose.

## Tier 2 — Digital correlation interferometer

True correlation requires the two receive channels to share one clock; independent USB SDRs cannot be made coherent off the shelf. The practical amateur route is a KrakenSDR (five RTL-SDR channels on a shared oscillator) with a noise source for phase calibration and a GNU Radio correlator, following the work published by Marcus Leech / CCERA. This is a substantial software project and is deliberately sequenced last.

## The alternative role for dish B

For long-duration spectroscopy, the second dish may contribute more as a **Dicke-switched comparison radiometer**: an RF switch alternates the receiver rapidly between the on-source dish and a cold-sky reference dish, canceling receiver gain drift. Year-scale measurements (the annual Doppler curve) are limited by calibration stability rather than raw sensitivity, so this configuration competes seriously with interferometry for the second dish's time.
