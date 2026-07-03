# Hydrogen Line Receiver

A small-dish receiver for the 21 cm (1420.4058 MHz) emission line of neutral hydrogen. The goal is detection of Galactic hydrogen, followed by Doppler measurements of the Milky Way's spiral arms and, over a year of operation, the annual modulation from Earth's orbital velocity.

The 1400–1427 MHz band is internationally protected for radio astronomy, which is why this experiment works even from an urban rooftop.

## Signal chain

```
20 dBi mesh parabolic (1.4 GHz center)
  → wideband LNA (Nooelec LaNA, powered via USB)
  → 1420 MHz bandpass filter (GPIO Labs)
  → LMR-400 coax
  → Airspy Mini SDR (bias tee disabled)
  → host (SDR++ for commissioning; Virgo/ezRA for science capture)
```

Two architectural rules drive this design:

**Amplify before you lose.** The system noise figure is set almost entirely by the first amplifier, so the LNA mounts directly at the dish feed. Coax losses after the LNA are benign; before it, they are unrecoverable.

**Power architecture.** The bandpass filter is not DC-tolerant — DC on either port destroys it. The LNA is therefore USB-powered and the SDR's bias tee must remain disabled. If a bias-tee-powered configuration is ever needed, a DC block belongs on the filter's receiver side.

An integrated alternative to the LNA + filter pair is the Nooelec SAWbird+ H1 (two LNAs cascaded around a 1420 MHz SAW filter, ~0.8 dB noise figure); the discrete chain was chosen for availability and because the parts are individually reusable.

## Why these parts

- **Airspy Mini** over an RTL-SDR: a 12-bit ADC (vs 8-bit) buys roughly 24 dB of dynamic range, and a 0.5 ppm TCXO keeps a spectral line from drifting with temperature — both matter when the signal is a small bump on the noise floor and the science is in its frequency.
- **Mesh dish** over the traditional DIY foam-board horn: roughly 10× the gain, weatherproof, and reusable for L-band satellite work (a swappable boom converts it to 1.7 GHz for GOES reception).
- **Filter order** is itself an experiment: LNA-first gives the best noise figure; filter-first protects the LNA from overload by strong out-of-band transmitters (urban cell sites). Both configurations will be measured.

## Mounting

The dish rides a ballasted tripod on a small flat roof. Wind, not weight, is the design load: every leg is ballasted, and a permanent installation will use a non-penetrating ballasted mount rather than membrane penetrations. A coax grounding block is installed at the building entry.

## Observation plan

1. **Commissioning:** waterfall at 1420.4 MHz, gain calibration to avoid front-end overload.
2. **First detection:** pointed integrations (60–300 s) toward the Galactic plane in Cygnus; the on-plane vs off-plane spectral difference is the detection.
3. **Rotation curve:** spectra at multiple galactic longitudes; Doppler offsets map spiral-arm velocities.
4. **Long duration:** fixed-pointing daily spectra for a year to trace Earth's ±30 km/s orbital velocity as a ±142 kHz annual sinusoid, plus continuous drift-scan mapping. See [operations](operations.md).
