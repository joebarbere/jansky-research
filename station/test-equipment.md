# Test Equipment

Instrument-grade results from hobby-grade parts require verification: every filter, amplifier, and antenna in the signal chain gets measured before first light.

## Bench inventory

| Instrument | Role |
|---|---|
| HackRF One | Signal generator: weak CW tones and sweeps (e.g. 1380–1460 MHz through the bandpass filter into a second SDR). Not used as the science receiver — 8-bit ADC and no TCXO make it the weakest receiver on the bench, but transmit capability makes it the only generator. |
| LiteVNA 64 | Vector network analyzer to 6 GHz. Filter S21 sweeps, antenna SWR at 1420 MHz, cable checks. Chosen over entry-level NanoVNA units, which operate in degraded harmonic mode above ~300 MHz and are unreliable at L-band. |
| tinySA Ultra+ | Spectrum analyzer: RFI surveys, LNA gain verification. An external attenuator stays on the input permanently — maximum safe input is about +6 dBm and front-end damage is not warrantied. |
| Passives | 50 Ω terminator, SMA attenuator set (3–30 dB), DC blocks, adapter kit. |

## Used-SDR health checks (HackRF example)

A second-hand SDR earns trust through a triage sequence, cheapest test first:

1. **Enumeration:** `hackrf_info`; flash current firmware; `hackrf_debug --si5351c` to confirm the clock synthesizer responds.
2. **Mechanical:** inspect the SMA center pin and USB jack under cable strain.
3. **Front-end amplifier:** tune a steady local FM broadcast carrier and toggle the RX amp in software — expect ~14 dB of change. No change indicates the input MMIC amplifier is blown (the classic ESD failure on used units; a ~$3 SMD repair).
4. **Frequency accuracy:** receive a known stable carrier (NOAA weather radio, 162.4–162.55 MHz) simultaneously on the device under test and on a TCXO-equipped reference SDR; the offset difference is the unit's ppm error, stored and applied as a software correction.
5. **Spur scan:** with a 50 Ω terminator on the input, sweep the full range and look for a flat noise floor.

## Counterfeit avoidance

The hobby SDR market is heavily cloned. Working rules adopted for this project: purchase only through manufacturer-listed sellers (rtl-sdr.com store links; tinysa.org's authorized-seller list; the Airspy US distributor's domestic-warehouse page), verify model-identifying details (genuine RTL-SDR Blog units use black metal enclosures and the R860 tuner; blue/green-cased "RTL-SDR Blog" listings are counterfeit), and treat any listing priced 30%+ below the known-good price as a clone regardless of its photos.
