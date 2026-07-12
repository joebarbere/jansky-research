# Idea note — the PWN X-ray/radio polarization mismatch (IXPE Lighthouse follow-up)

**Status: candidate idea (2026-07-12), not a plan.** The standing GATE-0 applies before any
plan is written: full-text novelty pass (read the IXPE paper's radio references — the
comparison may already be published) + a data-URL check.

## The trigger

Dinsmore, J. T., et al. 2026, *IXPE Polarizations of the Lighthouse Pulsar, Trail, and
Filament* (ApJ accepted, [arXiv:2604.22914](https://arxiv.org/abs/2604.22914)); NASA release
2026-07. A 1 Ms IXPE observation of PSR J1101−6101 / the Lighthouse Nebula:

- The escape **filament** is polarized at PD 55 ± 18% with the EVPA implying a magnetic
  field **parallel to the filament axis** — magnetic turbulence weaker than models assumed.
- The **trail's X-ray polarization is nearly orthogonal to its radio polarization**,
  implying the X-ray- and radio-emitting leptons occupy spatially distinct regions —
  i.e., more than one acceleration mechanism at work.

That second bullet is the radio-side hook: it is a claim *about radio data*, testable and
extensible from public archives.

## Candidate slice shapes (unranked, unvetted)

1. **A radio-polarimetry census of PWN trails/filaments with X-ray polarimetry** — collect
   the (small) set of PWNe with published or archival radio polarization (Lighthouse /
   IGR J11014−6103, the Guitar Nebula, J1509−5850, the Mouse, …) and ask whether radio-vs-
   X-ray B-field orthogonality is a one-off or a pattern. Sources: ATCA archive (ATOA),
   RACS/POSSUM Stokes Q/U via CASDA, published EVPA maps; IXPE data are public.
2. **Narrower reproduction leg**: RM-corrected radio EVPA of the Lighthouse trail from
   archival ATCA data, compared quantitatively against the IXPE map — recover-a-known
   shaped, honest about beam/RM systematics.

## Honest caveats (why this is a note, not a plan)

- Dec −61° and **imaging polarimetry**, not catalogue work — heavier than the usual slice;
  closer to the `rmstructure`/`stokesv` end of effort than the catalogue end.
- POSSUM coverage/DR status at this field must be checked (data-URL gate), and ATCA
  archival continuum polarimetry needs calibration work we have not done before.
- **Novelty risk is real**: the orthogonality claim comes *from* comparing to existing
  radio maps, so the census version must add something (uniform re-reduction, more
  objects, RM correction) beyond what the IXPE team already cited.

## Fit with existing muscles

Polarization/RM experience (`rmsky`, `rmstructure`, `rmdipole`), CASDA/cutout tooling
(`casda-cutout-fetch`, `radio-cutout` skills), and the house recover-a-known discipline.
