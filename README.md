# jansky-research

**Amateur radio-astronomy research, end to end.** A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky *teaches* radio
astronomy, this repo *does* original amateur research — surveys the landscape, builds tested
tooling, analyzes public data, automates it, and writes it up — with honesty as the first rule.

It began as a single vertical slice (survey → gap → tool → automated pipeline → reproducible paper)
and has grown into a **deep-research survey plus a set of self-contained research slices**, each:
a tested CPU-only tool reusing jansky's helpers, run on real public data, put through an
adversarial science-review gate, and written up — wins *and* negatives reported plainly.

## Status — twenty-three slices + a synthesis, honestly tallied

| Slice | Tool | Outcome |
|-------|------|---------|
| FRB burst-statistics | `jansky_research.frbstats` | ✅ reproduced the CHIME repeater **width** result |
| Ultra-steep-spectrum hunt | `jansky_research.spectra` | ➖ USS candidates **failed** the de Gasperin cross-check (negative) |
| FRB repeater periodicity | `jansky_research.frbperiod` | ✅ recovered FRB 20180916B's **16.35-day** period |
| SETI drift-search benchmark | `jansky_research.driftsearch` | ➖ benchmark built; the apparent "Voyager detection" was a **DC-spike artifact** caught in development — the tool reports an honest null |
| HI rotation curve | `jansky_research.hi` | ✅ recovered the **flat** (non-Keplerian) inner Milky Way curve |
| VLASS multi-epoch variability | `jansky_research.vlass` | ✅/➖ 703 deg² census: catalogue variability is **artifact-dominated**, but image-confirms **FK Comae Berenices** |
| Peaked-spectrum (GPS/CSS) selection | `jansky_research.peaked` | ✅ three-frequency curvature selector; **100% recovery** of a known HFP sample, high purity vs MHz-peaked |
| Southern peaked-spectrum (GLEAM-X×RACS) | `jansky_research.southern` | ✅ multi-band curvature that **measures** the turnover ν_pk; 90 candidates over a 3° cone, two systematic fixes |
| Radio–optical offsets (ICRF3 × Gaia × MOJAVE) | `jansky_research.offsets` | ✅ reproduces the AGN offset **excess tail** (≫ Rayleigh) **and its alignment with the parsec-scale jet** (KS p=3×10⁻²²) |
| Pulsar radio spectral indices (ATNF) | `jansky_research.pulsarspec` | ✅ reproduces the **steep** mean pulsar spectrum (α≈−1.8); MSPs ≈ normal pulsars |
| Sub-threshold radio stacking (SDSS quasars × VLASS-SE) | `jansky_research.stacking` | ✅ image-plane stacking with **injection-recovery** bias calibration; mean flux of undetected quasars |
| Multi-decade VLBI variability (Astrogeo) | `jansky_research.vlbi` | ✅ **control-floor** method recovers OJ 287 & BL Lac; blazars ~1.7× more variable than steady CSO controls |
| Solar type III exciter speed (e-Callisto) | `jansky_research.solarbursts` | ✅ drift→beam speed **~0.14 c** on a clean isolated burst (Newkirk inversion) |
| Galactic Faraday rotation sky (Taylor+2009) | `jansky_research.rmsky` | ✅ plane enhancement ratio **5.4**, inner-Galaxy RM **sign antisymmetry** recovered |
| Pulsar P–Ṗ diagram (ATNF) | `jansky_research.ppdot` | ✅ three classes over ~5 orders in B, 98% above the death line; Crab validates |
| Inner-heliosphere type III (Wind/WAVES) | `jansky_research.windwaves` | ✅ beam tracked to **~10 R⊙** (Alfvén surface); honest peak-time/R² caveats |
| Interplanetary type III (STEREO/WAVES HFR) | `jansky_research.swaves` | ✅ beam tracked to **0.38 AU**; honest R²-inflation (few independent time samples) caveat |
| 3D type III triangulation (STEREO-A+B DF) | `jansky_research.triangulate` | ✅ geometric vs plasma-frequency distance correlate at **r=0.989**; source localized in 3D |
| Euclidean source counts (NVSS) | `jansky_research.sourcecounts` | ✅ recovers the canonical **Hopkins 2003** 1.4 GHz counts; sub-Euclidean slope −1.91 |
| RACS Stokes-V coherent emitters (ASKAP via CASDA) | `jansky_research.stokesv` | ✅/➖ leakage-floor selection validated; forced photometry **recovers I**, but single-epoch **V is variability-limited** (honest) |
| Type III occurrence census vs the solar cycle (e-Callisto × SILSO) | `jansky_research.ecallisto_census` | ✅/➖ coverage-corrected census statistic **recovers an injected solar-cycle proportionality** (r=0.97); SILSO is real, event stream is synthetic — real multi-cycle ingest is future work (honest) |
| torch-fdmt: device-portable Fast DM Transform | `jansky_research.fdmt`+`singlepulse` | ✅ oracle-validated pure-PyTorch FDMT; **real Crab DM recovered to 0.3%** (giant pulse S/N 14); honest benchmark: FDMT-on-CPU beats brute-on-GPU 3.6× |
| RACS Stokes-V discovery: two-epoch forced photometry | `jansky_research.stokesv_discovery` | ✅/➖ recovers **GJ 65 (BL+UV Cet)** with a **10σ 4.2-yr V change**; all else quiescent (median 5σ limit 0.83 mJy); no new detections survived the novelty bar (honest) |
| Type III synthesis: corona → 0.4 AU (4 instruments) | `jansky_research.type3synthesis` | ✅ unified drift-to-distance ladder; **geometric check on the model distance** (same-event r=0.989) |

A long run of recover-a-known validations and methodology contributions, two honest negatives (the USS
candidates and the SETI "Voyager detection"), one mixed result (VLASS catalogue variability is
artifact-dominated but image-confirms a real variable star), and **zero overclaims that survived
review** — the science-reviewer caught the USS candidates evaporating against the authoritative
catalog, the SETI "detection" being an instrument artifact, and (every slice) at least one real
physics/citation/statistics fix before write-up. The negatives are arguably the most instructive part.
Each slice's honest assessment is in `survey/*-findings.md`.

New *data* slices are paused — the reliable no-auth sources (VizieR, SPDF, e-Callisto, Astrogeo,
CIRADA/VLASS) and the clean recover-a-known method space are largely used up. Work has shifted to
**synthesis and right-sized infrastructure**: the type III synthesis above (`plans/30`) and the
streaming **e-Callisto Airflow-on-Podman ingest pipeline** (`plans/31`) are built — Airflow now serves a
*frequently-updated* archive (daily schedule, catchup/backfill, per-station fan-out) instead of a static
CSV, and a server-less file-DAG runner (**Snakemake**, `workflow/Snakefile`, `plans/32`) drives the
static slices (`make figures`). Two larger new-domain efforts remain scoped in
`plans/28-breakthrough-listen-singlepulse.md` and
`plans/29-lotss-deep-144mhz-counts.md`. (The `stokesv` RACS Stokes-V slice, once CASDA-blocked, is now
complete — CASDA's cutout service recovered and the forced-photometry leg is wired.)

### Papers

Every slice is written up as its own **AASTeX paper** under `papers/<slice>/` (authored by Joseph
Barbere, with Claude credited via an AI-use disclosure + a `\software{}` citation — an AI/LLM is not
an eligible author). Each is built reproducibly with containerized tectonic, takes every headline
number from a pipeline-generated `generated/macros.tex` (no figure typed by hand), and is honest
about what it is — mostly recover-a-known validations and methodology, with two honest negatives:

| Paper | `papers/…` | Framing |
|-------|-----------|---------|
| FRB burst statistics, validated on CHIME/FRB Cat 1 | `frbstats/` | validation (tested, reproducible tool) |
| Recovering FRB 20180916B's 16.35-day period | `frbperiod/` | validation |
| The flat inner Milky Way rotation curve from LAB HI | `hi/` | validation |
| A CPU-only SETI drift-search benchmark + Voyager-1 null | `driftsearch/` | benchmark + honest negative |
| TGSS×NVSS USS selection is dominated by the flux scale | `spectra/` | cautionary negative |
| VLASS multi-epoch variability: a 703 deg² census + FK Com | `vlass/` | methodology + validation (recovers FK Com) |
| Three-frequency curvature selection of peaked-spectrum sources | `peaked/` | methodology + two recover-a-known validations |
| Measuring the turnover: southern peaked sources from GLEAM-X + RACS | `southern/` | methodology + measured-turnover candidate list |
| The AGN radio–optical offset excess and its alignment with the jet (ICRF3 × Gaia × MOJAVE) | `offsets/` | reproduction (excess + jet alignment) + reproducible catalogue |
| The steep radio spectra of pulsars from the ATNF catalogue | `pulsarspec/` | reproduction + MSP/normal comparison |
| Image-plane stacking with injection-recovery: SDSS quasars in VLASS-SE | `stacking/` | methodology + calibrated population-mean flux |
| Multi-decade parsec-scale VLBI variability from Astrogeo | `vlbi/` | control-floor method + recover-a-known (OJ 287, BL Lac) |
| A solar type III exciter speed from an e-Callisto dynamic spectrum | `solarbursts/` | method + recover-a-known |
| The Galactic Faraday rotation sky from the Taylor et al. 2009 RM catalogue | `rmsky/` | reproduction (plane enhancement + sign antisymmetry) |
| The pulsar P–Ṗ diagram from the ATNF catalogue | `ppdot/` | reproduction + population classes |
| Tracking an inner-heliosphere type III beam with Wind/WAVES | `windwaves/` | method + recover-a-known (to the Alfvén surface) |
| Tracking a type III beam to 0.4 AU with STEREO/WAVES | `swaves/` | method + recover-a-known (genuinely interplanetary) |
| 3D triangulation of a type III source with STEREO-A+B direction-finding | `triangulate/` | method + independent geometric-vs-plasma distance cross-check |
| Recovering the canonical 1.4 GHz Euclidean source counts from NVSS | `sourcecounts/` | reproduction (Hopkins 2003) |
| A type III beam from the corona to 0.4 AU, geometrically validated | `type3synthesis/` | synthesis + same-event geometric check on the model distance |
| A reproducible RACS Stokes-V coherent-emitter pipeline (and single-epoch limits) | `stokesv/` | tooling + honest single-epoch/variability result |
| A streaming e-Callisto burst-ingest pipeline with cross-station coincidence QC | `ecallisto_pipeline/` | automation pattern + coincidence-vetted burst events (rejects single-station RFI) |
| A coverage-corrected type III occurrence census (method + recover-a-known) | `ecallisto_census/` | census statistic + recover-a-known validation toward a multi-cycle census |
| A pure-PyTorch Fast DM Transform (device-portable dedispersion) | `torchfdmt/` | tool + oracle validation + real Crab recover-a-known + honest CPU/GPU benchmark |
| Two-epoch RACS Stokes-V forced photometry of nearby M dwarfs | `stokesv_discovery/` | method + GJ 65 variability recovery + upper-limit census |

`make paper` builds every slice's PDF; `make arxiv` runs the bundled **`arxiv-submit` skill**
(`.claude/skills/arxiv-submit/`) to assemble and validate an upload package per paper
(`papers/<slice>/arxiv-submission/`: the LaTeX-source tarball with its `.bbl`, plus a `metadata.yaml`
capturing every arXiv submission property and a `CHECKLIST.md`). The orchestration is right-sized: the
static slices build through a server-less **Snakemake** file-DAG (`workflow/Snakefile`, run by
`make figures`), while a **streaming Apache Airflow pipeline on rootless Podman** (`airflow/`) ingests
the frequently-updated e-Callisto archive.

### Where to publish (and where not to)

These papers are mostly **reproductions and honest negatives**, so the venue is matched to the actual
contribution — the *tooling and reproducibility*, not a novelty claim:

- **Software / citable archive:** the toolkit is meant for [JOSS](https://joss.theoj.org) (see
  `joss/paper.md`) and a [Zenodo](https://zenodo.org) DOI on release (`.zenodo.json`, `CITATION.cff`).
- **A short note in the literature:** the frbstats validation is condensed to a
  [Research Note of the AAS](https://journals.aas.org/research-notes/) (`papers/frbstats/rnaas.tex`,
  built by `make paper`).
- **arXiv:** reserved for the papers with a genuinely fresh angle — `type3synthesis/` (the lead: a
  type III beam tracked corona→0.4 AU across four instruments with an independent geometric check on the
  density-model distance), `vlass/` (a 703 deg² census with a real recovery, FK Comae Berenices, plus
  the QL-systematics methodology), `peaked/` (a three-frequency curvature selector with two
  recover-a-known validations and the TGSS-upper-limit + resolution-floor method), `triangulate/` (the
  3D triangulation it builds on), and `ecallisto_pipeline/` (the streaming Airflow-on-Podman ingest
  pattern, astro-ph.IM). `frbstats/` is **not** an arXiv target — its only fresh angle was the Airflow
  pattern, now re-homed on the streaming archive; it stays a JOSS/RNAAS tool. The pure
  reproductions/negatives are **not** posted as a preprint batch — arXiv moderation expects a
  contribution, and "I reproduced a known result" belongs in the repo + Zenodo.

## How it relates to `jansky`

This repo **depends on `jansky` as a library** and reuses its tested helpers (`jansky.transients`,
`jansky.rfi`, `jansky.timing`, `jansky.seti`, `jansky.sourcecounts`, `jansky.formats`,
`jansky.data`, …) rather than reimplementing them. It mirrors jansky's conventions: `uv`-managed,
ruff + mypy + pytest with an 85% coverage floor, Podman containers, and a `plans/NN-slug.md`
workflow. The `jansky` dependency is a local path source (`../jansky`) for cross-repo dev, switching
to the pinned git tag `jansky@v0.1.0` for clean-checkout CI. See `pyproject.toml`.

## Quickstart

```bash
# Requires the jansky repo checked out next to this one (../jansky) for local dev.
uv sync                                   # env + jansky (from ../jansky)
make test                                 # unit tests (offline, on synthetic fixtures)
make cov                                  # tests + 85% coverage floor

# Run a slice on real public data (each writes results/ + figures + macros into papers/<slice>/):
uv run python -m jansky_research.pipeline     # FRB burst statistics (CHIME catalog)
uv run python -m jansky_research.frbperiod    # FRB repeater periodicity
uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3   # USS hunt
uv run python -m jansky_research.driftsearch  # SETI injection-recovery benchmark
uv run python -m jansky_research.hi           # Milky Way HI rotation curve
uv run python -m jansky_research.vlass --ra 190 --dec 20 --radius 15  # VLASS variability census (needs --extra vlass)
# (append --offline to run any slice on its synthetic fixture, no network)

# The papers + orchestration:
make figures                              # build every static slice via the Snakemake DAG (offline; needs --extra workflow)
make paper                                # tectonic -> all papers/<slice>/main.pdf (in a container)
make arxiv                                # assemble + validate an arXiv package per paper
make reproduce                            # fetch -> figures -> papers -> arXiv packages, end to end

# Streaming ingest (the e-Callisto archive) runs on Airflow + Podman:
make airflow-up COMPOSE="uvx podman-compose" && make dag-test DATE=2011-09-14
make ecallisto-day DATE=20110914          # the same day's scan WITHOUT Airflow (the shared worker)
```

See `REPRODUCING.md` for the full reproduction, the right-sized-orchestration notes
(Snakemake static / Airflow streaming), and offline mode.

## Layout

```
jansky-research/
  src/jansky_research/   # the tooling package (tested-helper pattern, 85% floor) — one module per slice
    data.py              # dataset registry + offline synthetic fallback
    frbstats.py spectra.py frbperiod.py driftsearch.py hi.py vlass.py peaked.py southern.py
    offsets.py pulsarspec.py stacking.py vlbi.py solarbursts.py rmsky.py ppdot.py
    windwaves.py swaves.py triangulate.py sourcecounts.py type3synthesis.py
    ecallisto_catalog.py stokesv.py   # e-Callisto streaming worker; RACS Stokes-V coherent emitters
    pipeline.py          # the FRB pipeline (shared by Make / notebook / Snakemake)
    report.py            # figure/macro emitters -> paper inputs
  survey/                # PERMANENT: literature.md, github-landscape.md, gap-analysis.md,
                         #   candidate-gaps.md + *-scan.md (backlog), and each slice's *-findings.md
  workflow/              # Snakefile: the server-less file-DAG that builds the static slices' inputs
  airflow/               # Airflow-on-Podman stack + the streaming e-Callisto ingest DAG
  papers/<slice>/        # one AASTeX paper per slice (main.tex + refs.bib tracked;
                         #   figures/, generated/, arxiv-submission/ are produced by make)
                         #   frbstats/ also has rnaas.tex (a Research Note of the AAS)
  joss/                  # JOSS software paper (paper.md + paper.bib)
  CITATION.cff           # "Cite this repository"; .zenodo.json drives Zenodo archival
  containers/            # tectonic paper-build image
  .claude/skills/        # arxiv-submit (assemble + validate an arXiv upload package)
  .claude/agents/        # dataset-analyst, pipeline-runner, results-interpreter (+ reused jansky)
  plans/                 # numbered slice specs (00-29); the lasting record is each slice's
                         #   survey/*-findings.md + papers/<slice>/. #28-29 are planned, not started
```

## Method & gates

Each slice follows the same discipline, recorded as a `plans/NN-slug.md` spec: build a tested tool →
run it on real public data → **science-reviewer gate** (an adversarial pass that has caught real
blockers every time) → honest write-up. Findings are framed as *candidates / validations / limits*,
never dressed-up discoveries; every reported number is reproducible. The non-chosen survey gaps are
preserved as a backlog in `survey/candidate-gaps.md`.

## License

MIT — see [LICENSE](LICENSE).
