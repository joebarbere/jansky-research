# Versioning policy

`jansky-research` follows [Semantic Versioning 2.0.0](https://semver.org/). This file defines
**what the public API is** and gives a deterministic recipe for the next version number, so any
contributor (or Claude session) can decide a bump without a judgement call.

## Single source of truth

The version lives in two tracked files and **must match**:

- `pyproject.toml` → `[project] version`
- `CITATION.cff` → `version`

Zenodo takes the version from the **git release tag** (`vX.Y.Z`), not from a file, so the tag must
agree with those two. `.zenodo.json` deliberately carries **no** `version` field (Zenodo fills it
from the tag); leave it that way unless a version field appears there.

`CHANGELOG.md` (Keep a Changelog format) is the human-readable record; its `## [Unreleased]`
section is the input to the next-version recipe below.

## What the public API is

A change is **breaking** only if it breaks one of these:

- the `python -m jansky_research.<slice>` command-line entry points (their names, required
  arguments, and documented flags);
- each slice module's `__all__` and its `run()` signature;
- the public functions of `data.py`, `pipeline.py`, and `report.py`;
- the schema of any results file consumed downstream.

Papers, notebooks, findings, figures, macros, tests, CI config, and internal helpers are **not**
part of the public API — changing them is at most a PATCH.

## Bump rules

| Bump | When | Examples |
|------|------|----------|
| **MAJOR** `X.0.0` | Backward-incompatible public-API change | remove or rename a slice module or CLI; change a `run()` / public-function signature incompatibly; drop a supported Python version; change a downstream-consumed results-file schema |
| **MINOR** `1.X.0` | Backward-compatible addition | a new slice / module; a new public function or `--flag`; a new optional extra; additive behaviour |
| **PATCH** `1.0.X` | No public-API change | bug fix; paper / doc / test-only change; a dependency pin with no API effect; figure / macro regeneration |

## The next-version recipe

`scripts/next_version.py` automates this; the logic is:

1. `git describe --tags --abbrev=0` → the last released tag (`v1.0.0` after the initial release).
2. Read `CHANGELOG.md`'s `## [Unreleased]` section (cross-check against
   `git log <last-tag>..HEAD` for anything missed).
3. Classify each entry by its `###` subsection: **Added / Changed / Deprecated / Removed / Fixed /
   Security**. An entry is breaking if it is under **Removed**, or under **Changed** and its text
   begins with `**BREAKING**`.
4. **Highest severity present wins:**
   - any breaking entry → **MAJOR**;
   - else any **Added** entry → **MINOR**;
   - else → **PATCH**.

Run it any time:

```console
$ python scripts/next_version.py
next: 1.1.0 (MINOR — Unreleased has Added entries and no breaking changes)
```

It exits non-zero if `## [Unreleased]` is empty (nothing to release), so it doubles as a
"did this PR forget its changelog entry?" check.

## Initial release

There is no prior tag before `v1.0.0`, so the recipe's step 1 has nothing to read. The first
release is a **deliberate `1.0.0`** (not `0.x`): publishing the JOSS paper commits to a public API
and asserts a real, citable scholarly tool, which is exactly the SemVer 1.0.0 moment. When no tag
exists, `scripts/next_version.py` says so and recommends `1.0.0`. Every bump after that follows the
recipe mechanically.

## Per-PR obligation

Every PR that changes anything a user could observe **adds an entry** to the `## [Unreleased]`
section of `CHANGELOG.md`, under the right subsection. Paper/doc-only or test-only PRs still add a
`Fixed`/`Changed` line (they are PATCH-level, but the changelog stays complete). This keeps the
next-version recipe accurate.
