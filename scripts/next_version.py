#!/usr/bin/env python3
"""Recommend the next release version from CHANGELOG.md's ``## [Unreleased]`` section.

Implements the recipe in ``VERSIONING.md``: read the last git tag, classify the pending
``Unreleased`` entries (Added / Changed / Removed / ...), and apply "highest severity wins" —
any breaking entry -> MAJOR, else any Added -> MINOR, else PATCH. Prints the recommendation and
its reasoning. Exits non-zero when there is nothing to release, so it doubles as a
"did this PR forget its changelog entry?" check.

Pure stdlib. The parsing / decision logic lives in importable functions
(``parse_unreleased`` / ``recommend``) so it is unit-tested without touching git.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CHANGELOG = REPO_ROOT / "CHANGELOG.md"

# Keep a Changelog subsections, in the order we report them.
SECTIONS = ("Added", "Changed", "Deprecated", "Removed", "Fixed", "Security")

_SECTION_RE = re.compile(r"^###\s+(\w+)\s*$")
_BULLET_RE = re.compile(r"^\s*[-*]\s+(.+)$")


def parse_unreleased(text: str) -> dict[str, list[str]]:
    """Return ``{section: [bullet, ...]}`` for the ``## [Unreleased]`` block of *text*.

    Only the six Keep a Changelog subsections are collected; empty sections are omitted.
    Sub-bullets (indented continuation lines) are ignored — only top-level entries count.
    """
    lines = text.splitlines()
    # Find the Unreleased heading and the next top-level (## ) heading after it.
    start = None
    for i, line in enumerate(lines):
        if re.match(r"^##\s+\[Unreleased\]", line, re.IGNORECASE):
            start = i + 1
            break
    if start is None:
        return {}
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("## ") and not lines[i].startswith("###"):
            end = i
            break

    found: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines[start:end]:
        sec_match = _SECTION_RE.match(line)
        if sec_match:
            name = sec_match.group(1).capitalize()
            current = name if name in SECTIONS else None
            continue
        if current is None:
            continue
        # Only un-indented bullets are top-level entries.
        if line and not line[0].isspace():
            bullet = _BULLET_RE.match(line)
            if bullet:
                found.setdefault(current, []).append(bullet.group(1).strip())
    return {k: v for k, v in found.items() if v}


def _is_breaking(sections: dict[str, list[str]]) -> bool:
    """A change is breaking if it is Removed, or a Changed entry marked ``**BREAKING**``."""
    if sections.get("Removed"):
        return True
    return any(entry.lstrip().startswith("**BREAKING**") for entry in sections.get("Changed", []))


def recommend(last_version: str | None, sections: dict[str, list[str]]) -> tuple[str, str]:
    """Return ``(next_version, reason)`` given the last released version and pending sections.

    ``last_version`` is ``None`` before the first tag; the initial release is a deliberate 1.0.0.
    Raises ``ValueError`` when there is nothing to release.
    """
    if not sections:
        raise ValueError("no entries under ## [Unreleased] — nothing to release")

    if last_version is None:
        return "1.0.0", "initial release — no prior tag (see VERSIONING.md)"

    major, minor, patch = (int(p) for p in last_version.lstrip("v").split("."))
    present = ", ".join(s for s in SECTIONS if s in sections)
    if _is_breaking(sections):
        return f"{major + 1}.0.0", f"MAJOR — Unreleased has a breaking change ({present})"
    if sections.get("Added"):
        return (
            f"{major}.{minor + 1}.0",
            f"MINOR — Unreleased has Added entries, no breaking changes ({present})",
        )
    return f"{major}.{minor}.{patch + 1}", f"PATCH — no public-API change ({present})"


def last_tag() -> str | None:
    """The most recent ``vX.Y.Z`` git tag, or ``None`` if the repo has no tags yet."""
    try:
        out = subprocess.run(
            ["git", "describe", "--tags", "--abbrev=0"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    tag = out.stdout.strip()
    return tag or None


def main() -> int:
    if not CHANGELOG.exists():
        print(f"error: {CHANGELOG} not found", file=sys.stderr)
        return 2
    sections = parse_unreleased(CHANGELOG.read_text(encoding="utf-8"))
    try:
        version, reason = recommend(last_tag(), sections)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"next: {version} ({reason})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
