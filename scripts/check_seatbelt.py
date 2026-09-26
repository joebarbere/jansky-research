#!/usr/bin/env python
"""Check that the AWS cost seatbelt has not drifted (plan 96).

``infra/seatbelt.json`` is the source of truth. Two things can drift from it:

1. the copy in the sibling ``aws-ai`` repo (``../aws-ai/iam/cost-seatbelt.json``), and
2. the inline policy actually attached to the ``AdministratorAccess`` permission set.

Both are compared as parsed JSON, so whitespace and key order do not matter. The live check
needs an Identity Center session (``aws sso login --profile joebarbere-admin``); without one it
is reported as skipped, not passed. Exit status is non-zero on any drift.

    make seatbelt-check
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TRUTH = REPO / "infra" / "seatbelt.json"
COPY = REPO.parent / "aws-ai" / "iam" / "cost-seatbelt.json"
PERMISSION_SET = "AdministratorAccess"


def _aws(*args: str) -> str:
    env = {**os.environ, "AWS_PROFILE": os.environ.get("AWS_PROFILE", "joebarbere-admin")}
    out = subprocess.run(
        ["aws", *args, "--output", "text"], capture_output=True, text=True, env=env, check=True
    )
    return out.stdout.strip()


def live_policy() -> dict | None:
    """The inline policy on the permission set, or None if AWS is unreachable."""
    try:
        inst = _aws("sso-admin", "list-instances", "--query", "Instances[0].InstanceArn")
        for ps in _aws(
            "sso-admin",
            "list-permission-sets",
            "--instance-arn",
            inst,
            "--query",
            "PermissionSets[]",
        ).split():
            name = _aws(
                "sso-admin",
                "describe-permission-set",
                "--instance-arn",
                inst,
                "--permission-set-arn",
                ps,
                "--query",
                "PermissionSet.Name",
            )
            if name == PERMISSION_SET:
                raw = _aws(
                    "sso-admin",
                    "get-inline-policy-for-permission-set",
                    "--instance-arn",
                    inst,
                    "--permission-set-arn",
                    ps,
                    "--query",
                    "InlinePolicy",
                )
                return json.loads(raw) if raw and raw != "None" else {}
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"  live: SKIPPED (no AWS session?): {exc}", file=sys.stderr)
        return None
    return {}


def main() -> int:
    truth = json.loads(TRUTH.read_text())
    drift = False
    if COPY.exists():
        ok = json.loads(COPY.read_text()) == truth
        drift |= not ok
        print(f"  aws-ai copy: {'OK' if ok else 'DRIFTED'} ({COPY})")
    else:
        print(f"  aws-ai copy: SKIPPED (not checked out at {COPY})")
    live = live_policy()
    if live is not None:
        ok = live == truth
        drift |= not ok
        print(f"  live {PERMISSION_SET} inline policy: {'OK' if ok else 'DRIFTED'}")
    return 1 if drift else 0


if __name__ == "__main__":
    raise SystemExit(main())
