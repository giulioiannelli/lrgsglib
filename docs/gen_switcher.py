#!/usr/bin/env python3
"""Generate the pydata-sphinx-theme version switcher JSON from git tags.

The version switcher dropdown in the docs navbar is populated from a
``switcher.json`` file. This script builds that file from the repository's
release tags (``vX.Y.Z``), so the version list stays in sync with GitHub with
no manual editing.

Layout produced (one entry per release tag, newest first, plus a ``dev``
entry for the in-development docs built from ``main``)::

    [
      {"name": "dev",          "version": "dev",   "url": ".../dev/"},
      {"name": "0.2.0 (stable)","version": "0.2.0", "url": ".../0.2.0/", "preferred": true},
      {"name": "0.1.0",        "version": "0.1.0", "url": ".../0.1.0/"}
    ]

Usage::

    python docs/gen_switcher.py --base-url https://user.github.io/lrgsglib \
        --output gh-pages/switcher.json

The ``url`` of each entry is ``<base-url>/<version>/`` (and ``<base-url>/dev/``
for the dev entry). The newest release is tagged ``(stable)`` and marked
``preferred`` so the theme highlights it as the default.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys

# Match annotated/lightweight release tags like "v0.2.0" or "v1.10.3".
TAG_RE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")


def git_release_tags() -> list[tuple[int, int, int]]:
    """Return sorted (newest-first) release versions parsed from git tags."""
    try:
        out = subprocess.run(
            ["git", "tag", "--list", "v*"],
            check=True, capture_output=True, text=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    versions = []
    for line in out.splitlines():
        m = TAG_RE.match(line.strip())
        if m:
            versions.append(tuple(int(g) for g in m.groups()))
    return sorted(set(versions), reverse=True)


def build_switcher(base_url: str) -> list[dict]:
    base = base_url.rstrip("/")
    entries: list[dict] = [
        {"name": "dev", "version": "dev", "url": f"{base}/dev/"},
    ]
    for i, (major, minor, patch) in enumerate(git_release_tags()):
        v = f"{major}.{minor}.{patch}"
        entry: dict = {"name": v, "version": v, "url": f"{base}/{v}/"}
        if i == 0:  # newest release = stable / preferred
            entry["name"] = f"{v} (stable)"
            entry["preferred"] = True
        entries.append(entry)
    return entries


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-url", required=True,
                   help="Deployed site root, e.g. https://user.github.io/lrgsglib")
    p.add_argument("--output", default="-",
                   help="Output path for switcher.json ('-' for stdout)")
    args = p.parse_args()

    data = build_switcher(args.base_url)
    text = json.dumps(data, indent=2)
    if args.output == "-":
        print(text)
    else:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
        print(f"wrote {len(data)} entries to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
