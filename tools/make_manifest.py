#!/usr/bin/env python3
"""Write MANIFEST.sha256: one SHA-256 line per git-tracked file.

The manifest identifies the released files byte for byte.  Regenerate
after any change to the repository contents, and before tagging a
release:

    python3 tools/make_manifest.py

Verify a checkout against the manifest:

    sha256sum -c MANIFEST.sha256
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "MANIFEST.sha256"


def main() -> int:
    tracked = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True,
        text=True, check=True).stdout.split("\0")
    lines = []
    for name in sorted(filter(None, tracked)):
        if name == MANIFEST.name:
            continue
        digest = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        lines.append(f"{digest}  {name}")
    MANIFEST.write_text("\n".join(lines) + "\n")
    print(f"{MANIFEST.name}: {len(lines)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
