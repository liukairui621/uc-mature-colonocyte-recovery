"""Create a relocatable working copy without changing the audited source snapshot.

The executed project retained its original absolute analysis root in scripts and
provenance records. This helper copies text assets to a user-selected directory
and rewrites that root in the copy. Frozen hashes in the portable copy will no
longer match the audit snapshot; use the original repository for provenance.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


OLD_ROOT = "/root/projects/UC_Treatment_Recovery"
TEXT_SUFFIXES = {".R", ".r", ".py", ".json", ".md", ".txt", ".tsv", ".csv", ".yml", ".yaml"}
SKIP_DIRS = {".git", "inputs", "runs", "logs", "renv", "submission_bmc_genomics", "submission_work"}


def ignored(_path: str, names: list[str]) -> set[str]:
    return {name for name in names if name in SKIP_DIRS}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="New empty directory for the portable working copy")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    target = args.target.expanduser().resolve()
    if target.exists():
        raise SystemExit(f"Target already exists: {target}")
    shutil.copytree(source, target, ignore=ignored)
    new_root = target.as_posix()
    changed = 0
    for path in target.rglob("*"):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if OLD_ROOT in text:
            path.write_text(text.replace(OLD_ROOT, new_root), encoding="utf-8")
            changed += 1
    print(f"Portable copy created at {target}; rewrote {changed} text files.")


if __name__ == "__main__":
    main()
