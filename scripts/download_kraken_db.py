#!/usr/bin/env python3
"""Resolve a ``KRAKEN_DB`` config value to a usable on-disk path.

Resolution rules (in order):
  1. If the value already exists on disk, return it verbatim.
  2. Otherwise, if it looks like an HTTP(S) URL ending in ``.tar.gz``,
     download the archive with ``wget`` (or ``curl`` as a fallback),
     extract it into a directory named after the tarball stem, and
     return that directory.
  3. Otherwise, raise a clear error explaining the two accepted forms.

Re-running with the same URL is idempotent: if the destination directory
already contains a built Kraken2 database (``hash.k2d``), the download
and extract steps are skipped.

Usage:
    # As a library (preferred — called from the Snakefile):
    from download_kraken_db import resolve_kraken_db
    db_path = resolve_kraken_db(config.get("KRAKEN_DB"))

    # As a CLI (handy for debugging):
    python scripts/download_kraken_db.py <path-or-URL>
        # → prints the resolved local path on stdout
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

_TAR_GZ_URL_RE = re.compile(r"^https?://.+\.tar\.gz$", re.IGNORECASE)
_BUILT_DB_MARKER = "hash.k2d"  # canonical "kraken2-build finished" sentinel


def _looks_like_url(value: str) -> bool:
    return bool(_TAR_GZ_URL_RE.match(value))


def _derive_dbname(url: str) -> str:
    """Return the segment between the final ``/`` and ``.tar.gz`` in *url*."""
    last_segment = url.rsplit("/", 1)[-1]
    if not last_segment.lower().endswith(".tar.gz"):
        return ""
    return last_segment[: -len(".tar.gz")]


def _download(url: str, dest: Path) -> None:
    """Fetch *url* into *dest*, preferring ``wget -c`` for resume support."""
    if shutil.which("wget"):
        cmd = ["wget", "-c", "-O", str(dest), url]
    elif shutil.which("curl"):
        cmd = ["curl", "-L", "--fail", "-C", "-", "-o", str(dest), url]
    else:
        raise RuntimeError(
            "Neither 'wget' nor 'curl' is available; cannot download "
            f"Kraken2 database from {url}"
        )
    print(f"[KRAKEN_DB] Downloading {url} → {dest}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def _extract(tarball: Path, dest_dir: Path) -> None:
    print(f"[KRAKEN_DB] Extracting {tarball} → {dest_dir}/", file=sys.stderr)
    dest_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["tar", "-xvzf", str(tarball), "-C", str(dest_dir)],
        check=True,
    )


def resolve_kraken_db(kraken_db: str | None) -> str:
    """Resolve *kraken_db* to a local directory path.

    See module docstring for the rules. Raises ``ValueError`` /
    ``FileNotFoundError`` / ``RuntimeError`` for the various failure modes.
    """
    if not kraken_db:
        raise ValueError(
            "KRAKEN_DB is required. Provide either:\n"
            "  - a path to an existing Kraken2 database directory, or\n"
            "  - an HTTP(S) URL to a Kraken2 .tar.gz archive (will be "
            "downloaded automatically), e.g.\n"
            "    https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz"
        )

    # 1. Already on disk → use as-is.
    if Path(kraken_db).exists():
        return kraken_db

    # 2. Looks like a downloadable archive → fetch + extract.
    if _looks_like_url(kraken_db):
        dbname = _derive_dbname(kraken_db)
        if not dbname:
            raise ValueError(
                f"Could not derive a directory name from URL: {kraken_db!r}"
            )

        db_dir = Path(dbname)
        # Idempotency: skip the download if the DB is already extracted.
        if (db_dir / _BUILT_DB_MARKER).is_file():
            return str(db_dir)

        tarball = Path(f"{dbname}.tar.gz")
        _download(kraken_db, tarball)
        _extract(tarball, db_dir)

        if not (db_dir / _BUILT_DB_MARKER).is_file():
            raise FileNotFoundError(
                f"Extraction completed but '{_BUILT_DB_MARKER}' was not found "
                f"in '{db_dir}/'. The archive at {kraken_db} does not appear "
                "to contain a built Kraken2 database."
            )

        # Free the disk now that the DB is extracted.
        try:
            tarball.unlink()
        except OSError:
            pass

        print(f"[KRAKEN_DB] Database ready at '{db_dir}/'", file=sys.stderr)
        return str(db_dir)

    # 3. Neither a real path nor a downloadable URL → friendly error.
    raise FileNotFoundError(
        f"KRAKEN_DB={kraken_db!r} does not exist on disk and is not a "
        "downloadable Kraken2 database URL.\n"
        "Provide either:\n"
        "  - a path to an existing Kraken2 database directory, or\n"
        "  - an HTTP(S) URL ending in '.tar.gz', e.g.\n"
        "    https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz"
    )


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print(
            "Usage: download_kraken_db.py <KRAKEN_DB path or URL>",
            file=sys.stderr,
        )
        return 2
    try:
        print(resolve_kraken_db(argv[0]))
        return 0
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed (exit {exc.returncode}): "
            f"{' '.join(map(str, exc.cmd))}",
            file=sys.stderr,
        )
        return exc.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
