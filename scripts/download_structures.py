#!/usr/bin/env python3
"""
Download kinase structures from the RCSB PDB REST API.

Writes mmCIF (default) or PDB-format files to an output directory.
Fails loudly on invalid or unresolvable PDB IDs — no silent skipping.

Usage examples
--------------
# Download two structures in mmCIF format:
python3 scripts/download_structures.py \\
    --pdb_ids 1ATP 2ITO \\
    --output_dir benchmark/ground_truth

# Download all pilot targets listed in metadata file:
python3 scripts/download_structures.py \\
    --metadata_file benchmark/metadata/kinase_pilot_targets.json \\
    --output_dir benchmark/ground_truth \\
    --format pdb

# Download without overwriting existing files:
python3 scripts/download_structures.py \\
    --pdb_ids 1ATP \\
    --output_dir benchmark/ground_truth \\
    --skip_existing
"""

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RCSB_MMCIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

VALID_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

RETRY_DELAYS = [2, 4, 8, 16]  # seconds between retries


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def validate_pdb_id(pdb_id: str) -> str:
    """Return the uppercase PDB ID after validating its format.

    Raises:
        ValueError: if the ID does not match the PDB format (digit + 3 alphanumeric).
    """
    pdb_id = pdb_id.strip().upper()
    if not VALID_PDB_ID.match(pdb_id):
        raise ValueError(
            f"Invalid PDB ID format: '{pdb_id}'. "
            "Expected one digit followed by three alphanumeric characters (e.g. '1ATP')."
        )
    return pdb_id


def download_url(url: str, dest: Path, retries: int = 3) -> None:
    """Download *url* to *dest*, retrying on transient network errors.

    Raises:
        urllib.error.HTTPError: on HTTP 404 (invalid PDB ID) — not retried.
        RuntimeError:           if all retries are exhausted.
    """
    for attempt, delay in enumerate([0] + RETRY_DELAYS[:retries]):
        if delay:
            print(f"  Retry {attempt}/{retries} after {delay}s…", file=sys.stderr)
            time.sleep(delay)
        try:
            urllib.request.urlretrieve(url, dest)
            return
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise  # Hard failure — PDB ID does not exist
            if attempt >= retries:
                raise RuntimeError(
                    f"Failed to download {url} after {retries + 1} attempts: {exc}"
                ) from exc
            # Other HTTP errors: retry
        except urllib.error.URLError as exc:
            if attempt >= retries:
                raise RuntimeError(
                    f"Network error downloading {url}: {exc}"
                ) from exc


def download_structure(
    pdb_id: str,
    output_dir: Path,
    fmt: str = "mmcif",
    skip_existing: bool = False,
) -> Path:
    """Download a single PDB structure.

    Args:
        pdb_id:       PDB accession code (e.g. '1ATP').
        output_dir:   Directory where the file will be saved.
        fmt:          'mmcif' (saves as <ID>.cif) or 'pdb' (saves as <ID>.pdb).
        skip_existing: If True and the target file already exists, skip download.

    Returns:
        Path to the downloaded file.

    Raises:
        ValueError:             on invalid PDB ID format.
        urllib.error.HTTPError: if the PDB ID is not found on RCSB (HTTP 404).
        RuntimeError:           on persistent network failure.
    """
    pdb_id = validate_pdb_id(pdb_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    if fmt == "mmcif":
        url = RCSB_MMCIF_URL.format(pdb_id=pdb_id)
        dest = output_dir / f"{pdb_id}.cif"
    elif fmt == "pdb":
        url = RCSB_PDB_URL.format(pdb_id=pdb_id)
        dest = output_dir / f"{pdb_id}.pdb"
    else:
        raise ValueError(f"Unknown format '{fmt}'. Choose 'mmcif' or 'pdb'.")

    if skip_existing and dest.exists():
        print(f"  {pdb_id}: already exists at {dest}, skipping.")
        return dest

    print(f"  Downloading {pdb_id} ({fmt}) → {dest} …")
    try:
        download_url(url, dest)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            sys.exit(
                f"ERROR: PDB ID '{pdb_id}' not found on RCSB (HTTP 404). "
                f"Check that the accession code is correct."
            )
        raise

    print(f"  {pdb_id}: OK ({dest.stat().st_size // 1024} KB)")
    return dest


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download kinase structures from RCSB PDB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--pdb_ids",
        nargs="+",
        metavar="PDB_ID",
        help="One or more PDB IDs to download (e.g. 1ATP 2ITO).",
    )
    source_group.add_argument(
        "--metadata_file",
        metavar="FILE",
        help=(
            "Path to kinase_pilot_targets.json. "
            "Downloads all entries with inclusion_status='included'."
        ),
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        metavar="DIR",
        help="Directory where downloaded files will be saved.",
    )
    parser.add_argument(
        "--format",
        choices=["mmcif", "pdb"],
        default="mmcif",
        help="File format to download (default: mmcif).",
    )
    parser.add_argument(
        "--skip_existing",
        action="store_true",
        help="Skip download if the output file already exists.",
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    # Collect PDB IDs
    if args.pdb_ids:
        pdb_ids = args.pdb_ids
    else:
        meta_path = Path(args.metadata_file)
        if not meta_path.exists():
            sys.exit(f"ERROR: Metadata file not found: {meta_path}")
        with open(meta_path) as fh:
            metadata = json.load(fh)
        targets = metadata.get("targets", [])
        pdb_ids = [
            t["pdb_id"]
            for t in targets
            if t.get("inclusion_status") == "included"
            and not t.get("engineering_test_fixture", False)
        ]
        if not pdb_ids:
            print(
                "No targets with inclusion_status='included' found in metadata. "
                "Nothing to download.",
                file=sys.stderr,
            )
            return

    print(f"Downloading {len(pdb_ids)} structure(s) → {output_dir}/")

    errors: list[tuple[str, str]] = []
    for pdb_id in pdb_ids:
        try:
            download_structure(
                pdb_id,
                output_dir,
                fmt=args.format,
                skip_existing=args.skip_existing,
            )
        except Exception as exc:
            print(f"  ERROR for {pdb_id}: {exc}", file=sys.stderr)
            errors.append((pdb_id, str(exc)))

    print(f"\nDone. {len(pdb_ids) - len(errors)}/{len(pdb_ids)} downloaded successfully.")

    if errors:
        print("Failures:", file=sys.stderr)
        for pdb_id, msg in errors:
            print(f"  {pdb_id}: {msg}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
