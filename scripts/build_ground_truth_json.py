#!/usr/bin/env python3
"""
Convert downloaded mmCIF / PDB structure files into the ground-truth JSON format
expected by generate_sar.py's `structure_only` metrics mode.

For each input structure file, this script:
  1. Parses the file with Biopython (MMCIFParser or PDBParser).
  2. Iterates the first MODEL only.
  3. Collects C-alpha atoms from STANDARD amino-acid residues across all
     (or specified) chains, in chain-then-sequence-number order.
  4. Writes `{PDB_ID}_ground_truth.json` containing the coordinate list
     and extraction provenance.

Output JSON schema (consumed by generate_sar.py::load_ground_truth):
  {
    "pdb_id":          str,            # uppercase PDB accession
    "source_file":     str,            # path to the input structure file
    "source_format":   "mmcif"|"pdb",
    "chains_included": [str, ...],     # chain IDs whose residues are included
    "n_residues":      int,            # number of C-alpha positions extracted
    "extraction_notes": str,
    "coordinates":     [[x, y, z], ...] # float, Å, ordered by chain then seq
  }

COORDINATE ORDERING AND CHAIN HANDLING
---------------------------------------
By default ALL chains are included, concatenated in the order Biopython
iterates them (which follows the order they appear in the structure file).
Use --chains to restrict to a comma-separated list of chain IDs.

For RMSD and contact-map comparison the prediction and ground-truth arrays
are trimmed to their common length (min(N_pred, N_gt)) before alignment.
This is safe for engineering validation but means residue identity is NOT
verified.  Phase 2 will require explicit sequence alignment.

Requirements: biopython >= 1.81, numpy >= 1.24

Usage examples
--------------
# Process all .cif / .pdb files in a directory:
python3 scripts/build_ground_truth_json.py \\
    --structure_dir benchmark/ground_truth \\
    --output_dir benchmark/ground_truth

# Process specific files, keep only chain A:
python3 scripts/build_ground_truth_json.py \\
    --structure_files benchmark/ground_truth/1ATP.cif \\
    --output_dir benchmark/ground_truth \\
    --chains A

# Overwrite existing JSON files:
python3 scripts/build_ground_truth_json.py \\
    --structure_dir benchmark/ground_truth \\
    --output_dir benchmark/ground_truth \\
    --overwrite
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SUPPORTED_EXTENSIONS = {".cif", ".mmcif", ".pdb", ".ent"}


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------


def _format_from_suffix(suffix: str) -> str:
    if suffix in (".cif", ".mmcif"):
        return "mmcif"
    return "pdb"


def _parse_structure(path: Path):
    """Parse a structure file and return a Biopython Structure object."""
    try:
        from Bio.PDB import MMCIFParser, PDBParser
    except ImportError as exc:
        raise ImportError(
            "Biopython is required.  Install with: pip install biopython>=1.81"
        ) from exc

    fmt = _format_from_suffix(path.suffix.lower())
    if fmt == "mmcif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    return parser.get_structure(path.stem.upper(), str(path)), fmt


def extract_ca_coords_ordered(
    structure,
    chains: list[str] | None = None,
) -> tuple[list[list[float]], list[str]]:
    """Extract C-alpha coordinates from the first model of a Biopython Structure.

    Only standard amino-acid residues (het-flag == " ") that contain a "CA"
    atom are included.  HETATM residues (ligands, water) are skipped.

    Args:
        structure: Bio.PDB.Structure.Structure (first model is used)
        chains:    If given, only chains whose ID is in this list are included.
                   If None, all chains are included.

    Returns:
        coords:         List of [x, y, z] floats, in chain-then-sequence order.
        chains_used:    Ordered list of chain IDs that contributed residues.

    Raises:
        ValueError: if no C-alpha atoms are found after applying chain filter.
    """
    import numpy as np  # already a hard dependency

    first_model = next(iter(structure))

    coords: list[list[float]] = []
    chains_with_residues: list[str] = []

    for chain in first_model:
        chain_id = chain.get_id()
        if chains is not None and chain_id not in chains:
            continue

        chain_coords: list[list[float]] = []
        for residue in chain:
            het_flag, _seq, _ins = residue.get_id()
            if het_flag.strip() != "":
                # HETATM — skip (ligands, solvent, modified residues)
                continue
            if "CA" not in residue:
                # Residue present but CA not resolved (e.g., poor density)
                continue
            vec = residue["CA"].get_vector().get_array()
            chain_coords.append([float(vec[0]), float(vec[1]), float(vec[2])])

        if chain_coords:
            coords.extend(chain_coords)
            chains_with_residues.append(chain_id)

    if not coords:
        chain_msg = (
            f"chains={chains}" if chains else "all chains"
        )
        raise ValueError(
            f"No C-alpha atoms found in structure '{structure.id}' "
            f"({chain_msg}).  "
            "Check that the structure contains standard amino-acid residues "
            "and that the requested chain IDs exist."
        )

    return coords, chains_with_residues


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------


def build_ground_truth_json(
    path: Path,
    output_dir: Path,
    chains: list[str] | None,
    overwrite: bool,
) -> Path:
    """Parse one structure file and write its ground-truth JSON.

    Args:
        path:       Path to the input structure file (.cif or .pdb).
        output_dir: Directory where the output JSON is written.
        chains:     Chain filter (None = all chains).
        overwrite:  If False, raise FileExistsError when output already exists.

    Returns:
        Path to the written JSON file.

    Raises:
        FileExistsError: if the output exists and overwrite=False.
        ValueError:      if no C-alpha atoms can be extracted.
        ImportError:     if Biopython is not installed.
        RuntimeError:    if the structure file cannot be parsed.
    """
    pdb_id = path.stem.upper()
    out_file = output_dir / f"{pdb_id}_ground_truth.json"

    if out_file.exists() and not overwrite:
        raise FileExistsError(
            f"{out_file} already exists. Pass --overwrite to replace it."
        )

    structure, fmt = _parse_structure(path)

    coords, chains_used = extract_ca_coords_ordered(structure, chains=chains)

    chain_note = (
        f"chains restricted to {chains}" if chains else "all chains included"
    )
    record = {
        "pdb_id": pdb_id,
        "source_file": str(path.resolve()),
        "source_format": fmt,
        "chains_included": chains_used,
        "n_residues": len(coords),
        "extraction_notes": (
            f"C-alpha atoms only; model 1; {chain_note}; "
            "HETATM and unresolved residues skipped; "
            "produced by scripts/build_ground_truth_json.py"
        ),
        "coordinates": coords,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fh:
        json.dump(record, fh, indent=2)

    return out_file


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert structure files (.cif/.pdb) to ground-truth JSON "
            "for generate_sar.py structure_only mode."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    src_group = parser.add_mutually_exclusive_group(required=True)
    src_group.add_argument(
        "--structure_dir",
        metavar="DIR",
        help="Directory containing .cif / .pdb files to process.",
    )
    src_group.add_argument(
        "--structure_files",
        nargs="+",
        metavar="FILE",
        help="One or more structure files to process.",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        metavar="DIR",
        help="Directory where {PDB_ID}_ground_truth.json files are written.",
    )
    parser.add_argument(
        "--chains",
        metavar="A,B",
        default=None,
        help=(
            "Comma-separated chain IDs to include (e.g. 'A' or 'A,B'). "
            "Default: all chains."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing ground-truth JSON files.",
    )

    args = parser.parse_args()

    # Resolve chain filter
    chains: list[str] | None = None
    if args.chains:
        chains = [c.strip() for c in args.chains.split(",") if c.strip()]
        if not chains:
            sys.exit("ERROR: --chains given but no valid chain IDs parsed.")

    # Collect input files
    if args.structure_dir:
        src_dir = Path(args.structure_dir)
        if not src_dir.is_dir():
            sys.exit(f"ERROR: Not a directory: {src_dir}")
        files = sorted(
            p for p in src_dir.iterdir()
            if p.suffix.lower() in SUPPORTED_EXTENSIONS
        )
        if not files:
            sys.exit(
                f"ERROR: No supported structure files found in {src_dir}. "
                f"Supported extensions: {sorted(SUPPORTED_EXTENSIONS)}. "
                "Run scripts/download_structures.py first."
            )
    else:
        files = [Path(f) for f in args.structure_files]
        missing = [f for f in files if not f.exists()]
        if missing:
            sys.exit(
                f"ERROR: File(s) not found: {', '.join(str(f) for f in missing)}"
            )

    output_dir = Path(args.output_dir)

    print(
        f"Building ground-truth JSON for {len(files)} structure(s) "
        f"→ {output_dir}/"
    )
    if chains:
        print(f"  Chain filter: {chains}")

    errors: list[tuple[str, str]] = []

    for path in files:
        pdb_id = path.stem.upper()
        print(f"  {path.name} … ", end="", flush=True)
        try:
            out = build_ground_truth_json(path, output_dir, chains, args.overwrite)
            # Read back n_residues for the summary line
            with open(out) as fh:
                info = json.load(fh)
            print(
                f"OK — {info['n_residues']} residues, "
                f"chains {info['chains_included']} → {out.name}"
            )
        except Exception as exc:
            print("ERROR", file=sys.stderr)
            print(f"    {exc}", file=sys.stderr)
            errors.append((pdb_id, str(exc)))

    n_ok = len(files) - len(errors)
    print(f"\nDone: {n_ok}/{len(files)} ground-truth files written.")

    if errors:
        print("Failures:", file=sys.stderr)
        for pdb_id, msg in errors:
            print(f"  {pdb_id}: {msg}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
