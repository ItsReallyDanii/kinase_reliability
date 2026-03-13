#!/usr/bin/env python3
"""
Extract ligand (HETATM) information from mmCIF or PDB structure files.

Requires Biopython >= 1.81.

Writes a JSON inventory listing each non-solvent, non-buffer small molecule
found in each structure.  Useful for confirming ligand presence before
running the SAR pipeline.

Usage examples
--------------
# Extract from all .cif files in a directory:
python3 scripts/extract_ligands.py \\
    --structure_dir benchmark/ground_truth \\
    --output_file benchmark/metadata/ligand_inventory.json

# Extract from specific files:
python3 scripts/extract_ligands.py \\
    --structure_files benchmark/ground_truth/1ATP.cif benchmark/ground_truth/2ITO.cif \\
    --output_file benchmark/metadata/ligand_inventory.json

# Include solvent and buffer molecules:
python3 scripts/extract_ligands.py \\
    --structure_dir benchmark/ground_truth \\
    --output_file benchmark/metadata/ligand_inventory.json \\
    --include_solvent
"""

import argparse
import json
import sys
from pathlib import Path

# Common solvent / buffer residue codes excluded by default
SOLVENT_CODES = frozenset(
    [
        "HOH",  # water
        "WAT",  # water (alternate)
        "DOD",  # deuterium water
        "H2O",  # water
        "EDO",  # ethylene glycol
        "PEG",  # polyethylene glycol
        "GOL",  # glycerol
        "SO4",  # sulphate
        "PO4",  # phosphate
        "ACT",  # acetate
        "MES",  # MES buffer
        "HEPES",
        "TRIS",
        "MPD",  # 2-methyl-2,4-pentanediol
        "BME",  # beta-mercaptoethanol
        "DTT",  # dithiothreitol
        "FMT",  # formate
        "ACE",  # acetyl group (N-terminal cap)
        "NH2",  # amino group (C-terminal cap)
        "UNX",  # unknown atom
        "UNL",  # unknown ligand
    ]
)

# Metal ions to report separately
METAL_CODES = frozenset(
    ["MG", "ZN", "CA", "FE", "MN", "CU", "NI", "CO", "NA", "K", "CL"]
)


def _parse_structure(path: Path):
    """Return a Biopython Structure, auto-detecting format from extension."""
    try:
        from Bio.PDB import MMCIFParser, PDBParser
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for extract_ligands. "
            "Install with: pip install biopython>=1.81"
        ) from exc

    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    elif suffix in (".pdb", ".ent"):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            f"Unrecognised file extension '{suffix}' for {path}. "
            "Expected .cif, .mmcif, .pdb, or .ent."
        )

    return parser.get_structure(path.stem, str(path))


def extract_ligands(
    structure,
    pdb_id: str,
    include_solvent: bool = False,
) -> list[dict]:
    """Return a list of ligand records from a Biopython Structure.

    Each record contains:
        pdb_id, chain_id, residue_name, residue_seq_id,
        n_atoms, is_metal, excluded_as_solvent

    Args:
        structure:       Biopython Structure object.
        pdb_id:          PDB accession code (used for labelling).
        include_solvent: If True, include water/buffer molecules.

    Returns:
        List of dicts, one per HETATM residue found.
    """
    records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                het_flag, _seq_id, _ins_code = residue.get_id()
                if het_flag.strip() == "":
                    continue  # standard amino acid or nucleotide

                res_name = residue.get_resname().strip()
                is_solvent = res_name in SOLVENT_CODES
                is_metal = res_name in METAL_CODES

                if is_solvent and not include_solvent:
                    continue

                records.append(
                    {
                        "pdb_id": pdb_id,
                        "chain_id": chain.get_id(),
                        "residue_name": res_name,
                        "residue_seq_id": str(residue.get_id()[1]),
                        "n_atoms": len(list(residue.get_atoms())),
                        "is_metal": is_metal,
                        "excluded_as_solvent": is_solvent,
                    }
                )
        break  # first model only

    return records


def process_structure_file(
    path: Path, include_solvent: bool = False
) -> tuple[str, list[dict]]:
    """Parse a single structure file and extract ligand records.

    Returns:
        (pdb_id, list_of_records)

    Raises:
        ValueError:  on unrecognised file format.
        ImportError: if Biopython is not installed.
        Exception:   re-raised with filename context on parse failure.
    """
    pdb_id = path.stem.upper()
    try:
        structure = _parse_structure(path)
        records = extract_ligands(structure, pdb_id, include_solvent=include_solvent)
        return pdb_id, records
    except Exception as exc:
        raise RuntimeError(f"Failed to process {path}: {exc}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract ligand metadata from mmCIF/PDB structure files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--structure_dir",
        metavar="DIR",
        help="Directory containing .cif / .pdb files.",
    )
    source_group.add_argument(
        "--structure_files",
        nargs="+",
        metavar="FILE",
        help="One or more structure files to process.",
    )

    parser.add_argument(
        "--output_file",
        required=True,
        metavar="FILE",
        help="Output JSON file for the ligand inventory.",
    )
    parser.add_argument(
        "--include_solvent",
        action="store_true",
        help="Include water / buffer molecules in output (excluded by default).",
    )

    args = parser.parse_args()

    # Collect input files
    if args.structure_dir:
        structure_dir = Path(args.structure_dir)
        if not structure_dir.is_dir():
            sys.exit(f"ERROR: Not a directory: {structure_dir}")
        files = sorted(
            p for p in structure_dir.iterdir()
            if p.suffix.lower() in (".cif", ".mmcif", ".pdb", ".ent")
        )
        if not files:
            sys.exit(
                f"ERROR: No .cif / .pdb files found in {structure_dir}. "
                "Download structures first with scripts/download_structures.py"
            )
    else:
        files = [Path(f) for f in args.structure_files]
        missing = [f for f in files if not f.exists()]
        if missing:
            sys.exit(f"ERROR: File(s) not found: {', '.join(str(f) for f in missing)}")

    print(f"Processing {len(files)} structure file(s)…")

    inventory: dict[str, list] = {}
    errors: list[tuple[str, str]] = []

    for path in files:
        print(f"  {path.name} … ", end="", flush=True)
        try:
            pdb_id, records = process_structure_file(
                path, include_solvent=args.include_solvent
            )
            inventory[pdb_id] = records
            ligand_names = sorted({r["residue_name"] for r in records if not r["is_metal"] and not r["excluded_as_solvent"]})
            print(f"OK — {len(records)} HETATM record(s); ligands: {ligand_names or 'none'}")
        except Exception as exc:
            print(f"ERROR", file=sys.stderr)
            print(f"    {exc}", file=sys.stderr)
            errors.append((path.name, str(exc)))

    # Write output
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "description": "Ligand inventory extracted from structure files",
        "source": "scripts/extract_ligands.py",
        "include_solvent": args.include_solvent,
        "n_structures": len(inventory),
        "structures": inventory,
    }

    with open(output_path, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"\nInventory written to {output_path}")
    print(f"Processed {len(files)} file(s): {len(files) - len(errors)} OK, {len(errors)} failed.")

    if errors:
        print("Failures:", file=sys.stderr)
        for name, msg in errors:
            print(f"  {name}: {msg}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
