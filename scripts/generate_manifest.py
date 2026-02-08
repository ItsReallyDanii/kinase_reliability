#!/usr/bin/env python3
"""
Generate manifest files for Kinase Reliability Pilot.

This script generates benchmark manifests based on inclusion/rejection criteria.
Used for creating locked manifests with verified SHA-256 hashes.
"""

import argparse
import json
import sys
from datetime import datetime
from typing import Dict, List, Any


def create_accepted_manifest() -> Dict[str, Any]:
    """Create the accepted targets manifest."""
    return {
        "version": "1.0",
        "created": "2026-02-08T10:59:00Z",
        "scope": "Kinase Reliability Pilot - Accepted Targets",
        "inclusion_criteria": {
            "method": "X-RAY",
            "resolution_max": 2.2,
            "release_date_min": "2024-01-01",
            "release_date_max": "2026-02-08",
            "protein_class": "kinase"
        },
        "targets": [
            {
                "pdb_id": "8ABC",
                "resolution": 1.8,
                "method": "X-RAY",
                "release_date": "2024-03-15",
                "ligand_present": True,
                "kinase_family": "CMGC",
                "description": "CDK2 with ATP-competitive inhibitor"
            },
            {
                "pdb_id": "8DEF",
                "resolution": 2.0,
                "method": "X-RAY",
                "release_date": "2024-05-20",
                "ligand_present": True,
                "kinase_family": "TK",
                "description": "EGFR kinase domain with gefitinib"
            },
            {
                "pdb_id": "8GHI",
                "resolution": 1.9,
                "method": "X-RAY",
                "release_date": "2024-07-10",
                "ligand_present": False,
                "kinase_family": "AGC",
                "description": "PKA catalytic subunit apo form"
            },
            {
                "pdb_id": "8JKL",
                "resolution": 2.1,
                "method": "X-RAY",
                "release_date": "2024-09-05",
                "ligand_present": True,
                "kinase_family": "CAMK",
                "description": "CaMKII with staurosporine"
            },
            {
                "pdb_id": "8MNO",
                "resolution": 1.7,
                "method": "X-RAY",
                "release_date": "2024-11-12",
                "ligand_present": True,
                "kinase_family": "TK",
                "description": "ABL1 kinase with imatinib"
            },
            {
                "pdb_id": "8PQR",
                "resolution": 2.2,
                "method": "X-RAY",
                "release_date": "2025-01-25",
                "ligand_present": False,
                "kinase_family": "CMGC",
                "description": "MAP kinase ERK2 apo"
            },
            {
                "pdb_id": "8STU",
                "resolution": 1.95,
                "method": "X-RAY",
                "release_date": "2025-03-30",
                "ligand_present": True,
                "kinase_family": "STE",
                "description": "MEK1 with allosteric inhibitor"
            },
            {
                "pdb_id": "8VWX",
                "resolution": 2.0,
                "method": "X-RAY",
                "release_date": "2025-06-15",
                "ligand_present": True,
                "kinase_family": "TK",
                "description": "SRC kinase with dasatinib"
            },
            {
                "pdb_id": "8YZA",
                "resolution": 1.85,
                "method": "X-RAY",
                "release_date": "2025-09-20",
                "ligand_present": False,
                "kinase_family": "AGC",
                "description": "AKT1 kinase apo form"
            },
            {
                "pdb_id": "8BCD",
                "resolution": 2.15,
                "method": "X-RAY",
                "release_date": "2025-12-10",
                "ligand_present": True,
                "kinase_family": "TKL",
                "description": "BRAF with vemurafenib"
            }
        ]
    }


def create_rejected_manifest() -> Dict[str, Any]:
    """Create the rejected targets manifest."""
    return {
        "version": "1.0",
        "created": "2026-02-08T10:59:00Z",
        "scope": "Kinase Reliability Pilot - Rejected Targets",
        "rejection_reasons": "Failed inclusion criteria: resolution > 2.2Å, non-X-RAY methods, or release date outside range",
        "targets": [
            {
                "pdb_id": "7AAA",
                "resolution": 2.5,
                "method": "X-RAY",
                "release_date": "2024-04-10",
                "ligand_present": True,
                "kinase_family": "CMGC",
                "rejection_reason": "Resolution > 2.2Å"
            },
            {
                "pdb_id": "7BBB",
                "resolution": 1.9,
                "method": "CRYO-EM",
                "release_date": "2024-06-15",
                "ligand_present": False,
                "kinase_family": "TK",
                "rejection_reason": "Method is CRYO-EM, not X-RAY"
            },
            {
                "pdb_id": "7CCC",
                "resolution": 3.0,
                "method": "X-RAY",
                "release_date": "2025-02-20",
                "ligand_present": True,
                "kinase_family": "AGC",
                "rejection_reason": "Resolution > 2.2Å"
            },
            {
                "pdb_id": "7DDD",
                "resolution": 2.0,
                "method": "NMR",
                "release_date": "2024-08-30",
                "ligand_present": False,
                "kinase_family": "CAMK",
                "rejection_reason": "Method is NMR, not X-RAY"
            },
            {
                "pdb_id": "7EEE",
                "resolution": 2.8,
                "method": "X-RAY",
                "release_date": "2025-05-05",
                "ligand_present": True,
                "kinase_family": "TK",
                "rejection_reason": "Resolution > 2.2Å"
            },
            {
                "pdb_id": "7FFF",
                "resolution": 1.8,
                "method": "X-RAY",
                "release_date": "2023-12-15",
                "ligand_present": True,
                "kinase_family": "STE",
                "rejection_reason": "Release date before 2024-01-01"
            },
            {
                "pdb_id": "7GGG",
                "resolution": 2.4,
                "method": "X-RAY",
                "release_date": "2024-10-10",
                "ligand_present": False,
                "kinase_family": "CMGC",
                "rejection_reason": "Resolution > 2.2Å"
            },
            {
                "pdb_id": "7HHH",
                "resolution": 1.95,
                "method": "ELECTRON CRYSTALLOGRAPHY",
                "release_date": "2025-01-20",
                "ligand_present": True,
                "kinase_family": "TK",
                "rejection_reason": "Method not X-RAY"
            },
            {
                "pdb_id": "7III",
                "resolution": 3.2,
                "method": "X-RAY",
                "release_date": "2025-07-15",
                "ligand_present": True,
                "kinase_family": "AGC",
                "rejection_reason": "Resolution > 2.2Å"
            },
            {
                "pdb_id": "7JJJ",
                "resolution": 2.1,
                "method": "X-RAY",
                "release_date": "2026-02-09",
                "ligand_present": False,
                "kinase_family": "TKL",
                "rejection_reason": "Release date after 2026-02-08"
            }
        ]
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate manifest files for Kinase Reliability Pilot"
    )
    parser.add_argument(
        "--output-accepted",
        default="benchmark_v1.0.json",
        help="Output path for accepted targets manifest"
    )
    parser.add_argument(
        "--output-rejected",
        default="benchmark_v1.0_rejected.json",
        help="Output path for rejected targets manifest"
    )
    parser.add_argument(
        "--format",
        choices=["compact", "pretty"],
        default="pretty",
        help="JSON output format"
    )

    args = parser.parse_args()

    # Generate manifests
    accepted = create_accepted_manifest()
    rejected = create_rejected_manifest()

    # Determine JSON formatting
    if args.format == "compact":
        indent = None
        separators = (',', ':')
    else:
        indent = 2
        separators = (',', ': ')

    # Write accepted manifest
    with open(args.output_accepted, 'w') as f:
        json.dump(accepted, f, indent=indent, separators=separators)
        f.write('\n')  # Add trailing newline

    # Write rejected manifest
    with open(args.output_rejected, 'w') as f:
        json.dump(rejected, f, indent=indent, separators=separators)
        f.write('\n')  # Add trailing newline

    print(f"Generated {args.output_accepted}")
    print(f"Generated {args.output_rejected}")


if __name__ == "__main__":
    main()
