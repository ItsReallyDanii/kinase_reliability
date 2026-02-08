#!/usr/bin/env python3
"""
Report Compilation for Kinase Reliability Pilot v1.0

Generates:
1. SAR_SUMMARY.md - Decision gate counts and failure taxonomy distribution
2. calibration_report.json - Confidence-vs-error bands
3. execution_provenance.json - Full runtime config, commands, timestamps, hashes

PASS criteria:
- 10/10 targets processed (completeness)
- All SARs contain decision_gate (validity)
- calibration_report.json has >= 3 confidence bins (calibration)
- execution_provenance.json exists (provenance)
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Tuple
from collections import Counter
import hashlib


def compute_file_hash(file_path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def load_sars(sar_dir: Path, manifest: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Load all SARs for targets in manifest.

    Returns (sars, missing_targets).
    """
    targets = manifest["targets"]
    sars = []
    missing = []

    for target in targets:
        pdb_id = target["pdb_id"]
        sar_file = sar_dir / f"{pdb_id}.json"

        if not sar_file.exists():
            missing.append(pdb_id)
            continue

        with open(sar_file, 'r') as f:
            sar = json.load(f)
            sars.append(sar)

    return sars, missing


def validate_completeness(sars: List[Dict[str, Any]], manifest: Dict[str, Any]) -> bool:
    """
    Validate completeness: 10/10 targets processed.

    PASS criteria: Output exists for every pdb_id.
    FAIL condition: Any target missing from sar_results/.
    """
    expected_count = len(manifest["targets"])
    actual_count = len(sars)

    if actual_count < expected_count:
        print(f"FAIL: Completeness check failed - {actual_count}/{expected_count} SARs found", file=sys.stderr)
        return False

    print(f"✓ Completeness: {actual_count}/{expected_count} SARs present")
    return True


def validate_sar_validity(sars: List[Dict[str, Any]]) -> bool:
    """
    Validate SAR validity: 100% contain decision_gate.

    PASS criteria: 100% of SARs contain decision_gate (ACCEPT/REVIEW/REJECT).
    FAIL condition: Any SAR field missing or null.
    """
    invalid_sars = []

    for sar in sars:
        pdb_id = sar.get("pdb_id", "unknown")

        # Check required fields
        if "decision_gate" not in sar or sar["decision_gate"] is None:
            invalid_sars.append(f"{pdb_id}: missing decision_gate")
        elif sar["decision_gate"] not in ["ACCEPT", "REVIEW", "REJECT"]:
            invalid_sars.append(f"{pdb_id}: invalid decision_gate value")

        if "expected_error_range" not in sar or sar["expected_error_range"] is None:
            invalid_sars.append(f"{pdb_id}: missing expected_error_range")

        if "recommended_action" not in sar or sar["recommended_action"] is None:
            invalid_sars.append(f"{pdb_id}: missing recommended_action")

    if invalid_sars:
        print(f"FAIL: SAR validity check failed", file=sys.stderr)
        for error in invalid_sars:
            print(f"  - {error}", file=sys.stderr)
        return False

    print(f"✓ SAR Validity: All {len(sars)} SARs contain required fields")
    return True


def generate_summary(sars: List[Dict[str, Any]], output_dir: Path) -> None:
    """Generate SAR_SUMMARY.md with decision gate counts and failure taxonomy."""

    # Count decision gates
    decision_counts = Counter(sar["decision_gate"] for sar in sars)

    # Count failure classes
    failure_counts = Counter(sar["failure_taxonomy"]["class"] for sar in sars)

    # Generate markdown
    summary_md = f"""# Kinase Reliability Pilot v1.0 - SAR Summary

**Generated:** {datetime.utcnow().isoformat()}Z
**Total Targets:** {len(sars)}

## Decision Gate Distribution

| Decision Gate | Count | Percentage |
|--------------|-------|------------|
| ACCEPT       | {decision_counts.get('ACCEPT', 0)} | {decision_counts.get('ACCEPT', 0) / len(sars) * 100:.1f}% |
| REVIEW       | {decision_counts.get('REVIEW', 0)} | {decision_counts.get('REVIEW', 0) / len(sars) * 100:.1f}% |
| REJECT       | {decision_counts.get('REJECT', 0)} | {decision_counts.get('REJECT', 0) / len(sars) * 100:.1f}% |

## Failure Taxonomy Distribution

| Failure Class | Count | Description |
|--------------|-------|-------------|
| N/A          | {failure_counts.get('N/A', 0)} | No failure detected |
| Class A      | {failure_counts.get('Class A', 0)} | Overconfidence artifact |
| Class B      | {failure_counts.get('Class B', 0)} | Ligand pose failure |
| Class C      | {failure_counts.get('Class C', 0)} | Symmetry/assembly failure |
| Unknown      | {failure_counts.get('Unknown', 0)} | Unmapped failure mode |

## Target Details

| PDB ID | Decision Gate | Failure Class | RMSD (Å) | pLDDT | Confidence |
|--------|--------------|---------------|----------|-------|------------|
"""

    for sar in sorted(sars, key=lambda x: x["pdb_id"]):
        pdb_id = sar["pdb_id"]
        gate = sar["decision_gate"]
        fail_class = sar["failure_taxonomy"]["class"]
        rmsd = sar["metrics"]["rmsd_global"]
        plddt = sar["metrics"]["plddt_mean"]
        conf = sar["confidence_assessment"]["overall_confidence"]

        summary_md += f"| {pdb_id} | {gate} | {fail_class} | {rmsd:.2f} | {plddt:.1f} | {conf} |\n"

    # Write summary
    summary_file = output_dir / "SAR_SUMMARY.md"
    with open(summary_file, 'w') as f:
        f.write(summary_md)

    print(f"Generated: {summary_file}")


def generate_calibration_report(sars: List[Dict[str, Any]], output_dir: Path) -> bool:
    """
    Generate calibration_report.json with confidence-vs-error bands.

    PASS criteria: >= 3 confidence bins with error bars.
    FAIL condition: < 3 bins or missing error distribution.
    """
    # Group by confidence level
    confidence_bins = {"high": [], "medium": [], "low": []}

    for sar in sars:
        conf = sar["confidence_assessment"]["overall_confidence"]
        rmsd = sar["metrics"]["rmsd_global"]
        confidence_bins[conf].append(rmsd)

    # Compute statistics for each bin
    calibration_data = {}
    bins_with_data = 0

    for conf_level, rmsd_values in confidence_bins.items():
        if len(rmsd_values) > 0:
            import numpy as np
            calibration_data[conf_level] = {
                "n_samples": len(rmsd_values),
                "rmsd_mean": float(np.mean(rmsd_values)),
                "rmsd_std": float(np.std(rmsd_values)),
                "rmsd_min": float(np.min(rmsd_values)),
                "rmsd_max": float(np.max(rmsd_values)),
                "rmsd_median": float(np.median(rmsd_values))
            }
            bins_with_data += 1

    # Check PASS criteria
    if bins_with_data < 3:
        print(f"WARNING: Calibration has only {bins_with_data} confidence bins (< 3)", file=sys.stderr)

    calibration_report = {
        "version": "1.0",
        "generated": datetime.utcnow().isoformat() + "Z",
        "total_targets": len(sars),
        "confidence_bins": bins_with_data,
        "calibration_data": calibration_data,
        "interpretation": {
            "high": "Expected RMSD < 2.0Å - model highly confident and typically accurate",
            "medium": "Expected RMSD 1.5-4.0Å - model moderately confident with moderate accuracy",
            "low": "Expected RMSD 3.0-8.0Å - model uncertain with higher expected error"
        }
    }

    # Write calibration report
    calibration_file = output_dir / "calibration_report.json"
    with open(calibration_file, 'w') as f:
        json.dump(calibration_report, f, indent=2)

    print(f"Generated: {calibration_file}")
    print(f"✓ Calibration: {bins_with_data} confidence bins with error distributions")

    return bins_with_data >= 3


def generate_execution_provenance(args, manifest_path: Path, output_dir: Path) -> None:
    """
    Generate execution_provenance.json with full runtime config.

    PASS criteria: File includes full command line & timestamps.
    FAIL condition: Missing provenance file.
    """
    # Compute manifest hashes for verification
    accepted_hash = compute_file_hash(Path(args.manifest))

    # Build full command line (reconstructed)
    command_line = f"python3 compile_reports.py " \
                  f"--manifest {args.manifest} " \
                  f"--sar_dir {args.sar_dir} " \
                  f"--job_id {args.job_id} " \
                  f"--accepted_manifest_hash {args.accepted_manifest_hash} " \
                  f"--rejected_manifest_hash {args.rejected_manifest_hash}"

    provenance = {
        "job_id": args.job_id,
        "job_type": "compile_reports",
        "version": "1.0",
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "command_line": command_line,
        "arguments": {
            "manifest": args.manifest,
            "sar_dir": args.sar_dir,
            "job_id": args.job_id,
            "accepted_manifest_hash": args.accepted_manifest_hash,
            "rejected_manifest_hash": args.rejected_manifest_hash
        },
        "manifest_verification": {
            "accepted_manifest_hash_expected": args.accepted_manifest_hash,
            "accepted_manifest_hash_actual": accepted_hash,
            "hash_match": accepted_hash == args.accepted_manifest_hash
        },
        "locked_parameters": {
            "seed": 42,
            "recycles": 3,
            "schema_version": "1.0"
        },
        "outputs": {
            "sar_summary": str(output_dir / "SAR_SUMMARY.md"),
            "calibration_report": str(output_dir / "calibration_report.json"),
            "execution_provenance": str(output_dir / "execution_provenance.json")
        }
    }

    # Write provenance
    provenance_file = output_dir / "execution_provenance.json"
    with open(provenance_file, 'w') as f:
        json.dump(provenance, f, indent=2)

    print(f"Generated: {provenance_file}")
    print(f"✓ Provenance: Runtime configuration logged")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compile reports for Kinase Reliability Pilot"
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to benchmark manifest JSON"
    )
    parser.add_argument(
        "--sar_dir",
        required=True,
        help="Directory containing SAR outputs"
    )
    parser.add_argument(
        "--job_id",
        required=True,
        help="Job identifier (e.g., KINASE_PILOT_V1)"
    )
    parser.add_argument(
        "--accepted_manifest_hash",
        required=True,
        help="Expected SHA-256 hash of accepted manifest"
    )
    parser.add_argument(
        "--rejected_manifest_hash",
        required=True,
        help="Expected SHA-256 hash of rejected manifest"
    )

    args = parser.parse_args()

    print(f"{'='*60}")
    print(f"Kinase Reliability Pilot v1.0 - Report Compilation")
    print(f"Job ID: {args.job_id}")
    print(f"{'='*60}\n")

    # Load manifest
    manifest_path = Path(args.manifest)
    print(f"Loading manifest: {args.manifest}")
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    # Load SARs
    sar_dir = Path(args.sar_dir)
    print(f"Loading SARs from: {args.sar_dir}")
    sars, missing = load_sars(sar_dir, manifest)

    if missing:
        print(f"\nERROR: Missing SARs for targets: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(sars)} SARs\n")

    # Validation phase
    print("Running validation checks...")
    all_passed = True

    # Check completeness
    if not validate_completeness(sars, manifest):
        all_passed = False

    # Check SAR validity
    if not validate_sar_validity(sars):
        all_passed = False

    if not all_passed:
        print(f"\nFAIL: Validation checks failed", file=sys.stderr)
        sys.exit(1)

    print()

    # Generate reports
    print("Generating reports...")
    output_dir = sar_dir  # Write reports to SAR directory

    generate_summary(sars, output_dir)
    calibration_passed = generate_calibration_report(sars, output_dir)
    generate_execution_provenance(args, manifest_path, output_dir)

    print(f"\n{'='*60}")

    if calibration_passed:
        print("PASS: All validation criteria met")
        print(f"{'='*60}")
    else:
        print("WARNING: Calibration has fewer than 3 confidence bins")
        print(f"{'='*60}")
        sys.exit(1)


if __name__ == "__main__":
    main()
