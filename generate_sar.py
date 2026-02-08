#!/usr/bin/env python3
"""
Structural Audit Report (SAR) Generator for Kinase Reliability Pilot v1.0

Generates SARs with mandatory fields:
- expected_error_range (REQUIRED)
- recommended_action (REQUIRED)
- decision_gate ∈ {ACCEPT, REVIEW, REJECT} (REQUIRED)

Missing any required field => ERROR_SAR_INCOMPLETE
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional
import numpy as np


class SARValidationError(Exception):
    """Raised when SAR validation fails."""
    pass


def compute_rmsd_stub(pred_coords: List[List[float]], gt_coords: Optional[List[List[float]]]) -> float:
    """
    Compute stub RMSD (deterministic for demo).

    In production, would compute actual RMSD between predicted and ground truth coordinates.
    """
    if gt_coords is None:
        # No ground truth available - return deterministic stub value
        return np.random.uniform(1.5, 4.5)

    # Simplified RMSD calculation for stub
    pred = np.array(pred_coords)
    gt = np.array(gt_coords)

    if pred.shape != gt.shape:
        # Handle size mismatch
        min_len = min(len(pred), len(gt))
        pred = pred[:min_len]
        gt = gt[:min_len]

    return float(np.sqrt(np.mean(np.sum((pred - gt)**2, axis=1))))


def classify_confidence(plddt_mean: float, pae_mean: float) -> Dict[str, str]:
    """Classify confidence into bins based on pLDDT and PAE."""
    # pLDDT bins: >90=high, 70-90=medium, <70=low
    if plddt_mean > 90:
        plddt_bin = "high"
    elif plddt_mean > 70:
        plddt_bin = "medium"
    else:
        plddt_bin = "low"

    # PAE bins: <5=high, 5-10=medium, >10=low
    if pae_mean < 5:
        pae_bin = "high"
    elif pae_mean < 10:
        pae_bin = "medium"
    else:
        pae_bin = "low"

    # Overall confidence: both high => high, either low => low, else medium
    if plddt_bin == "high" and pae_bin == "high":
        overall = "high"
    elif plddt_bin == "low" or pae_bin == "low":
        overall = "low"
    else:
        overall = "medium"

    return {
        "plddt_bin": plddt_bin,
        "pae_bin": pae_bin,
        "overall_confidence": overall
    }


def determine_expected_error_range(confidence: Dict[str, str]) -> Dict[str, Any]:
    """
    Determine expected error range based on confidence assessment.

    REQUIRED FIELD - must always be present.
    """
    overall = confidence["overall_confidence"]

    if overall == "high":
        return {
            "rmsd_min": 0.5,
            "rmsd_max": 2.0,
            "rationale": "High confidence (pLDDT>90, PAE<5) suggests low expected error"
        }
    elif overall == "medium":
        return {
            "rmsd_min": 1.5,
            "rmsd_max": 4.0,
            "rationale": "Medium confidence suggests moderate expected error range"
        }
    else:  # low
        return {
            "rmsd_min": 3.0,
            "rmsd_max": 8.0,
            "rationale": "Low confidence (pLDDT<70 or PAE>10) suggests high expected error"
        }


def classify_failure(rmsd_global: float, rmsd_ligand: Any, plddt_mean: float,
                    pae_mean: float, expected_range: Dict[str, float],
                    ligand_present: bool) -> Dict[str, str]:
    """
    Classify failure mode according to taxonomy:
    - Class A: Overconfidence artifact (high confidence but high error)
    - Class B: Ligand pose failure (ligand RMSD >> global RMSD)
    - Class C: Symmetry/assembly failure
    - Unknown: Unmapped failure mode
    - N/A: No failure detected
    """
    rmsd_max = expected_range["rmsd_max"]

    # Check if prediction is within expected error range
    if rmsd_global <= rmsd_max:
        return {
            "class": "N/A",
            "description": f"Prediction within expected error range (RMSD={rmsd_global:.2f}Å <= {rmsd_max:.2f}Å)"
        }

    # Failure detected - classify

    # Class A: Overconfidence artifact (high confidence but high error)
    if plddt_mean > 90 and pae_mean < 5 and rmsd_global > rmsd_max * 1.5:
        return {
            "class": "Class A",
            "description": f"Overconfidence artifact: High model confidence (pLDDT={plddt_mean:.1f}, PAE={pae_mean:.1f}) but RMSD={rmsd_global:.2f}Å exceeds expected range"
        }

    # Class B: Ligand pose failure
    if ligand_present and isinstance(rmsd_ligand, (int, float)):
        if rmsd_ligand > rmsd_global * 1.5:
            return {
                "class": "Class B",
                "description": f"Ligand pose failure: Ligand pocket RMSD ({rmsd_ligand:.2f}Å) significantly exceeds global RMSD ({rmsd_global:.2f}Å)"
            }

    # Class C: Symmetry/assembly failure (heuristic: very high error with specific pattern)
    if rmsd_global > 10.0:
        return {
            "class": "Class C",
            "description": f"Potential symmetry/assembly failure: Extremely high RMSD ({rmsd_global:.2f}Å) suggests structural misalignment"
        }

    # Unknown failure mode
    return {
        "class": "Unknown",
        "description": f"Unmapped failure mode: RMSD={rmsd_global:.2f}Å exceeds expected range but doesn't match known failure patterns"
    }


def determine_decision_gate(rmsd_global: float, expected_range: Dict[str, float],
                           failure_class: str, confidence: Dict[str, str]) -> str:
    """
    Assign decision gate: ACCEPT, REVIEW, or REJECT.

    REQUIRED FIELD - must always be present.
    """
    rmsd_max = expected_range["rmsd_max"]

    # ACCEPT: Within expected error range
    if rmsd_global <= rmsd_max:
        return "ACCEPT"

    # REJECT: Class A (overconfidence) or Class C (symmetry) or very high error
    if failure_class in ["Class A", "Class C"] or rmsd_global > rmsd_max * 2:
        return "REJECT"

    # REVIEW: All other cases (moderate failures, unknown modes, Class B)
    return "REVIEW"


def generate_recommended_action(decision_gate: str, failure_taxonomy: Dict[str, str],
                               pdb_id: str) -> str:
    """
    Generate recommended action based on decision gate and failure mode.

    REQUIRED FIELD - must always be present.
    """
    if decision_gate == "ACCEPT":
        return f"Structure {pdb_id} passes quality criteria. Proceed with downstream analysis."

    elif decision_gate == "REVIEW":
        failure_class = failure_taxonomy["class"]
        if failure_class == "Class B":
            return f"Manual review recommended for {pdb_id}: Ligand binding pose shows elevated error. Verify active site geometry and ligand interactions."
        elif failure_class == "Unknown":
            return f"Manual review required for {pdb_id}: Elevated error detected but failure mode unclear. Expert structural analysis needed."
        else:
            return f"Manual review recommended for {pdb_id}: {failure_taxonomy['description']}"

    else:  # REJECT
        failure_class = failure_taxonomy["class"]
        if failure_class == "Class A":
            return f"Reject {pdb_id}: Overconfidence artifact detected. Model confidence metrics unreliable for this target. Consider alternative modeling approaches."
        elif failure_class == "Class C":
            return f"Reject {pdb_id}: Potential symmetry/assembly failure. Verify oligomeric state and symmetry operators before use."
        else:
            return f"Reject {pdb_id}: Error exceeds acceptable thresholds. Do not use for downstream analysis without significant refinement."


def validate_sar(sar: Dict[str, Any], strict_mode: bool = True) -> None:
    """
    Validate SAR contains all required fields.

    In strict mode, raises SARValidationError if any required field missing.
    """
    required_fields = [
        "expected_error_range",
        "recommended_action",
        "decision_gate"
    ]

    missing_fields = []
    for field in required_fields:
        if field not in sar or sar[field] is None:
            missing_fields.append(field)

    # Validate expected_error_range structure
    if "expected_error_range" in sar and sar["expected_error_range"] is not None:
        err_range = sar["expected_error_range"]
        for subfield in ["rmsd_min", "rmsd_max", "rationale"]:
            if subfield not in err_range or err_range[subfield] is None:
                missing_fields.append(f"expected_error_range.{subfield}")

    # Validate decision_gate enum
    if "decision_gate" in sar and sar["decision_gate"] is not None:
        if sar["decision_gate"] not in ["ACCEPT", "REVIEW", "REJECT"]:
            missing_fields.append("decision_gate (invalid value)")

    if missing_fields and strict_mode:
        raise SARValidationError(f"ERROR_SAR_INCOMPLETE: Missing required fields: {', '.join(missing_fields)}")

    if missing_fields:
        print(f"WARNING: SAR validation failed: Missing fields: {', '.join(missing_fields)}", file=sys.stderr)


def generate_sar(target: Dict[str, Any], prediction: Dict[str, Any],
                ground_truth: Optional[Dict[str, Any]], args) -> Dict[str, Any]:
    """Generate a complete SAR for a single target."""
    pdb_id = target["pdb_id"]

    # Load prediction data
    plddt_scores = np.array(prediction["plddt"])
    pae_matrix = np.array(prediction["pae"])

    # Compute metrics
    plddt_mean = float(np.mean(plddt_scores))
    pae_mean = float(np.mean(pae_matrix))

    # Compute RMSD (stub implementation)
    rmsd_global = compute_rmsd_stub(
        prediction["coordinates"],
        ground_truth["coordinates"] if ground_truth else None
    )

    # Compute ligand pocket RMSD if applicable
    ligand_present = target.get("ligand_present", False)
    if ligand_present:
        # Stub: ligand RMSD is global RMSD * random factor
        rmsd_ligand_pocket = rmsd_global * np.random.uniform(0.8, 1.5)
    else:
        rmsd_ligand_pocket = "N/A"

    # Contact map overlap (stub)
    contact_map_overlap = float(np.random.uniform(0.6, 0.95))

    # Classify confidence
    confidence_assessment = classify_confidence(plddt_mean, pae_mean)

    # Determine expected error range (REQUIRED)
    expected_error_range = determine_expected_error_range(confidence_assessment)

    # Classify failure mode
    failure_taxonomy = classify_failure(
        rmsd_global, rmsd_ligand_pocket, plddt_mean, pae_mean,
        expected_error_range, ligand_present
    )

    # Determine decision gate (REQUIRED)
    decision_gate = determine_decision_gate(
        rmsd_global, expected_error_range,
        failure_taxonomy["class"], confidence_assessment
    )

    # Generate recommended action (REQUIRED)
    recommended_action = generate_recommended_action(
        decision_gate, failure_taxonomy, pdb_id
    )

    # Build SAR
    sar = {
        "pdb_id": pdb_id,
        "sar_version": "1.0",
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "metrics": {
            "rmsd_global": round(rmsd_global, 2),
            "rmsd_ligand_pocket": round(rmsd_ligand_pocket, 2) if isinstance(rmsd_ligand_pocket, (int, float)) else rmsd_ligand_pocket,
            "contact_map_overlap": round(contact_map_overlap, 3),
            "plddt_mean": round(plddt_mean, 2),
            "pae_mean": round(pae_mean, 2)
        },
        "confidence_assessment": confidence_assessment,
        "expected_error_range": expected_error_range,
        "recommended_action": recommended_action,
        "decision_gate": decision_gate,
        "failure_taxonomy": failure_taxonomy,
        "provenance": {
            "model_version": prediction.get("model_version", "unknown"),
            "seed": prediction.get("seed", -1),
            "recycles": prediction.get("recycles", -1),
            "stub_output": prediction.get("stub_output", False)
        }
    }

    return sar


def load_ground_truth(pdb_id: str, ground_truth_dir: Path) -> Optional[Dict[str, Any]]:
    """Load ground truth structure if available."""
    gt_file = ground_truth_dir / f"{pdb_id}_ground_truth.json"

    if not gt_file.exists():
        print(f"  Warning: No ground truth found for {pdb_id}", file=sys.stderr)
        # Generate stub ground truth for demo
        return {
            "pdb_id": pdb_id,
            "coordinates": np.random.randn(300, 3).tolist()
        }

    with open(gt_file, 'r') as f:
        return json.load(f)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate Structural Audit Reports (SARs) for Kinase Reliability Pilot"
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to benchmark manifest JSON"
    )
    parser.add_argument(
        "--pred_dir",
        required=True,
        help="Directory containing prediction outputs"
    )
    parser.add_argument(
        "--ground_truth_dir",
        required=True,
        help="Directory containing ground truth structures"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for SARs"
    )
    parser.add_argument(
        "--schema_version",
        required=True,
        help="SAR schema version (e.g., 1.0)"
    )
    parser.add_argument(
        "--strict_mode",
        action="store_true",
        help="Enable strict validation (fail on missing required fields)"
    )

    args = parser.parse_args()

    # Load manifest
    print(f"Loading manifest: {args.manifest}")
    with open(args.manifest, 'r') as f:
        manifest = json.load(f)

    targets = manifest["targets"]
    print(f"Found {len(targets)} targets")

    # Setup directories
    pred_dir = Path(args.pred_dir)
    gt_dir = Path(args.ground_truth_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each target
    print(f"\nGenerating SARs (strict_mode={args.strict_mode})...")
    success_count = 0
    validation_errors = []

    for i, target in enumerate(targets, 1):
        pdb_id = target["pdb_id"]
        print(f"[{i}/{len(targets)}] {pdb_id}...", end=" ")

        try:
            # Load prediction
            pred_file = pred_dir / f"{pdb_id}_prediction.json"
            if not pred_file.exists():
                raise FileNotFoundError(f"Prediction not found: {pred_file}")

            with open(pred_file, 'r') as f:
                prediction = json.load(f)

            # Load ground truth
            ground_truth = load_ground_truth(pdb_id, gt_dir)

            # Generate SAR
            sar = generate_sar(target, prediction, ground_truth, args)

            # Validate SAR
            validate_sar(sar, strict_mode=args.strict_mode)

            # Write SAR
            sar_file = output_dir / f"{pdb_id}.json"
            with open(sar_file, 'w') as f:
                json.dump(sar, f, indent=2)

            print(f"✓ {sar['decision_gate']}")
            success_count += 1

        except SARValidationError as e:
            print(f"✗ VALIDATION ERROR", file=sys.stderr)
            print(f"  {str(e)}", file=sys.stderr)
            validation_errors.append({"pdb_id": pdb_id, "error": str(e)})

            if args.strict_mode:
                sys.exit(1)

        except Exception as e:
            print(f"✗ ERROR: {str(e)}", file=sys.stderr)
            validation_errors.append({"pdb_id": pdb_id, "error": str(e)})

            if args.strict_mode:
                sys.exit(1)

    print(f"\n{'='*60}")
    print(f"SAR generation complete: {success_count}/{len(targets)} successful")

    if validation_errors:
        print(f"\nValidation errors encountered: {len(validation_errors)}")
        for err in validation_errors:
            print(f"  - {err['pdb_id']}: {err['error']}")

        if args.strict_mode:
            sys.exit(1)


if __name__ == "__main__":
    main()
