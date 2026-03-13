#!/usr/bin/env python3
"""
Structural Audit Report (SAR) Generator for Kinase Reliability Pilot v1.0

Generates SARs with mandatory fields:
- expected_error_range (REQUIRED)
- recommended_action (REQUIRED)
- decision_gate ∈ {ACCEPT, REVIEW, REJECT} (REQUIRED)

Missing any required field => ERROR_SAR_INCOMPLETE

Metrics mode (--metrics_mode):
  synthetic          : Deterministic stub metrics derived from prediction arrays.
                       No real structural comparison.  All outputs flagged
                       metrics_status="synthetic".  This is the default when
                       prediction["stub_output"] == True.
  structure_only     : Real Kabsch RMSD and real contact map computed from
                       coordinates in prediction JSON vs ground truth JSON.
                       Useful for engineering validation before real AF3 runs.
                       Flags metrics_status="structure_only".
  inference_backed   : Reserved for when real AF3 predictions are available.
                       Currently raises NotImplementedError.

IMPORTANT: Do NOT use synthetic outputs as scientific results.  Check
provenance["stub_output"] and provenance["metrics_status"] before
interpreting any SAR.
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


# ---------------------------------------------------------------------------
# STUB / SYNTHETIC metric helpers
# (only called when metrics_mode == "synthetic")
# All randomness is isolated here and behind an explicit guard.
# ---------------------------------------------------------------------------


def _stub_rmsd_global(pred_coords: List[List[float]],
                      gt_coords: Optional[List[List[float]]]) -> float:
    """Deterministic stub RMSD for synthetic pipeline runs.

    WARNING: not a real structural metric.  Returns a value derived from
    input array sizes only.  Used exclusively when metrics_mode="synthetic".
    """
    if gt_coords is None:
        # Cannot compute even a stub without ground truth; return sentinel
        return float(np.linalg.norm(np.array(pred_coords[:5]))) if pred_coords else 15.0

    pred = np.array(pred_coords)
    gt = np.array(gt_coords)
    min_len = min(len(pred), len(gt))
    if min_len == 0:
        return 0.0
    pred = pred[:min_len]
    gt = gt[:min_len]
    # Simple unaligned RMSD on truncated arrays — not Kabsch-aligned
    return float(np.sqrt(np.mean(np.sum((pred - gt) ** 2, axis=1))))


def _stub_rmsd_ligand(rmsd_global: float, rng: np.random.RandomState) -> float:
    """Synthetic ligand RMSD.  Uses seeded RNG, not np.random global."""
    return float(rmsd_global * rng.uniform(0.8, 1.5))


def _stub_contact_map_overlap(rng: np.random.RandomState) -> float:
    """Synthetic contact map overlap.  Uses seeded RNG."""
    return float(rng.uniform(0.6, 0.95))


def _make_stub_rng(pdb_id: str, seed: int) -> np.random.RandomState:
    """Seeded per-target RNG so synthetic values are deterministic."""
    pdb_hash = sum(ord(c) for c in pdb_id)
    return np.random.RandomState(seed + pdb_hash)


# ---------------------------------------------------------------------------
# REAL metric helpers
# (called when metrics_mode in {"structure_only", "inference_backed"})
# ---------------------------------------------------------------------------


def _real_rmsd_global(pred_coords: List[List[float]],
                      gt_coords: List[List[float]]) -> float:
    """Real Kabsch-aligned global RMSD using metrics.rmsd."""
    try:
        from metrics.rmsd import compute_rmsd, truncate_to_common_length
    except ImportError as exc:
        raise ImportError(
            "metrics.rmsd not importable. Ensure you are running from the "
            "repository root and numpy is installed."
        ) from exc

    pred = np.array(pred_coords, dtype=float)
    gt = np.array(gt_coords, dtype=float)
    if len(pred) != len(gt):
        pred, gt = truncate_to_common_length(pred, gt)
    return compute_rmsd(pred, gt, align=True)


def _real_contact_map_overlap(pred_coords: List[List[float]],
                               gt_coords: List[List[float]],
                               threshold: float = 8.0) -> float:
    """Real contact map overlap (Jaccard) using metrics.contact_map."""
    try:
        from metrics.contact_map import compute_contact_map, compute_contact_map_overlap
    except ImportError as exc:
        raise ImportError(
            "metrics.contact_map not importable. Ensure you are running from "
            "the repository root and numpy is installed."
        ) from exc

    n = min(len(pred_coords), len(gt_coords))
    pred = np.array(pred_coords[:n], dtype=float)
    gt = np.array(gt_coords[:n], dtype=float)
    map_pred = compute_contact_map(pred, threshold=threshold)
    map_gt = compute_contact_map(gt, threshold=threshold)
    return compute_contact_map_overlap(map_pred, map_gt)


# ---------------------------------------------------------------------------
# Confidence classification
# ---------------------------------------------------------------------------


def classify_confidence(plddt_mean: float, pae_mean: float) -> Dict[str, str]:
    """Classify confidence into bins based on pLDDT and PAE."""
    if plddt_mean > 90:
        plddt_bin = "high"
    elif plddt_mean > 70:
        plddt_bin = "medium"
    else:
        plddt_bin = "low"

    if pae_mean < 5:
        pae_bin = "high"
    elif pae_mean < 10:
        pae_bin = "medium"
    else:
        pae_bin = "low"

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


# ---------------------------------------------------------------------------
# Expected error range (REQUIRED SAR field)
# ---------------------------------------------------------------------------


def determine_expected_error_range(confidence: Dict[str, str]) -> Dict[str, Any]:
    """Determine expected error range based on confidence assessment.

    REQUIRED FIELD — must always be present.
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
    else:
        return {
            "rmsd_min": 3.0,
            "rmsd_max": 8.0,
            "rationale": "Low confidence (pLDDT<70 or PAE>10) suggests high expected error"
        }


# ---------------------------------------------------------------------------
# Failure taxonomy
# ---------------------------------------------------------------------------


def classify_failure(rmsd_global: float, rmsd_ligand: Any, plddt_mean: float,
                     pae_mean: float, expected_range: Dict[str, float],
                     ligand_present: bool) -> Dict[str, str]:
    """Classify failure mode.

    - Class A: Overconfidence artifact (high confidence but high error)
    - Class B: Ligand pose failure (ligand RMSD >> global RMSD)
    - Class C: Symmetry/assembly failure
    - Unknown: Unmapped failure mode
    - N/A: No failure detected
    """
    rmsd_max = expected_range["rmsd_max"]

    if rmsd_global <= rmsd_max:
        return {
            "class": "N/A",
            "description": (
                f"Prediction within expected error range "
                f"(RMSD={rmsd_global:.2f}Å <= {rmsd_max:.2f}Å)"
            )
        }

    if plddt_mean > 90 and pae_mean < 5 and rmsd_global > rmsd_max * 1.5:
        return {
            "class": "Class A",
            "description": (
                f"Overconfidence artifact: High model confidence "
                f"(pLDDT={plddt_mean:.1f}, PAE={pae_mean:.1f}) but "
                f"RMSD={rmsd_global:.2f}Å exceeds expected range"
            )
        }

    if ligand_present and isinstance(rmsd_ligand, (int, float)):
        if rmsd_ligand > rmsd_global * 1.5:
            return {
                "class": "Class B",
                "description": (
                    f"Ligand pose failure: Ligand pocket RMSD ({rmsd_ligand:.2f}Å) "
                    f"significantly exceeds global RMSD ({rmsd_global:.2f}Å)"
                )
            }

    if rmsd_global > 10.0:
        return {
            "class": "Class C",
            "description": (
                f"Potential symmetry/assembly failure: Extremely high RMSD "
                f"({rmsd_global:.2f}Å) suggests structural misalignment"
            )
        }

    return {
        "class": "Unknown",
        "description": (
            f"Unmapped failure mode: RMSD={rmsd_global:.2f}Å exceeds expected "
            f"range but doesn't match known failure patterns"
        )
    }


# ---------------------------------------------------------------------------
# Decision gate (REQUIRED SAR field)
# ---------------------------------------------------------------------------


def determine_decision_gate(rmsd_global: float, expected_range: Dict[str, float],
                             failure_class: str, confidence: Dict[str, str]) -> str:
    """Assign ACCEPT, REVIEW, or REJECT.

    REQUIRED FIELD — must always be present.
    """
    rmsd_max = expected_range["rmsd_max"]

    if rmsd_global <= rmsd_max:
        return "ACCEPT"

    if failure_class in ["Class A", "Class C"] or rmsd_global > rmsd_max * 2:
        return "REJECT"

    return "REVIEW"


# ---------------------------------------------------------------------------
# Recommended action (REQUIRED SAR field)
# ---------------------------------------------------------------------------


def generate_recommended_action(decision_gate: str, failure_taxonomy: Dict[str, str],
                                 pdb_id: str) -> str:
    """Generate recommended action based on decision gate and failure mode.

    REQUIRED FIELD — must always be present.
    """
    if decision_gate == "ACCEPT":
        return (
            f"Structure {pdb_id} passes quality criteria. "
            f"Proceed with downstream analysis."
        )

    elif decision_gate == "REVIEW":
        failure_class = failure_taxonomy["class"]
        if failure_class == "Class B":
            return (
                f"Manual review recommended for {pdb_id}: Ligand binding pose "
                f"shows elevated error. Verify active site geometry and ligand interactions."
            )
        elif failure_class == "Unknown":
            return (
                f"Manual review required for {pdb_id}: Elevated error detected "
                f"but failure mode unclear. Expert structural analysis needed."
            )
        else:
            return (
                f"Manual review recommended for {pdb_id}: "
                f"{failure_taxonomy['description']}"
            )

    else:  # REJECT
        failure_class = failure_taxonomy["class"]
        if failure_class == "Class A":
            return (
                f"Reject {pdb_id}: Overconfidence artifact detected. Model confidence "
                f"metrics unreliable for this target. Consider alternative modeling approaches."
            )
        elif failure_class == "Class C":
            return (
                f"Reject {pdb_id}: Potential symmetry/assembly failure. Verify "
                f"oligomeric state and symmetry operators before use."
            )
        else:
            return (
                f"Reject {pdb_id}: Error exceeds acceptable thresholds. Do not use "
                f"for downstream analysis without significant refinement."
            )


# ---------------------------------------------------------------------------
# SAR validation
# ---------------------------------------------------------------------------


def validate_sar(sar: Dict[str, Any], strict_mode: bool = True) -> None:
    """Validate SAR contains all required fields.

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

    if "expected_error_range" in sar and sar["expected_error_range"] is not None:
        err_range = sar["expected_error_range"]
        for subfield in ["rmsd_min", "rmsd_max", "rationale"]:
            if subfield not in err_range or err_range[subfield] is None:
                missing_fields.append(f"expected_error_range.{subfield}")

    if "decision_gate" in sar and sar["decision_gate"] is not None:
        if sar["decision_gate"] not in ["ACCEPT", "REVIEW", "REJECT"]:
            missing_fields.append("decision_gate (invalid value)")

    if missing_fields and strict_mode:
        raise SARValidationError(
            f"ERROR_SAR_INCOMPLETE: Missing required fields: {', '.join(missing_fields)}"
        )

    if missing_fields:
        print(
            f"WARNING: SAR validation failed: Missing fields: {', '.join(missing_fields)}",
            file=sys.stderr
        )


# ---------------------------------------------------------------------------
# Ground truth loading
# ---------------------------------------------------------------------------


def load_ground_truth(pdb_id: str, ground_truth_dir: Path,
                      allow_missing: bool = True) -> Optional[Dict[str, Any]]:
    """Load ground truth structure JSON if available.

    Args:
        pdb_id:          Target ID.
        ground_truth_dir: Directory to search.
        allow_missing:   If False, raise FileNotFoundError when not found.

    Returns:
        Parsed JSON dict or None if missing (and allow_missing=True).

    Note: Ground truth JSON files must contain a "coordinates" key with
    an (N, 3) array of C-alpha positions.  They are produced by running
    a real structure through the ingestion pipeline, NOT by random generation.
    """
    gt_file = ground_truth_dir / f"{pdb_id}_ground_truth.json"

    if not gt_file.exists():
        msg = f"No ground truth found for {pdb_id} at {gt_file}"
        if not allow_missing:
            raise FileNotFoundError(msg)
        print(f"  Warning: {msg}", file=sys.stderr)
        return None

    with open(gt_file, 'r') as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Core SAR generation
# ---------------------------------------------------------------------------


def generate_sar(target: Dict[str, Any], prediction: Dict[str, Any],
                 ground_truth: Optional[Dict[str, Any]],
                 metrics_mode: str, seed: int) -> Dict[str, Any]:
    """Generate a complete SAR for a single target.

    Args:
        target:       Target entry from manifest.
        prediction:   Prediction JSON (from run_inference.py output).
        ground_truth: Ground truth JSON (coordinates from real structure) or None.
        metrics_mode: "synthetic" | "structure_only" | "inference_backed"
        seed:         Seed used for inference (for seeded stub RNG).

    Returns:
        SAR dict with all required fields populated.
    """
    pdb_id = target["pdb_id"]

    # --- Confidence from prediction arrays ---
    plddt_scores = np.array(prediction["plddt"])
    pae_matrix = np.array(prediction["pae"])
    plddt_mean = float(np.mean(plddt_scores))
    pae_mean = float(np.mean(pae_matrix))

    confidence_assessment = classify_confidence(plddt_mean, pae_mean)
    expected_error_range = determine_expected_error_range(confidence_assessment)

    ligand_present = target.get("ligand_present", False)

    # --- Metrics: strict separation between stub and real ---
    if metrics_mode == "synthetic":
        # All synthetic values — clearly labelled in provenance
        rng = _make_stub_rng(pdb_id, seed)
        rmsd_global = _stub_rmsd_global(
            prediction["coordinates"],
            ground_truth["coordinates"] if ground_truth else None
        )
        rmsd_ligand_pocket: Any = (
            _stub_rmsd_ligand(rmsd_global, rng) if ligand_present else "N/A"
        )
        contact_map_overlap = _stub_contact_map_overlap(rng)
        metrics_status = "synthetic"

    elif metrics_mode == "structure_only":
        # Real metrics from coordinate arrays (no AF3 inference yet)
        if ground_truth is None:
            raise ValueError(
                f"metrics_mode='structure_only' requires ground truth for {pdb_id}. "
                "Provide --ground_truth_dir with real structure JSON files."
            )
        rmsd_global = _real_rmsd_global(
            prediction["coordinates"],
            ground_truth["coordinates"]
        )
        # Ligand pocket RMSD: not yet implemented for structure_only mode
        rmsd_ligand_pocket = "N/A"
        contact_map_overlap = _real_contact_map_overlap(
            prediction["coordinates"],
            ground_truth["coordinates"]
        )
        metrics_status = "structure_only"

    elif metrics_mode == "inference_backed":
        raise NotImplementedError(
            "metrics_mode='inference_backed' is reserved for real AF3 predictions. "
            "AF3 integration is not yet implemented. "
            "See docs/CURRENT_STATUS.md for Phase 2 prerequisites."
        )
    else:
        raise ValueError(
            f"Unknown metrics_mode '{metrics_mode}'. "
            "Choose: synthetic | structure_only | inference_backed"
        )

    # --- Failure taxonomy and decision gate ---
    failure_taxonomy = classify_failure(
        rmsd_global, rmsd_ligand_pocket, plddt_mean, pae_mean,
        expected_error_range, ligand_present
    )

    decision_gate = determine_decision_gate(
        rmsd_global, expected_error_range,
        failure_taxonomy["class"], confidence_assessment
    )

    recommended_action = generate_recommended_action(
        decision_gate, failure_taxonomy, pdb_id
    )

    # --- Build SAR ---
    sar = {
        "pdb_id": pdb_id,
        "sar_version": "1.0",
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "metrics": {
            "rmsd_global": round(rmsd_global, 2),
            "rmsd_ligand_pocket": (
                round(rmsd_ligand_pocket, 2)
                if isinstance(rmsd_ligand_pocket, (int, float))
                else rmsd_ligand_pocket
            ),
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
            "stub_output": prediction.get("stub_output", False),
            "metrics_status": metrics_status,
        }
    }

    return sar


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate Structural Audit Reports (SARs) for Kinase Reliability Pilot"
    )
    parser.add_argument("--manifest", required=True,
                        help="Path to benchmark manifest JSON")
    parser.add_argument("--pred_dir", required=True,
                        help="Directory containing prediction outputs")
    parser.add_argument("--ground_truth_dir", required=True,
                        help="Directory containing ground truth structure JSONs")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for SARs")
    parser.add_argument("--schema_version", required=True,
                        help="SAR schema version (e.g., 1.0)")
    parser.add_argument("--strict_mode", action="store_true",
                        help="Enable strict validation (fail on missing required fields)")
    parser.add_argument(
        "--metrics_mode",
        choices=["synthetic", "structure_only", "inference_backed"],
        default=None,
        help=(
            "Metric computation mode.  Default: auto-detect from prediction "
            "provenance (stub_output=True => synthetic, else structure_only)."
        ),
    )

    args = parser.parse_args()

    with open(args.manifest, 'r') as f:
        manifest = json.load(f)

    targets = manifest["targets"]
    print(f"Loading manifest: {args.manifest} ({len(targets)} targets)")

    pred_dir = Path(args.pred_dir)
    gt_dir = Path(args.ground_truth_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating SARs (strict_mode={args.strict_mode})…")
    success_count = 0
    validation_errors = []

    for i, target in enumerate(targets, 1):
        pdb_id = target["pdb_id"]
        print(f"[{i}/{len(targets)}] {pdb_id}…", end=" ")

        try:
            pred_file = pred_dir / f"{pdb_id}_prediction.json"
            if not pred_file.exists():
                raise FileNotFoundError(f"Prediction not found: {pred_file}")

            with open(pred_file, 'r') as f:
                prediction = json.load(f)

            # Auto-detect metrics mode unless explicitly specified
            if args.metrics_mode is not None:
                metrics_mode = args.metrics_mode
            elif prediction.get("stub_output", False):
                metrics_mode = "synthetic"
            else:
                metrics_mode = "structure_only"

            # Load ground truth (required for structure_only mode)
            allow_missing = (metrics_mode == "synthetic")
            ground_truth = load_ground_truth(pdb_id, gt_dir, allow_missing=allow_missing)

            seed = prediction.get("seed", 42)
            sar = generate_sar(target, prediction, ground_truth, metrics_mode, seed)
            validate_sar(sar, strict_mode=args.strict_mode)

            sar_file = output_dir / f"{pdb_id}.json"
            with open(sar_file, 'w') as f:
                json.dump(sar, f, indent=2)

            print(f"OK [{sar['decision_gate']}] metrics={metrics_mode}")
            success_count += 1

        except SARValidationError as e:
            print(f"VALIDATION ERROR", file=sys.stderr)
            print(f"  {str(e)}", file=sys.stderr)
            validation_errors.append({"pdb_id": pdb_id, "error": str(e)})
            if args.strict_mode:
                sys.exit(1)

        except Exception as e:
            print(f"ERROR: {str(e)}", file=sys.stderr)
            validation_errors.append({"pdb_id": pdb_id, "error": str(e)})
            if args.strict_mode:
                sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"SAR generation complete: {success_count}/{len(targets)} successful")

    if validation_errors:
        print(f"\nErrors: {len(validation_errors)}")
        for err in validation_errors:
            print(f"  - {err['pdb_id']}: {err['error']}")
        if args.strict_mode:
            sys.exit(1)


if __name__ == "__main__":
    main()
