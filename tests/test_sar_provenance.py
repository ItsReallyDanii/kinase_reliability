"""
Tests for SAR provenance and synthetic-vs-real status labelling.

Validates that:
- SAR generation correctly labels metrics_status in provenance
- Stub outputs are always flagged metrics_status="synthetic"
- Random values are never emitted on the structure_only path
- Required SAR fields are always present
- The generate_sar function raises on bad metrics_mode, not silently wrong outputs
"""

import json
from unittest.mock import patch

import numpy as np
import pytest

from generate_sar import (
    SARValidationError,
    classify_confidence,
    determine_expected_error_range,
    generate_sar,
    validate_sar,
    _make_stub_rng,
    _stub_contact_map_overlap,
    _stub_rmsd_ligand,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_fake_prediction(n: int = 100, plddt_mean: float = 75.0,
                          pae_mean: float = 8.0, stub: bool = True) -> dict:
    """Create a minimal prediction dict for testing."""
    rng = np.random.RandomState(0)
    coords = rng.randn(n, 3).tolist()
    plddt = (np.ones(n) * plddt_mean).tolist()
    pae = (np.ones((n, n)) * pae_mean).tolist()
    return {
        "pdb_id": "FAKE",
        "coordinates": coords,
        "plddt": plddt,
        "pae": pae,
        "model_version": "af3_stub" if stub else "af3_prod",
        "seed": 42,
        "recycles": 3,
        "stub_output": stub,
    }


def _make_fake_ground_truth(n: int = 100) -> dict:
    rng = np.random.RandomState(1)
    return {
        "pdb_id": "FAKE",
        "coordinates": rng.randn(n, 3).tolist(),
    }


def _make_target(pdb_id: str = "FAKE", ligand_present: bool = False) -> dict:
    return {"pdb_id": pdb_id, "ligand_present": ligand_present}


# ---------------------------------------------------------------------------
# Provenance: stub path always labelled "synthetic"
# ---------------------------------------------------------------------------


def test_stub_sar_has_synthetic_metrics_status():
    pred = _make_fake_prediction(stub=True)
    gt = _make_fake_ground_truth()
    sar = generate_sar(_make_target(), pred, gt, metrics_mode="synthetic", seed=42)
    assert sar["provenance"]["metrics_status"] == "synthetic"


def test_stub_sar_has_stub_output_true():
    pred = _make_fake_prediction(stub=True)
    gt = _make_fake_ground_truth()
    sar = generate_sar(_make_target(), pred, gt, metrics_mode="synthetic", seed=42)
    assert sar["provenance"]["stub_output"] is True


# ---------------------------------------------------------------------------
# Provenance: structure_only path labelled correctly
# ---------------------------------------------------------------------------


def test_structure_only_sar_has_correct_status():
    pred = _make_fake_prediction(stub=False)
    gt = _make_fake_ground_truth()
    sar = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=42)
    assert sar["provenance"]["metrics_status"] == "structure_only"


def test_structure_only_requires_ground_truth():
    pred = _make_fake_prediction(stub=False)
    with pytest.raises(ValueError, match="ground truth"):
        generate_sar(_make_target(), pred, None, metrics_mode="structure_only", seed=42)


# ---------------------------------------------------------------------------
# Real metrics: no random values on structure_only path
# ---------------------------------------------------------------------------


def test_structure_only_rmsd_is_deterministic():
    pred = _make_fake_prediction(n=50, stub=False)
    gt = _make_fake_ground_truth(n=50)

    sar1 = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=42)
    sar2 = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=42)

    assert sar1["metrics"]["rmsd_global"] == sar2["metrics"]["rmsd_global"]
    assert sar1["metrics"]["contact_map_overlap"] == sar2["metrics"]["contact_map_overlap"]


def test_structure_only_rmsd_not_random_across_seeds():
    """RMSD on structure_only path must not change with different seeds (real metric)."""
    pred = _make_fake_prediction(n=50, stub=False)
    gt = _make_fake_ground_truth(n=50)

    sar_s42 = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=42)
    sar_s99 = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=99)

    assert sar_s42["metrics"]["rmsd_global"] == sar_s99["metrics"]["rmsd_global"], (
        "RMSD on structure_only path must be independent of seed (no randomness)"
    )


# ---------------------------------------------------------------------------
# Inference-backed raises NotImplementedError
# ---------------------------------------------------------------------------


def test_inference_backed_raises():
    pred = _make_fake_prediction(stub=False)
    gt = _make_fake_ground_truth()
    with pytest.raises(NotImplementedError, match="AF3"):
        generate_sar(_make_target(), pred, gt, metrics_mode="inference_backed", seed=42)


# ---------------------------------------------------------------------------
# Unknown metrics_mode raises ValueError
# ---------------------------------------------------------------------------


def test_bad_metrics_mode_raises():
    pred = _make_fake_prediction()
    gt = _make_fake_ground_truth()
    with pytest.raises(ValueError, match="metrics_mode"):
        generate_sar(_make_target(), pred, gt, metrics_mode="bogus_mode", seed=42)


# ---------------------------------------------------------------------------
# Required SAR fields always present
# ---------------------------------------------------------------------------


def test_required_fields_present_synthetic():
    pred = _make_fake_prediction(stub=True)
    gt = _make_fake_ground_truth()
    sar = generate_sar(_make_target(), pred, gt, metrics_mode="synthetic", seed=42)

    assert "expected_error_range" in sar
    assert "recommended_action" in sar
    assert "decision_gate" in sar

    err_range = sar["expected_error_range"]
    assert "rmsd_min" in err_range
    assert "rmsd_max" in err_range
    assert "rationale" in err_range

    assert sar["decision_gate"] in {"ACCEPT", "REVIEW", "REJECT"}
    assert len(sar["recommended_action"]) > 0


def test_required_fields_present_structure_only():
    pred = _make_fake_prediction(stub=False)
    gt = _make_fake_ground_truth()
    sar = generate_sar(_make_target(), pred, gt, metrics_mode="structure_only", seed=42)

    validate_sar(sar, strict_mode=True)  # should not raise


def test_validate_sar_raises_on_missing_field():
    broken_sar = {
        "pdb_id": "FAKE",
        "sar_version": "1.0",
        "recommended_action": "do something",
        "decision_gate": "ACCEPT",
        # expected_error_range deliberately missing
    }
    with pytest.raises(SARValidationError, match="expected_error_range"):
        validate_sar(broken_sar, strict_mode=True)


def test_validate_sar_raises_on_invalid_decision_gate():
    broken_sar = {
        "pdb_id": "FAKE",
        "sar_version": "1.0",
        "recommended_action": "do something",
        "decision_gate": "MAYBE",  # invalid
        "expected_error_range": {"rmsd_min": 0.5, "rmsd_max": 2.0, "rationale": "test"},
    }
    with pytest.raises(SARValidationError):
        validate_sar(broken_sar, strict_mode=True)


# ---------------------------------------------------------------------------
# Stub RNG is seeded per-target (determinism)
# ---------------------------------------------------------------------------


def test_stub_rng_is_deterministic():
    rng1 = _make_stub_rng("8ABC", 42)
    rng2 = _make_stub_rng("8ABC", 42)
    assert rng1.uniform() == rng2.uniform()


def test_stub_rng_differs_across_targets():
    rng_abc = _make_stub_rng("8ABC", 42)
    rng_def = _make_stub_rng("8DEF", 42)
    val_abc = rng_abc.uniform()
    val_def = rng_def.uniform()
    # Different targets get different RNG streams
    assert val_abc != val_def, "Stub RNG should differ per target ID"


# ---------------------------------------------------------------------------
# Confidence classification sanity
# ---------------------------------------------------------------------------


def test_high_confidence_classification():
    conf = classify_confidence(plddt_mean=92.0, pae_mean=3.0)
    assert conf["plddt_bin"] == "high"
    assert conf["pae_bin"] == "high"
    assert conf["overall_confidence"] == "high"


def test_low_confidence_classification():
    conf = classify_confidence(plddt_mean=50.0, pae_mean=15.0)
    assert conf["plddt_bin"] == "low"
    assert conf["pae_bin"] == "low"
    assert conf["overall_confidence"] == "low"


def test_expected_error_range_values():
    high_range = determine_expected_error_range({"overall_confidence": "high"})
    low_range = determine_expected_error_range({"overall_confidence": "low"})
    assert high_range["rmsd_max"] < low_range["rmsd_max"], (
        "High confidence should have tighter error bound than low confidence"
    )
