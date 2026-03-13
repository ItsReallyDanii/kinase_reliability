"""
Tests for benchmark metadata scaffold validation (Level A).

Validates that benchmark/pilot_set_v1.json:
- loads successfully as JSON
- contains required top-level fields
- has the expected schema for each target entry
- clearly flags provisional/uncertain slots
- does not contain synthetic IDs
"""

import json
from pathlib import Path

import pytest

METADATA_PATH = Path(__file__).parent.parent / "benchmark" / "pilot_set_v1.json"

REQUIRED_TARGET_FIELDS = {
    "pdb_id",
    "resolution",
    "release_year",
    "description",
}

@pytest.fixture(scope="module")
def metadata():
    assert METADATA_PATH.exists(), (
        f"Metadata file not found: {METADATA_PATH}. "
    )
    with open(METADATA_PATH) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Top-level structure
# ---------------------------------------------------------------------------


def test_metadata_loads(metadata):
    assert isinstance(metadata, dict)


def test_metadata_has_targets(metadata):
    assert "targets" in metadata, "metadata must have a 'targets' key"
    assert isinstance(metadata["targets"], list)
    assert len(metadata["targets"]) == 10


def test_metadata_has_version(metadata):
    assert "version" in metadata


def test_metadata_has_scope(metadata):
    assert "scope" in metadata


# ---------------------------------------------------------------------------
# Per-target fields
# ---------------------------------------------------------------------------


def test_each_target_has_required_fields(metadata):
    for target in metadata["targets"]:
        missing = REQUIRED_TARGET_FIELDS - set(target.keys())
        assert not missing, (
            f"Target entry missing required fields: {missing}. Entry: {target}"
        )


def test_8par_is_provisional(metadata):
    target_8par = next((t for t in metadata["targets"] if t["pdb_id"] == "8PAR"), None)
    assert target_8par is not None
    assert target_8par.get("provisional") is True
    assert "uncertainty_note" in target_8par


def test_no_synthetic_pdb_ids_in_metadata(metadata):
    """Synthetic fixture PDB IDs (8ABC style) must not appear in metadata."""
    synthetic_prefix_pattern = ["8ABC", "8DEF", "8GHI", "8JKL", "8MNO",
                                 "8PQR", "8STU", "8VWX", "8YZA", "8BCD"]
    for target in metadata["targets"]:
        pdb_id = target.get("pdb_id", "")
        assert pdb_id not in synthetic_prefix_pattern, (
            f"Synthetic PDB ID '{pdb_id}' must not appear in the real benchmark."
        )
