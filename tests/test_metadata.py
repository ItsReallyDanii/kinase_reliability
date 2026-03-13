"""
Tests for benchmark metadata scaffold validation.

Validates that benchmark/metadata/kinase_pilot_targets.json:
- loads successfully as JSON
- contains required top-level fields
- has the expected schema for each target entry
- clearly flags engineering test fixtures as non-benchmark
- does not claim any TBD entries are real targets
"""

import json
from pathlib import Path

import pytest

METADATA_PATH = Path(__file__).parent.parent / "benchmark" / "metadata" / "kinase_pilot_targets.json"

REQUIRED_TARGET_FIELDS = {
    "pdb_id",
    "kinase_name",
    "kinase_family",
    "experimental_method",
    "resolution_A",
    "ligand_name",
    "ligand_code",
    "release_date",
    "inclusion_status",
    "split",
    "notes",
    "engineering_test_fixture",
}

VALID_INCLUSION_STATUSES = {"included", "excluded", "pending_review"}
VALID_SPLITS = {"pilot_engineering_test", "pilot", "holdout", "excluded", "tbd"}
VALID_METHODS = {"X-RAY DIFFRACTION", "CRYO-EM", "NMR"}  # only first expected in pilot


@pytest.fixture(scope="module")
def metadata():
    assert METADATA_PATH.exists(), (
        f"Metadata file not found: {METADATA_PATH}. "
        "Create it with the benchmark metadata scaffold."
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
    assert len(metadata["targets"]) > 0


def test_metadata_has_schema_version(metadata):
    assert "schema_version" in metadata


def test_metadata_has_status_field(metadata):
    """Status field must be present to communicate draft state."""
    assert "status" in metadata, (
        "metadata must have a 'status' field clarifying whether targets are finalised"
    )


def test_metadata_has_inclusion_criteria(metadata):
    assert "inclusion_criteria" in metadata


# ---------------------------------------------------------------------------
# Per-target fields
# ---------------------------------------------------------------------------


def test_each_target_has_required_fields(metadata):
    for target in metadata["targets"]:
        # Skip comment-only entries
        if set(target.keys()) == {"_comment"}:
            continue
        missing = REQUIRED_TARGET_FIELDS - set(target.keys())
        assert not missing, (
            f"Target entry missing required fields: {missing}. Entry: {target}"
        )


def test_inclusion_status_valid(metadata):
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        status = target.get("inclusion_status")
        assert status in VALID_INCLUSION_STATUSES, (
            f"Invalid inclusion_status '{status}' for {target.get('pdb_id')}. "
            f"Expected one of {VALID_INCLUSION_STATUSES}"
        )


def test_split_valid(metadata):
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        split = target.get("split")
        assert split in VALID_SPLITS, (
            f"Invalid split '{split}' for {target.get('pdb_id')}. "
            f"Expected one of {VALID_SPLITS}"
        )


def test_engineering_fixture_flag_is_bool(metadata):
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        flag = target.get("engineering_test_fixture")
        assert isinstance(flag, bool), (
            f"engineering_test_fixture must be bool, got {type(flag)} for "
            f"{target.get('pdb_id')}"
        )


# ---------------------------------------------------------------------------
# Truth checks: engineering fixtures must be labelled correctly
# ---------------------------------------------------------------------------


def test_engineering_fixtures_are_excluded(metadata):
    """Engineering test fixtures must have inclusion_status='excluded'."""
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        if target.get("engineering_test_fixture"):
            assert target["inclusion_status"] == "excluded", (
                f"Engineering test fixture '{target['pdb_id']}' must have "
                f"inclusion_status='excluded', got '{target['inclusion_status']}'"
            )


def test_engineering_fixtures_not_marked_as_included(metadata):
    """No engineering fixture should be in the 'included' state."""
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        if target.get("engineering_test_fixture"):
            assert target["inclusion_status"] != "included", (
                f"Engineering fixture '{target['pdb_id']}' must NOT be "
                f"marked included — it is not a real benchmark target."
            )


def test_tbd_entries_have_pending_review_status(metadata):
    """TBD split entries must not be marked 'included'."""
    for target in metadata["targets"]:
        if "_comment" in target and "pdb_id" not in target:
            continue
        if target.get("split") == "tbd":
            assert target.get("inclusion_status") != "included", (
                f"TBD entry '{target.get('pdb_id')}' must not claim 'included' status "
                "until real target selection is complete."
            )


def test_no_synthetic_pdb_ids_in_metadata(metadata):
    """Synthetic fixture PDB IDs (8ABC style) must not appear in metadata."""
    synthetic_prefix_pattern = ["8ABC", "8DEF", "8GHI", "8JKL", "8MNO",
                                 "8PQR", "8STU", "8VWX", "8YZA", "8BCD"]
    for target in metadata["targets"]:
        pdb_id = target.get("pdb_id", "")
        assert pdb_id not in synthetic_prefix_pattern, (
            f"Synthetic PDB ID '{pdb_id}' from benchmark_v1.0.json must not appear "
            "in the real benchmark metadata file."
        )
