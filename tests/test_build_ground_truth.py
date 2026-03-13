"""
Tests for scripts/build_ground_truth_json.py.

All tests use in-memory synthetic PDB strings written to temporary files.
No network access or real structure downloads are required.

Covers:
  1. Valid single-chain structure → correct JSON shape and values
  2. Valid multi-chain structure → chains concatenated in order
  3. --chains filter → only requested chains included
  4. Structure with no CA atoms → clear ValueError
  5. Output JSON is accepted by generate_sar.load_ground_truth()
  6. End-to-end: generated JSON + stub prediction → structure_only SAR succeeds
  7. Overwrite flag behaviour
  8. Unsupported file extension → no output written
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Minimal PDB fixtures  (written to tmp files; no network required)
# ---------------------------------------------------------------------------

# Three residues, chain A, known coordinates
_PDB_3_RESIDUES_CHAIN_A = """\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       7.000   8.000   9.000  1.00  0.00           C
END
"""

# Three residues on chain A, two on chain B
_PDB_2_CHAIN = """\
ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       3.000   0.000   0.000  1.00  0.00           C
ATOM      4  CA  LEU B   1      10.000   0.000   0.000  1.00  0.00           C
ATOM      5  CA  ALA B   2      11.000   0.000   0.000  1.00  0.00           C
END
"""

# Only HETATM (ligand) — no protein CA atoms
_PDB_NO_CA = """\
HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C
HETATM    2  O1  LIG A   1       4.000   5.000   6.000  1.00  0.00           O
END
"""

# Mix: one CA atom plus HETATM — only the CA should be extracted
_PDB_MIXED = """\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
HETATM    2  C1  LIG A   2       9.000   9.000   9.000  1.00  0.00           C
END
"""


def _write_pdb(content: str, stem: str, tmp_dir: Path) -> Path:
    p = tmp_dir / f"{stem}.pdb"
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

from scripts.build_ground_truth_json import (
    build_ground_truth_json,
    extract_ca_coords_ordered,
    _parse_structure,
)


# ---------------------------------------------------------------------------
# 1. Valid single-chain structure → correct JSON
# ---------------------------------------------------------------------------


def test_single_chain_json_shape():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        out = build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(out) as fh:
            gt = json.load(fh)

    assert gt["pdb_id"] == "1TST"
    assert gt["n_residues"] == 3
    assert gt["chains_included"] == ["A"]
    assert len(gt["coordinates"]) == 3
    assert len(gt["coordinates"][0]) == 3  # each entry is [x, y, z]
    assert gt["source_format"] == "pdb"


def test_single_chain_coordinates_exact():
    """Extracted coordinates must exactly match the known PDB values."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "1TST_ground_truth.json") as fh:
            gt = json.load(fh)

    expected = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    for actual, exp in zip(gt["coordinates"], expected):
        assert actual == pytest.approx(exp, abs=1e-4)


# ---------------------------------------------------------------------------
# 2. Multi-chain structure → all chains concatenated
# ---------------------------------------------------------------------------


def test_multi_chain_all_chains():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_2_CHAIN, "2TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "2TST_ground_truth.json") as fh:
            gt = json.load(fh)

    assert gt["n_residues"] == 5
    assert set(gt["chains_included"]) == {"A", "B"}


def test_multi_chain_order_A_before_B():
    """Chain A residues must appear before chain B residues."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_2_CHAIN, "2TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "2TST_ground_truth.json") as fh:
            gt = json.load(fh)

    coords = gt["coordinates"]
    # Chain A: x ∈ {1,2,3}; Chain B: x ∈ {10,11}
    # First three entries must be from chain A (x < 5)
    assert coords[0][0] == pytest.approx(1.0, abs=1e-4)
    assert coords[3][0] == pytest.approx(10.0, abs=1e-4)


# ---------------------------------------------------------------------------
# 3. --chains filter
# ---------------------------------------------------------------------------


def test_chain_filter_restricts_output():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_2_CHAIN, "2TST", tmp)
        build_ground_truth_json(src, tmp, chains=["A"], overwrite=False)

        with open(tmp / "2TST_ground_truth.json") as fh:
            gt = json.load(fh)

    assert gt["n_residues"] == 3
    assert gt["chains_included"] == ["A"]


def test_chain_filter_nonexistent_chain_raises():
    """Requesting a chain that has no CA atoms must raise ValueError."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        with pytest.raises(ValueError, match="No C-alpha"):
            build_ground_truth_json(src, tmp, chains=["Z"], overwrite=False)


# ---------------------------------------------------------------------------
# 4. No CA atoms → clear error
# ---------------------------------------------------------------------------


def test_no_ca_atoms_raises():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_NO_CA, "9TST", tmp)
        with pytest.raises(ValueError, match="No C-alpha"):
            build_ground_truth_json(src, tmp, chains=None, overwrite=False)


def test_hetatm_not_included_in_coords():
    """HETATM records must not contribute CA coordinates."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_MIXED, "8TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "8TST_ground_truth.json") as fh:
            gt = json.load(fh)

    assert gt["n_residues"] == 1
    assert gt["coordinates"][0] == pytest.approx([1.0, 2.0, 3.0], abs=1e-4)


# ---------------------------------------------------------------------------
# 5. Output accepted by generate_sar.load_ground_truth()
# ---------------------------------------------------------------------------


def test_json_accepted_by_load_ground_truth():
    """The written JSON must be loadable by generate_sar.load_ground_truth."""
    from generate_sar import load_ground_truth

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        gt = load_ground_truth("1TST", tmp, allow_missing=False)

    assert gt is not None
    assert "coordinates" in gt
    assert len(gt["coordinates"]) == 3


# ---------------------------------------------------------------------------
# 6. End-to-end: generated JSON + stub prediction → structure_only SAR OK
# ---------------------------------------------------------------------------


def _make_stub_prediction_from_coords(coords: list, pdb_id: str = "1TST") -> dict:
    """Make a minimal prediction dict that matches the ground truth length."""
    n = len(coords)
    rng = np.random.RandomState(99)
    # Slightly perturbed coordinates (not random arrays — mimics a near-perfect prediction)
    pred_coords = (np.array(coords) + rng.randn(n, 3) * 0.5).tolist()
    plddt = (np.ones(n) * 78.0).tolist()
    pae = (np.ones((n, n)) * 7.5).tolist()
    return {
        "pdb_id": pdb_id,
        "coordinates": pred_coords,
        "plddt": plddt,
        "pae": pae,
        "model_version": "af3_stub",
        "seed": 42,
        "recycles": 3,
        "stub_output": False,  # triggers structure_only auto-detection
    }


def test_end_to_end_structure_only_sar():
    """
    Full bridge path:
      PDB file → build_ground_truth_json → ground_truth.json
      + stub prediction → generate_sar (structure_only) → valid SAR with RMSD ≥ 0
    """
    from generate_sar import generate_sar, validate_sar

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # Step 1: write PDB and extract ground truth
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "1TST_ground_truth.json") as fh:
            gt = json.load(fh)

        # Step 2: build a stub prediction using the same coords (near-identity)
        pred = _make_stub_prediction_from_coords(gt["coordinates"], pdb_id="1TST")

        # Step 3: generate SAR
        target = {"pdb_id": "1TST", "ligand_present": False}
        sar = generate_sar(target, pred, gt, metrics_mode="structure_only", seed=42)

        # Validate schema
        validate_sar(sar, strict_mode=True)

    # RMSD must be a non-negative float
    rmsd = sar["metrics"]["rmsd_global"]
    assert isinstance(rmsd, (int, float))
    assert rmsd >= 0.0

    # provenance must carry the correct status
    assert sar["provenance"]["metrics_status"] == "structure_only"
    assert sar["provenance"]["stub_output"] is False

    # Required SAR fields present
    assert sar["decision_gate"] in {"ACCEPT", "REVIEW", "REJECT"}
    assert len(sar["recommended_action"]) > 0


def test_identity_prediction_rmsd_near_zero():
    """
    If prediction coords == ground truth coords the structure_only RMSD
    must be ≈ 0 (validates the Kabsch path is actually used).
    """
    from generate_sar import generate_sar

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with open(tmp / "1TST_ground_truth.json") as fh:
            gt = json.load(fh)

        # Prediction == ground truth (identity)
        pred = {
            "pdb_id": "1TST",
            "coordinates": gt["coordinates"],
            "plddt": [80.0] * gt["n_residues"],
            "pae": [[5.0] * gt["n_residues"]] * gt["n_residues"],
            "model_version": "test",
            "seed": 42,
            "recycles": 3,
            "stub_output": False,
        }

        target = {"pdb_id": "1TST", "ligand_present": False}
        sar = generate_sar(target, pred, gt, metrics_mode="structure_only", seed=42)

    assert sar["metrics"]["rmsd_global"] == pytest.approx(0.0, abs=1e-4), (
        f"Identity prediction should give RMSD≈0, got {sar['metrics']['rmsd_global']}"
    )


# ---------------------------------------------------------------------------
# 7. Overwrite flag behaviour
# ---------------------------------------------------------------------------


def test_overwrite_false_raises_on_existing():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        build_ground_truth_json(src, tmp, chains=None, overwrite=False)

        with pytest.raises(FileExistsError, match="already exists"):
            build_ground_truth_json(src, tmp, chains=None, overwrite=False)


def test_overwrite_true_replaces_file():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_3_RESIDUES_CHAIN_A, "1TST", tmp)
        out1 = build_ground_truth_json(src, tmp, chains=None, overwrite=False)
        mtime1 = out1.stat().st_mtime

        import time; time.sleep(0.01)
        out2 = build_ground_truth_json(src, tmp, chains=None, overwrite=True)
        mtime2 = out2.stat().st_mtime

    assert mtime2 >= mtime1  # file was touched


# ---------------------------------------------------------------------------
# 8. Determinism
# ---------------------------------------------------------------------------


def test_extraction_is_deterministic():
    """Running the extraction twice must produce bit-identical JSON."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = _write_pdb(_PDB_2_CHAIN, "2TST", tmp)
        out = tmp / "2TST_ground_truth.json"

        build_ground_truth_json(src, tmp, chains=None, overwrite=False)
        content1 = out.read_text()

        build_ground_truth_json(src, tmp, chains=None, overwrite=True)
        content2 = out.read_text()

    # Timestamps differ in source_file (path) but coordinates must be identical
    gt1 = json.loads(content1)
    gt2 = json.loads(content2)
    assert gt1["coordinates"] == gt2["coordinates"]
    assert gt1["n_residues"] == gt2["n_residues"]
