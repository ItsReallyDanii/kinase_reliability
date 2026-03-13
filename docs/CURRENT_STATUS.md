# Current Status — Kinase Reliability Pilot

**Document type:** Truth-status claims ledger
**Date:** 2026-03-13
**Phase:** Level A — Engineering Complete

This document is the authoritative reference for what is real, what is
scaffolded, and what remains for future phases. It should be updated whenever
the repository's capabilities change materially.

---

## Top-Level Status

| Statement | True? |
|---|---|
| AlphaFold 3 has been run on real kinase targets | **NO** |
| The benchmark contains real PDB IDs | **YES** — `benchmark/pilot_set_v1.json` contains a Level A provisional pilot set |
| RMSD values in `sar_results/` reflect real structural comparisons | **NO** — stub outputs only. `structure_only` runs generate valid proxy metrics. |
| The SAR schema and governance scaffold are real | **YES** |
| Real RMSD code (Kabsch/SVD) is implemented | **YES** (`metrics/rmsd.py`) |
| Real contact map code is implemented | **YES** (`metrics/contact_map.py`) |
| structure_only pipeline path runnable end-to-end | **YES** — download → build_ground_truth_json → generate_sar |
| Real structures have been downloaded | **YES** — for the provisional pilot set |
| Real kinase benchmark targets have been selected | **YES** — Level A Provisional Pilot Set (10 targets, `8PAR` provisional) |

---

## Detailed Claims Ledger

### ✅ Implemented (Real)

#### Schema and Governance
- `schemas/sar_schema_v1.json`: SAR JSON schema v1.0, locked
- Required-field enforcement: `expected_error_range`, `recommended_action`, `decision_gate`
- Strict-mode validation raises `ERROR_SAR_INCOMPLETE` on missing fields
- SHA-256 integrity verification for manifest files
- Provenance tracking: full command line, timestamps, locked parameters logged
- Failure taxonomy: Class A/B/C/Unknown/N/A classification logic
- Decision gate logic: ACCEPT/REVIEW/REJECT assignment

#### Real Metric Implementations (Phase 1)
- `metrics/rmsd.py`: Kabsch/SVD aligned RMSD — deterministic, no random values
  - `compute_rmsd(coords_a, coords_b, align=True)` — numpy-only path
  - `compute_ca_rmsd_biopython(struct_a, struct_b)` — Biopython path (requires Bio installed)
- `metrics/contact_map.py`: Distance-based contact maps — deterministic, no random values
  - `compute_contact_map(coords, threshold=8.0)` — binary C-alpha contact matrix
  - `compute_contact_map_overlap(map_a, map_b)` — Jaccard overlap

#### Structure Ingestion Scripts (Phase 1/Level A)
- `scripts/download_structures.py`: Downloads mmCIF/PDB files from RCSB REST API
  - Fetched the 10 pilot structures to `pdb_ground_truth/`
  - Fails loudly on invalid IDs (HTTP 404)
  - Validates PDB ID format before requesting
  - Writes to predictable output paths
- `scripts/extract_ligands.py`: Extracts HETATM ligand records from mmCIF files
  - Uses Biopython MMCIF parser
  - Outputs ligand inventory JSON
  - Excludes solvent (HOH, WAT) and common buffers by default
- `scripts/build_ground_truth_json.py`: Bridges downloaded structures → SAR pipeline
  - Parses .cif / .pdb with Biopython; first model only
  - Extracts C-alpha coordinates for all (or specified) chains, in chain-then-sequence order
  - Writes `{PDB_ID}_ground_truth.json` with `coordinates` key consumed by `generate_sar.py`
  - No randomness; fails loudly on missing/malformed CA atoms
  - Completes the `structure_only` end-to-end path

#### Benchmark Metadata Format (Level A)
- `benchmark/pilot_set_v1.json`: Format defined with fields:
  - `pdb_id`, `resolution`, `release_year`, `description`
  - Contains precisely 10 Level A provisional IDs
  - Explicit provisional support (`provisional`: `true`) preserved for 8PAR.

---

### 🔧 Scaffolded (Infrastructure Real, Outputs Synthetic)

#### Pipeline Orchestration
- `run_inference.py`: Pipeline wrapper is real; inference is stub only
  - `af3_available = False` — always falls through to `run_alphafold3_stub()`
  - Stub generates deterministic coordinates from `pdb_id` hash + seed
  - All stub outputs flagged with `stub_output: true` in provenance

- `generate_sar.py`: SAR generation logic is real; metrics currently stub by default
  - Stub path: active when `prediction["stub_output"] == True` (default)
  - Real metric path: active when `--metrics_mode structure_only` and ground truth JSON available
  - Random values (`np.random.uniform`) are **isolated behind explicit stub guard**
  - Provenance records `metrics_status` field

- `compile_reports.py`: Report structure real; calibration stats derived from stub SARs

#### Benchmark Manifests (Synthetic Fixtures)
- `benchmark_v1.0.json`: Contains 10 synthetic PDB IDs (`8ABC`–`8BCD`)
  - These are **not real PDB accessions**
  - Used only for pipeline mechanics testing
  - Retained in place because locked hash constrains file modification
  - Quarantined with documentation in `benchmark/synthetic/README.md`

- `benchmark_v1.0_rejected.json`: Contains 10 synthetic rejected IDs (`7AAA`–`7JJJ`)
  - Same status as accepted manifest

#### Existing SAR Results (`sar_results/`)
- Generated from a synthetic fixture run (2026-02-08)
- **Not real scientific outputs**
- All RMSD values (~17Å) are artifacts of stub random generation
- 100% REJECT rate reflects stub behavior, not real AF3 performance
- Retained for pipeline integrity testing but must not be cited as results

---

### ❌ Not Yet Real (Future Phases)

#### Phase 2 Prerequisites
- Resolving the provisional status of real kinase pilot target `8PAR`
- AF3 integration (`af3_available = True` in `run_inference.py`)
- Real inference outputs in `sar_results_raw/`
- SAR generation from real AF3 coordinates vs real crystal structures

#### Phase 2+ Scope
- Real calibration analysis (confidence vs actual error)
- Scientific interpretation of RMSD/decision gate distributions
- Comparison across kinase families
- Any paper or preprint claims

---

## Stub vs Real — Decision Tree

```
Is stub_output=True in provenance?
├── YES → Synthetic outputs; do NOT interpret as science
│         Metrics: deterministic stubs (not structural measurements)
└── NO  → Check metrics_status in provenance:
          ├── "synthetic"               → Same as above
          ├── "structure_only"          → Real metrics, but no AF3 inference yet
          │                               (engineering validation only)
          └── "inference_backed"        → Real AF3 coordinates used
                                          (not yet implemented)
```

---

## Provenance Fields to Trust

When reading SARs, check:

```json
{
  "provenance": {
    "stub_output": true,          // true = synthetic; false = real
    "metrics_status": "synthetic" // synthetic | structure_only | inference_backed
  }
}
```

Any SAR without both fields should be treated as synthetic until verified.

---

## What Is Needed to Reach Phase 2

1. **Target selection**: Solidify provisional targets (i.e. resolve `8PAR`).
2. **Ground truth preparation**: Ground truth for all targets exists in `pdb_ground_truth/`.
   ```bash
   python3 scripts/download_structures.py --pdb_ids <IDS> --output_dir benchmark/ground_truth
   python3 scripts/build_ground_truth_json.py --structure_dir benchmark/ground_truth --output_dir benchmark/ground_truth
   python3 scripts/extract_ligands.py --structure_dir benchmark/ground_truth --output_file benchmark/metadata/ligand_inventory.json
   ```

3. **AF3 integration**: Set `af3_available = True` in `run_inference.py`
   - Implement `run_alphafold3()` with real AF3 API/CLI calls
   - Ensure output format matches stub structure

4. **Re-authorization**: Any manifest change requires new hash authorization

---

## History

| Date | Event |
|---|---|
| 2026-02-08 | Initial scaffold created; synthetic manifests locked |
| 2026-03-13 | Phase 0 + Phase 1: truth cleanup; real metrics implemented; synthetic fixtures quarantined |
