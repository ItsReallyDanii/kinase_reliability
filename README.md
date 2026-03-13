# Kinase Reliability Pilot

Governance, schema, and provenance scaffold for a planned computational audit
of AlphaFold 3 predictions on kinase structures.

---

## ⚠ Current Status: Level A — Engineering Complete

**AlphaFold 3 has NOT been run on any real targets yet.**
This repository currently provides:

- Governance and schema infrastructure (real)
- Pipeline scaffolding with stub/synthetic inference (real scaffold, synthetic outputs)
- SAR schema enforcement and provenance tracking (real)
- Real metric implementations (RMSD, contact map, ligand RMSD proxy) for engineering validation
- Real structure download and ground truth extraction scripts
- Level A Provisional Pilot Set (`benchmark/pilot_set_v1.json`) with downloaded structures and deterministic `structure_only` pipeline execution. Note that `8PAR` remains provisional.

See [docs/CURRENT_STATUS.md](docs/CURRENT_STATUS.md) for the full claims ledger.

---

## Claims Ledger

| Capability | Status | Notes |
|---|---|---|
| Governance / schema / provenance scaffold | **Implemented** | SAR schema v1.0, locked params, audit trail |
| SAR required-field enforcement (`decision_gate`, etc.) | **Implemented** | strict_mode validation |
| Pipeline orchestration (integrity → inference → SAR → report) | **Scaffolded** | Uses stub inference; real AF3 not integrated |
| Stub/synthetic inference runner | **Scaffolded** | Deterministic stubs only; `stub_output: true` in provenance |
| Real RMSD (Kabsch/SVD) | **Implemented** | `metrics/rmsd.py`; requires Biopython/numpy |
| Real contact map (distance-based) | **Implemented** | `metrics/contact_map.py` |
| Real structure download (RCSB) | **Implemented** | `scripts/download_structures.py` |
| Structure → ground-truth JSON bridge | **Implemented** | `scripts/build_ground_truth_json.py` |
| `structure_only` pipeline path (end-to-end) | **Implemented** | download → build_ground_truth_json → generate_sar |
| Real ligand extraction | **Implemented** | `scripts/extract_ligands.py`; requires Biopython |
| Benchmark target metadata scaffold | **Implemented** | `benchmark/pilot_set_v1.json` |
| Real kinase benchmark targets (2024+) | **Implemented** | Level A Provisional Pilot Set (10 targets; `8PAR` is provisional) |
| AlphaFold 3 inference on real targets | **Not yet real** | AF3 integration not implemented |
| Real calibration against AF3 outputs | **Not yet real** | Requires real inference first |
| Scientific evaluation / paper claims | **Not yet real** | Phase 2+ scope |

---

## Project Overview

The Kinase Reliability Pilot intends to execute a locked, reproducible benchmark
evaluating AlphaFold 3 performance on high-quality kinase crystal structures.
The pipeline generates Structural Audit Reports (SARs) with mandatory decision
gates and failure taxonomy classification.

**Current Pipeline State:** Stub inference only (AF3 not integrated)
**Schema State:** Locked v1.0
**Benchmark Targets:** Level A Provisional Pilot Set (`8PAR` provisional)

---

## Repository Structure

```
kinase_reliability/
├── README.md                              # This file
├── CLAUDE.md                              # AI assistant development guide
├── EXECUTION_CONTROL_DOCUMENT_LOCKED.md  # Formal execution specification (locked)
│
├── run_inference.py                       # AF3 inference runner (stub mode only)
├── generate_sar.py                        # SAR generator (real metrics + stub paths)
├── compile_reports.py                     # Report compiler
│
├── benchmark_v1.0.json                    # [SYNTHETIC FIXTURE] Pipeline test manifest
├── benchmark_v1.0_rejected.json           # [SYNTHETIC FIXTURE] Rejection examples
│
├── benchmark/
│   ├── synthetic/                         # Clearly labeled synthetic fixtures
│   │   └── README.md                      # Explains synthetic nature
│   ├── metadata/
│   │   └── kinase_pilot_targets.json      # Real benchmark metadata scaffold (TBD)
│   ├── targets/                           # Future: real target structures
│   └── ground_truth/                      # Future: downloaded ground truth
│
├── metrics/
│   ├── rmsd.py                            # Real Kabsch/SVD RMSD (no randomness)
│   └── contact_map.py                     # Real distance-based contact maps
│
├── scripts/
│   ├── generate_manifest.py               # Manifest generation utility
│   ├── download_structures.py             # Download real structures from RCSB
│   ├── build_ground_truth_json.py         # Convert .cif/.pdb → ground_truth JSON (bridge)
│   └── extract_ligands.py                 # Extract ligand metadata from structures
│
├── tests/
│   ├── test_rmsd_controls.py              # Identity / perturbation / mismatch controls
│   ├── test_contact_map.py                # Contact map determinism controls
│   ├── test_metadata.py                   # Benchmark metadata schema validation
│   ├── test_sar_provenance.py             # SAR synthetic-vs-real status labeling
│   └── test_build_ground_truth.py         # Bridge script end-to-end controls
│
├── schemas/
│   └── sar_schema_v1.json                 # SAR JSON schema (locked)
├── internal/
│   └── manifest_checksums.sha256          # Integrity checksums
├── docs/
│   ├── EXECUTION_CONTROL_DOCUMENT.md      # Execution specification
│   └── CURRENT_STATUS.md                  # Truth-status claims ledger
├── pdb_ground_truth/                      # Legacy ground truth dir (use benchmark/ground_truth/)
└── sar_results/                           # Generated outputs (synthetic fixture run)
```

---

## Quick Start

### Installation

```bash
pip install -r requirements.txt
```

### Run Metric Validation Controls

Validates that real RMSD and contact map implementations are correct before
any model inference:

```bash
pytest tests/ -v
```

### Download Real Structures and Run in structure_only Mode

```bash
# 1. Download known kinase structures
python3 scripts/download_structures.py \
  --pdb_ids 1ATP 2ITO 1IEP \
  --output_dir benchmark/ground_truth \
  --format mmcif

# 2. Convert to ground-truth JSON (the bridge step)
python3 scripts/build_ground_truth_json.py \
  --structure_dir benchmark/ground_truth \
  --output_dir benchmark/ground_truth

# 3. Extract ligand metadata
python3 scripts/extract_ligands.py \
  --structure_dir benchmark/ground_truth \
  --output_file benchmark/metadata/ligand_inventory.json

# 4. Run stub inference against the real PDB IDs
#    (requires a manifest listing those IDs — use benchmark/metadata/ entries
#    or create a local manifest; see docs/CURRENT_STATUS.md)

# 5. Generate SARs with real structural metrics
python3 generate_sar.py \
  --manifest <your_manifest.json> \
  --pred_dir ./sar_results_raw \
  --ground_truth_dir benchmark/ground_truth \
  --output_dir ./sar_results \
  --schema_version 1.0 \
  --metrics_mode structure_only \
  --strict_mode
```

**Note:** Step 4 currently uses stub inference (`stub_output: false` is not
available until AF3 is integrated). In the meantime, `structure_only` mode
validates the metric pipeline with real coordinates against themselves or
a deliberately perturbed copy — useful for confirming RMSD ≈ 0 on identity
comparisons before any model is run.

### Run Full Pipeline on Synthetic Fixtures

**Note:** This runs on synthetic PDB IDs with stub inference only.
Outputs are NOT real scientific results.

```bash
# 1. Verify manifest integrity
sha256sum -c internal/manifest_checksums.sha256

# 2. Run stub inference (NOT real AF3)
python3 run_inference.py \
  --manifest benchmark_v1.0.json \
  --model_version af3_stub \
  --seed 42 \
  --recycles 3 \
  --output_dir ./sar_results_raw

# 3. Generate SARs (stub metrics path)
python3 generate_sar.py \
  --manifest benchmark_v1.0.json \
  --pred_dir ./sar_results_raw \
  --ground_truth_dir ./pdb_ground_truth \
  --output_dir ./sar_results \
  --schema_version 1.0 \
  --strict_mode

# 4. Compile reports
python3 compile_reports.py \
  --manifest benchmark_v1.0.json \
  --sar_dir ./sar_results \
  --job_id SYNTHETIC_TEST_RUN \
  --accepted_manifest_hash f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218 \
  --rejected_manifest_hash 724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c
```

---

## Synthetic Fixture Warning

The default manifests (`benchmark_v1.0.json`, `benchmark_v1.0_rejected.json`)
use **synthetic PDB IDs** (e.g., `8ABC`, `8DEF`) that do not correspond to
real deposited structures. They exist solely to test pipeline mechanics.

- Do NOT cite outputs from these fixtures as scientific results
- Do NOT treat decision gates or RMSD values from stub runs as real AF3 evaluations
- The `sar_results/` directory in this repository contains outputs of a synthetic fixture run

See `benchmark/synthetic/README.md` for details.

---

## SAR Schema

Every SAR must contain three required fields:

| Field | Type | Constraint |
|---|---|---|
| `expected_error_range` | object | `rmsd_min`, `rmsd_max`, `rationale` all required |
| `recommended_action` | string | Non-empty |
| `decision_gate` | enum | `ACCEPT` \| `REVIEW` \| `REJECT` |

Missing any field → `ERROR_SAR_INCOMPLETE` (hard stop in strict mode).

### Decision Gates

- **ACCEPT**: RMSD within expected range → proceed with downstream analysis
- **REVIEW**: Moderate failures or unclear modes → manual expert review required
- **REJECT**: Critical failures (overconfidence, symmetry) → do not use without refinement

### Failure Taxonomy

- **Class A**: Overconfidence artifact (high confidence, high error)
- **Class B**: Ligand pose failure (ligand RMSD >> global RMSD)
- **Class C**: Symmetry/assembly failure (extreme misalignment)
- **Unknown**: Unmapped failure mode
- **N/A**: No failure detected (within expected error range)

---

## Locked Parameters

These values are frozen for the benchmark and must not be modified without
explicit re-authorization:

- **Seed**: 42
- **Recycles**: 3
- **Schema Version**: 1.0

See `EXECUTION_CONTROL_DOCUMENT_LOCKED.md` for the full governance specification.

---

## What Remains for Phase 2+

1. Real kinase target selection (2024+ PDB structures meeting inclusion criteria)
2. AlphaFold 3 integration (`af3_available = True` in `run_inference.py`)
3. Real ground truth structure ingestion into `benchmark/ground_truth/`
4. End-to-end pipeline run with real inference outputs
5. Calibration analysis against real AF3 confidence metrics
6. Scientific interpretation of results

---

## Governance

All changes to manifests, schemas, and locked parameters require explicit
authorization. See `docs/EXECUTION_CONTROL_DOCUMENT.md`.

## License

[Specify license]

## Contact

[Specify contact information]
