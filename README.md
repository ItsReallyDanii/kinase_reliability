# Kinase Reliability Pilot v1.0

Computational audit framework for AlphaFold 3 predictions on kinase structures.

## Overview

The Kinase Reliability Pilot executes a locked, reproducible benchmark to evaluate AlphaFold 3 performance on 10 high-quality kinase crystal structures. The pipeline generates Structural Audit Reports (SARs) with mandatory decision gates and failure taxonomy classification.

**State:** EXECUTE (locked)
**Version:** 1.0
**Timestamp:** 2026-02-08T10:59:00Z

## Key Features

- **Locked Configuration**: Seed=42, Recycles=3 (immutable)
- **Strict SAR Schema**: Enforces `expected_error_range`, `recommended_action`, `decision_gate ∈ {ACCEPT, REVIEW, REJECT}`
- **Failure Taxonomy**: Classes A/B/C for systematic error characterization
- **Integrity Verification**: SHA-256 checksums for all manifests and artifacts
- **Provenance Tracking**: Full runtime configuration and command history

## Scope

- **n=10 targets**: X-RAY structures, resolution ≤ 2.2Å, released 2024-01-01 to 2026-02-08
- **Protein class**: Kinases across CMGC, TK, AGC, CAMK, STE, TKL families
- **Constraints**: No wet-lab efficacy claims; results bounded to computational audit scope

## Quick Start

### Installation

```bash
# Clone repository
git clone <repository-url>
cd kinase_reliability

# Install dependencies
pip install -r requirements.txt
```

### Execution (4-Step Pipeline)

```bash
# 1. Verify manifest integrity
sha256sum -c internal/manifest_checksums.sha256

# 2. Run AF3 inference (locked config)
python3 run_inference.py \
  --manifest benchmark_v1.0.json \
  --model_version af3_prod \
  --seed 42 \
  --recycles 3 \
  --output_dir ./sar_results_raw

# 3. Generate SARs (strict schema)
python3 generate_sar.py \
  --manifest benchmark_v1.0.json \
  --pred_dir ./sar_results_raw \
  --ground_truth_dir ./pdb_ground_truth \
  --output_dir ./sar_results \
  --schema_version 1.0 \
  --strict_mode

# 4. Compile reports (calibration + provenance)
python3 compile_reports.py \
  --manifest benchmark_v1.0.json \
  --sar_dir ./sar_results \
  --job_id KINASE_PILOT_V1 \
  --accepted_manifest_hash f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218 \
  --rejected_manifest_hash 724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c
```

## Directory Structure

```
kinase_reliability/
├── benchmark_v1.0.json              # Accepted targets manifest (n=10)
├── benchmark_v1.0_rejected.json     # Rejected targets manifest
├── run_inference.py                 # AF3 inference runner
├── generate_sar.py                  # SAR generator
├── compile_reports.py               # Report compiler
├── requirements.txt                 # Python dependencies
├── .gitignore                       # Git ignore patterns
├── README.md                        # This file
├── CLAUDE.md                        # Development guide
├── internal/
│   └── manifest_checksums.sha256    # Integrity checksums
├── scripts/
│   └── generate_manifest.py         # Manifest generation utility
├── schemas/
│   └── sar_schema_v1.json          # SAR JSON schema
├── sar_results/                     # SAR outputs (generated)
├── pdb_ground_truth/                # Ground truth structures (user-provided)
└── docs/
    └── EXECUTION_CONTROL_DOCUMENT.md # Locked execution specification
```

## Outputs

### Per-Target SARs (`sar_results/<pdb_id>.json`)

Each SAR contains:
- **Metrics**: RMSD (global, ligand pocket), contact map overlap, pLDDT, PAE
- **Confidence Assessment**: High/medium/low bins for pLDDT and PAE
- **Expected Error Range**: RMSD bounds based on confidence (REQUIRED)
- **Recommended Action**: Actionable guidance (REQUIRED)
- **Decision Gate**: ACCEPT/REVIEW/REJECT (REQUIRED)
- **Failure Taxonomy**: Class A/B/C/Unknown/N/A classification
- **Provenance**: Model version, seed, recycles, stub flag

### Summary Reports (`sar_results/`)

- **SAR_SUMMARY.md**: Decision gate counts, failure taxonomy distribution, per-target table
- **calibration_report.json**: Confidence-vs-error bands (≥3 bins required)
- **execution_provenance.json**: Full command line, timestamps, manifest hashes

## PASS/FAIL Criteria

| Deliverable    | PASS Criteria                                      | FAIL Condition                     |
|----------------|----------------------------------------------------|------------------------------------|
| Integrity      | All 3 artifact hashes match authorization record   | Any hash mismatch                  |
| Completeness   | 10/10 targets processed                            | Any target missing from sar_results/ |
| SAR Validity   | 100% of SARs contain decision_gate                | Any SAR field missing or null      |
| Calibration    | calibration_report.json has ≥3 confidence bins    | <3 bins or missing error distribution |
| Provenance     | execution_provenance.json includes full command   | Missing provenance file            |

## Failure Taxonomy

- **Class A**: Overconfidence artifact (high confidence, high error)
- **Class B**: Ligand pose failure (ligand RMSD >> global RMSD)
- **Class C**: Symmetry/assembly failure (extreme misalignment)
- **Unknown**: Unmapped failure mode
- **N/A**: No failure detected (within expected error range)

## Decision Gates

- **ACCEPT**: Prediction within expected error range → proceed with downstream analysis
- **REVIEW**: Moderate failures or unclear modes → manual expert review required
- **REJECT**: Critical failures (overconfidence, symmetry) → do not use without refinement

## Locked Parameters

- **Seed**: 42 (fixed for reproducibility)
- **Recycles**: 3 (fixed for fair comparison)
- **Schema Version**: 1.0
- **Manifest Hashes**:
  - Accepted: `f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218`
  - Rejected: `724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c`

## Stub Mode (AF3 Unavailable)

If AlphaFold 3 is unavailable, the pipeline generates deterministic stub outputs:
- Reproducible pseudo-predictions based on `pdb_id` hash and `seed`
- SAR provenance includes `stub_output: true` flag
- Enables end-to-end pipeline testing without AF3 installation

## Governance

All changes to manifests, schemas, and locked parameters require explicit authorization. See `docs/EXECUTION_CONTROL_DOCUMENT.md` for the formal specification.

## Citation

If you use this pilot in your work, please cite:

```
Kinase Reliability Pilot v1.0 (2026)
Computational audit framework for AlphaFold 3 kinase predictions
https://github.com/<repository-url>
```

## License

[Specify license]

## Contact

[Specify contact information]
