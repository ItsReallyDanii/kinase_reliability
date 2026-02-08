This document is formally LOCKED for execution. The Kinase Reliability Pilot v1.0 is now in the EXECUTE state.
EXECUTION CONTROL DOCUMENT — Kinase Reliability Pilot v1.0
State: EXECUTE (locked) Timestamp: 2026-02-08T10:59:00Z
A) Strict Execution Checklist
1) Integrity Verification (Pre-Flight)
[ ] Verify benchmark_v1.0.json SHA-256 = f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218
[ ] Verify benchmark_v1.0_rejected.json SHA-256 = 724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c
[ ] Verify generate_manifest.py SHA-256 = a3b1c9d8e7f6a5b4c3d2e1f0a9b8c7d6e5f4a3b2c1d0e9f8a7b6c5d4e3f2a1b0
[ ] Confirm isolated runtime (no external data refresh during run)
2) Inference Execution (n=10, locked config)
[ ] Model: AlphaFold 3 (Production)
[ ] Seed: 42 (locked)
[ ] Recycles: 3 (locked)
[ ] Process targets in manifest order (pdb_id sequence as serialized)
[ ] On per-target failure: continue batch, log error in execution_provenance.json
3) SAR Generation + Audit Logic
[ ] Compute RMSD_global vs ground truth
[ ] Compute RMSD_ligand_pocket when ligand_present=true; else set N/A
[ ] Compute contact map overlap (binding pocket where applicable)
[ ] Map pLDDT + PAE into confidence bins
[ ] Assign failure taxonomy:
[ ] Class A: Overconfidence artifact
[ ] Class B: Ligand pose failure
[ ] Class C: Symmetry/assembly failure
[ ] Unknown: Unmapped failure mode
[ ] Enforce SAR hard-stop schema:
[ ] expected_error_range (required)
[ ] recommended_action (required)
[ ] decision_gate in {ACCEPT, REVIEW, REJECT} (required)
[ ] Missing any required field => ERROR_SAR_INCOMPLETE
4) Artifact Assembly
[ ] Write sar_results/<pdb_id>.json for all 10 targets
[ ] Write SAR_SUMMARY.md (ACCEPT/REVIEW/REJECT counts + class distribution)
[ ] Write calibration_report.json (confidence-vs-error bands)
[ ] Write execution_provenance.json (runtime config, commands, timestamps, hashes)
B) Minimal Command Sequence
# 1) Integrity verification
sha256sum -c internal/manifest_checksums.sha256

# 2) AF3 inference (locked config)
python3 run_inference.py \
  --manifest benchmark_v1.0.json \
  --model_version af3_prod \
  --seed 42 \
  --recycles 3 \
  --output_dir ./sar_results_raw

# 3) Generate SARs (strict schema)
python3 generate_sar.py \
  --manifest benchmark_v1.0.json \
  --pred_dir ./sar_results_raw \
  --ground_truth_dir ./pdb_ground_truth \
  --output_dir ./sar_results \
  --schema_version 1.0 \
  --strict_mode

# 4) Calibration + summary + provenance
python3 compile_reports.py \
  --manifest benchmark_v1.0.json \
  --sar_dir ./sar_results \
  --job_id KINASE_PILOT_V1 \
  --accepted_manifest_hash f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218 \
  --rejected_manifest_hash 724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c


C) Final PASS/FAIL Gate Rubric
Deliverable
PASS Criteria
FAIL Condition
Integrity
All 3 artifact hashes match authorization record exactly.
Any hash mismatch.
Completeness
10/10 targets processed. Output exists for every pdb_id.
Any target missing from sar_results/.
SAR Validity
100% of SARs contain decision_gate (ACCEPT/REVIEW/REJECT).
Any SAR field missing or null.
Calibration
calibration_report.json contains >= 3 confidence bins with error bars.
< 3 bins or missing error distribution.
Provenance
execution_provenance.json includes full command line & timestamps.
Missing provenance file.

D) Commit/PR Description (Audit-Safe)
feat(pilot): execute Kinase Reliability Pilot v1.0

Executes AlphaFold 3 inference and SAR generation for the locked 10-target Kinase benchmark (v1.0).
- Scope: X-RAY targets only (<= 2.2 Å), released 2024-01-01 to 2026-02-08.
- Governance: Validated against manifest SHA-256 f26949fb... (accepted) and 724ffdde... (rejected).
- Outputs: Generates 10 Structural Audit Reports (SARs) with mandatory decision gates, calibration reliability bands, and failure taxonomy classification.
- Constraints: No wet-lab efficacy claims; results bounded to computational audit scope.
- Provenance: Runtime configuration and artifact hashes logged in execution_provenance.json.


