# Claude Code Development Guide - Kinase Reliability Pilot

This guide helps Claude Code (and other AI assistants) work effectively with the Kinase Reliability Pilot codebase.

## Project Overview

The Kinase Reliability Pilot is a computational audit framework for AlphaFold 3 predictions. It processes 10 locked kinase targets through a 4-stage pipeline with strict integrity and validation requirements.

**State:** EXECUTE (locked)
**Critical:** Do NOT modify locked parameters, manifests, or schemas without explicit authorization

## Architecture

### Core Pipeline (4 Stages)

1. **Integrity Verification** (`sha256sum -c`)
   - Validates manifest checksums before execution
   - FAIL condition: Any hash mismatch

2. **Inference** (`run_inference.py`)
   - Runs AF3 with locked config (seed=42, recycles=3)
   - Processes targets in manifest serialized order
   - On failure: logs error, continues batch

3. **SAR Generation** (`generate_sar.py`)
   - Computes metrics (RMSD, pLDDT, PAE, contact maps)
   - Enforces strict schema (3 REQUIRED fields)
   - Assigns decision gates and failure taxonomy

4. **Report Compilation** (`compile_reports.py`)
   - Generates summary, calibration, provenance
   - Validates PASS/FAIL criteria
   - Writes audit-safe outputs

### Key Constraints

#### Locked Values (DO NOT MODIFY)

- `seed = 42`
- `recycles = 3`
- `schema_version = "1.0"`
- Manifest hashes (see `internal/manifest_checksums.sha256`)

#### Required SAR Fields

Every SAR MUST contain:
1. `expected_error_range` (object with rmsd_min, rmsd_max, rationale)
2. `recommended_action` (string, non-empty)
3. `decision_gate` (enum: ACCEPT/REVIEW/REJECT)

Missing any field → `ERROR_SAR_INCOMPLETE`

#### Processing Order

Targets MUST be processed in manifest serialized order (as listed in `benchmark_v1.0.json`). Do NOT sort, shuffle, or parallelize without preserving order.

## Common Development Tasks

### Adding New Metrics

When adding metrics to SARs:

1. Update `schemas/sar_schema_v1.json` (requires version bump)
2. Modify `generate_sar.py` to compute new metric
3. Update calibration logic in `compile_reports.py` if needed
4. Update `docs/EXECUTION_CONTROL_DOCUMENT.md` with justification

**Warning:** Schema changes require re-authorization and new locked hash.

### Debugging SAR Validation Errors

If you encounter `ERROR_SAR_INCOMPLETE`:

1. Check `generate_sar.py`: Ensure all 3 required fields are populated
2. Verify `decision_gate` is exactly one of: `ACCEPT`, `REVIEW`, `REJECT` (case-sensitive)
3. Validate `expected_error_range` contains `rmsd_min`, `rmsd_max`, `rationale`
4. Run with `--strict_mode` to get detailed error messages

### Testing Changes

```bash
# Test manifest generation
python3 scripts/generate_manifest.py --format pretty

# Verify checksums
sha256sum -c internal/manifest_checksums.sha256

# Test full pipeline (stub mode)
python3 run_inference.py --manifest benchmark_v1.0.json --model_version af3_stub --seed 42 --recycles 3 --output_dir ./test_output
python3 generate_sar.py --manifest benchmark_v1.0.json --pred_dir ./test_output --ground_truth_dir ./pdb_ground_truth --output_dir ./test_sars --schema_version 1.0 --strict_mode
python3 compile_reports.py --manifest benchmark_v1.0.json --sar_dir ./test_sars --job_id TEST_RUN --accepted_manifest_hash f26949fb43663df0c2d8d2f4f9b05f4f9edc58c53a0d7389ba49aa7e4b182218 --rejected_manifest_hash 724ffdde4439c8bb99e5a81f770971f69b61197e0af56bfea4c7509c0a28fd7c
```

### Understanding Failure Taxonomy

When debugging classification logic (`classify_failure` in `generate_sar.py`):

- **Class A (Overconfidence)**: `plddt_mean > 90 AND pae_mean < 5 AND rmsd > expected_max * 1.5`
- **Class B (Ligand Pose)**: `ligand_present AND rmsd_ligand > rmsd_global * 1.5`
- **Class C (Symmetry)**: `rmsd_global > 10.0Å`
- **Unknown**: Error exceeds expected range but no pattern match
- **N/A**: `rmsd_global <= expected_max`

### Modifying Decision Gate Logic

Decision gates in `determine_decision_gate`:

- **ACCEPT**: `rmsd <= expected_max`
- **REJECT**: `failure_class in [A, C] OR rmsd > expected_max * 2`
- **REVIEW**: All other cases

Changes to this logic require updating `docs/EXECUTION_CONTROL_DOCUMENT.md` with rationale.

## File Reference

### Critical Files (Do Not Modify Without Authorization)

- `EXECUTION_CONTROL_DOCUMENT_LOCKED.md` - Formal specification
- `benchmark_v1.0.json` - Accepted targets (locked hash)
- `benchmark_v1.0_rejected.json` - Rejected targets (locked hash)
- `internal/manifest_checksums.sha256` - Integrity hashes

### Modifiable Implementation Files

- `run_inference.py` - Inference logic (maintain interface contract)
- `generate_sar.py` - SAR generation (respect schema requirements)
- `compile_reports.py` - Report generation (maintain PASS/FAIL criteria)

### Configuration Files

- `schemas/sar_schema_v1.json` - SAR JSON schema (version locked)
- `requirements.txt` - Python dependencies

## Error Handling Patterns

### In `run_inference.py`

```python
try:
    # Process target
    result = run_alphafold3(...)
except Exception as e:
    # Log error, continue batch (DO NOT halt)
    execution_log.append({"pdb_id": pdb_id, "status": "error", "error": str(e)})
    continue
```

### In `generate_sar.py`

```python
# Validate SAR completeness
validate_sar(sar, strict_mode=args.strict_mode)
# If strict_mode: raises SARValidationError
# If not strict: warns but continues
```

### In `compile_reports.py`

```python
# Check PASS/FAIL criteria
if not validate_completeness(sars, manifest):
    print("FAIL: Completeness check failed", file=sys.stderr)
    sys.exit(1)  # Hard failure
```

## Integration with AlphaFold 3

### Current State (Stub Mode)

The scaffold uses deterministic stubs when AF3 is unavailable:

```python
af3_available = False  # Set to True when AF3 integrated

if af3_available:
    result = run_alphafold3(...)  # Production implementation
else:
    result = run_alphafold3_stub(...)  # Deterministic stub
```

### Production Integration

To integrate real AF3:

1. Set `af3_available = True` in `run_inference.py`
2. Implement `run_alphafold3()` function with AF3 API calls
3. Ensure output format matches stub structure
4. Remove `stub_output: true` flag from provenance
5. Update `requirements.txt` with AF3 dependencies

## Governance and Audit Trail

### Before Making Changes

1. Read `docs/EXECUTION_CONTROL_DOCUMENT.md`
2. Verify change doesn't affect locked parameters
3. Check if schema version bump needed
4. Plan for re-authorization if manifests change

### After Making Changes

1. Run full pipeline end-to-end
2. Verify all PASS criteria met
3. Update provenance with change rationale
4. Document changes in commit message
5. Update hashes if manifests modified

## Best Practices for Claude Code

1. **Always check locked state first**: Read `EXECUTION_CONTROL_DOCUMENT_LOCKED.md` before modifications
2. **Preserve manifest order**: Never shuffle or reorder targets
3. **Respect strict mode**: Don't bypass SAR validation in production
4. **Maintain provenance**: Log all configuration changes
5. **Test incrementally**: Verify each pipeline stage independently
6. **Document decisions**: Explain rationale for decision gate changes
7. **Version carefully**: Bump schema versions when adding required fields

## Useful Commands

```bash
# Quick validation
sha256sum benchmark_v1.0.json benchmark_v1.0_rejected.json scripts/generate_manifest.py

# Count decision gates in SARs
grep -h "decision_gate" sar_results/*.json | sort | uniq -c

# Check for missing required fields
python3 -c "import json, sys; sars = [json.load(open(f'sar_results/{p}.json')) for p in ['8ABC', '8DEF']]; missing = [p for p, s in zip(['8ABC', '8DEF'], sars) if 'decision_gate' not in s]; print('Missing:', missing)"

# Verify stub flag
grep -h "stub_output" sar_results_raw/*_prediction.json | sort | uniq -c
```

## Questions?

For governance questions, consult `docs/EXECUTION_CONTROL_DOCUMENT.md`.
For implementation questions, check function docstrings in Python scripts.
For schema questions, see `schemas/sar_schema_v1.json`.

**Remember:** When in doubt about locked parameters or schema changes, DO NOT PROCEED without explicit authorization.
