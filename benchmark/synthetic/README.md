# Synthetic Benchmark Fixtures

## WARNING: These are not real scientific targets

This directory documents the synthetic fixture manifests used for pipeline
mechanics testing. The PDB IDs in these manifests do **not** correspond to
real structures deposited in the RCSB Protein Data Bank.

---

## Synthetic Manifests

The following files at the repository root are **synthetic fixtures**:

| File | Description | Status |
|---|---|---|
| `../../benchmark_v1.0.json` | 10 synthetic accepted targets (`8ABC`–`8BCD`) | Synthetic fixture, locked hash |
| `../../benchmark_v1.0_rejected.json` | 10 synthetic rejected examples (`7AAA`–`7JJJ`) | Synthetic fixture, locked hash |

These files were locked in place at `2026-02-08T10:59:00Z` for pipeline
integrity testing. They **cannot be modified** without re-authorization (the
locked hashes in `internal/manifest_checksums.sha256` would be invalidated).

---

## Why are these retained at the root?

1. The pipeline uses them for end-to-end mechanics testing
2. Locked SHA-256 hashes bind their exact content to the governance trail
3. Moving or renaming them would invalidate the `sha256sum -c` integrity check

They are kept for pipeline testing purposes only.

---

## Synthetic PDB IDs Used

| Synthetic ID | Claimed kinase | Reality |
|---|---|---|
| `8ABC` | CDK2 | Not a real PDB accession |
| `8DEF` | EGFR | Not a real PDB accession |
| `8GHI` | PKA | Not a real PDB accession |
| `8JKL` | CaMKII | Not a real PDB accession |
| `8MNO` | ABL1 | Not a real PDB accession |
| `8PQR` | ERK2 | Not a real PDB accession |
| `8STU` | MEK1 | Not a real PDB accession |
| `8VWX` | SRC | Not a real PDB accession |
| `8YZA` | AKT1 | Not a real PDB accession |
| `8BCD` | BRAF | Not a real PDB accession |

---

## Real Targets

Real kinase benchmark targets will be added to:
```
benchmark/metadata/kinase_pilot_targets.json
```

That file contains the format specification and placeholder entries for the
actual 2024+ benchmark targets once selection is finalized.

---

## Outputs Generated From Synthetic Fixtures

Any outputs in `sar_results/` produced by running the pipeline against
`benchmark_v1.0.json` are **not scientific results**. They reflect:

- Deterministic stub coordinates (not AF3 predictions)
- RMSD values computed from random arrays (not structural comparisons)
- Decision gates and failure classes driven by stub metric ranges

Do not cite these outputs as reflecting AlphaFold 3 performance.
