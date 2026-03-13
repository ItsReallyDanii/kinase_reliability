"""
Real RMSD computation using the Kabsch algorithm (SVD-based optimal superposition).

All functions are deterministic given their inputs. No random values.

Two entry points:
  - numpy-only:   compute_rmsd(coords_a, coords_b)
  - Biopython:    compute_ca_rmsd_biopython(struct_a, struct_b)

The numpy path is always available. The Biopython path requires `biopython>=1.81`.
"""

from __future__ import annotations

import numpy as np
from typing import Tuple


# ---------------------------------------------------------------------------
# Core Kabsch algorithm (numpy only)
# ---------------------------------------------------------------------------


def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return the optimal rotation matrix R that minimises RMSD(P @ R, Q).

    Both P and Q must already be **centred** (zero mean).

    Uses SVD with a reflection-correction to guarantee a proper rotation
    (det = +1).

    Args:
        P: (N, 3) centred coordinates to be rotated
        Q: (N, 3) centred reference coordinates

    Returns:
        R: (3, 3) rotation matrix
    """
    H = P.T @ Q                             # cross-covariance (3x3)
    U, _S, Vt = np.linalg.svd(H)

    # Correct for reflection: ensure det(V U^T) = +1
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])

    R = Vt.T @ D @ U.T
    return R


def compute_rmsd(
    coords_a: np.ndarray | list,
    coords_b: np.ndarray | list,
    align: bool = True,
) -> float:
    """Compute RMSD between two sets of coordinates.

    Args:
        coords_a: (N, 3) array of coordinates (e.g. C-alpha positions, Å)
        coords_b: (N, 3) array of coordinates for the reference structure
        align:    If True (default), apply Kabsch optimal superposition before
                  computing RMSD.  If False, compute raw (unaligned) RMSD.

    Returns:
        RMSD in the same units as the input coordinates (typically Å).

    Raises:
        ValueError: if arrays have incompatible shapes or are empty.
    """
    coords_a = np.asarray(coords_a, dtype=float)
    coords_b = np.asarray(coords_b, dtype=float)

    if coords_a.ndim != 2 or coords_a.shape[1] != 3:
        raise ValueError(
            f"coords_a must be (N, 3), got {coords_a.shape}"
        )
    if coords_b.ndim != 2 or coords_b.shape[1] != 3:
        raise ValueError(
            f"coords_b must be (N, 3), got {coords_b.shape}"
        )
    if coords_a.shape != coords_b.shape:
        raise ValueError(
            f"Shape mismatch: {coords_a.shape} vs {coords_b.shape}. "
            "Trim to the same length before calling."
        )
    if len(coords_a) == 0:
        raise ValueError("Cannot compute RMSD of empty coordinate arrays.")

    if align:
        cen_a = coords_a.mean(axis=0)
        cen_b = coords_b.mean(axis=0)
        P = coords_a - cen_a
        Q = coords_b - cen_b
        R = kabsch_rotation(P, Q)
        P_rot = P @ R.T
        diff = P_rot - Q
    else:
        diff = coords_a - coords_b

    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def truncate_to_common_length(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Trim both coordinate arrays to their common (shorter) length.

    This is a safe fallback for comparing structures with different numbers
    of resolved residues.  A provenance note should always be recorded when
    trimming is applied.

    Args:
        coords_a: (N, 3) array
        coords_b: (M, 3) array

    Returns:
        Tuple of two (min(N,M), 3) arrays.
    """
    n = min(len(coords_a), len(coords_b))
    return coords_a[:n], coords_b[:n]


# ---------------------------------------------------------------------------
# Biopython path (optional; raises ImportError with a clear message if absent)
# ---------------------------------------------------------------------------


def extract_ca_coords(structure) -> np.ndarray:
    """Extract C-alpha coordinates from a Biopython Structure object.

    Args:
        structure: Bio.PDB.Structure.Structure instance

    Returns:
        (N, 3) float64 array of C-alpha positions in Å.

    Raises:
        ImportError: if Biopython is not installed.
        ValueError:  if no C-alpha atoms are found.
    """
    try:
        from Bio.PDB import PPBuilder  # noqa: F401 (import check)
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for extract_ca_coords. "
            "Install it with: pip install biopython>=1.81"
        ) from exc

    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca_coords.append(residue["CA"].get_vector().get_array())
        break  # use first model only

    if not ca_coords:
        raise ValueError(
            f"No C-alpha atoms found in structure '{structure.id}'. "
            "Check that the structure contains standard amino acid residues."
        )

    return np.array(ca_coords, dtype=float)


def compute_ca_rmsd_biopython(struct_a, struct_b, align: bool = True) -> float:
    """Compute global C-alpha RMSD between two Biopython Structure objects.

    Extracts C-alpha coordinates from both structures, trims to common length
    if necessary, then calls compute_rmsd().

    Args:
        struct_a: Bio.PDB.Structure.Structure (predicted)
        struct_b: Bio.PDB.Structure.Structure (reference / ground truth)
        align:    Apply Kabsch superposition (default True).

    Returns:
        C-alpha RMSD in Å.
    """
    coords_a = extract_ca_coords(struct_a)
    coords_b = extract_ca_coords(struct_b)

    if len(coords_a) != len(coords_b):
        coords_a, coords_b = truncate_to_common_length(coords_a, coords_b)

    return compute_rmsd(coords_a, coords_b, align=align)
