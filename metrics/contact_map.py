"""
Real contact map computation based on inter-residue C-alpha distances.

All functions are deterministic given their inputs.  No random values.

A contact is defined as a pair of residues (i, j) with |i - j| >= 2 whose
C-alpha atoms are within `threshold` Å of each other (default 8.0 Å).
"""

from __future__ import annotations

import numpy as np


def compute_contact_map(
    coords: np.ndarray | list,
    threshold: float = 8.0,
    min_seq_separation: int = 2,
) -> np.ndarray:
    """Compute a binary C-alpha contact map from coordinate array.

    Args:
        coords:              (N, 3) array of C-alpha positions in Å.
        threshold:           Distance cutoff in Å.  Pairs within this distance
                             are considered in contact (default 8.0 Å).
        min_seq_separation:  Minimum sequence separation |i - j| for a contact
                             to be counted (default 2; excludes immediate
                             neighbours and self).

    Returns:
        (N, N) uint8 numpy array.  Element [i, j] = 1 if residues i and j are
        in contact (distance < threshold AND |i-j| >= min_seq_separation),
        0 otherwise.  The matrix is symmetric.

    Raises:
        ValueError: if coords is not (N, 3) or N < 1.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Expected (N, 3) array, got {coords.shape}")
    if len(coords) == 0:
        raise ValueError("Cannot compute contact map for empty coordinate array.")

    n = len(coords)

    # Pairwise squared distances via broadcasting — O(N^2) but fine for typical
    # kinase domain sizes (~250–350 residues).
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]   # (N, N, 3)
    dist2 = np.sum(diff ** 2, axis=2)                              # (N, N)

    contact = (dist2 < threshold ** 2).astype(np.uint8)

    # Zero out diagonal and near-diagonal (sequence separation < min_seq_separation)
    for k in range(min_seq_separation):
        idx = np.arange(n - k)
        contact[idx, idx + k] = 0
        contact[idx + k, idx] = 0

    return contact


def compute_contact_map_overlap(
    map_a: np.ndarray,
    map_b: np.ndarray,
) -> float:
    """Compute Jaccard overlap between two binary contact maps.

    Jaccard = |intersection| / |union|

    Args:
        map_a: (N, N) binary contact map (0/1 or bool)
        map_b: (N, N) binary contact map (0/1 or bool), same shape as map_a

    Returns:
        Float in [0.0, 1.0].  Returns 0.0 if both maps have no contacts
        (union is empty).

    Raises:
        ValueError: if map shapes differ.
    """
    map_a = np.asarray(map_a, dtype=bool)
    map_b = np.asarray(map_b, dtype=bool)

    if map_a.shape != map_b.shape:
        raise ValueError(
            f"Contact map shape mismatch: {map_a.shape} vs {map_b.shape}"
        )

    intersection = int(np.logical_and(map_a, map_b).sum())
    union = int(np.logical_or(map_a, map_b).sum())

    if union == 0:
        return 0.0

    return float(intersection) / float(union)


def contact_map_from_structure(structure, threshold: float = 8.0) -> np.ndarray:
    """Compute C-alpha contact map from a Biopython Structure object.

    Args:
        structure: Bio.PDB.Structure.Structure instance
        threshold: Distance cutoff in Å (default 8.0 Å)

    Returns:
        (N, N) uint8 contact map.

    Raises:
        ImportError: if Biopython is not installed.
    """
    try:
        from metrics.rmsd import extract_ca_coords
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for contact_map_from_structure. "
            "Install it with: pip install biopython>=1.81"
        ) from exc

    coords = extract_ca_coords(structure)
    return compute_contact_map(coords, threshold=threshold)
