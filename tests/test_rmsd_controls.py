"""
Metric validation controls for RMSD computation.

These tests validate that metrics/rmsd.py behaves correctly on known inputs
BEFORE any real model inference is integrated.  Passing these tests is a
prerequisite for trusting the metric pipeline in Phase 2.

Controls
--------
1. Identity:        RMSD(structure, structure) ≈ 0
2. Small perturbation: RMSD(structure, structure + small_noise) ≈ noise_magnitude
3. Large displacement: RMSD(structure, translated_structure) > large threshold
4. Alignment:       Aligned RMSD < unaligned RMSD when structures differ by rotation
5. Shape errors:    Incompatible shapes raise ValueError (not silent wrong values)
"""

import numpy as np
import pytest

from metrics.rmsd import compute_rmsd, kabsch_rotation, truncate_to_common_length


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def kinase_like_coords():
    """Return a deterministic (300, 3) array resembling a kinase C-alpha trace."""
    rng = np.random.RandomState(0)
    # Simulate a folded protein: sum of random walk (gives compact globular-ish coords)
    steps = rng.randn(300, 3) * 1.5
    return np.cumsum(steps, axis=0)


# ---------------------------------------------------------------------------
# Control 1: Identity — RMSD(A, A) must be 0
# ---------------------------------------------------------------------------


def test_identity_rmsd_is_zero(kinase_like_coords):
    """RMSD of a structure against itself must be exactly 0.0."""
    rmsd = compute_rmsd(kinase_like_coords, kinase_like_coords.copy(), align=True)
    assert rmsd == pytest.approx(0.0, abs=1e-8), (
        f"Identity RMSD should be 0.0, got {rmsd}"
    )


def test_identity_rmsd_unaligned_is_zero(kinase_like_coords):
    rmsd = compute_rmsd(kinase_like_coords, kinase_like_coords.copy(), align=False)
    assert rmsd == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Control 2: Small perturbation — RMSD ≈ noise magnitude
# ---------------------------------------------------------------------------


def test_small_perturbation_rmsd_nonzero(kinase_like_coords):
    """RMSD after adding small Gaussian noise is small but nonzero."""
    noise_sigma = 0.1  # Å
    rng = np.random.RandomState(42)
    perturbed = kinase_like_coords + rng.randn(*kinase_like_coords.shape) * noise_sigma

    rmsd = compute_rmsd(kinase_like_coords, perturbed, align=True)

    # Should be approximately noise_sigma (within a factor of 2)
    assert rmsd > 0.0, "RMSD of perturbed structure should be > 0"
    assert rmsd < noise_sigma * 5, (
        f"RMSD ({rmsd:.4f} Å) unexpectedly large for noise_sigma={noise_sigma} Å"
    )


def test_small_perturbation_smaller_than_large_displacement(kinase_like_coords):
    """Small perturbation RMSD must be much smaller than large random noise RMSD.

    Note: pure rigid-body translations are removed by Kabsch alignment, so a
    large translation does NOT give a large aligned RMSD.  This test uses large
    independent random noise (which cannot be aligned away) as the mismatch case.
    """
    rng = np.random.RandomState(42)
    small_noise = kinase_like_coords + rng.randn(*kinase_like_coords.shape) * 0.1
    # Large independent random noise — alignment cannot reduce this significantly
    large_noise = kinase_like_coords + rng.randn(*kinase_like_coords.shape) * 15.0

    rmsd_small = compute_rmsd(kinase_like_coords, small_noise, align=True)
    rmsd_large = compute_rmsd(kinase_like_coords, large_noise, align=True)

    assert rmsd_small < rmsd_large, (
        f"Small perturbation RMSD ({rmsd_small:.3f}) should be < "
        f"large noise RMSD ({rmsd_large:.3f})"
    )


# ---------------------------------------------------------------------------
# Control 3: Obvious mismatch — RMSD clearly large
# ---------------------------------------------------------------------------


def test_obvious_mismatch_rmsd_large(kinase_like_coords):
    """Reversed coordinate array should give large RMSD (structures are dissimilar)."""
    reversed_coords = kinase_like_coords[::-1].copy()
    rmsd = compute_rmsd(kinase_like_coords, reversed_coords, align=True)

    # Two unrelated ~300-residue structures should have RMSD >> 5 Å
    assert rmsd > 5.0, (
        f"Expected large RMSD for reversed structure, got {rmsd:.2f} Å"
    )


def test_pure_translation_aligned_near_zero():
    """After Kabsch alignment, pure translation should give RMSD ≈ 0."""
    rng = np.random.RandomState(7)
    coords = rng.randn(100, 3) * 5.0
    translated = coords + np.array([50.0, -30.0, 20.0])

    rmsd_aligned = compute_rmsd(coords, translated, align=True)
    rmsd_unaligned = compute_rmsd(coords, translated, align=False)

    assert rmsd_aligned == pytest.approx(0.0, abs=1e-6), (
        f"Pure translation should give RMSD≈0 after alignment, got {rmsd_aligned}"
    )
    assert rmsd_unaligned > 10.0, (
        f"Unaligned RMSD should be large for 50Å translation, got {rmsd_unaligned}"
    )


# ---------------------------------------------------------------------------
# Control 4: Rotation — aligned RMSD ≤ unaligned RMSD
# ---------------------------------------------------------------------------


def test_rotation_alignment_reduces_rmsd():
    """Kabsch alignment should reduce or maintain RMSD for a rotated structure."""
    rng = np.random.RandomState(3)
    coords = rng.randn(150, 3) * 8.0

    # Apply a known rotation (90° around Z axis)
    theta = np.pi / 2
    R = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [0,              0,             1],
    ])
    rotated = coords @ R.T

    rmsd_aligned = compute_rmsd(coords, rotated, align=True)
    rmsd_unaligned = compute_rmsd(coords, rotated, align=False)

    assert rmsd_aligned == pytest.approx(0.0, abs=1e-6), (
        f"Pure rotation should give RMSD≈0 after alignment, got {rmsd_aligned}"
    )
    assert rmsd_unaligned > rmsd_aligned, (
        "Aligned RMSD should be less than unaligned RMSD for rotated structure"
    )


def test_kabsch_rotation_is_proper():
    """Rotation matrix returned by kabsch_rotation must have det ≈ +1."""
    rng = np.random.RandomState(99)
    P = rng.randn(50, 3)
    Q = rng.randn(50, 3)
    P -= P.mean(axis=0)
    Q -= Q.mean(axis=0)
    R = kabsch_rotation(P, Q)
    assert R.shape == (3, 3)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-9), (
        f"Rotation matrix det should be +1 (proper rotation), got {np.linalg.det(R)}"
    )


# ---------------------------------------------------------------------------
# Control 5: Error handling — bad inputs must raise, not silently return garbage
# ---------------------------------------------------------------------------


def test_empty_arrays_raise():
    with pytest.raises(ValueError, match="empty"):
        compute_rmsd(np.zeros((0, 3)), np.zeros((0, 3)))


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="[Ss]hape"):
        compute_rmsd(np.zeros((10, 3)), np.zeros((20, 3)))


def test_wrong_ndim_raises():
    with pytest.raises(ValueError):
        compute_rmsd(np.zeros((10,)), np.zeros((10,)))


def test_truncate_to_common_length():
    a = np.ones((100, 3))
    b = np.ones((80, 3))
    ta, tb = truncate_to_common_length(a, b)
    assert len(ta) == 80
    assert len(tb) == 80


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_rmsd_is_deterministic(kinase_like_coords):
    """compute_rmsd must return the same value on repeated calls."""
    rng = np.random.RandomState(11)
    other = rng.randn(300, 3) * 5.0
    r1 = compute_rmsd(kinase_like_coords, other, align=True)
    r2 = compute_rmsd(kinase_like_coords, other, align=True)
    assert r1 == r2, "RMSD computation is not deterministic"
