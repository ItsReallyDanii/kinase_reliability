"""
Validation controls for contact map computation.

Tests confirm that metrics/contact_map.py is deterministic and behaves
correctly before any real structure data is loaded.
"""

import numpy as np
import pytest

from metrics.contact_map import compute_contact_map, compute_contact_map_overlap


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def linear_chain():
    """50 residues equally spaced 3.8 Å apart along the X axis (idealised beta-strand)."""
    n = 50
    coords = np.zeros((n, 3))
    coords[:, 0] = np.arange(n) * 3.8  # 3.8 Å C-alpha spacing
    return coords


@pytest.fixture
def compact_glob():
    """100 residues packed tightly: many contacts guaranteed."""
    rng = np.random.RandomState(0)
    return rng.randn(100, 3) * 3.0  # small spread → many pairs within 8 Å


# ---------------------------------------------------------------------------
# Basic properties
# ---------------------------------------------------------------------------


def test_contact_map_shape(linear_chain):
    cmap = compute_contact_map(linear_chain)
    assert cmap.shape == (50, 50), f"Expected (50,50), got {cmap.shape}"


def test_contact_map_binary(compact_glob):
    cmap = compute_contact_map(compact_glob)
    assert set(np.unique(cmap)).issubset({0, 1}), "Contact map should be binary (0/1)"


def test_contact_map_symmetric(compact_glob):
    cmap = compute_contact_map(compact_glob)
    np.testing.assert_array_equal(cmap, cmap.T, err_msg="Contact map must be symmetric")


def test_contact_map_diagonal_zero(compact_glob):
    cmap = compute_contact_map(compact_glob)
    assert np.all(np.diag(cmap) == 0), "Diagonal (self-contacts) must be zero"


# ---------------------------------------------------------------------------
# Threshold behaviour
# ---------------------------------------------------------------------------


def test_wider_threshold_more_contacts(linear_chain):
    """Increasing the threshold should increase or maintain the number of contacts."""
    cmap_8 = compute_contact_map(linear_chain, threshold=8.0)
    cmap_12 = compute_contact_map(linear_chain, threshold=12.0)
    assert cmap_12.sum() >= cmap_8.sum(), (
        "Wider threshold should produce ≥ contacts"
    )


def test_zero_threshold_no_contacts(linear_chain):
    """With threshold=0, nothing is in contact."""
    cmap = compute_contact_map(linear_chain, threshold=0.0)
    assert cmap.sum() == 0, "With threshold=0 there should be no contacts"


def test_linear_chain_near_contacts_only(linear_chain):
    """In a 3.8Å-spaced linear chain, only nearest/next-nearest should contact at 8Å."""
    cmap = compute_contact_map(linear_chain, threshold=8.0, min_seq_separation=2)
    # Check that residues > 2 apart along chain are NOT in contact
    # (spacing 3.8Å × k; at k=3 distance = 11.4 Å > 8 Å)
    n = len(linear_chain)
    for i in range(n - 3):
        assert cmap[i, i + 3] == 0, (
            f"Residues {i} and {i+3} at distance {3 * 3.8:.1f}Å should not contact at 8Å"
        )


# ---------------------------------------------------------------------------
# Overlap (Jaccard) controls
# ---------------------------------------------------------------------------


def test_overlap_identity(compact_glob):
    """Overlap of a map with itself should be 1.0."""
    cmap = compute_contact_map(compact_glob)
    overlap = compute_contact_map_overlap(cmap, cmap)
    assert overlap == pytest.approx(1.0, abs=1e-9), (
        f"Identity overlap should be 1.0, got {overlap}"
    )


def test_overlap_empty_maps():
    """Overlap of two all-zero maps should be 0.0 (not NaN or error)."""
    empty = np.zeros((20, 20), dtype=np.uint8)
    overlap = compute_contact_map_overlap(empty, empty)
    assert overlap == 0.0


def test_overlap_disjoint_maps():
    """Two non-overlapping contact maps should have overlap 0.0."""
    map_a = np.zeros((10, 10), dtype=np.uint8)
    map_b = np.zeros((10, 10), dtype=np.uint8)
    map_a[0, 5] = map_a[5, 0] = 1
    map_b[1, 6] = map_b[6, 1] = 1
    overlap = compute_contact_map_overlap(map_a, map_b)
    assert overlap == pytest.approx(0.0, abs=1e-9)


def test_overlap_partial():
    """Overlap of maps sharing exactly half the contacts should be 1/3 (Jaccard)."""
    map_a = np.zeros((10, 10), dtype=np.uint8)
    map_b = np.zeros((10, 10), dtype=np.uint8)
    # 2 contacts in A, 2 in B, 1 shared => Jaccard = 1/3
    map_a[0, 5] = map_a[5, 0] = 1  # contact #1 (shared)
    map_a[1, 7] = map_a[7, 1] = 1  # contact #2 (A only)
    map_b[0, 5] = map_b[5, 0] = 1  # contact #1 (shared)
    map_b[2, 8] = map_b[8, 2] = 1  # contact #3 (B only)
    # intersection = 2 cells (symmetric), union = 6 cells → Jaccard = 2/6 = 1/3
    overlap = compute_contact_map_overlap(map_a, map_b)
    assert overlap == pytest.approx(1.0 / 3.0, abs=1e-9)


def test_overlap_shape_mismatch_raises():
    map_a = np.zeros((10, 10), dtype=np.uint8)
    map_b = np.zeros((8, 8), dtype=np.uint8)
    with pytest.raises(ValueError, match="shape"):
        compute_contact_map_overlap(map_a, map_b)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_contact_map_is_deterministic(compact_glob):
    """compute_contact_map must return identical results on repeated calls."""
    cmap1 = compute_contact_map(compact_glob)
    cmap2 = compute_contact_map(compact_glob)
    np.testing.assert_array_equal(cmap1, cmap2, err_msg="Contact map is not deterministic")


# ---------------------------------------------------------------------------
# Coherence between RMSD and contact map across controls
# ---------------------------------------------------------------------------


def test_contact_overlap_coherent_with_perturbation():
    """
    After a small perturbation, contact map overlap should remain high.
    After a large perturbation, it should drop noticeably.
    """
    rng = np.random.RandomState(5)
    coords = rng.randn(100, 3) * 4.0

    small_noise = coords + rng.randn(100, 3) * 0.05
    large_noise = coords + rng.randn(100, 3) * 10.0

    cmap_orig = compute_contact_map(coords)
    cmap_small = compute_contact_map(small_noise)
    cmap_large = compute_contact_map(large_noise)

    overlap_small = compute_contact_map_overlap(cmap_orig, cmap_small)
    overlap_large = compute_contact_map_overlap(cmap_orig, cmap_large)

    assert overlap_small > overlap_large, (
        f"Small perturbation overlap ({overlap_small:.3f}) should be > "
        f"large perturbation overlap ({overlap_large:.3f})"
    )
    assert overlap_small > 0.5, (
        f"Small perturbation (0.05Å) should preserve most contacts; "
        f"got overlap={overlap_small:.3f}"
    )


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_empty_coords_raises():
    with pytest.raises(ValueError, match="empty"):
        compute_contact_map(np.zeros((0, 3)))


def test_wrong_shape_raises():
    with pytest.raises(ValueError):
        compute_contact_map(np.zeros((10, 2)))
