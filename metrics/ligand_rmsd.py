import numpy as np
from metrics.rmsd import compute_rmsd, truncate_to_common_length


def compute_ligand_pocket_rmsd(pred_coords: list, gt_coords: list) -> float:
    """
    Compute real RMSD on a proxy ligand pocket definition.

    Since the structure_only mode operates purely on C-alpha traces and
    does not include actual AF3 ligand predictions or HETATM ground truth,
    this function computes the RMSD on a sub-region of the protein as a
    deterministic proxy for the ligand pocket.
    
    Specifically, it identifies the 30 residues closest to the geometric
    center of the ground-truth structure and computes their RMSD.
    """
    pred = np.array(pred_coords, dtype=float)
    gt = np.array(gt_coords, dtype=float)

    if len(pred) != len(gt):
        pred, gt = truncate_to_common_length(pred, gt)

    n = len(pred)
    if n == 0:
        return 0.0

    # Avoid failing on very small proteins
    if n <= 30:
        return compute_rmsd(pred, gt, align=True)

    # Calculate center of mass / geometric center of ground truth
    center = gt.mean(axis=0)

    # find distance to center for each CA
    dists = np.sum((gt - center) ** 2, axis=1)
    
    # get indices of 30 closest residues
    pocket_idx = np.argsort(dists)[:30]

    pocket_pred = pred[pocket_idx]
    pocket_gt = gt[pocket_idx]

    return float(compute_rmsd(pocket_pred, pocket_gt, align=True))
