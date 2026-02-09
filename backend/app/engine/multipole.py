"""
engine.py - Core multipole logic moved from the original single-file HTML/JS app.

This module is intentionally "pure python" (no web framework) so it can be imported
by any backend (FastAPI/Flask/etc).

Data model assumptions (mirrors the original mspgr.json used by the HTML):
- group: dict with keys:
    - "name": str
    - "members": list (optional)
    - "operators": list of operators, each operator is a dict with keys:
        - "R": 3x3 list (integers, crystal/fractional basis)
        - "t": length-3 list (fractions as floats)
        - "tr": optional, time-reversal factor (+1 or -1)
- wyckoff: dict with keys:
    - "label": str
    - "coord": str (e.g. "x,y,z" or "0,0,1/2")
    - "R": 3x3 list
    - "t": length-3 list

If your JSON differs slightly, adapt the field names here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union, cast
import math
import re
import numpy as np

from ..core.alignment_store import get_alignment


# --------------------------
# Basic linear algebra utils
# --------------------------

Vec3 = Tuple[float, float, float]
Mat3 = List[List[float]]

def vec_add(a: Sequence[float], b: Sequence[float]) -> List[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]

def vec_sub(a: Sequence[float], b: Sequence[float]) -> List[float]:
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

def vec_dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def matvec(M: Sequence[Sequence[float]], v: Sequence[float]) -> List[float]:
    return [
        M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
        M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
        M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2],
    ]

def matmul(A: Sequence[Sequence[float]], B: Sequence[Sequence[float]]) -> Mat3:
    out = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for k in range(3):
            aik = A[i][k]
            if aik == 0:
                continue
            for j in range(3):
                out[i][j] += aik * B[k][j]
    return out

def matrix_equal(A: Sequence[Sequence[float]], B: Sequence[Sequence[float]], eps: float = 1e-9) -> bool:
    for i in range(3):
        for j in range(3):
            if abs(A[i][j] - B[i][j]) > eps:
                return False
    return True

# --------------------------
# Tensor invariance solver (NumPy Accelerated)
# --------------------------

@dataclass(frozen=True)
class Op:
    R: Mat3
    tr: int = 1  # time reversal factor (±1). Defaults to +1.

def compute_multipoles(
    rank: int,
    ops: Sequence[Op],
    multipole_mode: str = "magnetic",
) -> List[List[float]]:
    """
    Returns a basis of invariant tensors in the *crystal/fractional* basis using SVD.
    Now NumPy accelerated for efficiency.
    Output: list of basis vectors of length 3**rank (flattened in base-3 index order).
    """
    if rank < 1:
        return []

    dim = 3 ** rank
    
    constraints = []
    
    # 1. Intrinsic symmetry (for rank>=3, enforce symmetry of last rank-1 indices)
    # This matches the legacy behaviour.
    if rank >= 3:
        # We need to enforce T_{...i j ...} = T_{...j i ...} 
        # specifically for the last (rank-1) indices as per original logic.
        # Original logic: "Only add swap-equations when the adjacent pair is out of order" 
        # over the last (rank-1) dimensions.
        
        # We can construct a projector or just direct constraints. 
        # Direct constraints are simpler for SVD.
        
        # We iterate over all indices of the LAST (rank-1) part.
        sub_rank = rank - 1
        sub_dim = 3 ** sub_rank
        
        # We need to link the full tensor index to the sub_rank index
        # Full index u = i0 * sub_dim + flat_sub
        
        # Let's generate the permutations for the sub-tensor
        # Vectorized approach:
        # Create an array of indices [0, ..., sub_dim-1]
        # Convert to digits
        # Check if sorted? No, we need explicit swap constraints.
        
        # For efficiency, we only generate constraints for adjacent swaps.
        # But for full symmetry of size K, we need K-1 adjacent swaps.
        # Here we need full symmetry of the sub-tensor? 
        # The original code checked specific adjacent pairs.
        
        # Let's re-implement the exact logic but using numpy construction
        for k in range(sub_rank - 1):
            # We want to enforce symmetry between index k and k+1 (0-based relative to sub-tensor)
            # These correspond to indices rank-sub_rank+k and rank-sub_rank+k+1 in the full tensor
            # i.e., indices 1+k and 1+k+1 if sub_rank = rank-1.
            
            # We build a permutation map.
            # Grid of indices (3, 3, ..., 3) [rank times]
            # We want P where P applies swap to axes.
            
            # The original code loops over i0 (first index) 
            # and then enforces symmetry on the rest.
            
            # Let's construct a sparse constraint matrix for this.
            # Or simpler:
            # Identity map
            indices = np.arange(dim).reshape([3] * rank)
            
            # Swapped map: swap axes (1+k) and (1+k+1)
            # The axes of `indices` are 0, 1, ..., rank-1.
            # The sub-tensor starts at axis 1.
            # So we swap axis 1+k and 1+k+1.
            
            axis_a = 1 + k
            axis_b = 1 + k + 1
            
            swapped_indices = np.swapaxes(indices, axis_a, axis_b).flatten()
            base_indices = indices.flatten()
            
            # Constraint: v[i] - v[swap(i)] = 0
            # We only need to do this for i < swap(i) to avoid duplicates
            mask = base_indices < swapped_indices
            
            if np.any(mask):
                # Rows: v[idx] - v[swapped_idx] = 0
                rows = base_indices[mask]
                cols_sub = swapped_indices[mask]
                
                # We can't easily append generic rows to a list if we want speed.
                # But building discrete rows for SVD is what we need.
                # Actually, SVD on (P - I) is better? No, P v = v.
                # So (P - I) v = 0.
                
                # Construct linear operator for (P - I)
                # It is a permutation matrix minus identity.
                # We don't need the full matrix.
                # Just the rows where they differ.
                
                # Or simply:
                # For each relevant swap, add P - I to constraints?
                # That is dim x dim. Too big if we stack many.
                # But we only need independent constraints.
                
                # Let's just accumulate the operator (P - I)
                # P_op = np.eye(dim)[swapped_indices] # This is huge (243x243 for rank 5)
                # Actually 243x243 is small for SVD. (243^2 = 59000 floats).
                # Even rank 5 (243 dim) is small.
                # Rank 6 (729) is still okay.
                
                # So we can just construct the full matrix.
                I = np.eye(dim)
                P_mat = np.eye(dim)[swapped_indices]
                constraints.append(P_mat - I)

    # 2. Symmetry constraints: (factor * R_kron - I) v = 0
    I = np.eye(dim)
    for op in ops:
        detR = np.linalg.det(op.R)
        factor = (op.tr * detR) if multipole_mode == "magnetic" else 1.0
        
        # Build Kronecker product R ⊗ ... ⊗ R
        # Use recursive kron or reduce
        M_op = np.array(op.R)
        for _ in range(rank - 1):
            M_op = np.kron(M_op, np.array(op.R))
            
        # Add constraint (factor * M_op - I) v = 0
        C = (factor * M_op) - I
        constraints.append(C)
        
    if not constraints:
        return np.eye(dim).tolist()
        
    # Stack all constraints
    A = np.vstack(constraints)
    
    # Solve A v = 0 using SVD
    # A is (N_constraints * dim) x dim
    U, S, Vh = np.linalg.svd(A)
    
    # Null space corresponds to singular values ~ 0
    tol = 1e-9
    basis = []
    
    # Vh rows are the eigenvectors
    # The last rows correspond to smallest singular values
    for i in range(dim):
        if S[i] < tol:
            basis.append(Vh[i].tolist())
            
    # Sort or normalize?
    # Original logic normalizes.
    # The vectors from SVD are already normalized.
    return basis


def _apply_op_to_tensor(vec: Sequence[float], rank: int, R: Sequence[Sequence[float]]) -> np.ndarray:
    if rank <= 0:
        return np.array(vec, dtype=float)
    kron = np.array(R, dtype=float)
    for _ in range(rank - 1):
        kron = np.kron(kron, R)
    return kron @ np.array(vec, dtype=float)


def _coord_close(a: Sequence[float], b: Sequence[float], tol: float = 1e-4) -> bool:
    return all(abs(x - y) < tol or abs(x - y - 1) < tol or abs(x - y + 1) < tol for x, y in zip(a, b))

def _build_orbit_mapping(
    group: Dict[str, Any],
    wyckoff: Dict[str, Any],
    tol: float = 1e-4,
) -> Tuple[List[List[float]], Dict[int, int]]:
    ops = group.get("operators", [])
    coord_str = wyckoff.get("coord", "")
    base = _parse_coord_triplet(coord_str, 0.123, 0.234, 0.345)
    base = [_frac01(v) for v in base]
    coords: List[List[float]] = []
    op_to_index: Dict[int, int] = {}
    for op_idx, op in enumerate(ops):
        coord = vec_add(matvec(op["R"], base), op.get("t", [0, 0, 0]))
        coord = [_frac01(v) for v in coord]
        found = None
        for i, existing in enumerate(coords):
            if _coord_close(coord, existing, tol):
                found = i
                break
        if found is None:
            coords.append(coord)
            found = len(coords) - 1
        op_to_index[op_idx] = found
    return coords, op_to_index


def _alignment_signs(
    group_index: Optional[int],
    wyckoff_index: Optional[Union[int, str]],
    rank: int,
) -> Dict[int, int]:
    if group_index is None or wyckoff_index is None:
        return {}
    data = get_alignment(group_index, str(wyckoff_index), rank)
    if not data:
        return {}
    signs: Dict[int, int] = {}
    for entry in data.get("alignment", []) or []:
        idx = entry.get("atom_index")
        rel = entry.get("relation_to_ref")
        if idx is None:
            continue
        if rel == "Same":
            signs[int(idx)] = 1
        elif rel == "Opposite":
            signs[int(idx)] = -1
        else:
            signs[int(idx)] = 0
    return signs


def _infer_alignment_signs_from_ops(
    dipole_basis: Sequence[Sequence[float]],
    orbit_ops: Sequence[Dict[str, Any]],
    tol: float = 1e-4,
) -> Dict[int, int]:
    """Infer Same/Opposite alignment signs from dipole basis and orbit ops."""
    if not dipole_basis or not orbit_ops:
        return {}

    ref_vecs = []
    for b in dipole_basis:
        ref_vecs.append(np.array(b, dtype=float))

    signs: Dict[int, int] = {}
    for idx, op in enumerate(orbit_ops):
        R = np.array(op["R"], dtype=float)
        tr = int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1
        detR = float(np.linalg.det(R))
        factor = tr * detR
        same = True
        opp = True
        for bvec in ref_vecs:
            transformed = factor * (R @ bvec)
            if np.any(np.abs(transformed - bvec) > tol):
                same = False
            if np.any(np.abs(transformed + bvec) > tol):
                opp = False
            if not same and not opp:
                break
        if same:
            signs[idx] = 1
        elif opp:
            signs[idx] = -1
        else:
            signs[idx] = 0
    return signs



# --------------------------
# No-SOC magnetic algorithm
# --------------------------

def _compute_no_soc_electric_basis(
    db_all: Optional[Sequence[Dict[str, Any]]],
    group: Dict[str, Any],
    wyckoff: Optional[Dict[str, Any]],
    rank: int,
) -> List[List[float]]:
    if rank <= 1:
        return [[1.0]]

    og_int = _parse_og_integer(group.get("name", ""))
    nuclear_group = _find_nuclear_group(db_all, og_int)
    if not nuclear_group or not wyckoff:
        return []

    nuclear_wyck = _find_matching_wyckoff(nuclear_group, wyckoff)
    if not nuclear_wyck:
        return []

    site_ops = get_site_ops(nuclear_group, nuclear_wyck)
    ops = [Op(R=op["R"], tr=1) for op in site_ops]
    return compute_multipoles(rank=rank - 1, ops=ops, multipole_mode="electric")


def _compute_no_soc_basis(
    db_all: Optional[Sequence[Dict[str, Any]]],
    group: Dict[str, Any],
    wyckoff: Optional[Dict[str, Any]],
    rank: int,
    ops_soc: Sequence[Op],
    group_index: Optional[int] = None,
    wyckoff_index: Optional[Union[int, str]] = None,
    lattice_matrix: Optional[Sequence[Sequence[float]]] = None,
) -> List[List[float]]:
    if wyckoff is None:
        return compute_multipoles(rank=rank, ops=ops_soc, multipole_mode="magnetic")

    dipole_basis = compute_multipoles(rank=1, ops=ops_soc, multipole_mode="magnetic")
    if not dipole_basis:
        return []

    electric_basis = _compute_no_soc_electric_basis(db_all, group, wyckoff, rank)
    if not electric_basis:
        return []

    basis_list: List[List[float]] = []
    for s in dipole_basis:
        spin_vec = np.array(s, dtype=float)
        for e in electric_basis:
            tensor = _tensor_product_spin(e, spin_vec.tolist())
            basis_list.append(tensor.tolist())

    return _combine_basis_vectors(basis_list)


def compute_no_soc_tensor_with_permutations(
    db_all: Optional[Sequence[Dict[str, Any]]],
    group: Dict[str, Any],
    wyckoff: Dict[str, Any],
    rank: int,
    spin_vector: Optional[Sequence[float]] = None,
    group_index: Optional[int] = None,
    wyckoff_index: Optional[Union[int, str]] = None,
    lattice_matrix: Optional[Sequence[Sequence[float]]] = None,
) -> List[List[float]]:
    """Construct no-SOC magnetic tensors by coupling electric (rank n-1) and dipoles."""
    if rank < 1:
        return []

    electric_basis = _compute_no_soc_electric_basis(db_all, group, wyckoff, rank)
    if not electric_basis:
        return []

    if spin_vector is not None:
        dipole_basis = [list(spin_vector)]
    else:
        site_ops_raw = get_site_ops(group, wyckoff)
        ops_soc = [Op(R=op["R"], tr=int(op.get("tr", 1)) or 1) for op in site_ops_raw]
        dipole_basis = compute_multipoles(rank=1, ops=ops_soc, multipole_mode="magnetic")
    if not dipole_basis:
        return []

    coords, op_to_index = _build_orbit_mapping(group, wyckoff)
    n_sites = len(coords)
    if n_sites == 0:
        return []

    op_for_site: List[Dict[str, Any]] = [None] * n_sites  # type: ignore[list-item]
    for op_idx, op in enumerate(group.get("operators", [])):
        site_idx = op_to_index.get(op_idx)
        if site_idx is None:
            continue
        if op_for_site[site_idx] is None:
            op_for_site[site_idx] = op
    for i, op in enumerate(op_for_site):
        if op is None:
            op_for_site[i] = {"R": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "t": [0.0, 0.0, 0.0], "tr": 1}

    alignment_signs = _alignment_signs(group_index, wyckoff_index, 1)
    if not alignment_signs:
        alignment_signs = _infer_alignment_signs_from_ops(dipole_basis, op_for_site)

    if alignment_signs:
        signs = [alignment_signs.get(i, 0) for i in range(n_sites)]
    else:
        signs = [1] * n_sites

    dim_site = 3 ** rank
    sub_rank = rank - 1

    basis_list: List[List[float]] = []
    for s in dipole_basis:
        spin_vec = np.array(s, dtype=float)
        for e in electric_basis:
            site_tensors = []
            for site_idx in range(n_sites):
                sign = signs[site_idx]
                if sign == 0:
                    site_tensors.append(np.zeros(dim_site, dtype=float))
                    continue
                op = op_for_site[site_idx]
                R = op["R"]
                if sub_rank <= 0:
                    e_trans = np.array(e, dtype=float)
                else:
                    e_trans = _apply_op_to_tensor(e, sub_rank, R)
                spin_site = spin_vec * sign
                tensor = _tensor_product_spin(e_trans.tolist(), spin_site.tolist())
                site_tensors.append(tensor)
            stacked = np.concatenate(site_tensors)
            basis_list.append(stacked.tolist())

    return _combine_basis_vectors(basis_list)

def _parse_og_integer(name: str) -> Optional[int]:
    if not name:
        return None
    match = re.search(r"OG:\s*(\d+)", name)
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def _find_nuclear_group(db_all: Optional[Sequence[Dict[str, Any]]], og_int: Optional[int]) -> Optional[Dict[str, Any]]:
    if og_int is None or not db_all:
        return None
    for g in db_all:
        g_og = _parse_og_integer(g.get("name", ""))
        if g_og == og_int:
            return g
    return None

def _normalize_label(label: Optional[str]) -> str:
    return (label or "").strip().lower()

def _normalize_coord(coord: Optional[str]) -> str:
    return (coord or "").replace(" ", "").lower()


def _find_matching_wyckoff(
    nuclear_group: Dict[str, Any],
    wyckoff: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    if not nuclear_group or not wyckoff:
        return None
    candidates = nuclear_group.get("wyckoff", [])
    label = wyckoff.get("label")
    coord = wyckoff.get("coord")
    label_norm = _normalize_label(label)
    coord_norm = _normalize_coord(coord)

    for w in candidates:
        if w.get("label") == label and w.get("coord") == coord:
            return w

    for w in candidates:
        if _normalize_label(w.get("label")) == label_norm and _normalize_coord(w.get("coord")) == coord_norm:
            return w

    for w in candidates:
        if w.get("label") == label or _normalize_label(w.get("label")) == label_norm:
            return w

    return None


def _normalize_frac_vec(v: Sequence[float]) -> List[float]:
    out: List[float] = []
    for x in v:
        y = x - math.floor(x + 1e-6)
        if abs(y - 1.0) < 1e-6:
            y = 0.0
        out.append(y)
    return out

def _coords_match(a: Sequence[float], b: Sequence[float], eps: float = 1e-4) -> bool:
    diff = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    return (
        abs(diff[0] - round(diff[0])) < eps and
        abs(diff[1] - round(diff[1])) < eps and
        abs(diff[2] - round(diff[2])) < eps
    )

def filter_nuclear_group_ops_by_alignment(
    db_all: Optional[Sequence[Dict[str, Any]]],
    group: Dict[str, Any],
    atoms: Sequence[Union[Sequence[float], Dict[str, Any]]],
    alignment_groups: Sequence[Union[Dict[str, Any], Sequence[int]]],
    eps: float = 1e-4,
) -> List[Dict[str, Any]]:
    """
    Apply all symmetry operations of the nuclear group (first group with same
    OG integer) and return only those ops that:
    1) Map atoms within the same alignment group, OR
    2) Induce a bijective permutation of the alignment groups.
    """
    og_int = _parse_og_integer(group.get("name", ""))
    nuclear_group = _find_nuclear_group(db_all, og_int)
    if not nuclear_group:
        return []

    ops = nuclear_group.get("operators", [])
    if not ops or not atoms or not alignment_groups:
        return []

    coords: List[List[float]] = []
    for a in atoms:
        if isinstance(a, dict):
            coord = a.get("coord")
        else:
            coord = a
        if not coord or len(coord) < 3:
            return []
        coords.append([float(coord[0]), float(coord[1]), float(coord[2])])

    groups_list: List[List[int]] = []
    for g in alignment_groups:
        if isinstance(g, dict):
            members = g.get("members") or []
        else:
            members = g
        groups_list.append([int(m) for m in members])

    index_to_group: Dict[int, int] = {}
    for gi, members in enumerate(groups_list):
        for m in members:
            index_to_group[m] = gi

    accepted_ops: List[Dict[str, Any]] = []

    for op in ops:
        mapped_indices: List[int] = []
        valid = True
        for i, coord in enumerate(coords):
            mapped = matvec(op["R"], coord)
            mapped = vec_add(mapped, op.get("t", [0.0, 0.0, 0.0]))
            mapped = _normalize_frac_vec(mapped)

            target_idx = None
            for j, target in enumerate(coords):
                if _coords_match(mapped, target, eps=eps):
                    target_idx = j
                    break
            if target_idx is None:
                valid = False
                break
            mapped_indices.append(target_idx)

        if not valid:
            continue

        # 1) Check same-group preservation
        preserves_groups = True
        for i, j in enumerate(mapped_indices):
            if index_to_group.get(i) != index_to_group.get(j):
                preserves_groups = False
                break
        if preserves_groups:
            accepted_ops.append(op)
            continue

        # 2) Check bijective permutation of groups
        group_map: Dict[int, int] = {}
        bijective = True
        for gi, members in enumerate(groups_list):
            target_group = None
            for m in members:
                mapped_idx = mapped_indices[m]
                tg = index_to_group.get(mapped_idx)
                if tg is None:
                    bijective = False
                    break
                if target_group is None:
                    target_group = tg
                elif target_group != tg:
                    bijective = False
                    break
            if not bijective:
                break
            group_map[gi] = target_group if target_group is not None else -1

        if not bijective:
            continue

        if len(set(group_map.values())) != len(groups_list):
            continue

        accepted_ops.append(op)

    return accepted_ops


def build_group_permutation_matrices(
    ops: Sequence[Dict[str, Any]],
    atoms: Sequence[Union[Sequence[float], Dict[str, Any]]],
    alignment_groups: Sequence[Union[Dict[str, Any], Sequence[int]]],
    eps: float = 1e-4,
) -> List[List[List[int]]]:
    """Return a permutation matrix (n_groups x n_groups) for each accepted op."""
    if not ops or not atoms or not alignment_groups:
        return []

    coords: List[List[float]] = []
    for a in atoms:
        if isinstance(a, dict):
            coord = a.get("coord")
        else:
            coord = a
        if not coord or len(coord) < 3:
            return []
        coords.append([float(coord[0]), float(coord[1]), float(coord[2])])

    groups_list: List[List[int]] = []
    for g in alignment_groups:
        if isinstance(g, dict):
            members = g.get("members") or []
        else:
            members = g
        groups_list.append([int(m) for m in members])

    index_to_group: Dict[int, int] = {}
    for gi, members in enumerate(groups_list):
        for m in members:
            index_to_group[m] = gi

    matrices: List[List[List[int]]] = []

    for op in ops:
        mapped_indices: List[int] = []
        valid = True
        for coord in coords:
            mapped = matvec(op["R"], coord)
            mapped = vec_add(mapped, op.get("t", [0.0, 0.0, 0.0]))
            mapped = _normalize_frac_vec(mapped)

            target_idx = None
            for j, target in enumerate(coords):
                if _coords_match(mapped, target, eps=eps):
                    target_idx = j
                    break
            if target_idx is None:
                valid = False
                break
            mapped_indices.append(target_idx)

        if not valid:
            continue

        group_map: Dict[int, int] = {}
        bijective = True
        for gi, members in enumerate(groups_list):
            target_group = None
            for m in members:
                mapped_idx = mapped_indices[m]
                tg = index_to_group.get(mapped_idx)
                if tg is None:
                    bijective = False
                    break
                if target_group is None:
                    target_group = tg
                elif target_group != tg:
                    bijective = False
                    break
            if not bijective:
                break
            group_map[gi] = target_group if target_group is not None else -1

        if not bijective:
            continue
        if len(set(group_map.values())) != len(groups_list):
            continue

        n = len(groups_list)
        P = [[0 for _ in range(n)] for _ in range(n)]
        for gi in range(n):
            gj = group_map.get(gi)
            if gj is None or gj < 0:
                continue
            P[gj][gi] = 1
        matrices.append(P)

    return matrices



def _tensor_product_spin(vec_spatial: Sequence[float], spin_vec: Sequence[float]) -> np.ndarray:
    """Return the tensor product Q (rank L-1) ⊗ S (rank 1) as a flattened rank-L vector.

    Both inputs are flattened vectors in base-3 ordering. The result is a
    length 3**L vector suitable for stacking per-site.
    """
    a = np.array(vec_spatial, dtype=float)
    b = np.array(spin_vec, dtype=float)
    # Important: order spin ⊗ spatial so the spin index is the first tensor index
    # (matching flattening used elsewhere: full_index = spin_index * sub_dim + spatial_index)
    return np.kron(b, a)


def _combine_basis_vectors(basis_list: Sequence[Sequence[float]], tol: float = 1e-8) -> List[List[float]]:
    if not basis_list:
        return []
    mat = np.array(basis_list, dtype=float)
    if mat.ndim != 2:
        return []
    try:
        _, s, vh = np.linalg.svd(mat, full_matrices=False)
    except np.linalg.LinAlgError:
        return [list(v) for v in basis_list]
    rank = int(np.sum(s > tol))
    if rank <= 0:
        return []
    return [row.tolist() for row in vh[:rank]]


# --------------------------
# Wyckoff helpers (site symmetry)
# --------------------------

def _parse_coord_expr(expr: str, x: float, y: float, z: float) -> float:
    """
    Very small parser for expressions like "1/2-x" used in the JSON.
    This mirrors the JS eval-based behaviour but avoids eval.
    Supported: +, -, *, /, parentheses, x,y,z, numbers, fractions like 1/2.
    """
    s = expr.strip().lower()
    s = s.replace("x", f"({x})").replace("y", f"({y})").replace("z", f"({z})")
    # allow only safe chars
    if not re.fullmatch(r"[0-9\.\+\-\*\/\(\)\s]+", s):
        # fallback: treat as 0
        return 0.0
    try:
        return float(eval(s, {"__builtins__": {}}, {}))
    except Exception:
        return 0.0

def _parse_coord_triplet(coord_str: str, x: float, y: float, z: float) -> List[float]:
    parts = [p.strip() for p in coord_str.split(",")]
    while len(parts) < 3:
        parts.append("0")
    return [_parse_coord_expr(parts[0], x, y, z),
            _parse_coord_expr(parts[1], x, y, z),
            _parse_coord_expr(parts[2], x, y, z)]

def _frac01(v: float) -> float:
    # normalize to [0,1)
    v = v - math.floor(v + 1e-6)
    if abs(v - 1.0) < 1e-6:
        v = 0.0
    return v

def get_site_ops(group: Dict[str, Any], wyckoff: Optional[Dict[str, Any]], eps: float = 1e-4) -> List[Dict[str, Any]]:
    """
    Mirrors the JS selection:
    - if wyckoff is None -> all group operators
    - else keep ops such that op.R * wR == wR and (op.R*wt + op.t - wt) is integer vector
    """
    ops = group.get("operators", [])
    if wyckoff is None:
        return ops[:]

    wR = wyckoff["R"]
    wt = wyckoff["t"]
    site_ops: List[Dict[str, Any]] = []
    for op in ops:
        RwR = matmul(op["R"], wR)
        if not matrix_equal(RwR, wR):
            continue
        transformed_t = vec_add(matvec(op["R"], wt), op.get("t", [0, 0, 0]))
        diff = vec_sub(transformed_t, wt)
        if (abs(diff[0] - round(diff[0])) < eps and
            abs(diff[1] - round(diff[1])) < eps and
            abs(diff[2] - round(diff[2])) < eps):
            site_ops.append(op)
    return site_ops


# --------------------------
# Public "one-shot" API used by backend
# --------------------------

def compute_results(
    group: Dict[str, Any],
    wyckoff: Optional[Dict[str, Any]],
    max_rank: int,
    mode: str = "magnetic",
    include_soc: bool = True,
    db_all: Optional[Sequence[Dict[str, Any]]] = None,
    magnetic_sites: Optional[Sequence[str]] = None,
    lattice_matrix: Optional[Sequence[Sequence[float]]] = None,
    frontend_atom_tensors_by_rank: Optional[Dict[int, List[List[List[float]]]]] = None,
    group_index: Optional[int] = None,
    wyckoff_index: Optional[Union[int, str]] = None,
) -> Dict[str, Any]:
    site_ops_raw = get_site_ops(group, wyckoff)
    ops: List[Op] = []
    for op in site_ops_raw:
        tr = int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1
        ops.append(Op(R=op["R"], tr=tr))

    ranks = []

    magnetic_labels = None
    if magnetic_sites:
        magnetic_labels = {str(site).lower() for site in magnetic_sites}

    lattice = np.eye(3, dtype=float)
    use_cartesian_ops = False
    if lattice_matrix is not None:
        try:
            candidate = np.array(lattice_matrix, dtype=float)
            if candidate.shape == (3, 3):
                det_candidate = float(np.linalg.det(candidate))
                if abs(det_candidate) > 1e-8:
                    lattice = candidate
                    use_cartesian_ops = True
        except Exception:
            use_cartesian_ops = False

    for r in range(1, max_rank + 1):
        if include_soc or mode != "magnetic":
            basis = compute_multipoles(rank=r, ops=ops, multipole_mode=mode)
        elif wyckoff is None:
            if magnetic_labels:
                summed_basis: List[List[float]] = []
                for w_idx, w in enumerate(group.get("wyckoff", [])):
                    label = str(w.get("label", "")).lower()
                    if label not in magnetic_labels:
                        continue
                    site_ops = get_site_ops(group, w)
                    ops_site = [Op(R=op["R"], tr=int(op.get("tr", 1)) or 1) for op in site_ops]
                    basis_site = _compute_no_soc_basis(
                        db_all=db_all,
                        group=group,
                        wyckoff=w,
                        rank=r,
                        ops_soc=ops_site,
                        group_index=group_index,
                        wyckoff_index=w_idx,
                        lattice_matrix=lattice_matrix,
                    )
                    if basis_site:
                        summed_basis.extend(basis_site)
                basis = _combine_basis_vectors(summed_basis)
            else:
                basis = []
        else:
            # Use new permutation-aware no-SOC implementation which builds global
            # stacked site tensors. To preserve the old API (per-site basis vectors
            # of length 3**rank), extract the reference-site block (site 0) from
            # each global basis vector.
            global_basis = compute_no_soc_tensor_with_permutations(
                db_all=db_all,
                group=group,
                wyckoff=wyckoff,
                rank=r,
                spin_vector=None,
                group_index=group_index,
                wyckoff_index=wyckoff_index,
                lattice_matrix=lattice_matrix,
            )
            if not global_basis:
                basis = []
            else:
                # Determine orbit size and per-site dim
                coords, _ = _build_orbit_mapping(group, wyckoff)
                n_sites = len(coords)
                dim_site = 3 ** r
                basis = []
                for gb in global_basis:
                    arr = np.array(gb, dtype=float)
                    if arr.size != n_sites * dim_site:
                        # unexpected shape; skip
                        continue
                    site0 = arr.reshape((n_sites, dim_site))[0]
                    basis.append(site0.tolist())
        ranks.append({
            "rank": r,
            "basis_count": len(basis),
            "basis": basis,  # list[list[float]] length 3**rank each
        })

    return {
        "group_name": group.get("name", ""),
        "wyckoff_label": wyckoff.get("label") if wyckoff else None,
        "site_ops_count": len(site_ops_raw),
        "site_ops": site_ops_raw,  # return raw so frontend can show them
        "ranks": ranks,
    }
