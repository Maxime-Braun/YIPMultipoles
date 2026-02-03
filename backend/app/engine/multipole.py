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
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import math
import re
import numpy as np


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
    # Note: The logic below creates constraint rows for permutations.
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



# --------------------------
# No-SOC magnetic algorithm
# --------------------------

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

def _find_matching_wyckoff(nuclear_group: Dict[str, Any], wyckoff: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    if not nuclear_group or not wyckoff:
        return None
    candidates = nuclear_group.get("wyckoff", [])
    label = wyckoff.get("label")
    coord = wyckoff.get("coord")
    label_norm = _normalize_label(label)
    coord_norm = _normalize_coord(coord)

    # Prefer exact label+coord match.
    for w in candidates:
        if w.get("label") == label and w.get("coord") == coord:
            return w

    # Try normalized match (spaces/case).
    for w in candidates:
        if _normalize_label(w.get("label")) == label_norm and _normalize_coord(w.get("coord")) == coord_norm:
            return w

    # Fallback: match by label only.
    for w in candidates:
        if w.get("label") == label or _normalize_label(w.get("label")) == label_norm:
            return w

    return None

def _tensor_product_spin(electric_basis: List[List[float]], spin_vec: Sequence[float]) -> List[List[float]]:
    if not electric_basis:
        return []
    block = len(electric_basis[0])
    out: List[List[float]] = []
    for b in electric_basis:
        vec = [0.0] * (3 * block)
        for s in range(3):
            if abs(spin_vec[s]) < 1e-12:
                continue
            offset = s * block
            for i, val in enumerate(b):
                vec[offset + i] = spin_vec[s] * val
        out.append(vec)
    return out

def _normalize_vec3(v: Sequence[float]) -> Optional[Tuple[float, float, float]]:
    if not v or len(v) < 3:
        return None
    norm = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    if norm < 1e-9:
        return None
    return (v[0] / norm, v[1] / norm, v[2] / norm)

def _soc_dipole_vector(basis_soc_rank1: List[List[float]]) -> Optional[Tuple[float, float, float]]:
    if not basis_soc_rank1:
        return None
    if len(basis_soc_rank1) == 1:
        return _normalize_vec3(basis_soc_rank1[0])

    # Build a representative dipole direction from the SOC basis.
    sx = sum(v[0] for v in basis_soc_rank1)
    sy = sum(v[1] for v in basis_soc_rank1)
    sz = sum(v[2] for v in basis_soc_rank1)
    return _normalize_vec3((sx, sy, sz))

def _compute_no_soc_basis(
    db_all: Optional[Sequence[Dict[str, Any]]],
    group: Dict[str, Any],
    wyckoff: Optional[Dict[str, Any]],
    rank: int,
    ops_soc: Sequence[Op],
) -> List[List[float]]:
    if wyckoff is None:
        return compute_multipoles(rank=rank, ops=ops_soc, multipole_mode="magnetic")

    og_int = _parse_og_integer(group.get("name", ""))
    nuclear_group = _find_nuclear_group(db_all, og_int)
    nuclear_wyck = _find_matching_wyckoff(nuclear_group, wyckoff) if nuclear_group else None
    if not nuclear_wyck:
        return compute_multipoles(rank=rank, ops=ops_soc, multipole_mode="magnetic")

    if rank <= 1:
        electric_basis = [[1.0]]
    else:
        if not nuclear_group:
            return compute_multipoles(rank=rank, ops=ops_soc, multipole_mode="magnetic")
        site_ops_nuclear = get_site_ops(nuclear_group, nuclear_wyck)
        ops_e = [Op(R=op["R"], tr=int(op.get("tr", 1)) or 1) for op in site_ops_nuclear]
        electric_basis = compute_multipoles(rank=rank - 1, ops=ops_e, multipole_mode="electric")

    if not electric_basis:
        return compute_multipoles(rank=rank, ops=ops_soc, multipole_mode="magnetic")

    spin_vec = _soc_dipole_vector(compute_multipoles(rank=1, ops=ops_soc, multipole_mode="magnetic"))
    return [] if spin_vec is None else _tensor_product_spin(electric_basis, spin_vec)


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
) -> Dict[str, Any]:
    site_ops_raw = get_site_ops(group, wyckoff)
    ops: List[Op] = []
    for op in site_ops_raw:
        tr = int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1
        ops.append(Op(R=op["R"], tr=tr))

    ranks = []
    force_forbidden_all = False

    if not include_soc and mode == "magnetic":
        soc_dipole = _soc_dipole_vector(
            compute_multipoles(rank=1, ops=ops, multipole_mode="magnetic")
        )
        if soc_dipole is None:
            force_forbidden_all = True
    for r in range(1, max_rank + 1):
        if force_forbidden_all:
            basis = []
        elif include_soc or mode != "magnetic":
            basis = compute_multipoles(rank=r, ops=ops, multipole_mode=mode)
        else:
            basis = _compute_no_soc_basis(
                db_all=db_all,
                group=group,
                wyckoff=wyckoff,
                rank=r,
                ops_soc=ops,
            )
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
