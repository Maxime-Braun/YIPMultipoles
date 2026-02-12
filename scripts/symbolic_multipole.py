#!/usr/bin/env python3
"""
Computes multipoles using a hybrid numerical-symbolic approach.

This script uses a fast numerical method (SVD) to find the nullspace of 
the symmetry constraints and then uses SymPy to generate a clean, symbolic 
representation of the basis vectors.

This avoids the performance issues of a pure-symbolic nullspace calculation.

Example Usage:
  ./scripts/symbolic_multipole.py --group 166.97 --wyckoff 9e --rank 3 --mode magnetic
"""
from __future__ import annotations

import json
from pathlib import Path
import math
import sys
import re
import argparse

import sympy as sp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / 'backend' / 'data' / 'mspgr.json'

# --- Utility functions from user's script ---

def load_db(path: Path):
    with open(path, 'r', encoding='utf-8') as fh:
        return json.load(fh)

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def matrix_equal(A, B, eps=1e-9):
    for i in range(3):
        for j in range(3):
            if abs(A[i][j] - B[i][j]) > eps:
                return False
    return True

def vec_add(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]

def vec_sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

def find_group_and_wyckoff(db, group_key_contains, wyck_label):
    for g in db:
        name = g.get('name','')
        if group_key_contains in name:
            for w in g.get('wyckoff', []):
                if w.get('label') == wyck_label:
                    return g, w
    return None, None

def get_site_ops(group, wyckoff, eps=1e-6):
    ops = group.get('operators', [])
    if wyckoff is None:
        return ops[:]
    wR = wyckoff.get('R')
    wt = wyckoff.get('t')
    site_ops = []
    for op in ops:
        RwR = matmul(op['R'], wR)
        if not matrix_equal(RwR, wR):
            continue
        transformed_t = vec_add([sum(op['R'][i][j]*wt[j] for j in range(3)) for i in range(3)], op.get('t', [0,0,0]))
        diff = vec_sub(transformed_t, wt)
        if (abs(diff[0]-round(diff[0])) < eps and
            abs(diff[1]-round(diff[1])) < eps and
            abs(diff[2]-round(diff[2])) < eps):
            site_ops.append(op)
    return site_ops

# --- Cartesian Transformation ---

def get_ideal_lattice_matrix(group):
    """
    Replicates the frontend heuristic for getting the transformation matrix
    from the crystal/fractional basis to a standard Cartesian basis.
    """
    name = group.get('name','')
    m = re.search(r'BNS:\s*(\d+)\.', name)
    n = int(m.group(1)) if m else 1
    a = 1.0
    b = 1.0
    c = 1.0
    alpha = 90.0
    beta = 90.0
    gamma = 90.0
    # Hexagonal/Trigonal systems
    if not (n >= 195) and (n >= 168 or n >= 143):
        gamma = 120.0
    
    u = math.radians(alpha)
    d = math.radians(beta)
    f = math.radians(gamma)
    p = math.cos(u)
    g = math.cos(d)
    x = math.cos(f)
    s = math.sin(f)
    mval = math.sqrt(max(0.0, 1 - p*p - g*g - x*x + 2*p*g*x))
    
    # matrix rows that transform fractional to cartesian
    return np.array([
        [a, b*x, c*g], 
        [0.0, b*s, c*(p - g*x)/s if abs(s) > 1e-12 else 0.0], 
        [0.0, 0.0, c*mval/(s if abs(s) > 1e-12 else 1.0)]
    ])

# --- Hybrid Calculation Engine ---

def get_symbolic_basis(rank: int, site_ops: list, multipole_mode: str, lattice_matrix: np.ndarray) -> list[sp.Matrix]:
    """
    Computes a simplified symbolic basis for the invariant tensors using a
    heuristic approach to clean up the basis from a floating-point rref.
    """
    dim = 3 ** rank
    constraints = []
    
    # 1. Intrinsic symmetry constraints (for rank>=3)
    if rank >= 3:
        indices = np.arange(dim).reshape([3] * rank)
        sub_rank = rank - 1
        for k in range(sub_rank - 1):
            axis_a = 1 + k
            axis_b = 1 + k + 1
            swapped_indices = np.swapaxes(indices, axis_a, axis_b)
            P_swap = np.eye(dim)[swapped_indices.flatten()]
            constraints.append(np.eye(dim) - P_swap)

    # 2. Symmetry operation constraints from the group
    I = np.eye(dim)
    for op in site_ops:
        R = np.array(op['R'])
        tr = int(op.get('tr', 1)) if op.get('tr', 1) is not None else 1
        detR = np.linalg.det(R)
        factor = (tr * detR) if multipole_mode == "magnetic" else 1.0
        R_kron = np.array(R)
        for _ in range(rank - 1):
            R_kron = np.kron(R_kron, R)
        C = I - (factor * R_kron)
        constraints.append(C)
        
    if not constraints:
        # If there are no constraints, the basis is the identity matrix
        crystal_basis = np.eye(dim)
    else:
        A = np.vstack(constraints)
        _u, s, vh = np.linalg.svd(A)
        tol = 1e-9
        crystal_basis = vh[s < tol].T
    
    if crystal_basis.shape[1] == 0:
        return []

    # 3. Transform the basis from crystal to Cartesian coordinates
    M_kron = np.array(lattice_matrix)
    for _ in range(rank - 1):
        M_kron = np.kron(M_kron, lattice_matrix)
    
    cartesian_basis = M_kron @ crystal_basis

    # 4. Simplify the Cartesian basis to a clean integer symbolic basis
    # Perform rref on the numerical basis first (fast but imprecise)
    M = sp.Matrix(cartesian_basis.T)
    rref_matrix, pivot_cols = M.rref()

    final_basis = []
    for row_idx in range(rref_matrix.rows):
        vec = rref_matrix.row(row_idx).T
        if not vec.is_zero:
            # Heuristic: normalize by the largest component to reveal ratios
            max_abs_val = max(abs(v) for v in vec)
            
            if max_abs_val < 1e-9:
                continue

            normalized_vec = vec / max_abs_val
            
            # nsimplify the normalized vector with a reasonable tolerance
            simplified_vec = normalized_vec.applyfunc(lambda x: sp.nsimplify(x, tolerance=1e-3, rational=True))
            
            # Clear denominators to get an integer vector
            lcm_denom = sp.lcm([sp.fraction(s)[1] for s in simplified_vec])
            int_vec = simplified_vec * lcm_denom

            # Simplify the integer vector by dividing by the GCD of its components
            g = sp.gcd([abs(n) for n in int_vec if n != 0])
            if g != 0:
                final_basis.append(int_vec / g)
            else:
                final_basis.append(int_vec)

    return final_basis


def analyze_and_print_tensor_constraints(basis_vectors: list[sp.Matrix], rank: int):
    """
    Analyzes the tensor basis to find which components are non-zero and which
    are equal or opposite to each other, then prints a clean summary.
    """
    if not basis_vectors:
        print('No non-trivial multipole allowed (all forbidden)')
        return

    # 1. Build the general symbolic tensor
    coeffs = sp.symbols(f'c0:{len(basis_vectors)}')
    dim = 3**rank
    general_solution = sum((c * v for c, v in zip(coeffs, basis_vectors)), sp.zeros(dim, 1))

    # 2. Identify all non-zero components
    def idx_to_name(idx, r):
        indices = []
        temp_idx = idx
        for _ in range(r):
            indices.insert(0, temp_idx % 3)
            temp_idx //= 3
        return f"T_{''.join(map(str, indices))}"

    non_zero_components = {}
    for i in range(dim):
        expr = sp.simplify(general_solution[i])
        if expr != 0:
            non_zero_components[i] = expr
    
    if not non_zero_components:
        print("All tensor components are zero.")
        return

    # 3. Group components into families based on equality
    print("\nAllowed Non-Zero Components and Their Relationships:")
    print("--------------------------------------------------")
    
    processed_indices = set()
    family_num = 1
    
    # Sort indices to have a consistent output order
    sorted_indices = sorted(non_zero_components.keys())

    for i in sorted_indices:
        if i in processed_indices:
            continue
        
        # Start a new family with the current component
        family_head_expr = non_zero_components[i]
        family_head_name = idx_to_name(i, rank)
        family_members = [family_head_name]
        processed_indices.add(i)

        # Find other components that are equal or opposite
        for j in sorted_indices:
            if i == j or j in processed_indices:
                continue
            
            other_expr = non_zero_components[j]
            
            # Check for equality
            if sp.simplify(family_head_expr - other_expr) == 0:
                family_members.append(f"{idx_to_name(j, rank)}")
                processed_indices.add(j)
            # Check for opposite
            elif sp.simplify(family_head_expr + other_expr) == 0:
                family_members.append(f"-{idx_to_name(j, rank)}")
                processed_indices.add(j)

        print(f"\nFamily {family_num}:")
        # Join the members with " = "
        print(f"  {' = '.join(family_members)}")
        family_num += 1

    # Print any remaining independent components
    independent_components = []
    for i in sorted_indices:
        if i not in processed_indices:
            independent_components.append(idx_to_name(i, rank))

    if independent_components:
        print(f"\nFamily {family_num} (Independent Components):")
        for comp in independent_components:
            print(f"  {comp}")


def main():
    parser = argparse.ArgumentParser(description="Hybrid symbolic multipole calculator.")
    parser.add_argument('--group', required=True, help="Key for the space group, e.g., '166.97'")
    parser.add_argument('--wyckoff', required=True, help="Wyckoff label, e.g., '9e'")
    parser.add_argument('--rank', type=int, required=True, help="Tensor rank")
    parser.add_argument('--mode', choices=['magnetic', 'electric'], default='magnetic', help="Multipole mode")

    args = parser.parse_args()

    db = load_db(DB_PATH)
    group, wyck = find_group_and_wyckoff(db, args.group, args.wyckoff)
    if group is None:
        print(f"Group containing '{args.group}' or Wyckoff '{args.wyckoff}' not found.")
        sys.exit(1)

    print(f"Found group: {group.get('name')}")
    print(f"Wyckoff: {wyck.get('label')} at coord={wyck.get('coord')}")

    site_ops = get_site_ops(group, wyck)
    print(f"Found {len(site_ops)} site symmetry operations.")

    # Get the transformation matrix to the Cartesian basis
    lattice_matrix = get_ideal_lattice_matrix(group)
    print("\nApplying Crystal-to-Cartesian transformation with matrix:")
    print(lattice_matrix)

    # Get the simplified symbolic basis using the new method
    print("\nComputing symbolic basis in Cartesian coordinates...")
    symbolic_basis_vectors = get_symbolic_basis(args.rank, site_ops, args.mode, lattice_matrix)
    
    if not symbolic_basis_vectors:
        print('No non-trivial multipole allowed (all forbidden).')
        return

    print(f"Nullspace dimension is: {len(symbolic_basis_vectors)}")

    # Analyze and print the tensor constraints in a simplified format
    print("\nResulting tensor is in the CARTESIAN basis.")
    analyze_and_print_tensor_constraints(symbolic_basis_vectors, args.rank)


if __name__ == '__main__':
    main()
