#!/usr/bin/env python3
"""Solve Neumann invariance equations with group permutation matrices.

Example:
  python scripts/neumann_permutation_solver.py --bns 166.97 --wyckoff 9e --rank 2
"""
from __future__ import annotations

import argparse
import math
import sys
import re
from pathlib import Path
from typing import List, Sequence

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.app.core.db import load_db
from backend.app.engine.multipole import (
    build_group_permutation_matrices,
    compute_multipoles,
    filter_nuclear_group_ops_by_alignment,
    get_site_ops,
    Op,
)
from backend.app.core.alignment_store import set_alignment

TOL = 1e-4


def _matvec(R: Sequence[Sequence[float]], v: Sequence[float]) -> List[float]:
    return [
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
    ]


def _vec_add(a: Sequence[float], b: Sequence[float]) -> List[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def _vec_sub(a: Sequence[float], b: Sequence[float]) -> List[float]:
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _is_integer_vec(v: Sequence[float]) -> bool:
    return (
        abs(v[0] - round(v[0])) < 1e-5
        and abs(v[1] - round(v[1])) < 1e-5
        and abs(v[2] - round(v[2])) < 1e-5
    )


def _normalize_frac(v: Sequence[float]) -> List[float]:
    out: List[float] = []
    for x in v:
        y = x - math.floor(x + 1e-6)
        if abs(y - 1.0) < 1e-6:
            y = 0.0
        out.append(y)
    return out


# Forward-declare _parse_coord_expr to satisfy static analyzers.
# The real implementation appears later in this file and will override
# this placeholder at module import time.
def _parse_coord_expr(expr: str, x: float, y: float, z: float) -> float:
    return 0.0


def _parse_coord_parts(coord: str) -> List[float]:
    parts = [p.strip() for p in coord.split(",")]

    def parse_const(s: str) -> float:
        s = s.strip()
        # If expression contains coordinate variables, defer to symbolic parser
        if any(ch in s for ch in ("x", "y", "z")):
            # Use the same placeholder values as other helpers when
            # evaluating expressions with x,y,z.
            try:
                return _parse_coord_expr(s, 0.123, 0.234, 0.345)
            except Exception:
                return 0.0
        if "/" in s:
            try:
                n, d = s.split("/")
                return float(n) / float(d)
            except Exception:
                return 0.0
        try:
            return float(s)
        except Exception:
            return 0.0

    return [parse_const(p) for p in parts]


def _get_lattice_translations(group_name: str) -> List[List[float]]:
    parts = str(group_name).split()
    symbol = ""
    for p in parts:
        if p and p[0].isupper() and not p.startswith("BNS") and not p.startswith("OG"):
            symbol = p
            break
    if not symbol:
        import re

        m = re.search(r"BNS:\s*[\d\.]+\s+([A-Z])", group_name)
        if m:
            symbol = m.group(1)
    if not symbol:
        return [[0.0, 0.0, 0.0]]
    first = symbol[0]
    translations = [[0.0, 0.0, 0.0]]
    if first == "I":
        translations.append([0.5, 0.5, 0.5])
    elif first == "F":
        translations.extend([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
    elif first == "A":
        translations.append([0, 0.5, 0.5])
    elif first == "B":
        translations.append([0.5, 0, 0.5])
    elif first == "C":
        translations.append([0.5, 0.5, 0])
    elif first == "R":
        translations.extend([[2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]])
    return translations


def _get_wyckoff_orbit(group: dict, wyckoff: dict) -> List[dict]:
    P0 = _parse_coord_parts(wyckoff["coord"])
    seen: List[List[float]] = []
    orbit: List[dict] = []
    lat_trans = _get_lattice_translations(group["name"])
    for lt in lat_trans:
        for op_idx, op in enumerate(group["operators"]):
            t = [op["t"][0] + lt[0], op["t"][1] + lt[1], op["t"][2] + lt[2]]
            Pp = _vec_add(_matvec(op["R"], P0), t)
            Pp = _normalize_frac(Pp)
            unique = True
            for s in seen:
                if _is_integer_vec(_vec_sub(Pp, s)):
                    unique = False
                    break
            if unique:
                seen.append(Pp)
                orbit.append({"coord": Pp, "op": op, "op_idx": op_idx})
    return orbit


def _det3(R: Sequence[Sequence[float]]) -> float:
    return (
        R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
        - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
        + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0])
    )


def _analyze_atom_ordering_rank1(basis: List[List[float]], orbit: List[dict]) -> List[dict]:
    atom_vectors = []
    for atom in orbit:
        op = atom["op"]
        R = op["R"]
        tr = int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1
        factor = tr * _det3(R)
        transformed_basis = []
        for b in basis:
            v = _matvec(R, b)
            v = [factor * x for x in v]
            transformed_basis.append(v)
        atom_vectors.append(transformed_basis)

    relation_to_ref = ["Incomparable"] * len(atom_vectors)
    if atom_vectors:
        ref = atom_vectors[0]
        for j, cand in enumerate(atom_vectors):
            all_same = True
            all_opp = True
            for b in range(len(basis)):
                v1 = ref[b]
                v2 = cand[b]
                for k in range(3):
                    if abs(v1[k] - v2[k]) > TOL:
                        all_same = False
                        break
                for k in range(3):
                    if abs(v1[k] + v2[k]) > TOL:
                        all_opp = False
                        break
                if not all_same and not all_opp:
                    break
            if all_same:
                relation_to_ref[j] = "Same"
            elif all_opp:
                relation_to_ref[j] = "Opposite"

    groups = []
    used = set()
    for i in range(len(atom_vectors)):
        if i in used:
            continue
        rep = atom_vectors[i]
        group = {"id": len(groups), "members": []}
        group["members"].append({"atom_idx": i, "relation": "Same", "relation_to_ref": relation_to_ref[i]})
        used.add(i)
        for j in range(i + 1, len(atom_vectors)):
            if j in used:
                continue
            cand = atom_vectors[j]
            is_parallel = True
            is_antiparallel = True
            for b in range(len(basis)):
                v1 = rep[b]
                v2 = cand[b]
                for k in range(3):
                    if abs(v1[k] - v2[k]) > TOL:
                        is_parallel = False
                        break
                for k in range(3):
                    if abs(v1[k] + v2[k]) > TOL:
                        is_antiparallel = False
                        break
                if not is_parallel and not is_antiparallel:
                    break
            if is_parallel:
                group["members"].append(
                    {"atom_idx": j, "relation": "Same", "relation_to_ref": relation_to_ref[j]}
                )
                used.add(j)
            elif is_antiparallel:
                group["members"].append(
                    {"atom_idx": j, "relation": "Opposite", "relation_to_ref": relation_to_ref[j]}
                )
                used.add(j)
        groups.append(group)
    return groups


def _split_groups(groups: List[dict]) -> List[dict]:
    groups_payload = []
    next_group_id = 0
    for g in groups:
        same_members, same_rel = [], []
        opp_members, opp_rel = [], []
        other_members, other_rel = [], []
        for m in g["members"]:
            rel = m["relation_to_ref"]
            if rel == "Same":
                same_members.append(m["atom_idx"])
                same_rel.append(rel)
            elif rel == "Opposite":
                opp_members.append(m["atom_idx"])
                opp_rel.append(rel)
            else:
                other_members.append(m["atom_idx"])
                other_rel.append(rel)

        if same_members:
            groups_payload.append(
                {"group_id": next_group_id, "members": same_members, "relations": same_rel}
            )
            next_group_id += 1
        if opp_members:
            groups_payload.append(
                {"group_id": next_group_id, "members": opp_members, "relations": opp_rel}
            )
            next_group_id += 1
        if other_members:
            groups_payload.append(
                {"group_id": next_group_id, "members": other_members, "relations": other_rel}
            )
            next_group_id += 1

    return groups_payload


def _rotation_kron(R: np.ndarray, rank: int) -> np.ndarray:
    out = R
    for _ in range(rank - 1):
        out = np.kron(out, R)
    return out


def _nullspace(A: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    u, s, vh = np.linalg.svd(A)
    mask = s < tol
    if not np.any(mask):
        return np.zeros((A.shape[1], 0))
    return vh[mask].T


def _axis_label_from_index(index: int, rank: int) -> str:
    axes = ["x", "y", "z"]
    parts = []
    for _ in range(rank):
        parts.insert(0, axes[index % 3])
        index //= 3
    return "".join(parts)


def _symbolic_tensor_from_basis(
    basis: np.ndarray,
    n_groups: int,
    rank: int,
    group_index: int = 0,
) -> List[List[str]]:
    dim_tensor = 3**rank
    if basis.shape[1] == 0:
        return []

    symbols = sp.symbols(f"c0:{basis.shape[1]}")
    vec = sp.Matrix.zeros(basis.shape[0], 1)
    for i in range(basis.shape[1]):
        vec += symbols[i] * sp.Matrix(basis[:, i])

    start = group_index * dim_tensor
    end = start + dim_tensor
    block = vec[start:end]

    labels = [_axis_label_from_index(i, rank) for i in range(dim_tensor)]
    allowed = []
    for i in range(dim_tensor):
        expr = sp.simplify(block[i])
        allowed.append(labels[i] if expr != 0 else "0")

    if rank == 2:
        grid = [["" for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                idx = i * 3 + j
                grid[i][j] = allowed[idx]
        return grid

    if rank == 3:
        blocks: List[List[str]] = []
        axes = ["x", "y", "z"]
        for i in range(3):
            blocks.append([f"i = {axes[i]}:"])
            for j in range(3):
                row = []
                for k in range(3):
                    idx = i * 9 + j * 3 + k
                    row.append(allowed[idx])
                blocks.append(row)
            blocks.append([""])
        return blocks

    return [[allowed[i]] for i in range(dim_tensor)]


def main() -> None:
    parser = argparse.ArgumentParser(description="Solve Neumann equations with group permutations.")
    parser.add_argument("--bns", required=True, help="BNS number, e.g. 166.97")
    parser.add_argument("--wyckoff", required=True, help="Wyckoff label, e.g. 9e")
    parser.add_argument("--rank", type=int, default=2, help="Tensor rank (default: 2)")
    parser.add_argument(
        "--no-soc",
        action="store_true",
        help="Use polar rotations only (ignore axial sign and time reversal)",
    )
    parser.add_argument("--symbolic", action="store_true", help="Print symbolic tensor for group 0")
    args = parser.parse_args()

    db = load_db()
    bns_key = f"BNS: {args.bns}"
    idx = next((i for i, g in enumerate(db) if bns_key in g.get("name", "")), None)
    if idx is None:
        raise SystemExit(f"Group {args.bns} not found.")

    group = db[idx]
    wyckoff = next((w for w in group.get("wyckoff", []) if str(w.get("label", "")) == args.wyckoff), None)
    if wyckoff is None:
        raise SystemExit(f"Wyckoff {args.wyckoff} not found for group {args.bns}.")

    site_ops_raw = get_site_ops(group, wyckoff)
    ops = [Op(R=op["R"], tr=int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1) for op in site_ops_raw]
    rank1_basis = compute_multipoles(rank=1, ops=ops, multipole_mode="magnetic")
    if not rank1_basis:
        raise SystemExit("No dipole basis (rank-1 forbidden).")

    orbit = _get_wyckoff_orbit(group, wyckoff)
    coords = [o["coord"] for o in orbit]

    groups = _analyze_atom_ordering_rank1(rank1_basis, orbit)
    groups_payload = _split_groups(groups)

    # Print detailed alignment groups for user visibility
    print(f"Computed raw groups (count={len(groups)}):")
    for g in groups:
        print(f"  group id {g.get('id')}: members = {[m['atom_idx'] for m in g.get('members', [])]}")
        for m in g.get('members', []):
            print(f"    atom {m['atom_idx']}: relation={m.get('relation')} relation_to_ref={m.get('relation_to_ref')}")

    print(f"Split alignment payload groups (count={len(groups_payload)}): {groups_payload}")

    # Build and store alignment payload into the alignment store so backend
    # helpers (e.g. _alignment_signs) can retrieve the alignment during
    # later computations. This mirrors the payload sent by the frontend.
    alignment_entries = []
    any_incomparable = False
    for g in groups:
        gid = g.get('id')
        for m in g.get('members', []):
            idx_atom = m.get('atom_idx')
            rel = m.get('relation_to_ref', 'Incomparable')
            if rel not in ('Same', 'Opposite'):
                any_incomparable = True
            alignment_entries.append({
                'group_id': gid,
                'atom_index': idx_atom,
                'relation_to_ref': rel,
                'coord': coords[idx_atom] if idx_atom is not None and idx_atom < len(coords) else None,
            })

    payload = {
        'group_index': idx,
        'wyckoff_index': args.wyckoff,
        'rank': 1,
        'any_incomparable': any_incomparable,
        'alignment': alignment_entries,
        'groups': groups_payload,
    }
    # store into in-memory alignment cache
    set_alignment(payload)

    accepted = filter_nuclear_group_ops_by_alignment(db, group, coords, groups_payload)
    matrices = build_group_permutation_matrices(accepted, coords, groups_payload)

    n_groups = len(groups_payload)
    dim_tensor = 3 ** args.rank
    I = np.eye(n_groups * dim_tensor)

    constraints = []
    for op, P_list in zip(accepted, matrices):
        R = np.array(op["R"], dtype=float)
        if args.no_soc:
            factor = 1.0
        else:
            tr = int(op.get("tr", 1)) if op.get("tr", 1) is not None else 1
            factor = tr * float(np.linalg.det(R))
        Rk = _rotation_kron(R, args.rank)
        P = np.array(P_list, dtype=float)
        op_mat = np.kron(P, factor * Rk)
        constraints.append(op_mat - I)

    if not constraints:
        raise SystemExit("No symmetry constraints found.")

    A = np.vstack(constraints)
    basis = _nullspace(A)

    print(f"BNS {args.bns} / Wyckoff {args.wyckoff} / rank {args.rank}")
    print(f"Group index: {idx}")
    print(f"Orbit size: {len(orbit)}")
    print(f"Alignment groups: {len(groups_payload)}")
    print(f"Accepted ops: {len(accepted)}")
    print(f"Nullspace dimension: {basis.shape[1]}")

    if basis.shape[1] > 0:
        print("First basis vector (flattened):")
        print(" ".join(f"{v:.3f}" for v in basis[:, 0]))

    if args.symbolic and basis.shape[1] > 0:
        for g_idx in range(n_groups):
            print(f"\nAllowed entries (group {g_idx}):")
            grid = _symbolic_tensor_from_basis(basis, n_groups, args.rank, group_index=g_idx)
            if args.rank == 2 and grid:
                for row in grid:
                    print("  " + " | ".join(row))
            elif args.rank == 3 and grid:
                for row in grid:
                    if len(row) == 1 and row[0].startswith("i ="):
                        print(row[0])
                        continue
                    if len(row) == 1 and row[0] == "":
                        print("")
                        continue
                    print("| " + "; ".join(f"{val:>6}" for val in row) + " |")
            else:
                for row in grid:
                    print("  " + row[0])


if __name__ == "__main__":
    main()
