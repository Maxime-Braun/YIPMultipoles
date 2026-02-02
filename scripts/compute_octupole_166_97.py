#!/usr/bin/env python3
"""Compute magnetic octupole (rank=3) for Wyckoff position 9e in group 166.97

Uses Neumann principle (T = factor * R⊗R⊗R * T) and SymPy for symbolic linear algebra.

Usage:
  ./scripts/compute_octupole_166_97.py
"""
from __future__ import annotations

import json
from pathlib import Path
import math
import sys
import textwrap

import sympy as sp
import numpy as np
import argparse

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / 'backend' / 'data' / 'mspgr.json'

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

def find_group_and_wyckoff(db, group_key_contains='166.97', wyck_label='9e'):
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

def flatten_index(i,j,k):
    return i*9 + j*3 + k

def build_symmetry_constraints(rank=3):
    # JS behaviour: for rank>=3 enforce symmetry of last (rank-1) indices
    # For rank=3: enforce symmetry j<->k
    eqs = []  # each eq is (u_idx, v_idx) meaning T_u - T_v = 0
    if rank >= 3:
        sub_rank = rank - 1
        sub_dim = 3 ** sub_rank
        for i0 in range(3):
            offset = i0 * sub_dim
            for flat in range(sub_dim):
                indices = []
                tmp = flat
                for _ in range(sub_rank):
                    indices.insert(0, tmp % 3)
                    tmp //= 3
                for k in range(sub_rank - 1):
                    if indices[k] <= indices[k+1]:
                        continue
                    swapped = indices[:]
                    swapped[k], swapped[k+1] = swapped[k+1], swapped[k]
                    flat_swapped = 0
                    for r in range(sub_rank):
                        flat_swapped = flat_swapped * 3 + swapped[r]
                    u = offset + flat
                    v = offset + flat_swapped
                    eqs.append((u, v))
    return eqs

def main():
    db = load_db(DB_PATH)
    group, wyck = find_group_and_wyckoff(db, '166.97', '9e')
    if group is None:
        print('Group 166.97 or Wyckoff 9e not found in DB')
        sys.exit(1)

    print(f"Found group: {group.get('name')}")
    print(f"Wyckoff: {wyck.get('label')} coord={wyck.get('coord')} t={wyck.get('t')}")

    site_ops = get_site_ops(group, wyck)
    print(f"Number of site symmetry ops: {len(site_ops)}")

    rank = 3
    dim = 3 ** rank

    # Build overall constraint matrix A such that A * vec(T) = 0
    rows = []

    I = sp.eye(dim)

    for op in site_ops:
        R = sp.Matrix(op['R'])
        tr = int(op.get('tr', 1)) if op.get('tr', 1) is not None else 1
        detR = int(round(sp.Matrix(op['R']).det()))
        factor = tr * detR  # magnetic mode

        # Kron product R ⊗ R ⊗ R
        K = sp.kronecker_product(sp.kronecker_product(R, R), R)
        C = (I - factor * K)
        # Append rows of C
        for r in range(C.rows):
            rows.append(list(C.row(r)))

    # Add intrinsic symmetry constraints
    for (u, v) in build_symmetry_constraints(rank):
        row = [sp.Integer(0)] * dim
        row[u] = sp.Integer(1)
        row[v] = sp.Integer(-1)
        rows.append(row)

    if len(rows) == 0:
        print('No constraints found; trivial full space')
        return

    A = sp.Matrix(rows)
    print(f'Constraint matrix size: {A.rows} x {A.cols}')

    null = A.nullspace()
    print(f'Number of independent octupole components (nullspace dim): {len(null)}')

    if len(null) == 0:
        print('No non-trivial magnetic octupole allowed (all forbidden)')
        return

    # Pretty-print basis vectors as tensor components
    def idx_to_comp(idx):
        i = idx // 9
        rem = idx % 9
        j = rem // 3
        k = rem % 3
        return f'T_{i}{j}{k}'

    for b_idx, vec in enumerate(null):
        print(f'\nBasis vector {b_idx+1}:')
        # normalize to rational simplified form
        v = sp.Matrix(vec)
        # Scale to smallest integer vector
        lcm_denom = 1
        for entry in v:
            if entry == 0:
                continue
            e = sp.Rational(entry)
            lcm_denom = sp.ilcm(lcm_denom, e.q)
        v_int = (v * lcm_denom).applyfunc(sp.Integer)
        # Simplify gcd
        nums = [abs(int(x)) for x in v_int if x != 0]
        if nums:
            g = nums[0]
            for n in nums[1:]:
                g = math.gcd(g, n)
            if g > 1:
                v_int = (v_int / g).applyfunc(sp.Integer)

        terms = []
        for idx in range(dim):
            val = v_int[idx]
            if val == 0:
                continue
            terms.append(f'{val}*{idx_to_comp(idx)}')
        if not terms:
            print('  (zero vector)')
        else:
            print('  ' + ' + '.join(terms))

    # --- Transform to Cartesian orthonormal basis ---
    def get_ideal_lattice_matrix(group):
        # replicate frontend heuristic for ideal lattice matrix
        name = group.get('name','')
        import re
        m = re.search(r'BNS:\s*(\d+)\.', name)
        n = int(m.group(1)) if m else 1
        a = 1.0
        b = 1.0
        c = 1.0
        alpha = 90.0
        beta = 90.0
        gamma = 90.0
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
        # matrix rows as in frontend: [[a, b*x, c*g],[0, b*s, c*(p-g*x)/s],[0,0,c*mval/s]]
        return [[a, b*x, c*g], [0.0, b*s, c*(p - g*x)/s if abs(s) > 1e-12 else 0.0], [0.0, 0.0, c*mval/(s if abs(s) > 1e-12 else 1.0)]]

    def transform_to_cartesian(vec, rank, M):
        dim = 3 ** rank
        out = [0.0] * dim
        # decode indices for output and input
        def decode(flat, r):
            idxs = []
            tmp = flat
            for _ in range(r):
                idxs.insert(0, tmp % 3)
                tmp //= 3
            return idxs

        for flat_out in range(dim):
            idxs_out = decode(flat_out, rank)
            ssum = 0.0
            for flat_in in range(dim):
                idxs_in = decode(flat_in, rank)
                prod = 1.0
                for d in range(rank):
                    prod *= M[idxs_out[d]][idxs_in[d]]
                # vec entries are sympy objects, convert to float
                ssum += float(vec[flat_in]) * prod
            out[flat_out] = ssum
        return out

    def gram_schmidt(vecs, tol=1e-12):
        orth = []
        for v in vecs:
            w = v[:]
            for u in orth:
                # proj of w on u
                dot = sum(wi * ui for wi, ui in zip(w, u))
                udot = sum(ui * ui for ui in u)
                if abs(udot) > 0:
                    factor = dot / udot
                    w = [wi - factor * ui for wi, ui in zip(w, u)]
            norm = math.sqrt(sum(x*x for x in w))
            if norm > tol:
                orth.append([x / norm for x in w])
        return orth

    M = get_ideal_lattice_matrix(group)
    print('\nIdeal lattice matrix (fractional -> Cartesian):')
    for row in M:
        print('  ', ['{:.6f}'.format(x) for x in row])

    cart_vecs = [transform_to_cartesian(vec, rank, M) for vec in null]
    orthonormal = gram_schmidt(cart_vecs)

    print(f'\nNumber of orthonormal Cartesian basis vectors: {len(orthonormal)}')
    for i, v in enumerate(orthonormal):
        print(f'\nCartesian orthonormal basis #{i+1}: (flattened length {len(v)})')
        # print as 3 slices
        for ia in range(3):
            rows = []
            for jb in range(3):
                row = []
                for kc in range(3):
                    idx = ia*9 + jb*3 + kc
                    row.append('{:+.6f}'.format(v[idx]))
                rows.append('[' + ', '.join(row) + ']')
            print(f'  T_{ia}:')
            for r in rows:
                print('    ', r)

    # --------------------
    # CLI for producing physical tensor
    # --------------------
    parser = argparse.ArgumentParser(description='Compute physical octupole tensor and optionally assemble crystal tensor')
    parser.add_argument('--coeffs', type=str, help='Comma-separated coefficients for Cartesian orthonormal basis or path to JSON/NPY file')
    parser.add_argument('--basis-index', type=int, help='Use a single orthonormal Cartesian basis vector by index (0-based)')
    parser.add_argument('--normalize', choices=['none','max','norm'], default='none', help='Normalize final tensor')
    parser.add_argument('--assemble-global', action='store_true', help='Assemble full crystal tensor by summing site tensors over Wyckoff orbit')
    parser.add_argument('--auto', action='store_true', help='Automatically pick a single physically measurable global tensor (maximizes global norm) when assembling')
    parser.add_argument('--signs', type=str, default=None, help='Comma-separated signs (+1,-1) for each orbit member when assembling global tensor')
    parser.add_argument('--save-json', type=str, default=None, help='Path to save final Cartesian tensor as JSON')
    parser.add_argument('--save-npy', type=str, default=None, help='Path to save final Cartesian tensor as .npy')
    parser.add_argument('--symbolic', action='store_true', help='Solve the symmetry system symbolically and output a symbolic tensor (crystal basis)')
    parser.add_argument('--symbolic-cart', action='store_true', help='Also convert the symbolic tensor to Cartesian coordinates (symbolic)')
    parser.add_argument('--save-symbolic', type=str, default=None, help='Path to save symbolic tensor JSON (entries as strings)')
    args = parser.parse_args()

    # If assembling global and no coefficients provided, enable auto selection
    if args.assemble_global and not args.coeffs and args.basis_index is None:
        args.auto = True

    # Build cart matrix (27 x n_basis) where columns are cart_vecs
    cart_mat = np.column_stack([np.array(c, dtype=float) for c in cart_vecs]) if len(cart_vecs)>0 else np.zeros((27,0))

    # If user requested symbolic solution, produce symbolic tensor and exit
    if args.symbolic:
        m = len(null)
        alphas = sp.symbols('a0:%d' % m)
        T_sym = sp.zeros(dim, 1)
        for k, vec in enumerate(null):
            T_sym += alphas[k] * sp.Matrix(vec)

        def idx_to_trip(idx):
            i = idx // 9
            rem = idx % 9
            j = rem // 3
            k2 = rem % 3
            return (i, j, k2)

        sym_tensor = [[[None for _ in range(3)] for _ in range(3)] for _ in range(3)]
        for idx in range(dim):
            i, j, k2 = idx_to_trip(idx)
            sym_tensor[i][j][k2] = sp.simplify(sp.expand(T_sym[idx, 0]))
        print('\nSymbolic tensor in crystal (fractional) basis:')
        for i in range(3):
            for j in range(3):
                row = []
                for k2 in range(3):
                    row.append(str(sym_tensor[i][j][k2]))
                print(f' T_{i}{j}: {row}')

        # Optionally convert to Cartesian symbolically
        if args.symbolic_cart:
            import re
            name = group.get('name','')
            mres = re.search(r'BNS:\s*(\d+)\.', name)
            n = int(mres.group(1)) if mres else 1
            a = sp.Rational(1)
            b = sp.Rational(1)
            c = sp.Rational(1)
            alpha_deg = sp.Rational(90)
            beta_deg = sp.Rational(90)
            gamma_deg = sp.Rational(90)
            if not (n >= 195) and (n >= 168 or n >= 143):
                gamma_deg = sp.Rational(120)
            u = sp.pi * alpha_deg / 180
            d = sp.pi * beta_deg / 180
            f = sp.pi * gamma_deg / 180
            p = sp.cos(u)
            g = sp.cos(d)
            x = sp.cos(f)
            s = sp.sin(f)
            mval = sp.sqrt(sp.Max(0, 1 - p**2 - g**2 - x**2 + 2*p*g*x))
            M_sym = sp.Matrix([[a, b*x, c*g], [0, b*s, c*(p - g*x)/s], [0, 0, c*mval/s]])

            Tcart_sym = [[[None for _ in range(3)] for _ in range(3)] for _ in range(3)]
            for out_idx in range(dim):
                i, j, k2 = idx_to_trip(out_idx)
                expr = sp.Integer(0)
                for in_idx in range(dim):
                    a_i, b_j, c_k = idx_to_trip(in_idx)
                    coeff = M_sym[i, a_i] * M_sym[j, b_j] * M_sym[k2, c_k]
                    expr += coeff * T_sym[in_idx, 0]
                Tcart_sym[i][j][k2] = sp.simplify(sp.expand(expr))

            print('\nSymbolic tensor in Cartesian basis:')
            for i in range(3):
                for j in range(3):
                    row = []
                    for k2 in range(3):
                        row.append(str(Tcart_sym[i][j][k2]))
                    print(f' T_cart_{i}{j}: {row}')

            if args.save_symbolic:
                out = {'crystal': [[ [str(sym_tensor[i][j][k]) for k in range(3)] for j in range(3)] for i in range(3)],
                       'cartesian': [[[str(Tcart_sym[i][j][k]) for k in range(3)] for j in range(3)] for i in range(3)]}
                with open(args.save_symbolic, 'w', encoding='utf-8') as fh:
                    json.dump(out, fh, indent=2)
                print(f'Saved symbolic tensor to {args.save_symbolic}')
            # (all-ones mode removed)
        else:
            if args.save_symbolic:
                out = {'crystal': [[ [str(sym_tensor[i][j][k]) for k in range(3)] for j in range(3)] for i in range(3)]}
                with open(args.save_symbolic, 'w', encoding='utf-8') as fh:
                    json.dump(out, fh, indent=2)
                print(f'Saved symbolic tensor to {args.save_symbolic}')
        return

    # If user requested automatic global tensor selection, do it now and exit.
    if args.assemble_global and args.auto:
        n_null = len(null)
        S = np.zeros((27, n_null), dtype=float)
        for k in range(n_null):
            vec_cr = np.array([float(x) for x in null[k]])
            G = np.zeros(27, dtype=float)
            seen_local = []
            for op in group.get('operators', []):
                Pp = [sum(op['R'][i][j] * wyck.get('t', [0,0,0])[j] for j in range(3)) + op.get('t', [0,0,0])[i] for i in range(3)]
                unique = True
                for s in seen_local:
                    if all(abs(Pp[i]-s[i]) < 1e-6 for i in range(3)):
                        unique = False
                        break
                if not unique: continue
                seen_local.append(Pp)

                detR = int(round(np.linalg.det(np.array(op['R']))))
                tr = int(op.get('tr', 1)) if op.get('tr', 1) is not None else 1
                factor = tr * detR
                Tprime = np.zeros(27, dtype=float)
                for out_idx in range(27):
                    i = out_idx // 9
                    rem = out_idx % 9
                    j = rem // 3
                    k2 = rem % 3
                    ssum = 0.0
                    for a in range(3):
                        for b in range(3):
                            for c in range(3):
                                in_idx = a*9 + b*3 + c
                                ssum += op['R'][i][a] * op['R'][j][b] * op['R'][k2][c] * vec_cr[in_idx]
                    Tprime[out_idx] = factor * ssum
                Tprime_cart = np.array(transform_to_cartesian(Tprime.tolist(), rank, M), dtype=float)
                G += Tprime_cart
            S[:, k] = G

        U, svals, Vt = np.linalg.svd(S, full_matrices=False)
        alpha = Vt[0, :]
        global_flat = S.dot(alpha)
        global_T = global_flat.reshape((3,3,3))
        if args.normalize == 'max':
            m = np.max(np.abs(global_T))
            if m > 0: global_T = global_T / m
        elif args.normalize == 'norm':
            nrm = np.linalg.norm(global_T)
            if nrm > 0: global_T = global_T / nrm
        if args.save_json:
            out = {'tensor_global': global_T.tolist()}
            with open(args.save_json, 'w', encoding='utf-8') as fh:
                json.dump(out, fh, indent=2)
            print(f'Saved AUTO GLOBAL Cartesian tensor to {args.save_json}')
        if args.save_npy:
            np.save(args.save_npy, global_T)
            print(f'Saved AUTO GLOBAL Cartesian tensor to {args.save_npy}')
        print('\nAUTO-selected alpha (right-singular vector):')
        print(alpha.tolist())
        return

    # Resolve desired cartesian flattened vector v_cart (length 27)
    v_cart = None
    if args.basis_index is not None:
        idx = args.basis_index
        if idx < 0 or idx >= len(orthonormal):
            print('Invalid basis index')
            return
        v_cart = np.array(orthonormal[idx], dtype=float)
    elif args.coeffs:
        s = args.coeffs
        # if path to file
        if s.endswith('.json') or s.endswith('.npy'):
            if s.endswith('.json'):
                with open(s, 'r', encoding='utf-8') as fh:
                    data = json.load(fh)
                coeffs = data.get('coeffs') if isinstance(data, dict) else data
            else:
                coeffs = np.load(s)
            coeffs = list(coeffs)
        else:
            coeffs = [float(x) for x in s.split(',') if x.strip()!='']

        if len(coeffs) != len(orthonormal):
            print(f'Coefficient length {len(coeffs)} does not match basis length {len(orthonormal)}')
            return
        v_cart = np.zeros(27, dtype=float)
        for a, b in zip(coeffs, orthonormal):
            v_cart += a * np.array(b, dtype=float)

    else:
        print('\nNo coefficients or basis-index provided; skipping physical tensor save. Use --coeffs or --basis-index to export a physical tensor.')

    if v_cart is not None or args.assemble_global:
        # normalization
        if args.normalize == 'max':
            m = np.max(np.abs(v_cart))
            if m > 0: v_cart = v_cart / m
        elif args.normalize == 'norm':
            nrm = np.linalg.norm(v_cart)
            if nrm > 0: v_cart = v_cart / nrm

        T_cart = v_cart.reshape((3,3,3))

        if args.save_json:
            out = {'tensor': T_cart.tolist()}
            with open(args.save_json, 'w', encoding='utf-8') as fh:
                json.dump(out, fh, indent=2)
            print(f'Saved Cartesian tensor to {args.save_json}')
        if args.save_npy:
            np.save(args.save_npy, T_cart)
            print(f'Saved Cartesian tensor to {args.save_npy}')

        # Assemble global if requested
        if args.assemble_global:
            # If no coefficients provided, optionally pick an automatic physical alpha
            if (args.coeffs is None and args.basis_index is None) and args.auto:
                # Build S matrix: column k is global Cartesian flattened tensor produced by null vector k
                n_null = len(null)
                S = np.zeros((27, n_null), dtype=float)
                for k in range(n_null):
                    vec_cr = np.array([float(x) for x in null[k]])
                    # sum contributions
                    G = np.zeros(27, dtype=float)
                    seen = []
                    for op in group.get('operators', []):
                        Pp = [sum(op['R'][i][j] * wyck.get('t', [0,0,0])[j] for j in range(3)) + op.get('t', [0,0,0])[i] for i in range(3)]
                        # uniqueness by fractional position
                        unique = True
                        for s in seen:
                            if all(abs(Pp[i]-s[i]) < 1e-6 for i in range(3)):
                                unique = False
                                break
                        if not unique: continue
                        seen.append(Pp)

                        # transform in crystal basis
                        detR = int(round(np.linalg.det(np.array(op['R']))))
                        tr = int(op.get('tr', 1)) if op.get('tr', 1) is not None else 1
                        factor = tr * detR
                        Tprime = np.zeros(27, dtype=float)
                        for out_idx in range(27):
                            i = out_idx // 9
                            rem = out_idx % 9
                            j = rem // 3
                            k2 = rem % 3
                            ssum = 0.0
                            for a in range(3):
                                for b in range(3):
                                    for c in range(3):
                                        in_idx = a*9 + b*3 + c
                                        ssum += op['R'][i][a] * op['R'][j][b] * op['R'][k2][c] * vec_cr[in_idx]
                            Tprime[out_idx] = factor * ssum
                        # convert to Cartesian
                        Tprime_cart = np.array(transform_to_cartesian(Tprime.tolist(), rank, M), dtype=float)
                        G += Tprime_cart
                    S[:, k] = G

                # SVD to get principal combination
                U, svals, Vt = np.linalg.svd(S, full_matrices=False)
                alpha = Vt[0, :]
                # assemble global_v
                global_flat = S.dot(alpha)
                global_T = global_flat.reshape((3,3,3))
                # write outputs if requested
                if args.normalize == 'max':
                    m = np.max(np.abs(global_T))
                    if m > 0: global_T = global_T / m
                elif args.normalize == 'norm':
                    nrm = np.linalg.norm(global_T)
                    if nrm > 0: global_T = global_T / nrm
                if args.save_json:
                    out = {'tensor_global': global_T.tolist()}
                    with open(args.save_json, 'w', encoding='utf-8') as fh:
                        json.dump(out, fh, indent=2)
                    print(f'Saved AUTO GLOBAL Cartesian tensor to {args.save_json}')
                if args.save_npy:
                    np.save(args.save_npy, global_T)
                    print(f'Saved AUTO GLOBAL Cartesian tensor to {args.save_npy}')
                print('\nAUTO-selected alpha (right-singular vector):')
                print(alpha.tolist())
                return
            # Need to express v_cart in terms of cart_vecs to get coefficients alpha in nullspace basis
            if cart_mat.shape[1] == 0:
                print('No basis cart vectors available to assemble global tensor')
                return
            # least squares / pseudo-inverse
            alpha, *_ = np.linalg.lstsq(cart_mat, v_cart, rcond=None)

            # Reconstruct crystal-basis vector
            null_mat = np.column_stack([np.array([float(x) for x in vec]) for vec in null])
            vec_crystal = null_mat.dot(alpha)

            # Create orbit by applying group operators
            def normalize_frac(v):
                return [v[i] - math.floor(v[i] + 1e-9) for i in range(3)]

            positions = []
            tensors = []
            seen = []
            for op in group.get('operators', []):
                # compute position
                Pp = [sum(op['R'][i][j] * wyck.get('t', [0,0,0])[j] for j in range(3)) + op.get('t', [0,0,0])[i] for i in range(3)]
                Pp = normalize_frac(Pp)
                # uniqueness
                unique = True
                for s in seen:
                    if all(abs(Pp[i]-s[i]) < 1e-6 for i in range(3)):
                        unique = False
                        break
                if not unique: continue
                seen.append(Pp)

                # transform crystal vec under op (apply R to indices and factor)
                detR = int(round(np.linalg.det(np.array(op['R']))))
                tr = int(op.get('tr', 1)) if op.get('tr', 1) is not None else 1
                factor = tr * detR
                # apply op to vec_crystal
                Tprime = np.zeros(27, dtype=float)
                for out_idx in range(27):
                    # decode out indices
                    i = out_idx // 9
                    rem = out_idx % 9
                    j = rem // 3
                    k = rem % 3
                    ssum = 0.0
                    for a in range(3):
                        for b in range(3):
                            for c in range(3):
                                in_idx = a*9 + b*3 + c
                                ssum += op['R'][i][a] * op['R'][j][b] * op['R'][k][c] * vec_crystal[in_idx]
                    Tprime[out_idx] = factor * ssum
                # Convert the transformed site tensor (crystal basis) to Cartesian
                Tprime_cart = np.array(transform_to_cartesian(Tprime.tolist(), rank, M), dtype=float)
                positions.append(Pp)
                tensors.append(Tprime_cart.reshape((3,3,3)))

            # apply signs if provided
            signs = None
            if args.signs:
                parts = [float(x) for x in args.signs.split(',') if x.strip()!='']
                if len(parts) != len(tensors):
                    print('Warning: number of signs does not match orbit size; ignoring signs')
                    signs = None
                else:
                    signs = parts

            # Sum to global
            global_T = np.zeros((3,3,3), dtype=float)
            for idx, Tt in enumerate(tensors):
                s = signs[idx] if signs is not None else 1.0
                global_T += s * Tt

            if args.normalize == 'max':
                m = np.max(np.abs(global_T))
                if m > 0: global_T = global_T / m
            elif args.normalize == 'norm':
                nrm = np.linalg.norm(global_T)
                if nrm > 0: global_T = global_T / nrm

            if args.save_json:
                out = {'tensor_global': global_T.tolist(), 'positions': positions}
                with open(args.save_json, 'w', encoding='utf-8') as fh:
                    json.dump(out, fh, indent=2)
                print(f'Saved GLOBAL Cartesian tensor and positions to {args.save_json}')
            if args.save_npy:
                np.save(args.save_npy, global_T)
                print(f'Saved GLOBAL Cartesian tensor to {args.save_npy}')

            print('\nOrbit positions (fractional) and per-site tensor norms:')
            for Pp, Tt in zip(positions, tensors):
                print(f'  pos={Pp} norm={np.linalg.norm(Tt)}')

if __name__ == '__main__':
    main()
