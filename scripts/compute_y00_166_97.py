import numpy as np
import sys

# --- PDF Helper ---
try:
    from fpdf import FPDF
    HAS_FPDF = True
except ImportError:
    HAS_FPDF = False
    print("WARNING: fpdf2 not installed. Use 'pdf.log(...)' to print to stdout only.")

class PDFLog:
    def __init__(self, filename="multipole_report.pdf"):
        self.filename = filename
        self.lines = []
        if HAS_FPDF:
            self.pdf = FPDF()
            self.pdf.add_page()
            # Use Helvetica which is standard. Courier is too.
            self.pdf.set_font("Courier", size=7) 
            self.pdf.set_auto_page_break(auto=True, margin=15)
        else:
            self.pdf = None

    def log(self, text=""):
        # Print to stdout
        print(text)
        # Add to PDF buffer
        if self.pdf:
            try:
                # Replace tabs with spaces
                safe_text = str(text).replace('\t', '    ')
                # Simple handling for very first line break
                if safe_text.strip() == "":
                    self.pdf.ln(4)
                    return

                self.pdf.multi_cell(0, 4, safe_text)
                self.pdf.ln(1) # Extra spacing
            except Exception as e:
                print(f"[PDF Error]: {e} on text: {text[:20]}...")

    def save(self):
        if self.pdf:
            print(f"\nSaving PDF report to {self.filename}...")
            self.pdf.output(self.filename)
            print("Done.")

# Global logger
logger = PDFLog("multipole_modes_report.pdf")

def log(text=""):
    logger.log(text)

# --- 1. Math Helpers (Spherical Harmonics & Tensors) ---

def factorial(n):
    res = 1
    for i in range(2, n + 1):
        res *= i
    return res

def assoc_legendre(l, m, x):
    # Compute P_l^m(x) using recursion
    # Simple explicit formulas for low ranks or standard recursion
    # Using scipy would be better but not in requirements.
    # Implementation of standard recurrence relation:
    # (l-m) P_l^m = x(2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    
    # Base cases
    pmm = 1.0
    if m > 0:
        somx2 = np.sqrt((1.0 - x) * (1.0 + x))
        fact = 1.0
        for i in range(1, m + 1):
            pmm *= -fact * somx2
            fact += 2.0
    
    if l == m:
        return pmm
    
    pmmp1 = x * (2.0 * m + 1.0) * pmm
    if l == m + 1:
        return pmmp1
    
    pll = 0.0
    for ll in range(m + 2, l + 1):
        pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
        
    return pll

def get_real_ylm_val(l, m, theta, phi):
    x = np.cos(theta)
    # Condon-Shortley phase is usually included in assoc_legendre defs, ensure consistency
    # My assoc_legendre above includes (-1)^m phase
    
    Plm = assoc_legendre(l, abs(m), x)
    
    norm = np.sqrt(((2.0 * l + 1.0) * factorial(l - abs(m))) / (4.0 * np.pi * factorial(l + abs(m))))
    
    Y = norm * Plm
    
    if m == 0:
        return Y
    elif m > 0:
        return np.sqrt(2.0) * Y * np.cos(m * phi)
    else:
        return np.sqrt(2.0) * Y * np.sin(abs(m) * phi)

def tensor_on_direction(tensor_flat, rank, nx, ny, nz):
    # Contract tensor T_{ijk...} with n_i n_j n_k ...
    # Tensor is flattened 3^rank
    # Ideally should be done efficiently.
    # T * (n (kron) n (kron) n ...)
    
    n = np.array([nx, ny, nz])
    
    # Recursive contraction
    # T(n) = T_{i...} n_i ...
    # View tensor as 3 x 3^(k-1)
    
    curr = tensor_flat
    
    for r in range(rank):
        # Contract last index
        # Reshape curr to (..., 3)
        shape = curr.shape
        prev_dim = shape[0] // 3
        
        # We need to contract the first index or last? Usually conventions vary.
        # Let's assume standard flattening: x, y, z varies fastest at the end? 
        # Actually usually last index varies fastest.
        # Let's contract one by one.
        
        # Reshape to (3, 3, ..., 3)
        t_reshaped = curr.reshape([3] * (rank - r))
        # Contract first dimension
        curr = np.tensordot(n, t_reshaped, axes=([0], [0]))
        
    return curr

def project_tensor_numerical(tensor_flat, rank):
    # Integrate over sphere to get coefficients for Ylm (l up to rank)
    # weights[l,m] = integral( T(Omega) * Ylm(Omega) )
    
    # Grid
    Nth = 20 + 2 * rank
    Nph = 40 + 2 * rank
    
    coeffs = [] # List of (l, m, weight)
    
    dTh = np.pi / Nth
    dPh = 2 * np.pi / Nph
    
    # We will compute weights for l from 0 to rank? 
    # Usually rank-k tensor projects to Y_l with l <= k and same parity.
    # But let's just sweep all l <= k.
    
    ylm_results = {} # (l,m) -> val
    
    for l in range(rank + 1):
        for m in range(-l, l + 1):
            ylm_results[(l,m)] = 0.0

    for i in range(Nth):
        theta = (i + 0.5) * dTh
        sin_t = np.sin(theta)
        w_grad = sin_t * dTh * dPh
        
        for j in range(Nph):
            phi = j * dPh
            nx = sin_t * np.cos(phi)
            ny = sin_t * np.sin(phi)
            nz = np.cos(theta)
            
            val = tensor_on_direction(tensor_flat, rank, nx, ny, nz)
            
            for l in range(rank + 1):
                # Parity check: Rank k tensor has parity (-1)^k (if polar) or (-1)^(k+1) (if axial)?
                # We are projecting geometry. T(n) is a function on sphere.
                # Just compute all.
                for m in range(-l, l + 1):
                    y = get_real_ylm_val(l, m, theta, phi)
                    ylm_results[(l,m)] += val * y * w_grad
                    
    # Format results
    res_list = []
    for l in range(rank + 1):
        for m in range(-l, l + 1):
            w = ylm_results[(l,m)]
            if abs(w) > 1e-4:
                res_list.append((l, m, w))
    return res_list


# --- 2. Symmetry Engines ---

def apply_op_matrix(R, rank):
    # Returns 3^k x 3^k matrix representing R acting on rank-k tensor
    # Flatten order: index 0 is x..x, etc.
    # T'_{i1..ik} = R_{i1 j1} ... R_{ik jk} T_{j1..jk}
    # This is calculating Kronecker product R (x) R (x) ... (x) R
    
    M = R
    for i in range(rank - 1):
        M = np.kron(M, R)
    return M

def transform_tensor(tensor_flat, rank, R, tr, det):
    # T' = det * tr * (R tensor R...) T
    # (Magnetic convention used in script)
    
    op_mat = apply_op_matrix(R, rank)
    # Assuming tr is time reversal (+1 or -1)
    # Assuming magnetic multipoles: transform with det(R) * tr
    # Wait, check frontend convention `applyOpToTensor`.
    # Frontend: 
    #   transformed = applyOpToTensor(baseVec, fullRank, R_frac)
    #   if magnetic: signFactor = tr * detR; localFrac *= signFactor
    
    factor = tr * det
    return factor * (op_mat @ tensor_flat)

def get_allowed_modes(site_ops, rank):
    # Solve system (Op - I) v = 0 for all Ops
    dim = 3**rank
    constraints = []
    
    I = np.eye(dim)
    
    for op in site_ops:
        R = op['R']
        tr = op['tr']
        det = np.linalg.det(R)
        
        M_op = apply_op_matrix(R, rank)
        factor = tr * det
        
        # Constraint: factor * M_op * v = v
        # => (factor * M_op - I) v = 0
        C = (factor * M_op) - I
        constraints.append(C)
        
    if not constraints:
        return np.eye(dim) # All allowed
        
    A = np.vstack(constraints)
    
    # SVD
    # A is (N_ops * dim) x dim
    # This might be big for rank 5 (243). 12 ops * 243 ~ 3000 rows. SVD is fine.
    
    U, S, Vh = np.linalg.svd(A)
    
    tol = 1e-9
    null_rows = []
    for i in range(dim):
        if S[i] < tol:
            null_rows.append(Vh[i])
            
    return null_rows


# --- 3. Main Runner ---

def run():
    log("--- Computing Allowed Magnetic Multipoles (Rank 1..5) for Site 9e / 166.97 ---")

    # --- Constants & Data ---
    ops_data = [
      {"R": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[0, -1, 0], [1, -1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[1, -1, 0], [0, -1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[0, 1, 0], [1, 0, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[-1, 0, 0], [-1, 1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[0, 1, 0], [-1, 1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[1, -1, 0], [1, 0, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[-1, 1, 0], [0, 1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[0, -1, 0], [-1, 0, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[1, 0, 0], [1, -1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}
    ]
    
    centerings = [
        np.array([0.0, 0.0, 0.0]),
        np.array([2.0/3.0, 1.0/3.0, 1.0/3.0]),
        np.array([1.0/3.0, 2.0/3.0, 2.0/3.0])
    ]

    gamma_rad = np.deg2rad(120.0)
    cg = np.cos(gamma_rad)
    sg = np.sin(gamma_rad)
    M = np.array([
        [1.0, -0.5, 0.0],
        [0.0,  sg,  0.0],
        [0.0,  0.0, 1.0]
    ])
    
    p0 = np.array([0.5, 0.0, 0.0])
    
    # 1. Build Orbit
    orbit = []
    seen = []
    
    for op in ops_data:
        R = np.array(op['R'])
        t = np.array(op['t'])
        tr = op.get('tr', 1)
        for tc in centerings:
            p_new = R @ p0 + t + tc
            p_wrap = p_new % 1.0
            
            is_uniq = True
            for s in seen:
                d = p_wrap - s
                d -= np.round(d)
                if np.linalg.norm(d) < 1e-4:
                    is_uniq = False
                    break
            
            if is_uniq:
                seen.append(p_wrap)
                orbit.append({
                    'pos': p_wrap,
                    'gen_R': R,
                    'gen_tr': tr
                })
    orbit.sort(key=lambda x: (round(x['pos'][0],2), round(x['pos'][1],2)))
    
    log(f"Orbit size: {len(orbit)}")
    
    # 2. Site Symmetry
    site_ops = []
    for op in ops_data:
        R = np.array(op['R'])
        t = np.array(op['t'])
        tr = op.get('tr', 1)
        for cent in centerings:
            p_new = R @ p0 + t + cent
            d = p_new - p0
            dist_int = np.linalg.norm(d - np.round(d))
            if dist_int < 1e-4:
                site_ops.append({'R': R, 'tr': tr})
                
    log(f"Site symmetry operations: {len(site_ops)}")
    
    # 3. Iterate Ranks
    for rank in range(1, 6):
        log(f"\n{'='*20} RANK {rank} {'='*20}")
        
        # A. Find Allowed Modes (in Fractional Basis)
        modes = get_allowed_modes(site_ops, rank)
        if len(modes) == 0:
            log("No allowed modes (Forbidden).")
            continue
            
        log(f"Allowed Modes: {len(modes)}")
        
        # B. For each mode, show projection
        for m_idx, mode_vec in enumerate(modes):
            # mode_vec is dim 3^rank fractional tensor
            log(f"\n  [Mode {m_idx + 1}/{len(modes)}]")
            
            # We will propagate to all atoms, but just showing the Generator (Atom 4 usually) first
            # Actually let's just show Generator (p0) and one other (e.g. Atom 0) to save space
            
            # Find generator index in orbit
            gen_idx = -1
            for i, at in enumerate(orbit):
                if np.allclose(at['pos'], p0, atol=1e-3):
                    gen_idx = i
                    break
            
            indices_to_show = [gen_idx, 0] # Show Generator and the first atom (usually rotated)
            
            for atom_idx in indices_to_show:
                if atom_idx < 0 or atom_idx >= len(orbit): continue
                
                atom = orbit[atom_idx]
                is_gen = (atom_idx == gen_idx)
                label = "Generator" if is_gen else f"Atom {atom_idx}"
                
                # Transform (Propagate)
                R = atom['gen_R']
                tr = atom['gen_tr']
                det = np.linalg.det(R)
                
                t_frac = transform_tensor(mode_vec, rank, R, tr, det)
                
                op_M = apply_op_matrix(M, rank)
                t_cart = op_M @ t_frac
                
                # Project to Ylm
                # Normalize?
                norm = np.linalg.norm(t_cart)
                if norm > 1e-9:
                    # t_cart /= norm # Keep scale relative to mode?
                    # Let's keep raw to see propagation signs
                    pass
                
                weights = project_tensor_numerical(t_cart, rank)
                
                # Print
                log(f"    --- {label} ({atom['pos'][0]:.2f}, {atom['pos'][1]:.2f}, {atom['pos'][2]:.2f}) ---")
                
                # Sort weights by l then m
                weights.sort(key=lambda x: (x[0], x[1]))
                
                if not weights:
                    log("      (No projection weights > 1e-4)")
                else:
                    line_items = []
                    for (l, m, w) in weights:
                        line_items.append(f"Y({l},{m}): {w:6.3f}")
                        if len(line_items) >= 4:
                            log("      " + ", ".join(line_items))
                            line_items = []
                    if line_items:
                        log("      " + ", ".join(line_items))
    
    logger.save()

if __name__ == "__main__":
    run()
