import numpy as np

def run():
    print("--- Verifying Frontend Logic Y00 Weights for Group 166.97, Wyckoff 9e ---")

    # --- 1. Constants & Data ---
    # Group 166.97 Operators (R-3m)
    ops_data = [
      {"R": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1},
      {"R": [[0, -1, 0], [1, -1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # C3+
      {"R": [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # C3-
      {"R": [[1, -1, 0], [0, -1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # C2'
      {"R": [[0, 1, 0], [1, 0, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # C2'
      {"R": [[-1, 0, 0], [-1, 1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # C2'
      {"R": [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # Inv
      {"R": [[0, 1, 0], [-1, 1, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # S6-
      {"R": [[1, -1, 0], [1, 0, 0], [0, 0, -1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # S6+
      {"R": [[-1, 1, 0], [0, 1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # m
      {"R": [[0, -1, 0], [-1, 0, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}, # m
      {"R": [[1, 0, 0], [1, -1, 0], [0, 0, 1]], "t": [0.0, 0.0, 0.0], "tr": 1}  # m
    ]
    
    # Lattice Centerings for R (Hexagonal axes)
    centerings = [
        np.array([0.0, 0.0, 0.0]),
        np.array([2.0/3.0, 1.0/3.0, 1.0/3.0]),
        np.array([1.0/3.0, 2.0/3.0, 2.0/3.0])
    ]

    # --- 2. Lattice Matrix Construction (Hexagonal) ---
    # a=b=1, c=1 (arbitrary for weight check), gamma=120
    # From getIdealLatticeMatrix logic:
    # M = [[a, b*cg, c*cb], [0, b*sg, ...], ...]
    gamma_deg = 120.0
    cg = np.cos(np.deg2rad(gamma_deg)) # -0.5
    sg = np.sin(np.deg2rad(gamma_deg)) # 0.866
    
    # M maps Fractional -> Cartesian
    # Columns are lattice vectors: a, b, c
    M = np.array([
        [1.0, -0.5, 0.0],
        [0.0,  sg,  0.0],
        [0.0,  0.0, 1.0]
    ])
    
    # Precompute Inverse for Cart -> Frac (if needed)
    M_inv = np.linalg.inv(M)

    # --- 3. Tensor Logic Replicas ---

    def transform_rank1_frac_to_cart(v_frac, M_mat):
        # v_cart = M * v_frac
        return M_mat @ v_frac

    def slice_spin_component(cart_tensor, spin_idx):
        # For rank 1, tensor is [x, y, z]
        return cart_tensor[spin_idx]

    # --- 4. Orbit Generation ---
    p0 = np.array([0.5, 0.0, 0.0]) # 9e
    orbit = []
    seen = []

    for i_op, op in enumerate(ops_data):
        R = np.array(op['R'])
        t = np.array(op['t'])
        tr = op['tr']
        
        # Site Symmetry Check: Does this op map p0 to p0?
        # Not strictly needed for orbit generation, but good for debugging.

        for i_cent, tc in enumerate(centerings):
            # p_new = R p0 + t + tc
            p_new = R @ p0 + t + tc
            p_wrap = p_new % 1.0
            
            # Check uniq
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
                    'op_R': R,
                    'op_tr': tr,
                    'is_inv': np.degrees(np.abs(np.linalg.det(R) + 1.0)) < 1e-3 # Det = -1?
                })
    
    print(f"Orbit Size: {len(orbit)}")

    # --- 5. Weight Calculation ---
    # We want to see the weight of Y00 for a generic Dipole input.
    # Input: Generic Rank 1 Tensor in Crystal Basis: v_gen_frac
    # Let's test 3 orthogonal directions in crystal basis to probe the space.
    
    test_inputs = [
        ("u (1,0,0)", np.array([1.0, 0.0, 0.0])),
        ("v (0,1,0)", np.array([0.0, 1.0, 0.0])),
        ("w (0,0,1)", np.array([0.0, 0.0, 1.0]))
    ]
    
    # Spin Indices: 0=X, 1=Y, 2=Z (Cartesian axes)
    spin_labels = ["Mx (Cart)", "My (Cart)", "Mz (Cart)"]

    print(f"\n{'Input Basis':<12} | {'Spin Comp':<10} | {'Atom 0 (Cart Val)':<20} | {'Y00 Weight':<15}")
    print("-" * 70)

    for label, v_seed in test_inputs:
        for spin_idx in range(3):
            # Calculate for first atom in orbit (Atom 0)
            # T_loc_frac = R * v_seed (Rank 1 transforms by R)
            # Note: For magnetic moments, is it pseudo-vector?
            # 166.97 is Colorless (tr=1).
            # Magnetic moment converts as Axial Vector: v' = det(R) * R * v.
            # Frontend Logic: 'magnetic' mode uses det(R) factor.
            
            # Let's check `getGenerators` logic in frontend.
            # It usually applies the factor det(R) for magnetic mode.
            
            atom0 = orbit[0]
            R0 = atom0['op_R']
            det0 = np.linalg.det(R0)
            
            # MAG MOMENT TRANSFORMATION: v_frac_local = det(R) * R * v_seed * tr
            v_loc_frac = det0 * (R0 @ v_seed) * atom0['op_tr']
            
            # Transform to Cartesian
            v_loc_cart = transform_rank1_frac_to_cart(v_loc_frac, M)
            
            # Slice Spin (Scalar X/Y/Z)
            val_scalar = list(v_loc_cart)[spin_idx]
            
            # Project to Y00
            # Y00 = 1/sqrt(4pi)
            # Integral = val_scalar * sqrt(4pi)
            weight_y00 = val_scalar * np.sqrt(4 * np.pi)
            
            # Output
            print(f"{label:<12} | {spin_labels[spin_idx]:<10} | {val_scalar:8.4f}             | {weight_y00:8.4f}")

    print("\nNote: Y00 Weight = CartesianComponent * sqrt(4pi)")
    print("If site symmetry forbids a component, it should be zero.")

if __name__ == "__main__":
    run()
