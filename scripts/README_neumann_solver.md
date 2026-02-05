# Neumann permutation solver

This script builds the Neumann invariance system

\[
M = (P \otimes R^{\otimes n}) M
\]

for all accepted symmetry operations and solves for the nullspace.

## Usage

```bash
python scripts/neumann_permutation_solver.py --bns 166.97 --wyckoff 9e --rank 2
```

## Output

- Orbit size and number of alignment groups
- Number of accepted operations
- Nullspace dimension
- First basis vector (flattened)

## Notes

- The rank is the spatial tensor rank (e.g., rank 2 for quadrupole).
- The script uses the permutation matrices from the alignment grouping logic.
