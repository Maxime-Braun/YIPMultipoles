# Group permutation analysis script

This script computes alignment groups and group-permutation matrices for a given
BNS number and Wyckoff label.

## Usage

```bash
python scripts/group_permutation_analysis.py --bns 166.97 --wyckoff 9e
```

## Output

The script prints:
- BNS/Wyckoff identifiers
- Orbit size and number of alignment groups
- Accepted nuclear-group symmetry operations
- A permutation matrix (one per operation)

## Troubleshooting

- **"Group not found"**: Confirm the BNS number matches the dataset in `backend/data/mspgr.json`.
- **"Wyckoff not found"**: Verify the label exists for the selected BNS group.
- **"No dipole basis"**: The rank-1 (dipole) basis is forbidden for this site.
