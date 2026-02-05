"""Compute system magnetic multipoles for BNS 198.9 with magnetic sites on 12b.

Usage:
  /path/to/python scripts/compute_system_multipoles_198_9.py

This uses the backend engine and the no-SOC system summation logic.
"""
from __future__ import annotations

from typing import Optional

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.app.core.db import load_db
from backend.app.engine.multipole import compute_results, get_site_ops, _compute_no_soc_basis, Op


DEFAULT_BNS = "198.9"
DEFAULT_MAGNETIC_SITES = ["12b"]
DEFAULT_MAX_RANK = 5


def find_group_index(db: list[dict], bns: str) -> Optional[int]:
    for idx, g in enumerate(db):
        name = g.get("name", "")
        if bns in name:
            return idx
    return None


def parse_args() -> tuple[str, list[str], int, bool, int]:
    import argparse

    parser = argparse.ArgumentParser(description="Compute system magnetic multipoles for a BNS group.")
    parser.add_argument("--bns", default=DEFAULT_BNS, help="BNS number (e.g. 198.9)")
    parser.add_argument("--sites", default=",".join(DEFAULT_MAGNETIC_SITES), help="Comma-separated magnetic site labels")
    parser.add_argument("--max-rank", type=int, default=DEFAULT_MAX_RANK, help="Maximum rank to compute")
    parser.add_argument("--debug-sites", action="store_true", help="Print per-site bases before summation")
    parser.add_argument("--debug-rank", type=int, default=3, help="Rank to print when --debug-sites is used")
    args = parser.parse_args()
    sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    return args.bns, sites, args.max_rank, args.debug_sites, args.debug_rank


def format_rank3(vec: list[float]) -> str:
    lines = []
    axes = ["x", "y", "z"]
    for i in range(3):
        lines.append(f"i = {axes[i]}:")
        for j in range(3):
            row = []
            for k in range(3):
                idx = i * 9 + j * 3 + k
                row.append(f"{vec[idx]: .6f}")
            lines.append("| " + "; ".join(row) + " |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    target_bns, magnetic_sites, max_rank, debug_sites, debug_rank = parse_args()
    db = load_db()
    idx = find_group_index(db, target_bns)
    if idx is None:
        raise SystemExit(f"Could not find BNS {target_bns} in DB.")

    group = db[idx]
    result = compute_results(
        group=group,
        wyckoff=None,
        max_rank=max_rank,
        mode="magnetic",
        include_soc=False,
        magnetic_sites=magnetic_sites,
        db_all=db,
    )

    print(f"Group: {result.get('group_name')}")
    print(f"Magnetic sites: {', '.join(magnetic_sites)}")
    for rank_info in result.get("ranks", []):
        rank = rank_info.get("rank")
        basis = rank_info.get("basis", [])
        print(f"\nRank {rank}:")
        if not basis:
            print("  (no allowed components)")
            continue
        if len(basis) != 1:
            print(f"  (unexpected {len(basis)} tensors; showing sum)")
        vec = basis[0]
        print("  tensor:")
        for i in range(0, len(vec), 3):
            row = vec[i:i+3]
            print("   ", " ".join(f"{v: .6f}" for v in row))

    if debug_sites:
        print("\nPer-site debug output:")
        for w in group.get("wyckoff", []):
            label = str(w.get("label", "")).lower()
            if label not in {s.lower() for s in magnetic_sites}:
                continue
            site_ops = get_site_ops(group, w)
            ops_site = [Op(R=op["R"], tr=int(op.get("tr", 1)) or 1) for op in site_ops]
            basis_site = _compute_no_soc_basis(
                db_all=db,
                group=group,
                wyckoff=w,
                rank=debug_rank,
                ops_soc=ops_site,
            )
            print(f"\nSite {w.get('label')}, rank {debug_rank} basis count: {len(basis_site)}")
            for idx_vec, vec in enumerate(basis_site):
                print(f"Basis {idx_vec + 1}:")
                print(format_rank3(vec))


if __name__ == "__main__":
    main()
