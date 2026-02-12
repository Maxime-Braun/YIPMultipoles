"""Compute no-SOC magnetic multipoles for BNS 136.499, Wyckoff 2a.

Usage:
  /path/to/python scripts/compute_no_soc_136_499_2a.py

Optional flags:
  --max-rank 3
  --label 2a
"""
from __future__ import annotations

from typing import Optional

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.app.core.db import load_db
from backend.app.engine.multipole import compute_results

DEFAULT_BNS = "136.499"
DEFAULT_WYCKOFF = "2a"
DEFAULT_MAX_RANK = 5


def find_group_index(db: list[dict], bns: str) -> Optional[int]:
    for idx, g in enumerate(db):
        name = g.get("name", "")
        if bns in name:
            return idx
    return None


def find_wyckoff(group: dict, label: str) -> Optional[dict]:
    target = label.strip().lower()
    for w in group.get("wyckoff", []):
        if str(w.get("label", "")).strip().lower() == target:
            return w
    return None


def parse_args() -> tuple[str, str, int]:
    parser = argparse.ArgumentParser(description="Compute no-SOC magnetic multipoles for a Wyckoff site.")
    parser.add_argument("--bns", default=DEFAULT_BNS, help="BNS number (e.g. 136.499)")
    parser.add_argument("--label", default=DEFAULT_WYCKOFF, help="Wyckoff label (e.g. 2a)")
    parser.add_argument("--max-rank", type=int, default=DEFAULT_MAX_RANK, help="Maximum rank to compute")
    args = parser.parse_args()
    return args.bns, args.label, args.max_rank


def main() -> None:
    bns, label, max_rank = parse_args()
    db = load_db()
    idx = find_group_index(db, bns)
    if idx is None:
        raise SystemExit(f"Could not find BNS {bns} in DB.")

    group = db[idx]
    wyckoff = find_wyckoff(group, label)
    if wyckoff is None:
        raise SystemExit(f"Could not find Wyckoff label {label} for BNS {bns}.")

    result = compute_results(
        group=group,
        wyckoff=wyckoff,
        max_rank=max_rank,
        mode="magnetic",
        include_soc=False,
        db_all=db,
    )

    print(f"Group: {result.get('group_name')}")
    print(f"Wyckoff: {label}")
    for rank_info in result.get("ranks", []):
        rank = rank_info.get("rank")
        basis = rank_info.get("basis", [])
        print(f"\nRank {rank}:")
        if not basis:
            print("  (no allowed components)")
            continue
        if len(basis) != 1:
            print(f"  (basis count: {len(basis)})")
        vec = basis[0]
        for i in range(0, len(vec), 3):
            row = vec[i:i + 3]
            print("   ", " ".join(f"{v: .6f}" for v in row))


if __name__ == "__main__":
    main()
