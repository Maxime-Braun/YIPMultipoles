"""Compute no-SOC magnetic multipoles for BNS 136.499, Wyckoff 2a using the new SOC→no-SOC filter.

Usage:
  /path/to/python scripts/compute_no_soc_136_499_2a_new.py
"""
from __future__ import annotations

from typing import Optional

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.app.core.db import load_db
from backend.app.engine.multipole import compute_results

TARGET_BNS = "136.499"
TARGET_WYCKOFF = "2a"
MAX_RANK = 5


def find_group_index(db: list[dict], bns: str) -> Optional[int]:
    for idx, g in enumerate(db):
        name = g.get("name", "")
        if bns in name:
            return idx
    return None


def find_wyckoff_index(group: dict, label: str) -> Optional[int]:
    target = label.strip().lower()
    for idx, w in enumerate(group.get("wyckoff", [])):
        if str(w.get("label", "")).strip().lower() == target:
            return idx
    return None


def main() -> None:
    db = load_db()
    idx = find_group_index(db, TARGET_BNS)
    if idx is None:
        raise SystemExit(f"Could not find BNS {TARGET_BNS} in DB.")

    group = db[idx]
    wyckoff_index = find_wyckoff_index(group, TARGET_WYCKOFF)
    if wyckoff_index is None:
        raise SystemExit(f"Could not find Wyckoff {TARGET_WYCKOFF} for BNS {TARGET_BNS}.")

    wyckoff = group.get("wyckoff", [])[wyckoff_index]

    result = compute_results(
        group=group,
        wyckoff=wyckoff,
        max_rank=MAX_RANK,
        mode="magnetic",
        include_soc=False,
        db_all=db,
        group_index=idx,
        wyckoff_index=wyckoff_index,
    )

    print(f"Group: {result.get('group_name')}")
    print(f"Wyckoff: {TARGET_WYCKOFF}")
    for rank_info in result.get("ranks", []):
        rank = rank_info.get("rank")
        basis = rank_info.get("basis", [])
        print(f"\nRank {rank}:")
        if not basis:
            print("  (no allowed components)")
            continue
        for b_idx, vec in enumerate(basis, start=1):
            print(f"  Basis {b_idx}:")
            for i in range(0, len(vec), 3):
                row = vec[i:i + 3]
                print("   ", " ".join(f"{v: .6f}" for v in row))


if __name__ == "__main__":
    main()
