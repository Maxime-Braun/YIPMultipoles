"""Compute no-SOC magnetic quadrupole (rank-2) for BNS 63.464, Wyckoff 8g.

Usage:
  /path/to/python scripts/compute_no_soc_63_464_8g_quadrupole.py
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

TARGET_BNS = "63.464"
TARGET_WYCKOFF = "8g"
TARGET_RANK = 2


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


def main() -> None:
    db = load_db()
    idx = find_group_index(db, TARGET_BNS)
    if idx is None:
        raise SystemExit(f"Could not find BNS {TARGET_BNS} in DB.")

    group = db[idx]
    wyckoff = find_wyckoff(group, TARGET_WYCKOFF)
    if wyckoff is None:
        raise SystemExit(f"Could not find Wyckoff {TARGET_WYCKOFF} for BNS {TARGET_BNS}.")

    result = compute_results(
        group=group,
        wyckoff=wyckoff,
        max_rank=TARGET_RANK,
        mode="magnetic",
        include_soc=False,
        db_all=db,
    )

    rank_info = next((r for r in result.get("ranks", []) if r.get("rank") == TARGET_RANK), None)
    basis = rank_info.get("basis", []) if rank_info else []

    print(f"Group: {result.get('group_name')}")
    print(f"Wyckoff: {TARGET_WYCKOFF}")
    print(f"Rank {TARGET_RANK} (quadrupole) basis count: {len(basis)}")
    if not basis:
        print("  (no allowed components)")
        return

    for basis_idx, vec in enumerate(basis, start=1):
        print(f"\nBasis {basis_idx} tensor:")
        for i in range(0, len(vec), 3):
            row = vec[i:i + 3]
            print("   ", " ".join(f"{v: .6f}" for v in row))


if __name__ == "__main__":
    main()
