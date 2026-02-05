from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

# Simple in-memory cache for alignment payloads
# Keyed by (group_index, wyckoff_index, rank)
_ALIGNMENT_STORE: Dict[Tuple[int, str, int], Dict[str, Any]] = {}


def _key(group_index: int, wyckoff_index: str, rank: int) -> Tuple[int, str, int]:
    return (group_index, str(wyckoff_index), rank)


def set_alignment(payload: Dict[str, Any]) -> None:
    group_index = int(payload.get("group_index", -1))
    wyckoff_index = str(payload.get("wyckoff_index", ""))
    rank = int(payload.get("rank", -1))
    _ALIGNMENT_STORE[_key(group_index, wyckoff_index, rank)] = payload


def get_alignment(
    group_index: int, wyckoff_index: str, rank: int
) -> Optional[Dict[str, Any]]:
    return _ALIGNMENT_STORE.get(_key(group_index, wyckoff_index, rank))
