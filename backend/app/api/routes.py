from fastapi import APIRouter
from starlette.concurrency import run_in_threadpool
from typing import Any, Dict, Optional

from ..core.alignment_store import get_alignment, set_alignment
from ..core.db import load_db
from ..domain.schemas import AlignmentRequest, ComputeRequest
from ..engine.multipole import compute_results

router = APIRouter(prefix="/api")

@router.get("/groups")
def groups() -> list[dict]:
    db = load_db()
    return [{"index": i, "name": g.get("name", f"group {i}")} for i, g in enumerate(db)]

@router.get("/wyckoff/{group_index}")
def wyckoff(group_index: int) -> dict:
    db = load_db()
    g = db[group_index]
    wlist = g.get("wyckoff", [])
    return {
        "group_name": g.get("name", ""),
        "wyckoff": [
            {"index": i, "label": w.get("label", ""), "coord": w.get("coord", "")}
            for i, w in enumerate(wlist)
        ],
    }

@router.get("/group/{group_index}")
def group_details(group_index: int):
    db = load_db()
    return db[group_index]

@router.post("/compute")
async def compute(req: ComputeRequest) -> dict:
    db = load_db()
    group = db[req.group_index]

    wyck: Optional[Dict[str, Any]]
    if req.wyckoff_index == "whole":
        wyck = None
    else:
        wyck = group.get("wyckoff", [])[int(req.wyckoff_index)]

    return await run_in_threadpool(
        compute_results,
        group=group,
        wyckoff=wyck,
        max_rank=req.max_rank,
        mode=req.mode,
        include_soc=req.include_soc,
        magnetic_sites=req.magnetic_sites,
        lattice_matrix=req.lattice_matrix,
        frontend_atom_tensors_by_rank=req.frontend_atom_tensors_by_rank,
        group_index=req.group_index,
        wyckoff_index=req.wyckoff_index,
        db_all=db,
    )


@router.post("/alignment")
def save_alignment(payload: AlignmentRequest) -> dict:
    data = payload.dict()

    groups = data.get("groups") or []
    alignment = data.get("alignment") or []

    if groups:
        new_groups = []
        group_id_map = {}
        next_id = 0

        for g in groups:
            members = g.get("members") or []
            relations = g.get("relations") or []
            same_members = []
            same_rel = []
            opp_members = []
            opp_rel = []
            other_members = []
            other_rel = []

            for idx, m in enumerate(members):
                rel = relations[idx] if idx < len(relations) else "Incomparable"
                if rel == "Same":
                    same_members.append(m)
                    same_rel.append(rel)
                elif rel == "Opposite":
                    opp_members.append(m)
                    opp_rel.append(rel)
                else:
                    other_members.append(m)
                    other_rel.append(rel)

            if same_members:
                new_groups.append({"group_id": next_id, "members": same_members, "relations": same_rel})
                group_id_map[(g.get("group_id"), "Same")] = next_id
                next_id += 1
            if opp_members:
                new_groups.append({"group_id": next_id, "members": opp_members, "relations": opp_rel})
                group_id_map[(g.get("group_id"), "Opposite")] = next_id
                next_id += 1
            if other_members:
                new_groups.append({"group_id": next_id, "members": other_members, "relations": other_rel})
                group_id_map[(g.get("group_id"), "Other")] = next_id
                next_id += 1

        if new_groups:
            data["groups"] = new_groups

            # Update alignment entries to new group ids
            for a in alignment:
                rel = a.get("relation_to_ref")
                key = (a.get("group_id"), rel)
                if key not in group_id_map:
                    key = (a.get("group_id"), "Other")
                if key in group_id_map:
                    a["group_id"] = group_id_map[key]
    set_alignment(data)
    return {
        "status": "ok",
        "group_index": payload.group_index,
        "wyckoff_index": payload.wyckoff_index,
        "rank": payload.rank,
        "count": len(payload.alignment),
    }


@router.get("/alignment/{group_index}/{wyckoff_index}/{rank}")
def fetch_alignment(group_index: int, wyckoff_index: str, rank: int) -> dict:
    data = get_alignment(group_index, wyckoff_index, rank)
    return {"status": "ok", "alignment": data}
