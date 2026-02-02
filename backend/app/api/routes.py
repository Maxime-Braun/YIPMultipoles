from fastapi import APIRouter
from starlette.concurrency import run_in_threadpool
from typing import Any, Dict, Optional

from ..core.db import load_db
from ..domain.schemas import ComputeRequest
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
        db_all=db,
    )
