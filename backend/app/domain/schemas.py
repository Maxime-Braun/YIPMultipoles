from pydantic import BaseModel, Field
from typing import List, Optional, Union

class ComputeRequest(BaseModel):
    group_index: int = Field(..., ge=0)
    wyckoff_index: Union[int, str] = "whole"  # "whole" or int
    max_rank: int = Field(3, ge=1, le=5)
    mode: str = Field("magnetic", pattern="^(magnetic|electric)$")
    include_soc: bool = True
    magnetic_sites: Optional[List[str]] = None


class AlignmentAtom(BaseModel):
    atom_index: int = Field(..., ge=0)
    relation_to_ref: str
    group_id: Optional[int] = None
    coord: List[float]
    coord_sym: Optional[str] = None
    op_str: Optional[str] = None


class AlignmentGroup(BaseModel):
    group_id: int = Field(..., ge=0)
    members: List[int]
    relations: List[str]


class AlignmentRequest(BaseModel):
    group_index: int = Field(..., ge=0)
    wyckoff_index: Union[int, str]
    rank: int = Field(..., ge=1, le=5)
    any_incomparable: bool = False
    alignment: List[AlignmentAtom]
    groups: Optional[List[AlignmentGroup]] = None
