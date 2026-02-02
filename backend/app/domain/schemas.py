from pydantic import BaseModel, Field
from typing import Union

class ComputeRequest(BaseModel):
    group_index: int = Field(..., ge=0)
    wyckoff_index: Union[int, str] = "whole"  # "whole" or int
    max_rank: int = Field(3, ge=1, le=5)
    mode: str = Field("magnetic", pattern="^(magnetic|electric)$")
    include_soc: bool = True
