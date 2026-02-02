import json
from .config import DB_PATH

_DB_CACHE = None

def load_db() -> list[dict]:
    global _DB_CACHE
    if _DB_CACHE is None:
        if not DB_PATH.exists():
            raise FileNotFoundError(f"Missing {DB_PATH}")
        _DB_CACHE = json.loads(DB_PATH.read_text(encoding="utf-8"))
    return _DB_CACHE
