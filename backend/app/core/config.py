from pathlib import Path

# This file is: backend/app/core/config.py
# parents[2] = backend/
BASE_DIR = Path(__file__).resolve().parents[2]

# Data
DB_PATH = BASE_DIR / "data" / "mspgr.json"

# Frontend build (Vite)
# dist is at project root: ../dist
DIST_DIR = BASE_DIR.parent / "dist"
