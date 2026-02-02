# Multipole Calculator — Local development

Quick steps to create a local environment and run the app.

Prerequisites:
- Python 3.10+ (python3)
- Node.js + npm (for the frontend)

1) Create virtualenv and install dependencies (recommended):

  ./scripts/setup_local_env.sh

This will create a Python virtual environment at `./.venv`, install Python packages
from `backend/requirements.txt`, and install frontend Node packages in `frontend/`.

2) Backend (after activating venv):

  source .venv/bin/activate
  uvicorn app.main:app --reload --host 127.0.0.1 --port 8000

The API will be available at http://127.0.0.1:8000 and the frontend (Vite) at
http://127.0.0.1:5173 if you run the frontend dev server.

3) Frontend:

  cd frontend
  npm run dev

Notes:
- The backend requirements are stored in `backend/requirements.txt`.
- If you'd rather use `make`:

  make setup
  make run-backend
  make run-frontend
