#!/usr/bin/env bash
set -euo pipefail

# Creates a local Python virtual environment, installs backend deps
# and installs frontend Node dependencies.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="$ROOT_DIR/.venv"

echo "Setting up local development environment..."

if [ ! -d "$VENV_DIR" ]; then
  echo "Creating virtual environment at $VENV_DIR"
  python3 -m venv "$VENV_DIR"
fi

echo "Activating virtual environment and upgrading pip..."
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip

echo "Installing backend Python packages from backend/requirements.txt"
pip install -r "$ROOT_DIR/backend/requirements.txt"

if [ -d "$ROOT_DIR/frontend" ]; then
  echo "Installing frontend dependencies (node/npm must be installed)"
  (cd "$ROOT_DIR/frontend" && npm install)
fi

echo "Setup complete. To activate the venv run: source $VENV_DIR/bin/activate"
echo "Start backend: uvicorn app.main:app --reload --host 127.0.0.1 --port 8000"
echo "Start frontend: (cd frontend && npm run dev)"
