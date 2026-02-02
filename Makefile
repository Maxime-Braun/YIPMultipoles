PYTHON?=python3
VENV=.venv

.PHONY: setup backend-deps frontend-deps run-backend run-frontend

setup: backend-deps frontend-deps

backend-deps:
	$(PYTHON) -m venv $(VENV)
	@echo "Activating venv and installing backend packages..."
	@. $(VENV)/bin/activate && python -m pip install --upgrade pip && pip install -r backend/requirements.txt

frontend-deps:
	@if [ -d frontend ]; then \
		echo "Installing frontend node dependencies..."; \
		cd frontend && npm install; \
	fi

run-backend:
	@. $(VENV)/bin/activate && uvicorn app.main:app --reload --host 127.0.0.1 --port 8000

run-frontend:
	cd frontend && npm run dev
