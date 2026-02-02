import json
from pathlib import Path
import sys

# Add the backend directory to the python path to allow imports
backend_dir = Path(__file__).resolve().parent.parent / 'backend'
sys.path.insert(0, str(backend_dir))

from app.engine.multipole import compute_results

def main():
    db_path = backend_dir / "data" / "mspgr.json"
    if not db_path.exists():
        print(f"Error: Database file not found at {db_path}")
        return

    with open(db_path, 'r', encoding='utf-8') as f:
        db = json.load(f)

    target_group_name = "166.97"
    target_wyckoff_label = "9e"

    group = None
    for g in db:
        if target_group_name in g.get("name", ""):
            group = g
            break

    if not group:
        print(f"Error: Group containing '{target_group_name}' not found.")
        return

    wyckoff = None
    for w in group.get("wyckoff", []):
        if w.get("label") == target_wyckoff_label:
            wyckoff = w
            break

    if not wyckoff:
        print(f"Error: Wyckoff position '{target_wyckoff_label}' not found in group '{group.get('name')}'.")
        return

    results = compute_results(
        group=group,
        wyckoff=wyckoff,
        max_rank=3,
        mode="magnetic",
        include_soc=True,
        db_all=db
    )

    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
