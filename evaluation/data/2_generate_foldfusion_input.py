#!/usr/bin/env python3
"""
Raw script to extract all unique UniProt IDs from
evaluation/data/evaluation_bins.json and write them one-per-line to
evaluation/foldfusion_input/uniprot_ids.txt.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def collect_uniprot_ids(obj: Any, acc: set[str]) -> None:
    """Recursively collect values for keys named 'uniprot_id'."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == "uniprot_id" and isinstance(value, str) and value.strip():
                acc.add(value.strip())
            if isinstance(value, dict | list):
                collect_uniprot_ids(value, acc)
    elif isinstance(obj, list):
        for item in obj:
            collect_uniprot_ids(item, acc)


def main() -> None:
    # Resolve repo root from this file's location
    repo_root = Path(__file__).resolve().parents[2]
    json_path = repo_root / "evaluation/data/evaluation_bins.json"
    out_path = repo_root / "evaluation/data/foldfusion_input/uniprot_ids.txt"

    if not json_path.is_file():
        raise SystemExit(f"Input JSON not found: {json_path}")

    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    ids: set[str] = set()
    collect_uniprot_ids(data, ids)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(sorted(ids)) + "\n", encoding="utf-8")
    print(f"Wrote {len(ids)} UniProt IDs to {out_path}")


if __name__ == "__main__":
    main()
