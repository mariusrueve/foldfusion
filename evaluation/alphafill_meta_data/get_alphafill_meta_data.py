#!/usr/bin/env python3
"""
Fetch AlphaFill metadata for UniProt IDs listed in evaluation/data/evaluation_bins.json
and save one JSON file per UniProt ID into evaluation/data/alphafill/{UNIPROT}.json.

API reference: https://alphafill.eu/man/alphafill-api/

Strategy:
- Try GET https://alphafill.eu/v1/aff/AF-{UNIPROT}-F1/json first (common AlphaFold model ID)
- If not found (404), query https://alphafill.eu/v1/aff/3d-beacon/{UNIPROT} and try to
  extract a valid AlphaFill ID (e.g., AF-{UNIPROT}-F{n}) from the response, then fetch its JSON.

This script uses only Python stdlib (urllib) to avoid extra dependencies.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import typing as t
from dataclasses import dataclass
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ALPHAFILL_BASE = "https://alphafill.eu/v1/aff"


@dataclass
class FetchResult:
    uniprot_id: str
    alphafill_id: t.Optional[str]
    status: str  # "ok" | "not_found" | "error"
    error: t.Optional[str] = None


def read_uniprot_ids(bins_path: str) -> list[str]:
    with open(bins_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    ids: set[str] = set()
    for bin_entry in data.get("ligand_bins", []):
        for e in bin_entry.get("entries", []):
            uid = e.get("uniprot_id")
            if isinstance(uid, str) and uid:
                ids.add(uid.strip())
    return sorted(ids)


def http_get_json(
    url: str, timeout: float = 30.0, headers: t.Optional[dict] = None
) -> t.Any:
    req = Request(
        url,
        headers={"User-Agent": "foldfusion/alphafill-meta-fetcher", **(headers or {})},
    )
    with urlopen(req, timeout=timeout) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        text = resp.read().decode(charset, errors="replace")
        return json.loads(text)


def try_fetch_alphafill_json_by_id(aff_id: str, timeout: float = 30.0) -> t.Any:
    url = f"{ALPHAFILL_BASE}/{aff_id}/json"
    return http_get_json(url, timeout=timeout)


def find_candidate_aff_ids_from_3dbeacon(
    uniprot_id: str, timeout: float = 30.0
) -> list[str]:
    """Query the 3d-beacon endpoint and extract strings that look like AF-{UID}-F{n}.

    Schema is not guaranteed, so we scan for plausible IDs in the JSON tree.
    """
    url = f"{ALPHAFILL_BASE}/3d-beacon/{uniprot_id}"
    try:
        obj = http_get_json(url, timeout=timeout)
    except Exception:
        return []

    pattern = re.compile(rf"\bAF-{re.escape(uniprot_id)}-F\d+\b")
    found: set[str] = set()

    def walk(x: t.Any):
        if isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, str):
            for m in pattern.findall(x):
                found.add(m)

    walk(obj)
    # Prefer F1 if available
    ordered = sorted(found, key=lambda s: (0 if s.endswith("-F1") else 1, s))
    return ordered


def fetch_alphafill_metadata(
    uniprot_id: str, *, timeout: float = 30.0, sleep_between: float = 0.2
) -> FetchResult:
    # First try the common canonical AlphaFold model id
    primary_id = f"AF-{uniprot_id}-F1"
    try:
        _ = try_fetch_alphafill_json_by_id(primary_id, timeout=timeout)
        return FetchResult(uniprot_id=uniprot_id, alphafill_id=primary_id, status="ok")
    except HTTPError as e:
        if e.code != 404:
            return FetchResult(
                uniprot_id=uniprot_id,
                alphafill_id=None,
                status="error",
                error=f"HTTP {e.code}: {e.reason}",
            )
        # 404 -> try to discover via 3d-beacon
    except URLError as e:
        return FetchResult(
            uniprot_id=uniprot_id,
            alphafill_id=None,
            status="error",
            error=f"URL error: {e.reason}",
        )
    except Exception as e:
        return FetchResult(
            uniprot_id=uniprot_id, alphafill_id=None, status="error", error=str(e)
        )

    # Fallback: discover IDs via 3d-beacon and try those
    time.sleep(sleep_between)
    candidates = find_candidate_aff_ids_from_3dbeacon(uniprot_id, timeout=timeout)
    for aff_id in candidates:
        try:
            _ = try_fetch_alphafill_json_by_id(aff_id, timeout=timeout)
            return FetchResult(uniprot_id=uniprot_id, alphafill_id=aff_id, status="ok")
        except HTTPError as e:
            if e.code == 404:
                continue
            return FetchResult(
                uniprot_id=uniprot_id,
                alphafill_id=None,
                status="error",
                error=f"HTTP {e.code}: {e.reason}",
            )
        except URLError as e:
            return FetchResult(
                uniprot_id=uniprot_id,
                alphafill_id=None,
                status="error",
                error=f"URL error: {e.reason}",
            )
        except Exception as e:
            return FetchResult(
                uniprot_id=uniprot_id, alphafill_id=None, status="error", error=str(e)
            )

    return FetchResult(uniprot_id=uniprot_id, alphafill_id=None, status="not_found")


def main(argv: list[str]) -> int:
    # Paths relative to repo root
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    bins_path = os.path.join(repo_root, "evaluation", "data", "evaluation_bins.json")
    out_dir = os.path.join(repo_root, "evaluation", "data", "alphafill")

    # Allow overriding via CLI args
    if len(argv) > 1:
        bins_path = os.path.abspath(argv[1])
    if len(argv) > 2:
        out_dir = os.path.abspath(argv[2])

    os.makedirs(out_dir, exist_ok=True)

    uniprot_ids = read_uniprot_ids(bins_path)
    if not uniprot_ids:
        print(f"No UniProt IDs found in {bins_path}")
        return 1

    print(f"Found {len(uniprot_ids)} UniProt IDs. Fetching AlphaFill metadata…")

    # Throttle to be polite
    sleep_between = 0.25
    successes = 0
    not_found = 0
    errors = 0

    for i, uid in enumerate(uniprot_ids, 1):
        result = fetch_alphafill_metadata(
            uid, timeout=45.0, sleep_between=sleep_between
        )
        if result.status == "ok" and result.alphafill_id:
            try:
                data = try_fetch_alphafill_json_by_id(result.alphafill_id, timeout=60.0)
                out_path = os.path.join(out_dir, f"{uid}.json")
                with open(out_path, "w", encoding="utf-8") as f:
                    json.dump(data, f, ensure_ascii=False, indent=2)
                successes += 1
                print(f"[{i}/{len(uniprot_ids)}] {uid}: saved ({result.alphafill_id})")
            except Exception as e:
                errors += 1
                print(
                    f"[{i}/{len(uniprot_ids)}] {uid}: error saving {result.alphafill_id}: {e}"
                )
        elif result.status == "not_found":
            not_found += 1
            print(f"[{i}/{len(uniprot_ids)}] {uid}: not found on AlphaFill")
        else:
            errors += 1
            print(f"[{i}/{len(uniprot_ids)}] {uid}: error: {result.error}")

        time.sleep(sleep_between)

    print(f"Done. OK={successes}, NotFound={not_found}, Errors={errors}")
    return 0 if successes > 0 else 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
