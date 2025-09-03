#!/usr/bin/env python3
"""
Build a benchmarking dataset for ligand- and protein-class–based bins.

Inputs (CLI):
  - items per bin (int)
  - output JSON path

Outputs:
  - Main JSON (pandas-friendly "rows" + detailed "ligand_bins" and "protein_bins")
  - TXT of unique UniProt IDs (one per line)
  - Per-UniProt metadata saved locally (UniProt JSON always; AlphaFill JSON if available)

Notes:
  - Uses RCSB Search/Data API to ensure real PDB validation entries containing requested ligands.
  - Uses UniProt REST API to fetch protein/gene/keywords/EC metadata.
  - Uses AlphaFill REST API to detect and (if present) save AlphaFill transplant metadata.
  - Skips any accession where required pieces are missing; keeps going until each bin reaches target size (if possible).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import random
import sys
import time
import unicodedata
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

try:
    import requests  # type: ignore
except Exception:
    print("This script requires 'requests' (pip install requests)", file=sys.stderr)
    raise

# ------------- Config (tweakables / sane defaults) -----------------

# Ligand bins: CCD ids per group (you can expand/tweak freely)
DEFAULT_BIN_DEFS: dict[str, list[str]] = {
    "ADP/ATP": ["ATP", "ADP", "ANP", "AGS"],  # ANP: AMP-PNP; AGS: ATP-γ-S
    "FAD": ["FAD"],
    "FMN": ["FMN"],
    "NAD(P)": ["NAD", "NAI", "NAP", "NDP"],  # NAD+, NADH, NADP+, NADPH
    "PLP": ["PLP"],
    "HEME": ["HEM", "HEC", "HEA", "HEB"],
    "Zn2+": ["ZN"],
    "Ca2+/Mg2+": ["CA", "MG"],
    "Sterol": ["CLR", "CHS"],
    "Misc.": ["SAM", "COA", "THF", "SF4", "FES", "F3S"],
}

# Bin synonyms / normalization helpers (so ADP goes into "ADP/ATP", etc.)
_BIN_SYNONYMS = {
    "adp/atp": "ADP/ATP",
    "adpatp": "ADP/ATP",
    "atp": "ADP/ATP",
    "adp": "ADP/ATP",
    "fad": "FAD",
    "fmn": "FMN",
    "nad": "NAD(P)",
    "nadp": "NAD(P)",
    "nad(p)": "NAD(P)",
    "nadh": "NAD(P)",
    "nadph": "NAD(P)",
    "plp": "PLP",
    "heme": "HEME",
    "haem": "HEME",
    "zn": "Zn2+",
    "zn2+": "Zn2+",
    "zinc": "Zn2+",
    "camg": "Ca2+/Mg2+",
    "ca2+mg2+": "Ca2+/Mg2+",
    "calciummagnesium": "Ca2+/Mg2+",
    "sterol": "Sterol",
    "cholesterol": "Sterol",
    "misc": "Misc.",
    "miscellaneous": "Misc.",
}

# APIs
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ALPHAFILL_JSON_URL = "https://alphafill.eu/v1/aff/{acc}/json"
ALPHAFILL_STATUS_URL = "https://alphafill.eu/v1/aff/{acc}/status"

USER_AGENT = "alphafill-benchmark-builder/2.0 (mailto:example@example.org)"
HTTP_TIMEOUT = 30
SLEEP_BETWEEN = 0.10  # Be gentle to public APIs
RETRIES = 3
BACKOFF = 0.6
RCSB_ROWS_PER_CHEM = 1200  # Pull ample results per CCD to fill bins reliably

# Defaults if not provided on CLI
DEFAULT_RESOLUTION_MAX = 2.8
DEFAULT_METHOD = "X-RAY DIFFRACTION"  # set to "ANY" to disable method filter
DEFAULT_SEED = 1337

# -------------------- Data classes --------------------


@dataclass
class PDBEvidence:
    resolution: Optional[float]
    method: Optional[str]
    matched_ligands: list[str] = field(default_factory=list)


@dataclass
class UniProtEntry:
    uniprot_id: str
    protein: Optional[str] = None
    gene: Optional[str] = None
    organism: Optional[str] = None
    length: Optional[int] = None
    ec_numbers: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    validation_pdb: list[str] = field(default_factory=list)
    pdb_evidence: dict[str, PDBEvidence] = field(default_factory=dict)
    all_supporting_pdbs: list[dict] = field(default_factory=list)
    alphafill_json_saved: Optional[str] = None
    uniprot_json_saved: Optional[str] = None


# -------------------- Utilities --------------------

_SESSION: Optional[requests.Session] = None


def _session() -> requests.Session:
    global _SESSION
    if _SESSION is None:
        _SESSION = requests.Session()
        _SESSION.headers.update({"User-Agent": USER_AGENT})
    return _SESSION


def _sleep_retry(attempt: int) -> None:
    time.sleep(SLEEP_BETWEEN + BACKOFF * attempt)


def _get_json(url: str) -> Optional[dict]:
    sess = _session()
    for attempt in range(1, RETRIES + 1):
        try:
            r = sess.get(url, timeout=HTTP_TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (404, 400):
                return None
        except Exception:
            pass
        _sleep_retry(attempt)
    return None


def _post_json(url: str, payload: dict) -> Optional[dict]:
    sess = _session()
    for attempt in range(1, RETRIES + 1):
        try:
            r = sess.post(url, json=payload, timeout=HTTP_TIMEOUT)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (404, 400):
                return None
        except Exception:
            pass
        _sleep_retry(attempt)
    return None


def _norm(s: str) -> str:
    t = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    t = t.lower()
    for ch in " _-/.(),+":
        t = t.replace(ch, "")
    return t


def _canonical_map() -> dict[str, str]:
    m: dict[str, str] = {_norm(k): k for k in DEFAULT_BIN_DEFS.keys()}
    m.update({_norm(k): v for k, v in _BIN_SYNONYMS.items()})
    return m


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def save_json(obj: dict, path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


# -------------------- RCSB helpers --------------------


def rcsb_search_entries_for_chem(
    chem_id: str,
    resolution_max: float,
    method: Optional[str] = DEFAULT_METHOD,
) -> set[str]:
    """
    Entries containing a chemical component ID with optional method and resolution filters.
    Search attribute: rcsb_chem_comp_container_identifiers.comp_id (exact match)
    """
    chem_id = (chem_id or "").upper()
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                "operator": "exact_match",
                "value": chem_id,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "value": resolution_max,
            },
        },
    ]
    if method and method.upper() != "ANY":
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": method,
                },
            }
        )

    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": RCSB_ROWS_PER_CHEM},
            "results_content_type": ["experimental"],
        },
    }
    data = _post_json(RCSB_SEARCH_URL, query) or {}
    time.sleep(SLEEP_BETWEEN)
    return {item["identifier"] for item in data.get("result_set", [])}


def rcsb_entry_details(pdb_id: str) -> tuple[Optional[float], Optional[str], list[str]]:
    url = RCSB_ENTRY_URL.format(pdb_id=pdb_id)
    data = _get_json(url) or {}
    time.sleep(SLEEP_BETWEEN)
    # resolution
    res = None
    try:
        rc = data.get("rcsb_entry_info", {}).get("resolution_combined")
        if isinstance(rc, list) and rc:
            # use minimal value
            res = min([float(x) for x in rc if isinstance(x, (int, float))])
        elif isinstance(rc, (int, float)):
            res = float(rc)
    except Exception:
        res = None
    # method
    method = None
    try:
        exptl = data.get("exptl")
        if isinstance(exptl, list) and exptl:
            method = exptl[0].get("method")
    except Exception:
        method = None
    # polymer entities
    eids = []
    try:
        eids = (
            data.get("rcsb_entry_container_identifiers", {}).get(
                "polymer_entity_ids", []
            )
            or []
        )
        eids = [str(x) for x in eids]
    except Exception:
        eids = []
    return res, method, eids


def rcsb_entry_uniprot_ids(pdb_id: str, entity_ids: Iterable[str]) -> set[str]:
    accs: set[str] = set()
    for eid in entity_ids:
        url = RCSB_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid)
        data = _get_json(url) or {}
        time.sleep(SLEEP_BETWEEN)
        ids = (data.get("rcsb_polymer_entity_container_identifiers", {}) or {}).get(
            "uniprot_ids"
        ) or []
        for acc in ids:
            if isinstance(acc, str):
                accs.add(acc)
    return accs


# -------------------- UniProt & AlphaFill --------------------

_UNIPROT_CACHE: dict[str, dict] = {}


def fetch_uniprot_json(acc: str) -> Optional[dict]:
    if acc in _UNIPROT_CACHE:
        return _UNIPROT_CACHE[acc]
    data = _get_json(UNIPROT_ENTRY_URL.format(acc=acc))
    _UNIPROT_CACHE[acc] = data or {}
    return _UNIPROT_CACHE[acc]


def extract_uniprot_fields(
    data: Optional[dict],
) -> tuple[
    Optional[str], Optional[str], Optional[str], Optional[int], list[str], list[str]
]:
    if not data:
        return None, None, None, None, [], []
    # protein/gene/organism/length/ECs/keywords
    pname = None
    try:
        pname = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value")
        )
    except Exception:
        pass
    gene = None
    try:
        genes = data.get("genes") or []
        if genes:
            gene = genes[0].get("geneName", {}).get("value")
    except Exception:
        pass
    org = None
    try:
        org = data.get("organism", {}).get("scientificName")
    except Exception:
        pass
    length = None
    try:
        length = int(data.get("sequence", {}).get("length"))
    except Exception:
        pass
    ecs: list[str] = []
    try:
        for e in (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("ecNumbers")
            or []
        ):
            v = e.get("value")
            if isinstance(v, str):
                ecs.append(v)
    except Exception:
        pass
    keywords: list[str] = []
    try:
        for kw in data.get("keywords", []) or []:
            v = kw.get("value")
            if isinstance(v, str):
                keywords.append(v)
        if len(keywords) > 25:
            keywords = keywords[:25]
    except Exception:
        pass
    return pname, gene, org, length, ecs, keywords


def alphafill_has_json(acc: str) -> bool:
    # Prefer /status first (cheap), fall back to trying /json
    status = _get_json(ALPHAFILL_STATUS_URL.format(acc=acc))
    if (
        status
        and isinstance(status, dict)
        and status.get("status") in {"finished", "success", "ok"}
    ):
        return True
    # Some servers may not expose /status; probe /json directly
    data = _get_json(ALPHAFILL_JSON_URL.format(acc=acc))
    return bool(data)


def download_alphafill_json(acc: str) -> Optional[dict]:
    return _get_json(ALPHAFILL_JSON_URL.format(acc=acc))


# -------------------- Protein class assignment --------------------


def classify_protein_classes(ue: UniProtEntry) -> set[str]:
    """Coarse classes from EC numbers & keywords."""
    classes: set[str] = set()
    ecs = ue.ec_numbers or []
    kws = [k.lower() for k in (ue.keywords or [])]

    for ec in ecs:
        if ec.startswith("1."):
            classes.add("Oxidoreductase")
        if ec.startswith("2."):
            classes.add("Transferase")
        if ec.startswith("2.7"):
            classes.add("Kinase")
        if ec.startswith("3."):
            classes.add("Hydrolase")
        if ec.startswith("4."):
            classes.add("Lyase")
        if ec.startswith("5."):
            classes.add("Isomerase")
        if ec.startswith("6."):
            classes.add("Ligase")

    if any("kinase" in k for k in kws):
        classes.add("Kinase")
    if any("oxidoreductase" in k for k in kws):
        classes.add("Oxidoreductase")
    if any("transferase" in k for k in kws):
        classes.add("Transferase")
    if any("hydrolase" in k for k in kws):
        classes.add("Hydrolase")
    if any("lyase" in k for k in kws):
        classes.add("Lyase")
    if any("isomerase" in k for k in kws):
        classes.add("Isomerase")
    if any("ligase" in k for k in kws):
        classes.add("Ligase")
    if any("transport" in k or "transmembrane" in k for k in kws):
        classes.add("Transporter")

    if not classes:
        classes.add("Unclassified")
    return classes


# -------------------- Core builders --------------------


def pick_best_validation_pdb(evid: dict[str, PDBEvidence]) -> Optional[str]:
    best_id, best_res = None, float("inf")
    for pdb_id, ev in evid.items():
        r = ev.resolution if isinstance(ev.resolution, (int, float)) else float("inf")
        if r < best_res:
            best_res, best_id = r, pdb_id
    return best_id


def build_ligand_bin(
    bin_name: str,
    comp_ids: list[str],
    items_needed: int,
    resolution_max: float,
    method: Optional[str],
    rng: random.Random,
    save_dirs: dict[str, Path],
    require_alphafill: bool = True,
) -> tuple[dict, list[UniProtEntry]]:
    """
    Build a single ligand bin filling up to items_needed UniProt entries
    with (a) a validating PDB containing the ligand and (b) (optionally) AlphaFill JSON present.
    """
    entry_to_mligs: dict[str, set[str]] = defaultdict(set)

    # 1) accumulate PDB entries per chem id
    for chem in rng.sample(comp_ids, k=len(comp_ids)):
        pdbs = rcsb_search_entries_for_chem(chem, resolution_max, method)
        for pdb_id in pdbs:
            entry_to_mligs[pdb_id].add(chem)

    if not entry_to_mligs:
        return (
            {"ligand_bin": bin_name, "ligand_comp_ids": comp_ids, "entries": []},
            [],
        )

    pdb_ids = list(entry_to_mligs.keys())
    rng.shuffle(pdb_ids)

    uniprot_to_evidence: dict[str, dict[str, PDBEvidence]] = defaultdict(dict)
    # walk entries and collect UniProt mappings
    for pdb_id in pdb_ids:
        matched = sorted(entry_to_mligs[pdb_id])
        res, meth, eids = rcsb_entry_details(pdb_id)
        if not eids:
            continue
        ups = rcsb_entry_uniprot_ids(pdb_id, eids)
        for acc in ups:
            uniprot_to_evidence[acc][pdb_id] = PDBEvidence(
                resolution=res, method=meth, matched_ligands=matched
            )
        # early stop if surplus
        if len(uniprot_to_evidence) >= items_needed * 3:
            break

    # 2) assemble UniProt entries, enforce AlphaFill availability if requested, and save metadata
    entries: list[UniProtEntry] = []
    for acc, evid in rng.sample(
        list(uniprot_to_evidence.items()), k=len(uniprot_to_evidence)
    ):
        if len(entries) >= items_needed:
            break

        # AlphaFill availability check
        has_af = alphafill_has_json(acc) if require_alphafill else True
        if not has_af:
            continue

        up_json = fetch_uniprot_json(acc)
        pname, gene, org, length, ecs, kws = extract_uniprot_fields(up_json)

        ue = UniProtEntry(
            uniprot_id=acc,
            protein=pname,
            gene=gene,
            organism=org,
            length=length,
            ec_numbers=ecs,
            keywords=kws,
            pdb_evidence=evid,
        )
        best_pdb = pick_best_validation_pdb(evid)
        if best_pdb:
            ue.validation_pdb = [best_pdb]
        # expand supporting PDBs
        for pid, ev in evid.items():
            ue.all_supporting_pdbs.append(
                {
                    "pdb_id": pid,
                    "resolution": ev.resolution,
                    "method": ev.method,
                    "matched_ligands": ev.matched_ligands,
                }
            )

        # save per-UniProt metadata
        # UniProt JSON
        if up_json:
            up_path = save_dirs["uniprot"] / f"{acc}.json"
            save_json(up_json, up_path)
            ue.uniprot_json_saved = str(up_path)

        # AlphaFill JSON (if available)
        af_json = download_alphafill_json(acc) if require_alphafill else None
        if af_json:
            af_path = save_dirs["alphafill"] / f"{acc}.json"
            save_json(af_json, af_path)
            ue.alphafill_json_saved = str(af_path)
        elif require_alphafill:
            # Race condition (status true but json fetch failed) — skip
            continue

        # must have a validating PDB to include
        if not ue.validation_pdb:
            continue

        entries.append(ue)

    # serialize minimal entries for the bin summary
    def _ser(e: UniProtEntry) -> dict:
        best_pdb = e.validation_pdb[0] if e.validation_pdb else None
        matched = (
            e.pdb_evidence[best_pdb].matched_ligands
            if (best_pdb and best_pdb in e.pdb_evidence)
            else []
        )
        return {
            "uniprot_id": e.uniprot_id,
            "protein": e.protein,
            "gene": e.gene,
            "validation_pdb": best_pdb,
            "matched_ligands": matched,
            "uniprot_json": e.uniprot_json_saved,
            "alphafill_json": e.alphafill_json_saved,
        }

    b = {
        "ligand_bin": bin_name,
        "ligand_comp_ids": comp_ids,
        "entries": [_ser(e) for e in entries],
    }
    return b, entries


# -------------------- CLI & orchestration --------------------


def resolve_bins(user_bins: list[str]) -> list[tuple[str, list[str]]]:
    if not user_bins:
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]
    if any(_norm(b) == "all" for b in user_bins):
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]
    cmap = _canonical_map()
    resolved: list[tuple[str, list[str]]] = []
    unknown: list[str] = []
    for b in user_bins:
        canon = cmap.get(_norm(b))
        if canon:
            resolved.append((canon, DEFAULT_BIN_DEFS[canon]))
        else:
            unknown.append(b)
    if unknown and not resolved:
        print(f"[WARN] Unknown bins {unknown}; using ALL defined bins", file=sys.stderr)
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]
    return resolved


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build ligand & protein-class benchmark bins against AlphaFill + PDB."
    )
    ap.add_argument(
        "items_per_bin",
        type=int,
        help="Target number of UniProt accessions per bin (e.g., 300).",
    )
    ap.add_argument(
        "outfile",
        type=str,
        help=(
            "Output path: either a JSON file (e.g., dataset.json) or a directory. "
            "If a directory is given, a file named "
            "'benchmark_dataset.json' will be created inside."
        ),
    )
    ap.add_argument(
        "--bins",
        nargs="*",
        default=[],
        help="Subset of ligand bins to build (names or synonyms). Default: all.",
    )
    ap.add_argument(
        "--resolution-max",
        type=float,
        default=DEFAULT_RESOLUTION_MAX,
        help=f"Max resolution cutoff (Å). Default: {DEFAULT_RESOLUTION_MAX}",
    )
    ap.add_argument(
        "--method",
        type=str,
        default=DEFAULT_METHOD,
        help='Experimental method filter (e.g., "X-RAY DIFFRACTION" or "ANY").',
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"RNG seed. Default: {DEFAULT_SEED}",
    )
    ap.add_argument(
        "--allow-missing-alphafill",
        action="store_true",
        help="If set, do NOT require AlphaFill JSON to include an entry.",
    )
    ap.add_argument("--dry-run", action="store_true", help="Parse inputs and exit.")
    args = ap.parse_args()

    items_per_bin: int = max(1, int(args.items_per_bin))

    # Coerce outfile: accept either a JSON file path or a directory
    raw_out = Path(args.outfile).absolute()

    def _coerce_outfile(p: Path) -> Path:
        # If an explicit .json file was provided, use as-is
        if p.suffix.lower() == ".json":
            return p
        # If an existing directory was provided, write inside it
        if p.exists() and p.is_dir():
            return p / "benchmark_dataset.json"
        # If no suffix, treat as directory (even if it doesn't exist yet)
        if p.suffix == "":
            return p / "benchmark_dataset.json"
        # If a non-.json suffix was provided, normalize to .json
        return p.with_suffix(".json")

    outfile = _coerce_outfile(raw_out)
    outbase = outfile.with_suffix("")
    meta_root = outbase.parent / f"{outbase.name}_meta"
    uniprot_dir = meta_root / "uniprot"
    alphafill_dir = meta_root / "alphafill"
    ids_txt = outbase.parent / f"{outbase.name}_uniprots.txt"

    if args.dry_run:
        print(f"[OK] Dry run. Would build {items_per_bin}/bin -> {outfile}")
        return 0

    ensure_dir(uniprot_dir)
    ensure_dir(alphafill_dir)

    rng = random.Random(args.seed)
    selected_bins = resolve_bins(args.bins)

    out: dict = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "params": {
            "items_per_bin": items_per_bin,
            "resolution_max": args.resolution_max,
            "method": args.method,
            "seed": args.seed,
            "bins_requested": [b for b, _ in selected_bins],
        },
        "ligand_bins": [],
        "protein_bins": [],
        "rows": [],  # flat table for pandas
    }

    all_entries: list[UniProtEntry] = []
    flat_rows: list[dict] = []

    for bin_name, comp_ids in selected_bins:
        print(f"[INFO] Building ligand bin '{bin_name}' with CCDs {comp_ids}")
        b, entries = build_ligand_bin(
            bin_name=bin_name,
            comp_ids=comp_ids,
            items_needed=items_per_bin,
            resolution_max=args.resolution_max,
            method=args.method,
            rng=rng,
            save_dirs={"uniprot": uniprot_dir, "alphafill": alphafill_dir},
            require_alphafill=(not args.allow_missing_alphafill),
        )
        if len(entries) < items_per_bin:
            print(
                (
                    f"[WARN] Bin '{bin_name}' filled {len(entries)}/{items_per_bin}. "
                    "Consider loosening filters or increasing RCSB_ROWS_PER_CHEM."
                ),
                file=sys.stderr,
            )
        out["ligand_bins"].append(b)
        all_entries.extend(entries)

        # flat rows for pandas
        for ue in entries:
            best_pdb = ue.validation_pdb[0] if ue.validation_pdb else None
            matched = (
                ue.pdb_evidence[best_pdb].matched_ligands
                if (best_pdb and best_pdb in ue.pdb_evidence)
                else []
            )
            flat_rows.append(
                {
                    "bin_type": "ligand",
                    "bin_name": bin_name,
                    "uniprot_id": ue.uniprot_id,
                    "protein": ue.protein,
                    "gene": ue.gene,
                    "validation_pdb": best_pdb,
                    "matched_ligand": "|".join(matched),
                    "uniprot_json": ue.uniprot_json_saved,
                    "alphafill_json": ue.alphafill_json_saved,
                }
            )

    # Derive protein class bins from the aggregated set
    if all_entries:
        print("[INFO] Building protein-class bins from collected entries")
        class_to_entries: dict[str, list[UniProtEntry]] = defaultdict(list)
        seen: dict[str, set[str]] = defaultdict(set)
        for ue in all_entries:
            for cls in classify_protein_classes(ue):
                if ue.uniprot_id not in seen[cls]:
                    seen[cls].add(ue.uniprot_id)
                    class_to_entries[cls].append(ue)

        protein_bins: list[dict] = []
        for cls, es in class_to_entries.items():
            take = es if len(es) <= items_per_bin else rng.sample(es, items_per_bin)
            protein_bins.append(
                {
                    "protein_class": cls,
                    "entries": [
                        {
                            "uniprot_id": e.uniprot_id,
                            "protein": e.protein,
                            "gene": e.gene,
                            "validation_pdb": (
                                e.validation_pdb[0] if e.validation_pdb else None
                            ),
                            "matched_ligands": (
                                e.pdb_evidence[e.validation_pdb[0]].matched_ligands
                                if (
                                    e.validation_pdb
                                    and e.validation_pdb[0] in e.pdb_evidence
                                )
                                else []
                            ),
                            "uniprot_json": e.uniprot_json_saved,
                            "alphafill_json": e.alphafill_json_saved,
                        }
                        for e in take
                    ],
                }
            )
            # flat rows
            for e in take:
                best_pdb = e.validation_pdb[0] if e.validation_pdb else None
                matched = (
                    e.pdb_evidence[best_pdb].matched_ligands
                    if (best_pdb and best_pdb in e.pdb_evidence)
                    else []
                )
                flat_rows.append(
                    {
                        "bin_type": "protein",
                        "bin_name": cls,
                        "uniprot_id": e.uniprot_id,
                        "protein": e.protein,
                        "gene": e.gene,
                        "validation_pdb": best_pdb,
                        "matched_ligand": "|".join(matched),
                        "uniprot_json": e.uniprot_json_saved,
                        "alphafill_json": e.alphafill_json_saved,
                    }
                )

        protein_bins.sort(key=lambda x: x["protein_class"])
        out["protein_bins"] = protein_bins

    # Attach flat rows for pandas
    out["rows"] = flat_rows

    # Write outputs
    ensure_dir(outfile.parent)
    with outfile.open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"[OK] Wrote dataset JSON: {outfile}")

    # Uniprot IDs TXT
    uniq_ids = sorted({row["uniprot_id"] for row in flat_rows})
    with ids_txt.open("w", encoding="utf-8") as fh:
        for acc in uniq_ids:
            fh.write(f"{acc}\n")
    print(f"[OK] Wrote UniProt ID list: {ids_txt}")

    # Final summary
    print(
        f"[SUM] Ligand bins: {len(out['ligand_bins'])}; "
        f"Protein-class bins: {len(out['protein_bins'])}; "
        f"Rows: {len(out['rows'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
