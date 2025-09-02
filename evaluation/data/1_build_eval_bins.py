#!/usr/bin/env python3

from __future__ import annotations

import difflib
import json
import random
import sys
import time
import unicodedata
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field

try:
    import requests  # type: ignore
except Exception:
    print(
        "This script requires the 'requests' package. "
        "Install with: pip install requests",
        file=sys.stderr,
    )
    raise

# -------------------------
# Configuration
# -------------------------

DEFAULT_BIN_DEFS: dict[str, list[str]] = {
    # Pragmatic chem_comp ID sets; tweak as needed
    "ADP/ATP": ["ATP", "ADP", "ANP", "AGS"],  # ANP: AMP-PNP; AGS: ATP-γ-S
    "FAD": ["FAD"],
    "FMN": ["FMN"],
    "NAD(P)": ["NAD", "NAI", "NAP", "NDP"],  # NAD+, NADH, NADP+, NADPH
    "PLP": ["PLP"],
    "HEME": ["HEM", "HEC", "HEA", "HEB"],  # b, c, a, b variants
    "Zn2+": ["ZN"],
    "Ca2+/Mg2+": ["CA", "MG"],
    "Sterol": ["CLR", "CHS"],  # cholesterol & cholesteryl sulfate
    "Misc.": ["SAM", "COA", "THF", "SF4", "FES", "F3S"],
}

# Synonyms and forgiving canonicalization for bin names
# (keys are **normalized** via _norm(); values are canonical bin names)
_BIN_SYNONYMS = {
    # ADP/ATP
    "adp/atp": "ADP/ATP",
    "adpatp": "ADP/ATP",
    "atp": "ADP/ATP",
    "adp": "ADP/ATP",
    # FAD
    "fad": "FAD",
    # FMN
    "fmn": "FMN",
    # NAD(P)
    "nadp": "NAD(P)",
    "nad": "NAD(P)",
    "nad(p)": "NAD(P)",
    "nadh": "NAD(P)",
    "nadph": "NAD(P)",
    # PLP
    "plp": "PLP",
    # HEME
    "heme": "HEME",
    "haem": "HEME",
    # Zn2+
    "zn": "Zn2+",
    "zn2+": "Zn2+",
    "zinc": "Zn2+",
    # Ca/Mg
    "camg": "Ca2+/Mg2+",
    "ca2+mg2+": "Ca2+/Mg2+",
    "calciummagnesium": "Ca2+/Mg2+",
    # Sterol
    "sterol": "Sterol",
    "cholesterol": "Sterol",
    # Misc
    "misc": "Misc.",
    "miscellaneous": "Misc.",
}

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"

USER_AGENT = "ligand-bin-builder/1.1 (+https://example.org)"
HTTP_TIMEOUT = 30
SLEEP_BETWEEN = 0.1  # be gentle
RETRIES = 3
BACKOFF = 0.6
RCSB_ROWS_PER_CHEM = (
    400  # cap number of hits we pull per chem id to keep runtime reasonable
)

# Reuse a single HTTP session to reduce connection overhead
_SESSION: requests.Session | None = None

# -------------------------
# Data Classes
# -------------------------


@dataclass
class PDBEvidence:
    resolution: float | None
    method: str | None
    matched_ligands: list[str] = field(default_factory=list)


@dataclass
class UniProtEntry:
    uniprot_id: str
    protein: str | None = None
    gene: str | None = None
    organism: str | None = None
    length: int | None = None
    ec_numbers: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    validation_pdb: list[str] = field(default_factory=list)
    pdb_evidence: dict[str, PDBEvidence] = field(default_factory=dict)
    all_supporting_pdbs: list[dict[str, object]] = field(default_factory=list)
    alphafold_url: str | None = None


# -------------------------
# Helpers
# -------------------------


def _norm(s: str) -> str:
    # Unicode fold (strip accents/superscripts), lower, drop spaces and
    # punctuation we don't need
    t = unicodedata.normalize("NFKD", s)
    t = t.encode("ascii", "ignore").decode("ascii")
    t = t.lower()
    for ch in [" ", "_", "-", "/", ".", ",", "(", ")"]:
        t = t.replace(ch, "")
    return t


def _canonical_map() -> dict[str, str]:
    m: dict[str, str] = {}
    for canon in DEFAULT_BIN_DEFS.keys():
        m[_norm(canon)] = canon
    # add synonyms
    for k, v in _BIN_SYNONYMS.items():
        m[_norm(k)] = v
    return m


# -------------------------
# HTTP helpers (with simple retries)
# -------------------------


def _get_json(url: str) -> dict | None:
    global _SESSION
    if _SESSION is None:
        try:
            _SESSION = requests.Session()  # type: ignore[attr-defined]
            _SESSION.headers.update({"User-Agent": USER_AGENT})
        except Exception:
            _SESSION = None
    for attempt in range(1, RETRIES + 1):
        try:
            if _SESSION is not None:
                resp = _SESSION.get(url, timeout=HTTP_TIMEOUT)
            else:
                resp = requests.get(
                    url, headers={"User-Agent": USER_AGENT}, timeout=HTTP_TIMEOUT
                )
            if resp.status_code == 200:
                return resp.json()
        except Exception:
            pass
        time.sleep(SLEEP_BETWEEN + BACKOFF * attempt)
    return None


def _post_json(url: str, payload: dict) -> dict | None:
    global _SESSION
    if _SESSION is None:
        try:
            _SESSION = requests.Session()  # type: ignore[attr-defined]
            _SESSION.headers.update({"User-Agent": USER_AGENT})
        except Exception:
            _SESSION = None
    for attempt in range(1, RETRIES + 1):
        try:
            if _SESSION is not None:
                resp = _SESSION.post(url, json=payload, timeout=HTTP_TIMEOUT)
            else:
                resp = requests.post(
                    url,
                    json=payload,
                    headers={"User-Agent": USER_AGENT},
                    timeout=HTTP_TIMEOUT,
                )
            if resp.status_code == 200:
                return resp.json()
        except Exception:
            pass
        time.sleep(SLEEP_BETWEEN + BACKOFF * attempt)
    return None


# -------------------------
# RCSB queries
# -------------------------


def rcsb_search_entries_for_chem(
    chem_id: str, resolution_max: float, method: str | None = "X-RAY DIFFRACTION"
) -> set[str]:
    """Return PDB IDs that contain chem_id, match method and resolution threshold.
    If method is None or 'ANY', don't filter by method.

    Uses a verified search attribute: rcsb_chem_comp_container_identifiers.comp_id
    with exact_match, combined with resolution and optional method filters.
    """

    chem_id = (chem_id or "").upper()

    nodes: list[dict] = [
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
            # Limit rows per chem id to keep the pipeline fast while still
            # finding many entries
            "paginate": {"start": 0, "rows": RCSB_ROWS_PER_CHEM},
            "results_content_type": ["experimental"],
        },
    }
    data = _post_json(RCSB_SEARCH_URL, query)
    time.sleep(SLEEP_BETWEEN)
    if not data or "result_set" not in data:
        return set()
    try:
        return {item["identifier"] for item in data.get("result_set", [])}
    except Exception:
        return set()


def rcsb_entry_details(pdb_id: str) -> tuple[float | None, str | None, list[str]]:
    """Return (best_resolution, method, polymer_entity_ids) for an entry."""
    url = RCSB_ENTRY_URL.format(pdb_id=pdb_id)
    data = _get_json(url)
    time.sleep(SLEEP_BETWEEN)
    if not data:
        return None, None, []
    # resolution_combined is a list; take min if present
    res = None
    try:
        rc = data.get("rcsb_entry_info", {}).get("resolution_combined")
        if isinstance(rc, list) and rc:
            res = min([x for x in rc if isinstance(x, int | float)])
        elif isinstance(rc, int | float):
            res = float(rc)
    except Exception:
        res = None
    method = None
    try:
        exptl = data.get("exptl")
        if isinstance(exptl, list) and exptl:
            method = exptl[0].get("method")
    except Exception:
        method = None
    entity_ids = []
    try:
        entity_ids = (
            data.get("rcsb_entry_container_identifiers", {}).get(
                "polymer_entity_ids", []
            )
            or []
        )
    except Exception:
        entity_ids = []
    return res, method, [str(eid) for eid in entity_ids]


def rcsb_entry_uniprot_ids(pdb_id: str, entity_ids: Iterable[str]) -> set[str]:
    """Return UniProt accessions mapped in the given polymer entities of this entry."""
    uniprots: set[str] = set()
    for eid in entity_ids:
        url = RCSB_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid)
        data = _get_json(url)
        time.sleep(SLEEP_BETWEEN)
        if not data:
            continue
        ids = (
            data.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids")
            or []
        )
        for acc in ids:
            if isinstance(acc, str):
                uniprots.add(acc)
    return uniprots


# -------------------------
# UniProt metadata
# -------------------------


_UNIPROT_CACHE: dict[
    str,
    tuple[str | None, str | None, str | None, int | None, list[str], list[str]],
] = {}


def fetch_uniprot_meta(
    acc: str,
) -> tuple[str | None, str | None, str | None, int | None, list[str], list[str]]:
    """Return (protein_name, gene, organism, length, ec_numbers, keywords)."""
    if acc in _UNIPROT_CACHE:
        return _UNIPROT_CACHE[acc]
    url = UNIPROT_URL.format(acc=acc)
    data = _get_json(url)
    time.sleep(SLEEP_BETWEEN)
    if not data:
        result = (None, None, None, None, [], [])
        _UNIPROT_CACHE[acc] = result
        return result
    protein_name = None
    try:
        protein_name = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value")
        )
    except Exception:
        protein_name = None
    gene = None
    try:
        genes = data.get("genes") or []
        if genes:
            gene = genes[0].get("geneName", {}).get("value")
    except Exception:
        gene = None
    organism = None
    try:
        organism = data.get("organism", {}).get("scientificName")
    except Exception:
        organism = None
    length = None
    try:
        length = int(data.get("sequence", {}).get("length"))
    except Exception:
        length = None
    ec_numbers: list[str] = []
    try:
        ec = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("ecNumbers", [])
        )
        for e in ec or []:
            v = e.get("value")
            if isinstance(v, str):
                ec_numbers.append(v)
    except Exception:
        pass
    keywords: list[str] = []
    try:
        for kw in data.get("keywords", []) or []:
            v = kw.get("value")
            if isinstance(v, str):
                keywords.append(v)
        # keep it compact
        if len(keywords) > 10:
            keywords = keywords[:10]
    except Exception:
        pass
    result = (protein_name, gene, organism, length, ec_numbers, keywords)
    _UNIPROT_CACHE[acc] = result
    return result


# -------------------------
# Assembly logic
# -------------------------


def build_bin(
    bin_name: str,
    comp_ids: list[str],
    resolution_max: float,
    method: str | None,
    max_per_bin: int | None,
    seed: int | None,
) -> tuple[dict, list[UniProtEntry]]:
    """Build a single bin as a dict matching the grouped JSON schema."""
    rng = random.Random(seed)

    # 1) Find all entries matching any comp_id
    entry_to_matched_ligs: dict[str, set[str]] = defaultdict(set)
    # Randomize order of comp_ids to diversify early hits
    for chem in rng.sample(comp_ids, k=len(comp_ids)):
        pdbs = rcsb_search_entries_for_chem(chem, resolution_max, method)
        for pdb_id in pdbs:
            entry_to_matched_ligs[pdb_id].add(chem)

    # 2) For each entry, get resolution/method and UniProt IDs
    uniprot_to_evidence: dict[str, dict[str, PDBEvidence]] = defaultdict(dict)
    pdb_meta_cache: dict[str, tuple[float | None, str | None, list[str]]] = {}

    # Early-stop threshold if user asked for a max per bin; gather a small cushion
    target_uniprots = (
        max_per_bin * 2 if isinstance(max_per_bin, int) and max_per_bin > 0 else None
    )

    # Iterate over a shuffled list of PDB IDs to get diverse early coverage
    pdb_ids = list(entry_to_matched_ligs.keys())
    rng.shuffle(pdb_ids)
    for pdb_id in pdb_ids:
        matched_ligs = entry_to_matched_ligs[pdb_id]
        res, meth, entity_ids = rcsb_entry_details(pdb_id)
        pdb_meta_cache[pdb_id] = (res, meth, entity_ids)
        if not entity_ids:
            continue
        uniprot_ids = rcsb_entry_uniprot_ids(pdb_id, entity_ids)
        for up in uniprot_ids:
            uniprot_to_evidence[up][pdb_id] = PDBEvidence(
                resolution=res, method=meth, matched_ligands=sorted(matched_ligs)
            )
        # If we already have enough UniProt entries, stop fetching more PDBs
        # for this bin
        if target_uniprots is not None and len(uniprot_to_evidence) >= target_uniprots:
            break

    # 3) For each UniProt, pick best validation PDB (lowest resolution)
    entries: list[UniProtEntry] = []
    for up, evid in uniprot_to_evidence.items():
        # Best PDB by lowest resolution (None treated as +inf)
        best_pdb = None
        best_res = float("inf")
        for pdb_id, ev in evid.items():
            r = (
                ev.resolution
                if isinstance(ev.resolution, int | float)
                else float("inf")
            )
            if r < best_res:
                best_res = r
                best_pdb = pdb_id
        # Fetch UniProt metadata
        pname, gene, org, length, ec_numbers, keywords = fetch_uniprot_meta(up)
        ue = UniProtEntry(
            uniprot_id=up,
            protein=pname,
            gene=gene,
            organism=org,
            length=length,
            ec_numbers=ec_numbers,
            keywords=keywords,
            validation_pdb=[best_pdb] if best_pdb else [],
            pdb_evidence=evid,
            alphafold_url=f"https://alphafold.ebi.ac.uk/entry/{up}",
        )
        # Expand all_supporting_pdbs for convenience
        for pdb_id, ev in evid.items():
            ue.all_supporting_pdbs.append(
                {
                    "pdb_id": pdb_id,
                    "resolution": ev.resolution,
                    "method": ev.method,
                    "matched_ligands": ev.matched_ligands,
                }
            )
        entries.append(ue)

    # 4) Randomly sample per-bin if requested
    if max_per_bin is not None and len(entries) > max_per_bin:
        entries = rng.sample(entries, max_per_bin)

    # 5) Serialize (minimal fields requested)
    def serialize_entry(ue: UniProtEntry) -> dict:
        best_pdb = ue.validation_pdb[0] if ue.validation_pdb else None
        matched_ligs = (
            ue.pdb_evidence.get(best_pdb).matched_ligands  # type: ignore[union-attr]
            if best_pdb and best_pdb in ue.pdb_evidence
            else []
        )
        return {
            "uniprot_id": ue.uniprot_id,
            "protein": ue.protein,
            "gene": ue.gene,
            "validation_pdb": best_pdb,
            "matched_ligands": matched_ligs,
        }

    return (
        {
            "ligand_bin": bin_name,
            "ligand_comp_ids": comp_ids,
            "entries": [serialize_entry(e) for e in entries],
        },
        entries,
    )


# -------------------------
# Protein class bins
# -------------------------


def classify_protein_classes(ue: UniProtEntry) -> set[str]:
    """Infer broad protein classes from EC numbers and UniProt keywords.
    Returns a set of class names such as 'Kinase', 'Oxidoreductase', etc.
    """
    classes: set[str] = set()
    ecs = ue.ec_numbers or []
    kws = [k.lower() for k in (ue.keywords or [])]

    # EC-based classes
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

    # Keyword-based hints
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

    return classes


# -------------------------
# CLI & utilities
# -------------------------


def _resolve_bins(user_bins: list[str], fuzzy: bool) -> list[tuple[str, list[str]]]:
    """Resolve user-provided bin names (case-insensitive, synonyms, optional fuzzy).
    Returns list of (canonical_name, comp_ids).
    """
    if not user_bins:
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]

    # special token 'all'
    if any(_norm(b) == "all" for b in user_bins):
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]

    cmap = _canonical_map()
    resolved: list[tuple[str, list[str]]] = []
    unknown: list[str] = []
    for b in user_bins:
        key = cmap.get(_norm(b))
        if key:
            resolved.append((key, DEFAULT_BIN_DEFS[key]))
        else:
            unknown.append(b)

    if not resolved and fuzzy and unknown:
        # Try difflib suggestions against canonical names
        canon_names = list(DEFAULT_BIN_DEFS.keys())
        for ub in unknown:
            guess = difflib.get_close_matches(ub, canon_names, n=1, cutoff=0.5)
            if guess:
                g = guess[0]
                resolved.append((g, DEFAULT_BIN_DEFS[g]))
                print(
                    f"[WARN] Unknown bin '{ub}'. Using closest match '{g}'.",
                    file=sys.stderr,
                )

    if not resolved:
        # Final graceful fallback: use all bins but warn prominently
        print(
            f"[WARN] None of the provided bins could be resolved: {unknown}. "
            "Using ALL bins.",
            file=sys.stderr,
        )
        return [(name, DEFAULT_BIN_DEFS[name]) for name in DEFAULT_BIN_DEFS.keys()]

    return resolved


def _write_output(obj: dict, outfile: str) -> None:
    s = json.dumps(obj, indent=2)
    if outfile == "-":
        sys.stdout.write(s + "\n")
    else:
        with open(outfile, "w", encoding="utf-8") as fh:
            fh.write(s)


def main() -> int:
    # Simple script configuration – edit here
    outfile = "evaluation/data/evaluation_bins.json"  # Use '-' for stdout
    user_bins: list[str] = []  # [] = all bins, or e.g. ["HEME", "ADP/ATP"]
    resolution_max: float = 2.5
    method: str | None = "X-RAY DIFFRACTION"  # or "ANY"/None
    max_per_bin: int | None = 30  # small for testing
    seed: int = 42
    fuzzy: bool = True
    dry_run: bool = False

    selected_bins = _resolve_bins(user_bins or [], fuzzy=fuzzy)

    out = {"ligand_bins": [], "protein_bins": []}
    ligand_rows: list[dict] = []
    all_entries: list[UniProtEntry] = []

    try:
        from tqdm import tqdm  # type: ignore

        iterator = tqdm(selected_bins, desc="Building bins")
    except Exception:
        iterator = selected_bins

    try:
        for bin_name, comp_ids in iterator:
            print(f"[INFO] Building bin '{bin_name}' with comp IDs {comp_ids}")
            if dry_run:
                b = {
                    "ligand_bin": bin_name,
                    "ligand_comp_ids": comp_ids,
                    "entries": [],
                }
                entries = []
            else:
                b, entries = build_bin(
                    bin_name=bin_name,
                    comp_ids=comp_ids,
                    resolution_max=resolution_max,
                    method=method,
                    max_per_bin=max_per_bin,
                    seed=seed,
                )
            all_entries.extend(entries)
            # Flat ligand rows for pandas
            for ue in entries:
                best_pdb = ue.validation_pdb[0] if ue.validation_pdb else None
                matched_ligs_list = (
                    ue.pdb_evidence.get(best_pdb).matched_ligands  # type: ignore[union-attr]
                    if best_pdb and best_pdb in ue.pdb_evidence
                    else []
                )
                ligand_rows.append(
                    {
                        "ligand_bin": bin_name,
                        "uniprot_id": ue.uniprot_id,
                        "protein": ue.protein,
                        "gene": ue.gene,
                        "validation_pdb": best_pdb,
                        "matched_ligands": "|".join(matched_ligs_list),
                    }
                )
            if not b["entries"]:
                print(
                    f"[WARN] Bin '{bin_name}' produced 0 entries "
                    f"(filters may be too strict).",
                    file=sys.stderr,
                )
            out["ligand_bins"].append(b)
    except KeyboardInterrupt:
        print("\n[WARN] Interrupted. Writing partial results...", file=sys.stderr)
        _write_output(out, outfile)
        print(f"[OK] Wrote partial {outfile if outfile != '-' else 'STDOUT'}")
        return 0

    # Build protein class bins from aggregated entries
    if all_entries:
        rng = random.Random(seed)
        class_to_entries: dict[str, list[UniProtEntry]] = defaultdict(list)
        class_seen: dict[str, set[str]] = defaultdict(set)
        for ue in all_entries:
            classes = classify_protein_classes(ue)
            if not classes:
                classes = {"Unclassified"}
            for cls in classes:
                if ue.uniprot_id not in class_seen[cls]:
                    class_seen[cls].add(ue.uniprot_id)
                    class_to_entries[cls].append(ue)

        def serialize_entry_min(ue: UniProtEntry) -> dict:
            best_pdb = ue.validation_pdb[0] if ue.validation_pdb else None
            matched_ligs = (
                ue.pdb_evidence.get(best_pdb).matched_ligands  # type: ignore[union-attr]
                if best_pdb and best_pdb in ue.pdb_evidence
                else []
            )
            return {
                "uniprot_id": ue.uniprot_id,
                "protein": ue.protein,
                "gene": ue.gene,
                "validation_pdb": best_pdb,
                "matched_ligands": matched_ligs,
            }

        protein_bins = []
        protein_rows: list[dict] = []
        for cls, entries_for_class in class_to_entries.items():
            pick = entries_for_class
            if (
                isinstance(max_per_bin, int)
                and max_per_bin > 0
                and len(pick) > max_per_bin
            ):
                pick = rng.sample(entries_for_class, max_per_bin)
            protein_bins.append(
                {
                    "protein_class": cls,
                    "entries": [serialize_entry_min(e) for e in pick],
                }
            )
            for ue in pick:
                best_pdb = ue.validation_pdb[0] if ue.validation_pdb else None
                matched_ligs_list = (
                    ue.pdb_evidence.get(best_pdb).matched_ligands  # type: ignore[union-attr]
                    if best_pdb and best_pdb in ue.pdb_evidence
                    else []
                )
                protein_rows.append(
                    {
                        "protein_class": cls,
                        "uniprot_id": ue.uniprot_id,
                        "protein": ue.protein,
                        "gene": ue.gene,
                        "validation_pdb": best_pdb,
                        "matched_ligands": "|".join(matched_ligs_list),
                    }
                )
        # stable order
        protein_bins.sort(key=lambda x: x["protein_class"])
        out["protein_bins"] = protein_bins
        # attach flat rows for pandas
        out["ligand_rows"] = ligand_rows
        out["protein_rows"] = protein_rows

    _write_output(out, outfile)

    print(f"[OK] Wrote {outfile if outfile != '-' else 'STDOUT'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
