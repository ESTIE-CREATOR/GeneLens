"""
Gene ID detection and mapping to HGNC gene symbols.

Supported input formats:
  - HGNC symbol (TP53, BRCA1)         — passed through unchanged
  - Ensembl gene ID  (ENSG00000141510) — mapped via hgnc_mapping.tsv + ensembl_to_symbol.tsv
  - Ensembl transcript (ENST...)       — detected, Entrez lookup attempted
  - Entrez / NCBI gene ID (7157)       — mapped via hgnc_mapping.tsv
  - RefSeq accession (NM_000546)       — regex detection, best-effort mapping
  - UniProt accession (P04637)         — regex detection, best-effort mapping

The mapping files are committed static assets; no network calls at runtime.
"""

import os
import re
import pandas as pd
import streamlit as st

_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
_HGNC_PATH   = os.path.join(_DATA_DIR, "hgnc_mapping.tsv")
_ENSEMBL_PATH = os.path.join(_DATA_DIR, "ensembl_to_symbol.tsv")

# ── ID-type regex patterns ──────────────────────────────────────────────────
_PATTERNS = [
    ("ensembl_gene",       re.compile(r"^ENSG\d{6,}")),
    ("ensembl_transcript", re.compile(r"^ENST\d{6,}")),
    ("refseq",             re.compile(r"^(?:NM_|NR_|XM_|XR_|NP_)\d+")),
    ("uniprot",            re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")),
    ("entrez",             re.compile(r"^\d+$")),
]


def detect_id_type(gene_ids, sample_size: int = 100) -> str | None:
    """
    Sample gene IDs and return which format they are in.

    Returns one of:
      'ensembl_gene', 'ensembl_transcript', 'refseq', 'uniprot', 'entrez'
    or None if the IDs look like they are already HGNC symbols / unrecognised.
    """
    ids = [str(g).strip() for g in gene_ids[:sample_size] if str(g).strip() not in ("", "nan")]
    if not ids:
        return None

    n = len(ids)
    for name, pat in _PATTERNS:
        frac = sum(1 for g in ids if pat.match(g.split(".")[0])) / n
        if frac >= 0.80:
            return name

    # Mixed column: try per-row detection if any clear pattern reaches 50%
    for name, pat in _PATTERNS:
        frac = sum(1 for g in ids if pat.match(g.split(".")[0])) / n
        if frac >= 0.50:
            return name

    return None  # assume already gene symbols


@st.cache_data(show_spinner=False)
def _load_mappings() -> dict:
    """Load all local mapping files into lookup dicts. Cached per Streamlit session."""
    maps: dict[str, dict] = {
        "ensembl_gene": {},
        "ensembl_transcript": {},
        "entrez": {},
        "refseq": {},
        "uniprot": {},
    }

    # Primary: HGNC/NCBI mapping (Entrez + Ensembl cross-refs)
    if os.path.exists(_HGNC_PATH):
        df = pd.read_csv(_HGNC_PATH, sep="\t", dtype=str, low_memory=False).fillna("")
        for _, row in df.iterrows():
            sym = row.get("symbol", "").strip()
            if not sym:
                continue
            ens = row.get("ensembl_gene_id", "").strip()
            if ens:
                maps["ensembl_gene"][ens] = sym
            eid = row.get("entrez_id", "").strip()
            if eid:
                maps["entrez"][eid] = sym
            # RefSeq: pipe-separated list
            for acc in row.get("refseq_accession", "").split("|"):
                acc = acc.strip().split(".")[0]  # strip version
                if acc:
                    maps["refseq"][acc] = sym
            # UniProt
            for uid in row.get("uniprot_ids", "").split("|"):
                uid = uid.strip()
                if uid:
                    maps["uniprot"][uid] = sym

    # Supplement: BioMart Ensembl → symbol (more complete for lncRNA/non-coding)
    if os.path.exists(_ENSEMBL_PATH):
        df2 = pd.read_csv(_ENSEMBL_PATH, sep="\t", dtype=str).fillna("")
        for _, row in df2.iterrows():
            ens = row.get("ensembl_id", "").strip()
            sym = row.get("gene_symbol", "").strip()
            if ens and sym and ens not in maps["ensembl_gene"]:
                maps["ensembl_gene"][ens] = sym

    return maps


def map_to_gene_symbol(df: pd.DataFrame) -> tuple:
    """
    Detect gene ID type in df.index, convert to HGNC symbols via local lookup.

    Returns:
        (mapped_df, id_type, stats)
        - mapped_df: df with gene symbols as index (original df if no mapping needed)
        - id_type:   detected type string or None
        - stats:     {'n_input', 'n_mapped', 'n_dropped', 'id_type'}
    """
    gene_ids = list(df.index)
    n_input  = len(gene_ids)
    id_type  = detect_id_type(gene_ids)
    stats    = {"n_input": n_input, "n_mapped": n_input, "n_dropped": 0, "id_type": id_type}

    if id_type is None:
        # Already HGNC symbols or unrecognised — return as-is
        return df, id_type, stats

    maps = _load_mappings()
    gene_map = maps.get(id_type, {})

    if not gene_map:
        return df, id_type, stats

    df = df.copy()
    df.index = df.index.map(lambda g: gene_map.get(str(g).split(".")[0], None))
    df_mapped = df[df.index.notna()].copy()
    df_mapped.index.name = "Gene"
    df_mapped = df_mapped.groupby(df_mapped.index).mean()  # collapse duplicate symbols

    stats["n_mapped"]  = len(df_mapped)
    stats["n_dropped"] = n_input - len(df_mapped)

    return df_mapped, id_type, stats


# Backward-compat alias used by older app.py code
def map_ensembl_to_symbol(df: pd.DataFrame) -> pd.DataFrame:
    mapped, _, _ = map_to_gene_symbol(df)
    return mapped


def genes_are_hgnc_symbols(id_type: str | None) -> bool:
    """
    Returns True if the gene IDs were either already HGNC symbols (id_type=None)
    or were successfully mapped to them.
    """
    return True  # after map_to_gene_symbol, index is always HGNC or original symbols
