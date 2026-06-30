"""
Gene ID detection and mapping to HGNC gene symbols.

Supported input formats:
  - HGNC symbol (TP53, BRCA1)          — passed through unchanged
  - Ensembl gene ID  (ENSG00000141510)  — mapped via hgnc_mapping.tsv + ensembl_to_symbol.tsv
  - Ensembl transcript (ENST...)        — detected, Entrez lookup attempted
  - Entrez / NCBI gene ID (7157)        — mapped via hgnc_mapping.tsv
  - RefSeq accession (NM_000546)        — regex detection, best-effort mapping
  - UniProt accession (P04637)          — regex detection, best-effort mapping

Mapping files are committed static assets — zero network calls at runtime.
"""

import os
import re
import pandas as pd
import streamlit as st

_DATA_DIR     = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
_HGNC_PATH    = os.path.join(_DATA_DIR, "hgnc_mapping.tsv")
_ENSEMBL_PATH = os.path.join(_DATA_DIR, "ensembl_to_symbol.tsv")

# Ordered: most-specific patterns first so Entrez (pure digits) doesn't shadow others
_PATTERNS = [
    ("ensembl_gene",       re.compile(r"^ENSG\d{6,}")),
    ("ensembl_transcript", re.compile(r"^ENST\d{6,}")),
    ("refseq",             re.compile(r"^(?:NM_|NR_|XM_|XR_|NP_)\d+")),
    ("uniprot",            re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])")),
    ("entrez",             re.compile(r"^\d+$")),
]


def detect_id_type(gene_ids, sample_size: int = 100):
    """
    Sample gene IDs and return the detected format, or None if already HGNC symbols.

    Returns one of: 'ensembl_gene', 'ensembl_transcript', 'refseq', 'uniprot', 'entrez', None
    """
    ids = [str(g).strip() for g in gene_ids[:sample_size] if str(g).strip() not in ("", "nan")]
    if not ids:
        return None

    n = len(ids)
    for name, pat in _PATTERNS:
        frac = sum(1 for g in ids if pat.match(g.split(".")[0])) / n
        if frac >= 0.80:
            return name

    # Mixed column: accept 50% majority
    for name, pat in _PATTERNS:
        frac = sum(1 for g in ids if pat.match(g.split(".")[0])) / n
        if frac >= 0.50:
            return name

    return None  # treat as HGNC symbols


@st.cache_data(show_spinner=False)
def _load_mappings() -> dict:
    """
    Load static mapping files into lookup dicts.
    Uses fully vectorised pandas — no iterrows over 193k rows.
    Result is cached once per server process, shared across all users.
    """
    maps: dict = {k: {} for k in ("ensembl_gene", "ensembl_transcript", "entrez", "refseq", "uniprot")}

    # ── Primary: NCBI gene_info (Entrez + Ensembl cross-refs, 193k genes) ──
    if os.path.exists(_HGNC_PATH):
        df = pd.read_csv(_HGNC_PATH, sep="\t", dtype=str, low_memory=False).fillna("")
        sym = df["symbol"]

        # Ensembl (single value per row)
        ens_mask = df["ensembl_gene_id"].str.len() > 0
        maps["ensembl_gene"].update(
            dict(zip(df.loc[ens_mask, "ensembl_gene_id"], sym[ens_mask]))
        )

        # Entrez (single value per row)
        eid_mask = df["entrez_id"].str.len() > 0
        maps["entrez"].update(
            dict(zip(df.loc[eid_mask, "entrez_id"], sym[eid_mask]))
        )

        # RefSeq (pipe-separated — explode)
        rs_mask = df["refseq_accession"].str.len() > 0
        if rs_mask.any():
            rs = df.loc[rs_mask, ["symbol", "refseq_accession"]].copy()
            rs = rs.assign(refseq_accession=rs["refseq_accession"].str.split("|")).explode("refseq_accession")
            rs["refseq_accession"] = rs["refseq_accession"].str.strip().str.split(".").str[0]
            rs = rs[rs["refseq_accession"].str.len() > 0]
            maps["refseq"].update(dict(zip(rs["refseq_accession"], rs["symbol"])))

        # UniProt (pipe-separated — explode)
        up_mask = df["uniprot_ids"].str.len() > 0
        if up_mask.any():
            up = df.loc[up_mask, ["symbol", "uniprot_ids"]].copy()
            up = up.assign(uniprot_ids=up["uniprot_ids"].str.split("|")).explode("uniprot_ids")
            up["uniprot_ids"] = up["uniprot_ids"].str.strip()
            up = up[up["uniprot_ids"].str.len() > 0]
            maps["uniprot"].update(dict(zip(up["uniprot_ids"], up["symbol"])))

    # ── Supplement: BioMart (49k entries, fills in lncRNA / non-coding not in NCBI) ──
    if os.path.exists(_ENSEMBL_PATH):
        df2 = pd.read_csv(_ENSEMBL_PATH, sep="\t", dtype=str).fillna("")
        existing = set(maps["ensembl_gene"])
        new_mask = (
            (df2["ensembl_id"].str.len() > 0) &
            (df2["gene_symbol"].str.len() > 0) &
            (~df2["ensembl_id"].isin(existing))
        )
        maps["ensembl_gene"].update(
            dict(zip(df2.loc[new_mask, "ensembl_id"], df2.loc[new_mask, "gene_symbol"]))
        )

    return maps


def map_to_gene_symbol(df: pd.DataFrame) -> tuple:
    """
    Detect gene ID type in df.index, convert to HGNC symbols via local lookup.

    Returns:
        (mapped_df, id_type, stats)
        - mapped_df: df with HGNC symbols as index (original df if no conversion needed)
        - id_type:   detected format string, or None
        - stats:     dict with n_input, n_mapped, n_dropped, id_type
    """
    gene_ids = list(df.index)
    n_input  = len(gene_ids)
    id_type  = detect_id_type(gene_ids)
    stats    = {"n_input": n_input, "n_mapped": n_input, "n_dropped": 0, "id_type": id_type}

    if id_type is None:
        return df, id_type, stats

    maps     = _load_mappings()
    gene_map = maps.get(id_type, {})

    if not gene_map:
        return df, id_type, stats

    df = df.copy()
    df.index = df.index.map(lambda g: gene_map.get(str(g).split(".")[0], None))
    df_mapped = df[df.index.notna()].copy()
    df_mapped.index.name = "Gene"
    df_mapped = df_mapped.groupby(df_mapped.index).mean()

    stats["n_mapped"]  = len(df_mapped)
    stats["n_dropped"] = n_input - len(df_mapped)

    return df_mapped, id_type, stats


def map_ensembl_to_symbol(df: pd.DataFrame) -> pd.DataFrame:
    """Backward-compat alias."""
    mapped, _, _ = map_to_gene_symbol(df)
    return mapped
