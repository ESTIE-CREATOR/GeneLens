import os
import pandas as pd
import streamlit as st

_MAPPING_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "ensembl_to_symbol.tsv")


@st.cache_data(show_spinner=False)
def _load_mapping() -> dict:
    path = os.path.abspath(_MAPPING_PATH)
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.dropna(subset=["ensembl_id", "gene_symbol"])
    df = df[df["gene_symbol"].str.strip() != ""]
    return dict(zip(df["ensembl_id"].str.strip(), df["gene_symbol"].str.strip()))


def is_ensembl_index(index: pd.Index) -> bool:
    sample = index[:20]
    return sum(str(g).startswith("ENSG") for g in sample) >= len(sample) // 2


def map_ensembl_to_symbol(df: pd.DataFrame) -> pd.DataFrame:
    """
    If df.index looks like Ensembl IDs (ENSG...), map them to gene symbols
    using the local BioMart flat file. Rows that cannot be mapped are dropped.
    Returns df with gene-symbol index, or the original df if no mapping needed.
    """
    if not is_ensembl_index(df.index):
        return df

    gene_map = _load_mapping()
    if not gene_map:
        return df

    df = df.copy()
    df.index = df.index.map(lambda g: gene_map.get(str(g).split(".")[0], None))
    df = df[df.index.notna()]
    df.index.name = "Gene"

    # Collapse duplicate symbols by taking the mean across samples
    df = df.groupby(df.index).mean()
    return df
