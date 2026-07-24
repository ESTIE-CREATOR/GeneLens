import os
import io
import gzip
import tempfile
import pandas as pd
import numpy as np
import requests

_GEO_CACHE = os.path.join(tempfile.gettempdir(), "genelens_geo_cache")

_MATRIX_KEYWORDS = ("count", "matrix", "fpkm", "tpm", "expression", "raw")
_SKIP_KEYWORDS = (
    "filelist", "readme", ".tar", ".pdf", ".xlsx", ".xls",
    ".bam", ".bai", ".bw", ".bigwig", ".idat", ".cel",
)

EXAMPLE_DATASETS = [
    {
        "accession": "GSE2565",
        "title": "Mouse Lung — LPS vs Saline",
        "description": "Affymetrix microarray · 24 samples · Mus musculus",
        "data_type": "microarray",
    },
    {
        "accession": "GSE4107",
        "title": "Colorectal Adenoma vs Normal Mucosa",
        "description": "Affymetrix microarray · 22 samples · Homo sapiens",
        "data_type": "microarray",
    },
    {
        "accession": "GSE96870",
        "title": "Mouse Airway Inflammation — RNA-seq Counts",
        "description": "RNA-seq integer counts · 22 samples · Mus musculus",
        "data_type": "rnaseq",
    },
]


def _candidate_supplementary_urls(gse) -> list:
    """Rank a GSE's supplementary files by how likely they are to be a counts/expression matrix."""
    urls = gse.metadata.get("supplementary_file", [])
    if isinstance(urls, str):
        urls = [urls]
    candidates = [u for u in urls if not any(s in u.lower() for s in _SKIP_KEYWORDS)]
    candidates.sort(key=lambda u: sum(kw in u.lower() for kw in _MATRIX_KEYWORDS), reverse=True)
    return candidates


def _download_and_parse_matrix(url: str):
    """Download one supplementary file and try to parse it as a genes × samples matrix."""
    https_url = url.replace("ftp://", "https://", 1)
    try:
        resp = requests.get(https_url, headers={"User-Agent": "Mozilla/5.0"}, timeout=120)
        resp.raise_for_status()
        raw = resp.content
    except Exception:
        return None

    if https_url.lower().endswith(".gz"):
        try:
            raw = gzip.decompress(raw)
        except OSError:
            return None

    text = raw.decode("utf-8", errors="replace")
    for sep in ("\t", ","):
        try:
            df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
        except Exception:
            continue
        if df.shape[1] >= 1 and df.shape[0] > 1:
            numeric = df.apply(pd.to_numeric, errors="coerce")
            if numeric.notna().mean().mean() > 0.9:
                return numeric
    return None


def _fetch_supplementary_matrix(gse):
    """Fall back to a GSE-level supplementary counts matrix when per-sample GSM tables are empty."""
    for url in _candidate_supplementary_urls(gse):
        df = _download_and_parse_matrix(url)
        if df is not None:
            return df
    return None


def fetch_geo_dataset(accession: str, cache_dir: str = _GEO_CACHE) -> tuple:
    """
    Download and parse a GEO dataset by accession number.

    Returns:
        (expr_df, sample_info_dict, is_counts)
        - expr_df: genes × samples DataFrame
        - sample_info_dict: {col_name: {'source': ..., 'characteristics': ...}}
        - is_counts: True if data looks like raw integer counts
    """
    try:
        import GEOparse
    except ImportError:
        raise ImportError(
            "GEOparse is not installed. Run: pip install GEOparse"
        )

    os.makedirs(cache_dir, exist_ok=True)
    accession = accession.strip().upper()
    if not accession.startswith("GSE"):
        raise ValueError(
            f"'{accession}' is not a valid GSE accession. "
            "Format: GSE followed by digits, e.g. GSE96870"
        )

    gse = GEOparse.get_GEO(geo=accession, destdir=cache_dir, silent=True)
    try:
        expr = gse.pivot_samples("VALUE")
    except (KeyError, ValueError):
        # RNA-seq datasets often have empty per-sample GSM tables (counts live
        # only in a GSE-level supplementary matrix file instead) — pivot_samples
        # crashes on these rather than returning an empty frame.
        expr = pd.DataFrame()

    if expr.empty:
        # RNA-seq series often store counts only in a GSE-level supplementary
        # matrix file rather than per-GSM tables — try downloading that instead.
        fallback = _fetch_supplementary_matrix(gse)
        if fallback is not None:
            expr = fallback

    if expr.empty:
        raise ValueError(
            f"No expression matrix found in {accession}. "
            "RNA-seq datasets often store counts in supplementary files — "
            "download the count matrix CSV from NCBI GEO and upload it directly."
        )

    expr = expr.dropna(how="all").fillna(0)

    # Map probe IDs → gene symbols when platform annotation is available
    platform_id = gse.metadata.get("platform_id", [None])[0]
    if platform_id and platform_id in gse.gpls:
        gpl = gse.gpls[platform_id]
        for col in ["Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol", "Gene_Symbol", "gene_assignment"]:
            if col in gpl.table.columns:
                symbol_series = (
                    gpl.table[col]
                    .fillna("")
                    .str.split("///")
                    .str[0]
                    .str.strip()
                )
                symbol_series.index = gpl.table.index
                valid = symbol_series[symbol_series != ""]
                expr.index = expr.index.map(lambda x: valid.get(x, x) if x in valid.index else x)
                expr = expr[~expr.index.duplicated(keep="first")]
                break

    # Build sample info and rename columns from GSM → sample title
    sample_info = {}
    rename_map = {}
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [gsm_name])[0]
        source = gsm.metadata.get("source_name_ch1", [""])[0]
        chars = "; ".join(gsm.metadata.get("characteristics_ch1", []))
        display = title
        sample_info[display] = {"source": source, "characteristics": chars}
        if gsm_name in expr.columns:
            rename_map[gsm_name] = display

    expr = expr.rename(columns=rename_map)

    # Heuristic: if median non-zero value > 50 and max > 1000, treat as raw counts
    flat = expr.values.flatten()
    nonzero = flat[flat > 0]
    is_counts = bool(len(nonzero) > 0 and np.percentile(nonzero, 50) > 50 and flat.max() > 1000)

    return expr, sample_info, is_counts
