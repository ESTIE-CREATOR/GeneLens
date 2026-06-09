import io
import pandas as pd


def parse_uploaded_file(file, already_log2: bool) -> pd.DataFrame:
    """Parse a gene expression matrix from CSV, TSV, or TXT.

    Detects separator from extension, falls back to the other common
    delimiter if the first attempt produces only one column. Skips lines
    starting with '#' (common in GEO supplementary text files).
    Converts Ensembl IDs to gene symbols when detected.
    """
    name = file.name.lower()
    raw = file.read()

    primary_sep = "\t" if (name.endswith(".tsv") or name.endswith(".txt")) else ","
    fallback_sep = "," if primary_sep == "\t" else "\t"

    df = None
    for sep in (primary_sep, fallback_sep):
        try:
            candidate = pd.read_csv(
                io.BytesIO(raw), sep=sep, index_col=0, comment="#"
            )
            if candidate.shape[1] >= 1:
                df = candidate
                break
        except Exception:
            continue

    if df is None:
        raise ValueError(
            "Could not parse the file as CSV or TSV. "
            "Ensure the first column contains gene names and remaining columns are samples."
        )

    df.index.name = "Gene"
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if not already_log2:
        df = df.astype(int)

    df = convert_ensembl_to_symbols(df)
    return df


def convert_ensembl_to_symbols(df: pd.DataFrame) -> pd.DataFrame:
    """
    If the DataFrame index looks like Ensembl IDs (ENSG...),
    convert them to gene symbols using MyGene.info API.
    Falls back gracefully if conversion fails.
    """
    index_vals = df.index.astype(str)
    is_ensembl = index_vals.str.startswith("ENSG").sum() > len(index_vals) * 0.5

    if not is_ensembl:
        return df  # already gene symbols, nothing to do

    try:
        import requests
        import json

        ensembl_ids = index_vals.tolist()
        # Query MyGene.info in batches of 1000
        symbol_map = {}
        batch_size = 1000
        for i in range(0, len(ensembl_ids), batch_size):
            batch = ensembl_ids[i:i + batch_size]
            response = requests.post(
                "https://mygene.info/v3/gene",
                data={
                    "ids": ",".join(batch),
                    "fields": "symbol",
                    "species": "human",
                },
                timeout=30,
            )
            if response.status_code == 200:
                for item in response.json():
                    if isinstance(item, dict) and "symbol" in item:
                        query = item.get("query", "")
                        symbol_map[query] = item["symbol"]

        if len(symbol_map) > 100:
            new_index = [symbol_map.get(g, g) for g in ensembl_ids]
            df = df.copy()
            df.index = new_index
            df.index.name = "Gene"
            # Drop rows where conversion failed (still look like ENSG...)
            df = df[~df.index.astype(str).str.startswith("ENSG")]

        return df

    except Exception:
        return df  # if anything fails, return original
