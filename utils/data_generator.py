import numpy as np
import pandas as pd


def generate_demo_data(n_genes: int = 500, n_samples: int = 12, seed: int = 42) -> pd.DataFrame:
    """
    Generate a realistic synthetic gene expression count matrix.
    Returns a DataFrame: rows = genes, columns = samples (6 control, 6 treated).
    """
    np.random.seed(seed)

    gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]

    # Add some real-sounding gene names for realism
    real_genes = [
        "TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "VEGFA", "TNF", "IL6",
        "CDKN2A", "PTEN", "RB1", "AKT1", "KRAS", "HRAS", "PIK3CA",
        "BCL2", "BAX", "CASP3", "MDM2", "CDK4", "CDK6", "CCND1",
        "STAT3", "NFKB1", "HIF1A", "SOD1", "SOD2", "CAT", "GPX1", "GSR"
    ]
    gene_names[:len(real_genes)] = real_genes

    # Base expression levels (log-normal distribution)
    base_expression = np.random.lognormal(mean=5, sigma=2, size=n_genes)

    # Control samples (6 replicates)
    control_data = np.random.negative_binomial(
        n=10,
        p=10 / (10 + base_expression[:, None]),
        size=(n_genes, 6)
    ).astype(float)

    # Treated samples — some genes are differentially expressed
    fold_changes = np.ones(n_genes)
    up_idx = np.random.choice(n_genes, size=40, replace=False)
    down_idx = np.random.choice(list(set(range(n_genes)) - set(up_idx)), size=40, replace=False)
    fold_changes[up_idx] = np.random.uniform(2.5, 8.0, size=40)
    fold_changes[down_idx] = np.random.uniform(0.1, 0.4, size=40)

    treated_expression = base_expression * fold_changes
    treated_data = np.random.negative_binomial(
        n=10,
        p=10 / (10 + treated_expression[:, None]),
        size=(n_genes, 6)
    ).astype(float)

    # Add small noise
    control_data += np.random.normal(0, 0.5, control_data.shape).clip(0)
    treated_data += np.random.normal(0, 0.5, treated_data.shape).clip(0)

    control_cols = [f"Control_{i+1}" for i in range(6)]
    treated_cols = [f"Treated_{i+1}" for i in range(6)]

    df = pd.DataFrame(
        np.hstack([control_data, treated_data]),
        index=gene_names,
        columns=control_cols + treated_cols
    )
    df = df.round(0).astype(int)
    df.index.name = "Gene"
    return df
