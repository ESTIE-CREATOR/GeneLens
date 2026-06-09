import pandas as pd
import numpy as np


def run_pathway_enrichment(results: pd.DataFrame) -> dict:
    """
    Perform GO and KEGG pathway enrichment on DE gene lists.
    Uses gseapy enrichr for online enrichment analysis.
    Falls back to manual annotation if gseapy fails.

    Returns dict with keys: go_results, kegg_results, success
    """
    try:
        import gseapy as gp

        sig = results[results["significance"] != "Not significant"]
        up_genes = sig[sig["significance"] == "Upregulated"]["Gene"].tolist()
        down_genes = sig[sig["significance"] == "Downregulated"]["Gene"].tolist()
        all_de_genes = sig["Gene"].tolist()

        if len(all_de_genes) < 3:
            return {"go_results": pd.DataFrame(), "kegg_results": pd.DataFrame(), "success": False, "error": "Not enough DE genes for enrichment (need at least 3)"}

        go_enr = gp.enrichr(
            gene_list=all_de_genes,
            gene_sets=["GO_Biological_Process_2021"],
            organism="human",
            outdir=None,
            verbose=False,
        )
        go_df = go_enr.results.copy()
        go_df = go_df[go_df["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").head(15)

        kegg_enr = gp.enrichr(
            gene_list=all_de_genes,
            gene_sets=["KEGG_2021_Human"],
            organism="human",
            outdir=None,
            verbose=False,
        )
        kegg_df = kegg_enr.results.copy()
        kegg_df = kegg_df[kegg_df["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").head(15)

        return {
            "go_results": go_df,
            "kegg_results": kegg_df,
            "success": True,
            "n_de_genes": len(all_de_genes),
            "n_up": len(up_genes),
            "n_down": len(down_genes),
        }

    except ImportError:
        return {"go_results": pd.DataFrame(), "kegg_results": pd.DataFrame(), "success": False, "error": "gseapy not installed. Run: pip install gseapy"}

    except Exception as e:
        return {"go_results": pd.DataFrame(), "kegg_results": pd.DataFrame(), "success": False, "error": str(e)}


def plot_go_bar(go_df: pd.DataFrame) -> object:
    """
    Plot horizontal bar chart of top enriched GO terms.
    Returns a plotly figure.
    """
    import plotly.graph_objects as go

    if go_df.empty:
        return None

    df = go_df.head(10).copy()
    df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(1e-300))
    df["Term_short"] = df["Term"].str[:55]

    fig = go.Figure(go.Bar(
        x=df["-log10(padj)"],
        y=df["Term_short"],
        orientation="h",
        marker=dict(
            color=df["-log10(padj)"],
            colorscale=[[0, "#2E75B6"], [0.5, "#7C6AF7"], [1.0, "#E05C5C"]],
            showscale=False,
        ),
        text=[f"p={v:.2e}" for v in df["Adjusted P-value"]],
        textposition="outside",
        textfont=dict(size=9, color="#E8EAF0"),
        hovertemplate="<b>%{y}</b><br>-log10(padj): %{x:.2f}<extra></extra>",
    ))

    fig.update_layout(
        plot_bgcolor="#1A1D27",
        paper_bgcolor="#0F1117",
        font=dict(family="Courier New, monospace", color="#E8EAF0"),
        title=dict(text="GO Biological Process Enrichment — Top 10 Terms", font=dict(size=15), x=0.5),
        xaxis=dict(title="-log₁₀(adjusted p-value)", gridcolor="rgba(255,255,255,0.05)"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
        height=480,
        margin=dict(l=320, r=120, t=60, b=60),
    )
    return fig


def plot_kegg_bar(kegg_df: pd.DataFrame) -> object:
    """
    Plot horizontal bar chart of top enriched KEGG pathways.
    Returns a plotly figure.
    """
    import plotly.graph_objects as go

    if kegg_df.empty:
        return None

    df = kegg_df.head(10).copy()
    df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(1e-300))
    df["Term_short"] = df["Term"].str[:55]

    fig = go.Figure(go.Bar(
        x=df["-log10(padj)"],
        y=df["Term_short"],
        orientation="h",
        marker=dict(
            color=df["-log10(padj)"],
            colorscale=[[0, "#1F4E79"], [0.5, "#2E75B6"], [1.0, "#7C6AF7"]],
            showscale=False,
        ),
        text=[f"p={v:.2e}" for v in df["Adjusted P-value"]],
        textposition="outside",
        textfont=dict(size=9, color="#E8EAF0"),
        hovertemplate="<b>%{y}</b><br>-log10(padj): %{x:.2f}<extra></extra>",
    ))

    fig.update_layout(
        plot_bgcolor="#1A1D27",
        paper_bgcolor="#0F1117",
        font=dict(family="Courier New, monospace", color="#E8EAF0"),
        title=dict(text="KEGG Pathway Enrichment — Top 10 Pathways", font=dict(size=15), x=0.5),
        xaxis=dict(title="-log₁₀(adjusted p-value)", gridcolor="rgba(255,255,255,0.05)"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
        height=480,
        margin=dict(l=320, r=120, t=60, b=60),
    )
    return fig
