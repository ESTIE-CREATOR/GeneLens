import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


COLORS = {
    "upregulated": "#E05C5C",
    "downregulated": "#5C8FE0",
    "not_significant": "#AAAAAA",
    "background": "#0F1117",
    "surface": "#1A1D27",
    "text": "#E8EAF0",
    "accent": "#7C6AF7",
}

LAYOUT_BASE = dict(
    plot_bgcolor=COLORS["surface"],
    paper_bgcolor=COLORS["background"],
    font=dict(family="Courier New, monospace", color=COLORS["text"], size=12),
    margin=dict(l=60, r=40, t=60, b=60),
)


def plot_volcano(
    results: pd.DataFrame,
    pval_threshold: float = 0.05,
    fc_threshold: float = 1.5,
    label_top_n: int = 15
) -> go.Figure:
    """Interactive volcano plot."""
    log2_fc_thresh = np.log2(fc_threshold)

    color_map = {
        "Upregulated": COLORS["upregulated"],
        "Downregulated": COLORS["downregulated"],
        "Not significant": COLORS["not_significant"],
    }

    fig = go.Figure()

    for group, color in color_map.items():
        mask = results["significance"] == group
        subset = results[mask]
        fig.add_trace(go.Scatter(
            x=subset["log2FoldChange"],
            y=subset["-log10pval"],
            mode="markers",
            name=group,
            marker=dict(
                color=color,
                size=5 if group == "Not significant" else 7,
                opacity=0.6 if group == "Not significant" else 0.85,
                line=dict(width=0.3, color="white") if group != "Not significant" else dict(width=0),
            ),
            text=subset["Gene"],
            hovertemplate="<b>%{text}</b><br>log2FC: %{x:.3f}<br>-log10(p): %{y:.3f}<extra></extra>",
        ))

    # Label top DE genes
    top = results[results["significance"] != "Not significant"].head(label_top_n)
    for _, row in top.iterrows():
        fig.add_annotation(
            x=row["log2FoldChange"], y=row["-log10pval"],
            text=row["Gene"], showarrow=True,
            arrowhead=2, arrowsize=0.8, arrowwidth=1,
            arrowcolor="rgba(255,255,255,0.3)",
            font=dict(size=9, color=COLORS["text"]),
            bgcolor="rgba(20,20,30,0.7)",
            bordercolor="rgba(255,255,255,0.1)",
        )

    # Threshold lines
    neg_log_pval = -np.log10(pval_threshold)
    fig.add_hline(y=neg_log_pval, line_dash="dash", line_color="rgba(255,255,255,0.3)", line_width=1)
    fig.add_vline(x=log2_fc_thresh, line_dash="dash", line_color="rgba(255,255,255,0.3)", line_width=1)
    fig.add_vline(x=-log2_fc_thresh, line_dash="dash", line_color="rgba(255,255,255,0.3)", line_width=1)

    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text="Volcano Plot — Differential Expression", font=dict(size=16), x=0.5),
        xaxis=dict(title="log₂ Fold Change", gridcolor="rgba(255,255,255,0.05)", zeroline=True, zerolinecolor="rgba(255,255,255,0.15)"),
        yaxis=dict(title="-log₁₀(p-value)", gridcolor="rgba(255,255,255,0.05)"),
        legend=dict(bgcolor="rgba(255,255,255,0.05)", bordercolor="rgba(255,255,255,0.1)", borderwidth=1),
        height=520,
    )
    return fig


def plot_heatmap(
    df: pd.DataFrame,
    results: pd.DataFrame,
    control_cols: list,
    treated_cols: list,
    top_n: int = 40
) -> go.Figure:
    """Z-scored heatmap of top DE genes."""
    top_genes = results[results["significance"] != "Not significant"]["Gene"].head(top_n).tolist()
    if len(top_genes) < 5:
        top_genes = results["Gene"].head(top_n).tolist()

    subset = df.loc[top_genes, control_cols + treated_cols].copy().astype(float) + 1
    log_subset = np.log2(subset)

    # Z-score across samples
    zscored = ((log_subset.T - log_subset.mean(axis=1)) / log_subset.std(axis=1).clip(0.01)).T

    # Annotation: significance marker
    sig_map = results.drop_duplicates(subset="Gene").set_index("Gene")["significance"]
    gene_labels = [
        f"▲ {g}" if sig_map.get(g) == "Upregulated"
        else f"▼ {g}" if sig_map.get(g) == "Downregulated"
        else g
        for g in top_genes
    ]

    fig = go.Figure(data=go.Heatmap(
        z=zscored.values,
        x=control_cols + treated_cols,
        y=gene_labels,
        colorscale=[
            [0.0, "#3A6BC4"],
            [0.5, "#1A1D27"],
            [1.0, "#C43A3A"],
        ],
        zmid=0,
        colorbar=dict(
            title=dict(text="Z-score", font=dict(color=COLORS["text"])),
            tickfont=dict(color=COLORS["text"]),
            len=0.8,
        ),
        hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>",
    ))

    # Add vertical line between control and treated
    n_ctrl = len(control_cols)
    fig.add_vline(x=n_ctrl - 0.5, line_color="rgba(255,255,255,0.5)", line_width=2)

    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text=f"Heatmap — Top {len(top_genes)} DE Genes (Z-scored)", font=dict(size=16), x=0.5),
        xaxis=dict(tickangle=-45),
        yaxis=dict(autorange="reversed", tickfont=dict(size=9)),
        height=max(500, len(top_genes) * 16 + 100),
    )
    return fig


def plot_pca(
    df: pd.DataFrame,
    control_cols: list,
    treated_cols: list
) -> go.Figure:
    """PCA plot of samples."""
    all_cols = control_cols + treated_cols
    data = np.log2(df[all_cols].values.T + 1)

    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    n_components = min(3, data_scaled.shape[0], data_scaled.shape[1])
    if n_components < 2:
        return None
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(data_scaled)
    var_exp = pca.explained_variance_ratio_ * 100

    groups = ["Control"] * len(control_cols) + ["Treated"] * len(treated_cols)
    colors_list = [COLORS["downregulated"]] * len(control_cols) + [COLORS["upregulated"]] * len(treated_cols)

    fig = go.Figure()
    for group, color in [("Control", COLORS["downregulated"]), ("Treated", COLORS["upregulated"])]:
        mask = [g == group for g in groups]
        idx = [i for i, m in enumerate(mask) if m]
        fig.add_trace(go.Scatter(
            x=coords[idx, 0],
            y=coords[idx, 1],
            mode="markers+text",
            name=group,
            marker=dict(color=color, size=14, opacity=0.85, line=dict(width=1.5, color="white")),
            text=all_cols[idx[0]:idx[-1]+1] if group == "Control" else [all_cols[i] for i in idx],
            textposition="top center",
            textfont=dict(size=9),
            hovertemplate=f"<b>%{{text}}</b><br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra></extra>",
        ))

    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text="PCA — Sample Clustering", font=dict(size=16), x=0.5),
        xaxis=dict(title=f"PC1 ({var_exp[0]:.1f}% variance)", gridcolor="rgba(255,255,255,0.05)"),
        yaxis=dict(title=f"PC2 ({var_exp[1]:.1f}% variance)" if len(var_exp) > 1 else "PC2", gridcolor="rgba(255,255,255,0.05)"),
        legend=dict(bgcolor="rgba(255,255,255,0.05)"),
        height=480,
    )
    return fig


def plot_deg_bar(results: pd.DataFrame) -> go.Figure:
    """Bar chart summary of up/down/not-sig genes."""
    counts = results["significance"].value_counts()
    categories = ["Upregulated", "Downregulated", "Not significant"]
    values = [counts.get(c, 0) for c in categories]
    bar_colors = [COLORS["upregulated"], COLORS["downregulated"], COLORS["not_significant"]]

    fig = go.Figure(go.Bar(
        x=categories,
        y=values,
        marker_color=bar_colors,
        marker_line=dict(width=0),
        text=values,
        textposition="outside",
        textfont=dict(color=COLORS["text"], size=14),
        hovertemplate="%{x}: %{y} genes<extra></extra>",
    ))

    fig.update_layout(
        **LAYOUT_BASE,
        title=dict(text="DE Gene Summary", font=dict(size=16), x=0.5),
        xaxis=dict(showgrid=False),
        yaxis=dict(title="Number of Genes", gridcolor="rgba(255,255,255,0.05)"),
        showlegend=False,
        height=380,
    )
    return fig
