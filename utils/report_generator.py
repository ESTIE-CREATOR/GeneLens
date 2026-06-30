import io
import base64
import datetime
import pandas as pd
import plotly.io as pio


def _fig_to_b64_png(fig, width: int = 900, height: int = 480) -> str:
    try:
        import kaleido  # noqa: F401 — just check it's available
        png = pio.to_image(fig, format="png", width=width, height=height, scale=1.5)
        return base64.b64encode(png).decode()
    except Exception:
        return ""


def _table_html(df: pd.DataFrame, fmt: dict | None = None) -> str:
    if fmt:
        df = df.copy()
        for col, fmtstr in fmt.items():
            if col in df.columns:
                df[col] = df[col].map(lambda x: fmtstr.format(x))
    return df.to_html(index=False, border=0, escape=True, classes="gentable")


def generate_html_report(
    dataset_label: str,
    summary: dict,
    results: pd.DataFrame,
    ml_results: dict,
    top_degs: pd.DataFrame,
    fig_volcano,
    fig_heatmap,
    fig_pca,
    fig_bar,
    interpretation: str = "",
) -> bytes:
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    # Try static PNG embeds (needs kaleido); fall back to inline HTML charts
    use_static = False
    try:
        import kaleido  # noqa: F401
        use_static = True
    except ImportError:
        pass

    def chart_block(fig, w=900, h=480):
        if use_static:
            b64 = _fig_to_b64_png(fig, w, h)
            if b64:
                return f'<img src="data:image/png;base64,{b64}" style="width:100%;max-width:{w}px;margin:12px 0;">'
        return pio.to_html(fig, full_html=False, include_plotlyjs="cdn", config={"displayModeBar": False})

    v_html = chart_block(fig_volcano, 900, 500)
    b_html = chart_block(fig_bar, 420, 380)
    hm_html = chart_block(fig_heatmap, 900, 600)
    pca_html = chart_block(fig_pca, 700, 480) if fig_pca else "<p style='color:#aaa'>PCA not available (need ≥2 samples per group).</p>"
    _placeholder = "<p style='color:#aaa;text-align:center'>Not available.</p>"
    roc_html = chart_block(ml_results["roc_fig"], 500, 420) if ml_results.get("roc_fig") else _placeholder
    imp_html = chart_block(ml_results["importance_fig"], 600, 480) if ml_results.get("importance_fig") else _placeholder

    deg_table = _table_html(
        top_degs[["Gene", "log2FoldChange", "pvalue", "padj", "significance"]].head(30),
        fmt={"log2FoldChange": "{:.3f}", "pvalue": "{:.2e}", "padj": "{:.2e}"},
    )

    interp_block = ""
    if interpretation:
        interp_block = f"""
        <div class="section"><h2>AI Biological Interpretation</h2>
        <div class="ai-box">{interpretation}</div></div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GeneLens Report — {dataset_label}</title>
<style>
  body {{font-family: 'Segoe UI', Arial, sans-serif; background:#0f1117; color:#e8eaf0; margin:0; padding:0;}}
  .page {{max-width:1100px; margin:0 auto; padding:40px 32px;}}
  h1 {{font-size:2rem; background:linear-gradient(90deg,#7C6AF7,#5C8FE0,#E05C5C);
       -webkit-background-clip:text; -webkit-text-fill-color:transparent; margin-bottom:4px;}}
  h2 {{color:#7C6AF7; font-size:1.1rem; letter-spacing:0.08em; text-transform:uppercase;
       border-bottom:1px solid rgba(124,106,247,0.3); padding-bottom:6px; margin-top:36px;}}
  .meta {{color:#6b7280; font-size:0.82rem; margin-bottom:32px;}}
  .cards {{display:flex; gap:16px; flex-wrap:wrap; margin:16px 0 24px;}}
  .card {{background:#1a1d27; border-radius:10px; padding:18px 26px; flex:1; min-width:130px; text-align:center;}}
  .card-val {{font-size:2rem; font-weight:700;}}
  .card-lbl {{font-size:0.78rem; color:#9ca3af; margin-top:4px;}}
  .section {{margin-top:40px;}}
  .two-col {{display:grid; grid-template-columns:2fr 1fr; gap:16px; align-items:start;}}
  .ai-box {{background:#1a1d27; border-left:3px solid #7C6AF7; padding:16px 20px;
            border-radius:0 8px 8px 0; line-height:1.7; font-size:0.9rem;}}
  table {{width:100%; border-collapse:collapse; font-size:0.83rem; margin-top:12px;}}
  th {{background:#1a1d27; color:#7C6AF7; padding:8px 12px; text-align:left; letter-spacing:0.04em;}}
  td {{padding:7px 12px; border-bottom:1px solid rgba(255,255,255,0.05); color:#d1d5db;}}
  tr:hover td {{background:rgba(124,106,247,0.06);}}
  .footer {{color:#4b5563; font-size:0.75rem; text-align:center; margin-top:60px; padding-top:20px;
            border-top:1px solid rgba(255,255,255,0.06);}}
  @media print {{body{{background:#fff;color:#111;}} .card,.ai-box{{background:#f5f5f5;}} h1{{color:#333;}} h2{{color:#333;}}}}
</style>
</head>
<body>
<div class="page">
  <h1>GeneLens Analysis Report</h1>
  <div class="meta">Dataset: {dataset_label} &nbsp;·&nbsp; Generated: {now} &nbsp;·&nbsp; <a href="https://genelens.streamlit.app" style="color:#5C8FE0">genelens.streamlit.app</a></div>

  <div class="section"><h2>Summary</h2>
  <div class="cards">
    <div class="card"><div class="card-val" style="color:#e8eaf0">{summary['total_genes']:,}</div><div class="card-lbl">Total Genes</div></div>
    <div class="card"><div class="card-val" style="color:#E05C5C">{summary['upregulated']:,}</div><div class="card-lbl">Upregulated</div></div>
    <div class="card"><div class="card-val" style="color:#5C8FE0">{summary['downregulated']:,}</div><div class="card-lbl">Downregulated</div></div>
    <div class="card"><div class="card-val" style="color:#6B7280">{summary['not_significant']:,}</div><div class="card-lbl">Not Significant</div></div>
    <div class="card"><div class="card-val" style="color:#7C6AF7">{ml_results['accuracy']:.3f}</div><div class="card-lbl">ML CV AUC</div></div>
  </div></div>

  <div class="section"><h2>Differential Expression</h2>
  <div class="two-col">{v_html}{b_html}</div></div>

  <div class="section"><h2>Heatmap — Top DE Genes</h2>{hm_html}</div>

  <div class="section"><h2>PCA — Sample Clustering</h2>{pca_html}</div>

  <div class="section"><h2>Machine Learning Classification</h2>
  <div class="two-col">{roc_html}{imp_html}</div>
  <p style="text-align:center;color:#9ca3af;font-size:0.85rem;">Cross-Validation AUC: <b style="color:#7C6AF7">{ml_results['accuracy']:.3f} ± {ml_results['std']:.3f}</b> &nbsp;|&nbsp; Features: <b>{ml_results['n_features']}</b> genes</p>
  </div>

  <div class="section"><h2>Top DE Genes</h2>{deg_table}</div>

  {interp_block}

  <div class="footer">
    Generated by <b>GeneLens</b> — Real Gene Expression Analysis<br>
    <a href="https://genelens.streamlit.app" style="color:#5C8FE0">genelens.streamlit.app</a> &nbsp;·&nbsp;
    To save as PDF: File → Print → Save as PDF
  </div>
</div>
</body>
</html>"""

    return html.encode("utf-8")
