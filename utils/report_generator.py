import io
import shutil
import datetime
from xml.sax.saxutils import escape
import pandas as pd
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from utils.pathway_enrichment import plot_go_bar, plot_kegg_bar

ACCENT = colors.HexColor("#7C6AF7")
DARK = colors.HexColor("#1a1d27")
MUTED = colors.HexColor("#6b7280")
FOOTER_COLOR = colors.HexColor("#374151")
ERROR_COLOR_HEX = "#c0392b"

_styles = getSampleStyleSheet()
_TITLE = ParagraphStyle("GLTitle", parent=_styles["Title"], textColor=ACCENT)
_H2 = ParagraphStyle("GLH2", parent=_styles["Heading2"], textColor=ACCENT, spaceBefore=14)
_H3 = ParagraphStyle("GLH3", parent=_styles["Heading3"], textColor=ACCENT, spaceBefore=10)
_META = ParagraphStyle("GLMeta", parent=_styles["Normal"], textColor=MUTED, fontSize=9)
_NOTE = ParagraphStyle("GLNote", parent=_styles["Normal"], alignment=TA_CENTER, fontSize=9, textColor=MUTED)
_FOOTER = ParagraphStyle("GLFooter", parent=_styles["Normal"], alignment=TA_CENTER, textColor=FOOTER_COLOR)
_CELL = ParagraphStyle("GLCell", parent=_styles["Normal"], fontSize=8, leading=10)

_chrome_checked = False

# Prefer a system-installed Chromium (e.g. `apt install chromium` via packages.txt)
# over kaleido's own downloaded "Chrome for Testing" copy. On Debian hosts like
# Streamlit Cloud, apt resolves the correct shared-library dependency names for
# whatever release is running; kaleido's standalone download has no such guarantee
# and can crash on launch ("browser closed immediately") if libs are missing.
_SYSTEM_CHROME_PATH = (
    shutil.which("chromium")
    or shutil.which("chromium-browser")
    or shutil.which("google-chrome")
    or shutil.which("google-chrome-stable")
)


def ensure_kaleido_chrome():
    """
    Fallback for hosts without a system Chromium: download kaleido's own
    "Chrome for Testing" once per process. Safe to call repeatedly.
    """
    global _chrome_checked
    if _chrome_checked or _SYSTEM_CHROME_PATH:
        return
    _chrome_checked = True
    try:
        import kaleido
        kaleido.get_chrome_sync()
    except Exception:
        pass


def render_chart_png(fig, width: int = 900, height: int = 500, scale: float = 2) -> bytes:
    """Render a plotly figure to PNG bytes, preferring the system Chromium if present."""
    import kaleido
    fig_dict = fig.to_dict() if hasattr(fig, "to_dict") else fig
    kopts = {"path": _SYSTEM_CHROME_PATH} if _SYSTEM_CHROME_PATH else {}
    return kaleido.calc_fig_sync(
        fig_dict,
        opts=dict(format="png", width=width, height=height, scale=scale),
        kopts=kopts,
    )


def _chart_flowable(fig, width_in: float, w_px: int = 900, h_px: int = 500, not_available_msg: str = "Not available."):
    """Return a single flowable for a chart: the image, or a visible error/reason note."""
    if fig is None:
        return Paragraph(not_available_msg, _NOTE)
    try:
        png = render_chart_png(fig, width=w_px, height=h_px, scale=2)
    except Exception as e:
        return Paragraph(
            f'<font color="{ERROR_COLOR_HEX}">Chart could not be rendered: {escape(str(e))[:200]}</font>',
            _styles["Normal"],
        )
    height_in = width_in * (h_px / w_px)
    img = Image(io.BytesIO(png), width=width_in * inch, height=height_in * inch)
    img.hAlign = "CENTER"
    return img


def _add_chart(story: list, fig, width_in: float, w_px: int = 900, h_px: int = 500, not_available_msg: str = "Not available."):
    """Append a chart image to the story (centered, with breathing room), or a visible error note."""
    story.append(_chart_flowable(fig, width_in, w_px, h_px, not_available_msg))
    story.append(Spacer(1, 14))


def _side_by_side(story: list, left, right, col_widths=(3.6 * inch, 3.1 * inch)):
    """Lay two chart flowables out side by side in a borderless table, with a spacer after."""
    t = Table([[left, right]], colWidths=list(col_widths), hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    story.append(t)
    story.append(Spacer(1, 14))


def _summary_table(summary: dict, ml_results: dict) -> Table:
    labels = ["Total Genes", "Upregulated", "Downregulated", "Not Significant", "ML CV AUC"]
    values = [
        f"{summary['total_genes']:,}", f"{summary['upregulated']:,}",
        f"{summary['downregulated']:,}", f"{summary['not_significant']:,}",
        f"{ml_results['accuracy']:.3f}",
    ]
    t = Table([labels, values], hAlign="LEFT", colWidths=[1.42 * inch] * 5)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), ACCENT),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#dddddd")),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return t


def _deg_table(top_degs: pd.DataFrame) -> Table:
    cols = ["Gene", "log2FoldChange", "pvalue", "padj", "significance"]
    deg = top_degs[cols].head(30)
    rows = [cols] + [
        [str(r["Gene"]), f"{r['log2FoldChange']:.3f}", f"{r['pvalue']:.2e}", f"{r['padj']:.2e}", str(r["significance"])]
        for _, r in deg.iterrows()
    ]
    t = Table(rows, hAlign="LEFT", repeatRows=1, colWidths=[1.3 * inch, 1.1 * inch, 1.1 * inch, 1.1 * inch, 1.6 * inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), ACCENT),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#dddddd")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f5f7")]),
    ]))
    return t


def _enrichment_table(df: pd.DataFrame) -> Table:
    rows = [["Term", "Overlap", "Adjusted P-value"]]
    for _, r in df.head(15).iterrows():
        rows.append([
            Paragraph(escape(str(r["Term"])), _CELL),
            str(r["Overlap"]),
            f"{r['Adjusted P-value']:.2e}",
        ])
    t = Table(rows, hAlign="LEFT", repeatRows=1, colWidths=[4.2 * inch, 1.1 * inch, 1.4 * inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), ACCENT),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#dddddd")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f5f7")]),
    ]))
    return t


def generate_pdf_report(
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
    enr_results: dict | None = None,
) -> bytes:
    ensure_kaleido_chrome()
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf, pagesize=LETTER,
        topMargin=0.6 * inch, bottomMargin=0.6 * inch,
        leftMargin=0.6 * inch, rightMargin=0.6 * inch,
    )

    story = [
        Paragraph("GeneLens Analysis Report", _TITLE),
        Paragraph(f"Dataset: {escape(dataset_label)} &nbsp;&nbsp;|&nbsp;&nbsp; Generated: {now}", _META),
        Spacer(1, 12),
        Paragraph("Summary", _H2),
        _summary_table(summary, ml_results),
        Spacer(1, 10),
        Paragraph("Differential Expression", _H2),
    ]

    _side_by_side(
        story,
        _chart_flowable(fig_volcano, width_in=4.0, w_px=900, h_px=500),
        _chart_flowable(fig_bar, width_in=2.7, w_px=420, h_px=380),
    )

    story.append(Paragraph("Heatmap — Top DE Genes", _H2))
    _add_chart(story, fig_heatmap, width_in=6.5, w_px=900, h_px=600)

    story.append(Paragraph("PCA — Sample Clustering", _H2))
    _add_chart(story, fig_pca, width_in=5.5, w_px=700, h_px=480,
               not_available_msg="PCA not available (need ≥2 samples per group).")

    story.append(Paragraph("Machine Learning Classification", _H2))
    if ml_results.get("roc_fig") is not None:
        _side_by_side(
            story,
            _chart_flowable(ml_results.get("roc_fig"), width_in=3.6, w_px=500, h_px=420),
            _chart_flowable(ml_results.get("importance_fig"), width_in=3.1, w_px=600, h_px=480),
        )
        story.append(Paragraph(
            f"Cross-Validation AUC: {ml_results['accuracy']:.3f} ± {ml_results['std']:.3f} "
            f"&nbsp;|&nbsp; Features: {ml_results['n_features']} genes",
            _NOTE,
        ))
    else:
        ml_error = ml_results.get("error")
        if ml_error:
            story.append(Paragraph(
                f'<font color="{ERROR_COLOR_HEX}">ML classification did not run: {escape(str(ml_error))[:300]}</font>',
                _styles["Normal"],
            ))
        else:
            story.append(Paragraph("ML classification not available.", _NOTE))
    story.append(Spacer(1, 10))

    story.append(PageBreak())
    story.append(Paragraph("Top DE Genes", _H2))
    story.append(_deg_table(top_degs))

    go_df = (enr_results or {}).get("go_results")
    kegg_df = (enr_results or {}).get("kegg_results")
    if (go_df is not None and not go_df.empty) or (kegg_df is not None and not kegg_df.empty):
        story.append(PageBreak())
        story.append(Paragraph("Pathway Enrichment — GO & KEGG", _H2))
        if go_df is not None and not go_df.empty:
            story.append(Paragraph("GO Biological Process", _H3))
            _add_chart(story, plot_go_bar(go_df), width_in=6.5, w_px=900, h_px=480)
            story.append(Spacer(1, 6))
            story.append(_enrichment_table(go_df))
        if kegg_df is not None and not kegg_df.empty:
            story.append(Paragraph("KEGG Pathways", _H3))
            _add_chart(story, plot_kegg_bar(kegg_df), width_in=6.5, w_px=900, h_px=480)
            story.append(Spacer(1, 6))
            story.append(_enrichment_table(kegg_df))

    if interpretation:
        story.append(Spacer(1, 14))
        story.append(Paragraph("AI Biological Interpretation", _H2))
        story.append(Paragraph(escape(interpretation).replace("\n", "<br/>"), _styles["Normal"]))

    story.append(Spacer(1, 30))
    story.append(Paragraph("<b>Generated by GeneLens — Real Gene Expression Analysis</b>", _FOOTER))
    story.append(Paragraph("<i>genelens.streamlit.app</i>", _FOOTER))

    doc.build(story)
    return buf.getvalue()
