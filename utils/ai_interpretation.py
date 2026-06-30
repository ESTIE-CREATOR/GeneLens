import anthropic


def generate_ai_interpretation(
    api_key: str,
    summary_stats: dict,
    top_up_genes: list,
    top_down_genes: list,
    ml_auc: float,
    top_features: list,
) -> str:
    """
    Call the Anthropic API to generate a biological interpretation
    of the differential expression results.

    Args:
        api_key: Anthropic API key
        summary_stats: Dict with total, upregulated, downregulated counts
        top_up_genes: List of top upregulated gene names
        top_down_genes: List of top downregulated gene names
        ml_auc: ROC AUC score from the ML classifier
        top_features: Top gene names from the RF feature importance

    Returns:
        Interpretation text string
    """
    prompt = f"""You are an expert computational biologist. Analyse these gene expression results and write a concise, scientifically accurate biological interpretation.

DIFFERENTIAL EXPRESSION SUMMARY:
- Total genes analysed: {summary_stats['total_genes']}
- Upregulated genes: {summary_stats['upregulated']}
- Downregulated genes: {summary_stats['downregulated']}

TOP UPREGULATED GENES: {', '.join(top_up_genes[:10])}
TOP DOWNREGULATED GENES: {', '.join(top_down_genes[:10])}

MACHINE LEARNING:
- Random Forest classifier AUC: {ml_auc:.3f}
- Top discriminating genes: {', '.join(top_features[:8])}

Write a biological interpretation covering:
1. What biological processes or pathways might these DE genes be involved in?
2. What do the upregulated and downregulated genes suggest about the cellular response?
3. What is the significance of the ML AUC score?
4. What follow-up experiments or analyses would you recommend?
5. Mention 2-3 relevant published databases or tools (e.g. GSEA, STRING, KEGG) researchers should use for pathway enrichment.

Keep the tone professional and scientific. Write in clear paragraphs, not bullet points. Limit to 350 words."""

    client = anthropic.Anthropic(api_key=api_key)

    message = client.messages.create(
        model="claude-sonnet-4-6-20251001",
        max_tokens=600,
        messages=[{"role": "user", "content": prompt}]
    )

    return message.content[0].text
