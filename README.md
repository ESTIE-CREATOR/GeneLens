# GeneLens 🧬
### Differential Expression · GO Enrichment · KEGG Pathways · ML Classification · AI Interpretation

A complete bioinformatics pipeline — from raw RNA-seq count matrices to pathway-level biological interpretation.

<img width="1919" height="949" alt="Image" src="https://github.com/user-attachments/assets/e5665890-5e5b-4216-9320-1a332603526f" />

---

## 🔬 Live App
👉 **[genelens.streamlit.app](https://genelens.streamlit.app)**

---

## What GeneLens Does

| Feature | Method |
|---------|--------|
| Differential Expression | Welch's t-test + Benjamini-Hochberg FDR correction |
| Volcano Plot | Interactive — labelled gene names, threshold lines |
| Heatmap | Z-scored expression across all samples |
| PCA | Sample clustering to verify group separation |
| GO Enrichment | GO_Biological_Process_2021 via Enrichr API |
| KEGG Pathways | KEGG_2021_Human via Enrichr API |
| ML Classification | Random Forest (200 trees, 5-fold CV, AUC + feature importance) |
| AI Interpretation | Anthropic Claude API — automated biological interpretation |

> GO and KEGG enrichment works for datasets with real gene names (e.g. SOD2, TP53).
> Large Ensembl ID datasets are gracefully skipped with a helpful message.

---

## Validated on Real Published Data

### Dataset 1 — COVID-19 PBMC RNA-seq (GSE152418, *Science* 2020)

**Scale:** 60,683 genes × 34 samples | 17 COVID-19 patients vs. 17 healthy controls

<!-- Add your COVID volcano screenshot URL here -->
<img width="959" height="474" alt="Image" src="https://github.com/user-attachments/assets/c3f17865-3161-43f2-b139-e6a2688458e9" />

*3,405 upregulated · 2,250 downregulated*

<img width="959" height="476" alt="Image" src="https://github.com/user-attachments/assets/fd47d0de-f499-4f62-9faf-d3c58d4ca24a" />

*Clean separation of COVID-19 patients from healthy controls on PC1*

**Key findings:**
- `CDC20`, `FOXM1`, `RRM2` downregulated → immune cell cycle arrest
- `CYP1B1` upregulated → oxidative stress response
- `TGFBR2`, `SETDB1` upregulated → immune evasion and epigenetic reprogramming

---

### Dataset 2 — Type 2 Diabetes Liver — Workflow Demonstration

**Purpose:** Demonstration of the full GeneLens analysis workflow and biological interpretation capabilities, connecting RNA-seq analysis to oxidative stress pathways studied in my undergraduate thesis.

**Structure:** 226 curated metabolic and oxidative stress genes × 13 samples | 5 Normal vs. 8 T2DM conditions

**Note:** This dataset was constructed as a biologically informed demonstration using gene expression directions and magnitudes derived from published T2DM transcriptomics literature (including Pihlajamäki et al. 2009, GSE15653). It is not a direct download from GEO — it is a curated demonstration dataset designed to illustrate the GeneLens pipeline on biologically relevant gene sets.

<!-- Add your T2DM screenshots here -->
<img width="1918" height="877" alt="Image" src="https://github.com/user-attachments/assets/d2528e19-c859-4b4f-a5ef-e8ab85319d89" />

*GPX1, GCLM, GSR downregulated — antioxidant enzymes from thesis research*

<img width="1919" height="943" alt="Image" src="https://github.com/user-attachments/assets/1f415041-4918-4d2a-932d-954a9172ff04" />

*Glutathione metabolic process — most significantly depleted GO term*

<img width="1919" height="948" alt="Image" src="https://github.com/user-attachments/assets/e79afb9c-c05f-4051-a77d-dbd09504edc8" />


**Key findings from demonstration:**
- `GPX1`, `GCLM`, `GSR` downregulated — antioxidant suppression
- `SOD2`, `CAT` strong downward trends — consistent with thesis biochemical findings
- `CYP2E1` upregulated → ROS generation
- **Top GO term:** Glucose homeostasis (GO:0042593, p=6.00e-07)
- **Top KEGG pathway:** Insulin signalling (p=1.86e-08)

> 📁 Full case study and methods: [genelens-t2dm-study](https://github.com/ESTIE-CREATOR/genelens-t2dm-study)

---

## Supported File Formats

| Format | Notes |
|--------|-------|
| `.csv` | Standard comma-separated count matrix |
| `.tsv` | Tab-separated count matrix |
| `.txt` | GEO Series Matrix files — auto-parsed |

**Expected structure:** rows = genes, columns = samples, first column = gene names.

---

## How to Run Locally

```bash
git clone https://github.com/ESTIE-CREATOR/GeneLens.git
cd GeneLens
python -m venv venv
source venv/bin/activate     # Mac/Linux
venv\Scripts\activate        # Windows
pip install -r requirements.txt
streamlit run app.py
```

Opens at `http://localhost:8501`

---

## Project Structure

```
GeneLens/
├── app.py                      ← Main Streamlit application
├── requirements.txt
├── .streamlit/config.toml      ← Dark theme
└── utils/
    ├── de_analysis.py          ← Welch t-test + BH FDR
    ├── visualisations.py       ← Volcano, heatmap, PCA, bar chart
    ├── ml_classifier.py        ← Random Forest + ROC curve
    ├── pathway_enrichment.py   ← GO and KEGG via Enrichr API
    ├── ai_interpretation.py    ← Claude API interpretation
    ├── file_parser.py          ← GEO Series Matrix + CSV/TSV parser
    └── data_generator.py       ← Demo dataset generator
```

---

## Tech Stack

| Layer | Tools |
|-------|-------|
| Frontend | Streamlit, Plotly |
| Analysis | Python, SciPy, pandas, NumPy |
| Pathway Enrichment | gseapy (Enrichr API) |
| Machine Learning | scikit-learn |
| AI Interpretation | Anthropic Claude API |
| Deployment | Streamlit Cloud |

---

## Related Research

📁 **T2DM Workflow Case Study:** [github.com/ESTIE-CREATOR/genelens-t2dm-study](https://github.com/ESTIE-CREATOR/genelens-t2dm-study)

---

## Built By

**Alabi Esther Oluwatimilehin**
BSc Biochemistry — First Class Honours, University of Medical Sciences, Nigeria
[github.com/ESTIE-CREATOR](https://github.com/ESTIE-CREATOR) · [linkedin.com/in/alabi-esther-essie](https://linkedin.com/in/alabi-esther-essie)

---
*MIT License*
