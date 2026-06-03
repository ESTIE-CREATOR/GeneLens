import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
from sklearn.pipeline import Pipeline


COLORS = {
    "upregulated": "#E05C5C",
    "downregulated": "#5C8FE0",
    "accent": "#7C6AF7",
    "background": "#0F1117",
    "surface": "#1A1D27",
    "text": "#E8EAF0",
}


def run_ml_classification(
    df: pd.DataFrame,
    results: pd.DataFrame,
    control_cols: list,
    treated_cols: list,
    top_n_features: int = 30,
    n_splits: int = 5,
    random_state: int = 42
) -> dict:
    """
    Train a Random Forest classifier to distinguish Control vs Treated.
    Uses top DE genes as features.

    Returns:
        dict with: accuracy, std, feature_importances, roc_fig, cv_scores
    """
    # Select top DE genes as features, filtering to genes present in df
    sig_genes = results[results["significance"] != "Not significant"]["Gene"].head(top_n_features).tolist()
    if len(sig_genes) < 5:
        sig_genes = results["Gene"].head(top_n_features).tolist()
    sig_genes = [g for g in sig_genes if g in df.index]

    all_cols = control_cols + treated_cols
    labels = [0] * len(control_cols) + [1] * len(treated_cols)

    # Feature matrix: log2-transformed expression
    feature_df = np.log2(df.loc[sig_genes, all_cols].values.T.astype(float) + 1)
    X = feature_df
    y = np.array(labels)

    # Pipeline: scale → RF
    pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("rf", RandomForestClassifier(
            n_estimators=200,
            max_depth=5,
            random_state=random_state,
            class_weight="balanced"
        ))
    ])

    # Cross-validation — requires at least 2 samples per class
    n_cls_min = min(len(control_cols), len(treated_cols))
    if n_cls_min >= 2:
        cv = StratifiedKFold(
            n_splits=min(n_splits, n_cls_min),
            shuffle=True,
            random_state=random_state,
        )
        cv_scores = cross_val_score(pipeline, X, y, cv=cv, scoring="roc_auc")
    else:
        cv_scores = np.array([float("nan")])

    # Fit on all data for feature importance + ROC
    pipeline.fit(X, y)
    rf_model = pipeline.named_steps["rf"]

    importances = rf_model.feature_importances_
    feature_importance_df = pd.DataFrame({
        "Gene": sig_genes,
        "Importance": importances
    }).sort_values("Importance", ascending=False).reset_index(drop=True)

    # ROC curve (on training data — for illustration with small n)
    X_scaled = pipeline.named_steps["scaler"].transform(X)
    probs = rf_model.predict_proba(X_scaled)[:, 1]
    fpr, tpr, _ = roc_curve(y, probs)
    roc_auc = auc(fpr, tpr)

    roc_fig = _plot_roc(fpr, tpr, roc_auc, cv_scores.mean())
    importance_fig = _plot_feature_importance(feature_importance_df.head(20))

    return {
        "accuracy": float(np.nanmean(cv_scores)),
        "std": float(np.nanstd(cv_scores)),
        "cv_scores": cv_scores,
        "roc_auc": roc_auc,
        "feature_importances": feature_importance_df,
        "roc_fig": roc_fig,
        "importance_fig": importance_fig,
        "n_features": len(sig_genes),
    }


def _plot_roc(fpr, tpr, roc_auc: float, cv_auc: float) -> go.Figure:
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=fpr, y=tpr,
        mode="lines",
        name=f"ROC (AUC = {roc_auc:.3f})",
        line=dict(color=COLORS["accent"], width=2.5),
        fill="tozeroy",
        fillcolor="rgba(124,106,247,0.15)",
    ))
    fig.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1],
        mode="lines",
        name="Random (AUC = 0.5)",
        line=dict(color=COLORS["text"], width=1, dash="dash"),
        opacity=0.4,
    ))
    fig.update_layout(
        plot_bgcolor=COLORS["surface"],
        paper_bgcolor=COLORS["background"],
        font=dict(family="Courier New, monospace", color=COLORS["text"]),
        title=dict(
            text=f"ROC Curve  |  CV AUC: {cv_auc:.3f}",
            font=dict(size=15), x=0.5
        ),
        xaxis=dict(title="False Positive Rate", gridcolor="rgba(255,255,255,0.05)"),
        yaxis=dict(title="True Positive Rate", gridcolor="rgba(255,255,255,0.05)"),
        legend=dict(bgcolor="rgba(255,255,255,0.05)"),
        height=420,
        margin=dict(l=60, r=40, t=60, b=60),
    )
    return fig


def _plot_feature_importance(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure(go.Bar(
        x=df["Importance"],
        y=df["Gene"],
        orientation="h",
        marker=dict(
            color=df["Importance"],
            colorscale=[[0, COLORS["downregulated"]], [0.5, COLORS["accent"]], [1, COLORS["upregulated"]]],
            showscale=False,
        ),
        hovertemplate="<b>%{y}</b><br>Importance: %{x:.4f}<extra></extra>",
    ))
    fig.update_layout(
        plot_bgcolor=COLORS["surface"],
        paper_bgcolor=COLORS["background"],
        font=dict(family="Courier New, monospace", color=COLORS["text"]),
        title=dict(text="Feature Importance — Top Classifier Genes", font=dict(size=15), x=0.5),
        xaxis=dict(title="Mean Decrease in Impurity", gridcolor="rgba(255,255,255,0.05)"),
        yaxis=dict(autorange="reversed"),
        height=480,
        margin=dict(l=120, r=40, t=60, b=60),
    )
    return fig
