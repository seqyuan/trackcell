"""Optional XGBoost enhancements for SCT integration (Seurat v4 RPCA remains default)."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import PCA
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import LabelEncoder

from .._sctransform import _as_genes_by_cells, get_residuals
from .prep_sct import _unpack_vst_out

try:
    from xgboost import XGBClassifier, XGBRegressor
except ImportError:  # pragma: no cover - optional at install time
    XGBClassifier = None
    XGBRegressor = None


def xgboost_available() -> bool:
    return XGBClassifier is not None


def _require_xgboost() -> None:
    if not xgboost_available():
        raise ImportError(
            "xgboost is required for XGBoost enhancement modes. "
            "Install with: pip install 'trackcell[ml]' or pip install xgboost"
        )


def _fit_batch_classifier(
    matrix: np.ndarray,
    y: np.ndarray,
    *,
    n_estimators: int,
    max_depth: int,
    learning_rate: float,
    subsample: float,
    colsample_bytree: float,
    seed: int,
):
    """Train an XGBoost batch classifier with a valid objective for n_classes."""
    n_classes = len(np.unique(y))
    common = dict(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        subsample=subsample,
        colsample_bytree=colsample_bytree,
        random_state=seed,
        n_jobs=1,
        verbosity=0,
    )
    if n_classes == 2:
        return XGBClassifier(objective="binary:logistic", **common)
    return XGBClassifier(objective="multi:softprob", num_class=n_classes, **common)


def _variance_candidate_genes(
    adata: AnnData,
    *,
    batch_key: str,
    sct_key: str,
    n_candidates: int,
) -> list[str]:
    uns = adata.uns.get(sct_key, {})
    batch_models = uns.get("batch_models", {})
    if not batch_models:
        raise ValueError(f"Missing {sct_key} batch_models; run sctransform with batch_key first.")
    rv_frames = []
    for batch_label, entry in batch_models.items():
        gene_attr = pd.DataFrame(entry.get("gene_attr", {})).T
        if "residual_variance" not in gene_attr.columns:
            continue
        rv_frames.append(gene_attr["residual_variance"].rename(batch_label))
    if not rv_frames:
        raise ValueError(f"No residual variance found in {sct_key} batch_models.")
    combined = pd.concat(rv_frames, axis=1)
    score = combined.max(axis=1, skipna=True).sort_values(ascending=False)
    return score.index[: min(n_candidates, len(score))].tolist()


def _residual_matrix(
    adata: AnnData,
    *,
    batch_key: str,
    sct_key: str,
    genes: Sequence[str],
    layer: Optional[str] = None,
) -> np.ndarray:
    """Cells x genes Pearson/deviance residuals for feature selection."""
    uns = adata.uns.get(sct_key, {})
    batch_models = uns["batch_models"]
    batches = adata.obs[batch_key].astype(str)
    genes = list(genes)
    out = np.zeros((adata.n_obs, len(genes)), dtype=np.float64)
    gene_to_col = {g: i for i, g in enumerate(genes)}

    for batch_label in batches.unique():
        mask = batches == batch_label
        entry = batch_models[batch_label]
        vst_out = _unpack_vst_out(entry)
        if layer is not None:
            counts = adata[mask].layers[layer]
        else:
            counts = adata[mask].X
        umi = _as_genes_by_cells(
            counts,
            adata.var_names.to_list(),
            adata.obs_names[mask].to_list(),
        )
        overlap = pd.Index(genes).intersection(vst_out["model_pars_fit"].index)
        if len(overlap) == 0:
            continue
        res = get_residuals(
            vst_out,
            umi,
            overlap,
            adata.var_names,
            cell_attr=vst_out["cell_attr"],
        )
        rows = np.where(mask.to_numpy())[0]
        for j, gene in enumerate(overlap):
            out[rows, gene_to_col[gene]] = res.iloc[j].to_numpy()
    return out


def select_integration_features_xgb(
    adata: AnnData,
    *,
    batch_key: str,
    sct_key: str = "sct",
    n_features: int = 2000,
    layer: Optional[str] = None,
    candidate_multiplier: int = 5,
    n_estimators: int = 200,
    max_depth: int = 6,
    learning_rate: float = 0.05,
    subsample: float = 0.8,
    colsample_bytree: float = 0.8,
    seed: int = 0,
) -> list[str]:
    """
    Optional enhancement: rank anchor genes by XGBoost batch-discrimination importance.

    Prefilters with the Seurat-style cross-batch residual variance screen, then
    trains a batch classifier on per-cell SCT residuals and returns the top
    ``n_features`` by gain importance. Does not replace reciprocal PCA anchor finding.
    """
    _require_xgboost()
    n_candidates = min(max(n_features * candidate_multiplier, n_features), 10000)
    candidates = _variance_candidate_genes(
        adata,
        batch_key=batch_key,
        sct_key=sct_key,
        n_candidates=n_candidates,
    )
    matrix = _residual_matrix(
        adata,
        batch_key=batch_key,
        sct_key=sct_key,
        genes=candidates,
        layer=layer,
    )
    labels = adata.obs[batch_key].astype(str).to_numpy()
    if len(np.unique(labels)) < 2:
        return candidates[:n_features]

    encoder = LabelEncoder()
    y = encoder.fit_transform(labels)
    clf = _fit_batch_classifier(
        matrix,
        y,
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        subsample=subsample,
        colsample_bytree=colsample_bytree,
        seed=seed,
    )
    clf.fit(matrix, y)
    importance = clf.feature_importances_
    order = np.argsort(importance)[::-1]
    selected = [candidates[i] for i in order[: min(n_features, len(candidates))]]
    return selected


def integrate_pair_xgb(
    data_ref: np.ndarray,
    data_query: np.ndarray,
    anchors: pd.DataFrame,
    *,
    n_delta_pcs: int = 20,
    n_estimators: int = 150,
    max_depth: int = 4,
    learning_rate: float = 0.08,
    subsample: float = 0.9,
    colsample_bytree: float = 0.8,
    seed: int = 0,
) -> np.ndarray:
    """
    Optional enhancement: integrate query toward reference with XGBoost delta prediction.

    Uses the same Seurat RPCA anchors; replaces only the IntegrateData-style
    anchor-weighted correction with a learned low-rank delta model.
    """
    _require_xgboost()
    if anchors.empty:
        return data_query.copy()

    cell1 = anchors["cell1"].astype(int).to_numpy()
    cell2 = anchors["cell2"].astype(int).to_numpy()
    integration_matrix = data_query[cell2] - data_ref[cell1]
    n_anchors = integration_matrix.shape[0]
    if n_anchors < 2:
        corrected = data_query.copy()
        for i in range(n_anchors):
            corrected[cell2[i]] += integration_matrix[i]
        return corrected

    x_train = data_query[cell2]
    y_train = integration_matrix
    n_comp = min(n_delta_pcs, n_anchors - 1, y_train.shape[1])
    n_comp = max(n_comp, 1)
    pca = PCA(n_components=n_comp, random_state=seed)
    y_scores = pca.fit_transform(y_train)

    reg = MultiOutputRegressor(
        XGBRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            subsample=subsample,
            colsample_bytree=colsample_bytree,
            random_state=seed,
            n_jobs=1,
            verbosity=0,
        )
    )
    reg.fit(x_train, y_scores)
    y_pred = reg.predict(data_query)
    delta = y_pred @ pca.components_
    return data_query + delta
