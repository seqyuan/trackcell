"""Gene-name utility helpers (Seurat/scanpy alias normalization)."""

from __future__ import annotations

from typing import Optional, Sequence

import pandas as pd


def _gene_suffix_aliases(name: str) -> list[str]:
    """Return a gene name plus Seurat ``.N`` / ``-N`` duplicate suffix variants."""
    name = str(name)
    aliases = [name]
    if "." in name:
        base, suffix = name.rsplit(".", 1)
        if suffix.isdigit():
            aliases.append(f"{base}-{suffix}")
    if "-" in name:
        base, suffix = name.rsplit("-", 1)
        if suffix.isdigit():
            aliases.append(f"{base}.{suffix}")
    return aliases


def build_gene_alias_lookup(target_genes: pd.Index) -> dict[str, str]:
    """Map gene names and suffix aliases to canonical names in ``target_genes``."""
    lookup: dict[str, str] = {}
    for gene in map(str, target_genes):
        for alias in _gene_suffix_aliases(gene):
            lookup.setdefault(alias, gene)
    return lookup


def _gene_name_alias(
    name: str, candidates: set[str] | dict[str, str]
) -> Optional[str]:
    """Map Seurat/scanpy duplicate suffix variants to a name in ``candidates``."""
    if isinstance(candidates, dict):
        for alias in _gene_suffix_aliases(name):
            hit = candidates.get(alias)
            if hit is not None:
                return hit
        return None
    if name in candidates:
        return name
    for alias in _gene_suffix_aliases(name)[1:]:
        if alias in candidates:
            return alias
    return None


def resolve_gene_alias(name: str, target_genes: pd.Index) -> Optional[str]:
    """Resolve ``name`` to a canonical entry in ``target_genes`` if possible."""
    return _gene_name_alias(str(name), build_gene_alias_lookup(target_genes))


def align_gene_list(genes: Sequence[str], target_genes: pd.Index) -> list[str]:
    """Align a list of gene names to standard names in ``target_genes``."""
    lookup = build_gene_alias_lookup(target_genes)
    out: list[str] = []
    for gene in genes:
        matched = _gene_name_alias(str(gene), lookup)
        if matched is not None and matched not in out:
            out.append(matched)
    return out


def normalize_model_pars(model_pars: pd.DataFrame) -> pd.DataFrame:
    """Normalize common R/Python intercept column spellings."""
    out = model_pars.copy()
    out.index = out.index.astype(str)
    rename = {
        col: col.strip("()")
        for col in out.columns
        if col.startswith("(") and col.endswith(")")
    }
    return out.rename(columns=rename) if rename else out


def align_gene_frame(frame: pd.DataFrame, target_genes: pd.Index) -> pd.DataFrame:
    """Align a gene-indexed frame to canonical names in ``target_genes``."""
    frame = normalize_model_pars(frame)
    lookup = build_gene_alias_lookup(target_genes)
    rows: list[pd.Series] = []
    index: list[str] = []
    for gene, row in frame.iterrows():
        matched = _gene_name_alias(str(gene), lookup)
        if matched is not None and matched not in index:
            rows.append(row)
            index.append(matched)
    if not rows:
        return frame.iloc[0:0]
    return pd.DataFrame(rows, index=pd.Index(index))
