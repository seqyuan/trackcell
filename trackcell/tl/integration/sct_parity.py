"""SCT clip-range helper for integration."""

from __future__ import annotations

import math


def seurat_sct_clip_range(n_cells: int) -> tuple[float, float]:
    """Seurat ``SCTransform`` / ``sct.clip.range``: ``±sqrt(ncol / 30)``."""
    bound = float(math.sqrt(max(n_cells, 1) / 30.0))
    return (-bound, bound)
