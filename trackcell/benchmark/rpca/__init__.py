"""SCT + RPCA integration benchmark vs Seurat reference."""

from trackcell.benchmark.rpca.data import (
    RPCA_BENCHMARK_ROOT,
    RPCA_REF_SEED,
    RPCA_SAMPLES,
    RPCA_SEEDS,
)
from trackcell.benchmark.rpca.report import write_rpca_benchmark_docs
from trackcell.benchmark.rpca.runner import run_rpca_benchmark_suite

__all__ = [
    "RPCA_BENCHMARK_ROOT",
    "RPCA_REF_SEED",
    "RPCA_SAMPLES",
    "RPCA_SEEDS",
    "run_rpca_benchmark_suite",
    "write_rpca_benchmark_docs",
]
