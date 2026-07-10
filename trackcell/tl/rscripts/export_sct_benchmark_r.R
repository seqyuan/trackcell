#!/usr/bin/env Rscript
# Export full SCT VST artifacts for trackcell SCT benchmark (per sample × seed).
#
# Filtering (CreateSeuratObject):
#   min.cells = 100  (gene expressed in >= 100 cells)
#   min.features = 50 (cell expresses >= 50 genes)
#
# Env: SCT_SAMPLE, SCT_SEED, SCT_BENCHMARK_ROOT (optional)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(sctransform)
  library(jsonlite)
})

benchmark_root <- Sys.getenv("SCT_BENCHMARK_ROOT", unset = "/Volumes/process/tmp/tcl_test")
data_dir <- file.path(benchmark_root, "data")
ref_dir <- file.path(benchmark_root, "ref")
sample <- Sys.getenv("SCT_SAMPLE", unset = "GSM8779707")
seed <- as.integer(Sys.getenv("SCT_SEED", unset = "1448145"))
ncells <- as.integer(Sys.getenv("SCT_NCELLS", unset = "2000"))
ngenes <- as.integer(Sys.getenv("SCT_NGENES", unset = "2000"))
n_hvg <- as.integer(Sys.getenv("SCT_N_HVG", unset = "3000"))

min_cells_gene <- 100L
min_features <- 50L

out_dir <- file.path(benchmark_root, "steps", "sct_benchmark", "r", paste0("seed_", seed), sample)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

cell_id <- function(s, bc) paste0(s, "_", bc)

counts <- Read10X(file.path(data_dir, sample), gene.column = 2)
obj <- CreateSeuratObject(
  counts = counts,
  project = sample,
  min.cells = min_cells_gene,
  min.features = min_features
)
obj <- RenameCells(obj, new.names = vapply(colnames(obj), function(bc) cell_id(sample, bc), character(1)))

umi <- GetAssayData(obj, assay = "RNA", slot = "counts")
writeLines(rownames(umi), file.path(ref_dir, paste0(sample, "_genes_min_cells.txt")))
writeLines(colnames(umi), file.path(out_dir, "cells_after_filter.txt"))

cell_attr <- data.frame(
  log_umi = log10(obj$nCount_RNA),
  row.names = colnames(obj)
)

set.seed(seed)
vst.out <- do.call(vst, list(
  umi = umi,
  cell_attr = cell_attr,
  n_cells = ncells,
  n_genes = ngenes,
  vst.flavor = "v2",
  return_corrected_umi = FALSE,
  verbosity = 0
))

write.csv(vst.out$model_pars, file.path(out_dir, "model_pars_step1.csv"))
write.csv(vst.out$model_pars_fit, file.path(out_dir, "model_pars_fit.csv"))
write.csv(vst.out$gene_attr, file.path(out_dir, "gene_attr.csv"))
writeLines(rownames(vst.out$gene_attr), file.path(out_dir, "genes_vst.txt"))
write.csv(data.frame(gmean = vst.out$genes_log_gmean), file.path(out_dir, "genes_log_gmean.csv"))
write.csv(
  data.frame(gmean = vst.out$genes_log_gmean_step1),
  file.path(out_dir, "genes_log_gmean_step1.csv")
)
writeLines(vst.out$cells_step1, file.path(out_dir, "cells_step1.txt"))
writeLines(names(vst.out$genes_log_gmean_step1), file.path(out_dir, "genes_step1.txt"))
write.csv(vst.out$y, file.path(out_dir, "residuals.csv"))

hvg <- rownames(vst.out$gene_attr)[order(vst.out$gene_attr$residual_variance, decreasing = TRUE)][seq_len(min(n_hvg, nrow(vst.out$gene_attr)))]
writeLines(hvg, file.path(out_dir, "hvg_top3000.txt"))

min_var <- vst.out$arguments$min_variance
if (is.null(min_var) || !is.numeric(min_var)) {
  min_var <- (median(umi@x) / 5)^2
}

write_json(
  list(
    sample = sample,
    seed = seed,
    filter = list(min_cells_gene = min_cells_gene, min_features = min_features),
    n_cells = ncol(umi),
    n_genes = nrow(umi),
    n_cells_step1 = length(vst.out$cells_step1),
    n_genes_step1 = length(vst.out$genes_log_gmean_step1),
    min_variance = min_var,
    bw_adjust = 3,
    gmean_eps = 1
  ),
  file.path(out_dir, "meta.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("Exported", sample, "seed", seed, "to", out_dir, "\n")
