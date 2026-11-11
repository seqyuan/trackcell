#!/usr/bin/env Rscript
# Regenerate the SCT golden vectors committed next to this script.
#
#   Rscript export.R tests/data/sct_golden
#
# The input matrix is committed as umi.csv (40 genes x 120 cells, deterministic
# integer counts). `n_cells` / `n_genes` cover the whole matrix, so sctransform's
# subsampling cannot introduce R-RNG-dependent output; the Python side can run
# the same matrix without receiving cells_step1 / genes_step1.
#
# Requires: sctransform (tested with 0.4.3), Matrix, jsonlite.
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[[1]] else "tests/data/sct_golden"
suppressPackageStartupMessages({
  library(Matrix)
  library(sctransform)
  library(jsonlite)
})

umi <- as.matrix(read.csv(file.path(dir, "umi.csv"), row.names = 1, check.names = FALSE))
umi <- Matrix(umi, sparse = TRUE)
cell_attr <- data.frame(log_umi = log10(colSums(umi)), row.names = colnames(umi))

set.seed(1448145)
v <- sctransform::vst(
  umi = umi,
  cell_attr = cell_attr,
  n_cells = ncol(umi),
  n_genes = nrow(umi),
  vst.flavor = "v2",
  return_corrected_umi = FALSE,
  verbosity = 0
)

write.csv(v$model_pars, file.path(dir, "model_pars.csv"))
write.csv(v$model_pars_fit, file.path(dir, "model_pars_fit.csv"))
write.csv(data.frame(gmean = v$genes_log_gmean), file.path(dir, "genes_log_gmean.csv"))
write.csv(v$gene_attr, file.path(dir, "gene_attr.csv"))
write.csv(as.matrix(v$y), file.path(dir, "residuals.csv"))
write_json(
  list(
    n_genes = nrow(umi),
    n_cells = ncol(umi),
    seed = 1448145,
    vst_flavor = "v2",
    sctransform_version = as.character(packageVersion("sctransform")),
    min_variance = v$arguments$min_variance,
    gmean_eps = 1
  ),
  file.path(dir, "meta.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat("Exported SCT golden vectors to", dir, "\n")
