#!/usr/bin/env Rscript
# Fit glmGamPoi_offset for trackcell (sctransform::fit_glmGamPoi_offset).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: fit_glmGamPoi_offset.R <umi.mtx> <genes.txt> <cells.txt> <cell_attr.csv> <out.csv> <allow_inf_theta>")
}

umi_path <- args[[1]]
genes_path <- args[[2]]
cells_path <- args[[3]]
attr_path <- args[[4]]
out_path <- args[[5]]
allow_inf_theta <- as.logical(args[[6]])

suppressPackageStartupMessages({
  library(Matrix)
  library(sctransform)
})

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  stop("glmGamPoi is not installed in this R environment.")
}

genes <- readLines(genes_path)
cells <- readLines(cells_path)
umi <- readMM(umi_path)
if (nrow(umi) != length(genes) || ncol(umi) != length(cells)) {
  stop(
    "UMI dimensions (", nrow(umi), " x ", ncol(umi),
    ") do not match gene/cell lists (", length(genes), " x ", length(cells), ")."
  )
}
rownames(umi) <- genes
colnames(umi) <- cells
umi <- as.matrix(umi)

data <- read.csv(attr_path, row.names = 1, check.names = FALSE)
data <- data[cells, , drop = FALSE]

fit_fun <- get("fit_glmGamPoi_offset", envir = asNamespace("sctransform"))
model_pars <- fit_fun(
  umi = umi,
  model_str = "y ~ log_umi",
  data = data,
  allow_inf_theta = allow_inf_theta
)
write.csv(model_pars, out_path, row.names = TRUE)
