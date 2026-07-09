#!/usr/bin/env Rscript
# Full sctransform::vst for trackcell R backend (Seurat 4.4.0 parity).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: run_vst.R <umi.mtx> <genes.txt> <cells.txt> <cell_attr.csv> <config.json> <out_dir>"
  )
}

umi_path <- args[[1]]
genes_path <- args[[2]]
cells_path <- args[[3]]
attr_path <- args[[4]]
config_path <- args[[5]]
out_dir <- args[[6]]

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
  library(sctransform)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cfg <- fromJSON(config_path)

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
umi <- as(umi, "CsparseMatrix")

cell_attr <- read.csv(attr_path, row.names = 1, check.names = FALSE)
cell_attr <- cell_attr[cells, , drop = FALSE]
if (!"log_umi" %in% colnames(cell_attr)) {
  cell_attr$log_umi <- log10(Matrix::colSums(umi) + 1)
}

seed <- cfg$seed %||% 1448145L
n_cells <- cfg$n_cells
if (!is.null(n_cells)) {
  n_cells <- min(as.integer(n_cells), ncol(umi))
}
n_genes <- cfg$n_genes %||% 2000L
vst_flavor <- cfg$vst_flavor
variable_features_n <- cfg$variable_features_n %||% 3000L
variable_features_rv_th <- cfg$variable_features_rv_th %||% 1.3
do_correct_umi <- isTRUE(cfg$do_correct_umi %||% TRUE)
residual_type <- cfg$residual_type %||% "pearson"
return_only_var_genes <- isTRUE(cfg$return_only_var_genes %||% TRUE)

clip_range <- cfg$clip_range
if (is.null(clip_range)) {
  clip_range <- c(-sqrt(ncol(umi) / 30), sqrt(ncol(umi) / 30))
}

set.seed(seed)
vst.args <- list(
  umi = umi,
  cell_attr = cell_attr,
  n_cells = n_cells,
  n_genes = as.integer(n_genes),
  return_corrected_umi = do_correct_umi && identical(residual_type, "pearson"),
  return_gene_attr = TRUE,
  return_cell_attr = TRUE,
  residual_type = residual_type,
  verbosity = 0
)
if (!is.null(vst_flavor)) {
  vst.args$vst.flavor <- vst_flavor
}

vst.out <- do.call(vst, vst.args)

feature.variance <- vst.out$gene_attr$residual_variance
if (is.null(names(feature.variance))) {
  names(feature.variance) <- rownames(vst.out$gene_attr)
}
feature.variance <- feature.variance[!is.na(feature.variance)]
if (length(feature.variance) == 0) {
  stop("No finite residual variances in vst output.")
}
feature.variance <- feature.variance[order(feature.variance, decreasing = TRUE)]
if (!is.null(variable_features_n)) {
  top.features <- names(feature.variance)[
    seq_len(min(as.integer(variable_features_n), length(feature.variance)))
  ]
} else {
  top.features <- names(feature.variance)[
    feature.variance >= variable_features_rv_th
  ]
}
top.features <- as.character(top.features)
top.features <- top.features[!is.na(top.features) & nzchar(top.features)]
residual.features <- top.features
if (length(residual.features) == 0) {
  stop("No variable features selected from vst output.")
}

y <- vst.out$y[residual.features, , drop = FALSE]
scale.data <- y
scale.data[scale.data < clip_range[1]] <- clip_range[1]
scale.data[scale.data > clip_range[2]] <- clip_range[2]
scale.data <- sweep(scale.data, 1, rowMeans(scale.data), "-")
scale.data <- as(scale.data, "CsparseMatrix")

writeLines(as.character(residual.features), file.path(out_dir, "variable_features.txt"))
writeMM(scale.data, file.path(out_dir, "scale_data.mtx"))
writeLines(as.character(rownames(scale.data)), file.path(out_dir, "scale_genes.txt"))
writeLines(as.character(colnames(scale.data)), file.path(out_dir, "scale_cells.txt"))
write.csv(vst.out$gene_attr, file.path(out_dir, "gene_attr.csv"))
write.csv(vst.out$model_pars_fit, file.path(out_dir, "model_pars_fit.csv"))

if (!is.null(vst.out$umi_corrected)) {
  umi_corrected <- as(vst.out$umi_corrected, "CsparseMatrix")
  writeMM(umi_corrected, file.path(out_dir, "counts.mtx"))
  writeLines(as.character(rownames(vst.out$umi_corrected)), file.path(out_dir, "counts_genes.txt"))
  writeLines(as.character(colnames(vst.out$umi_corrected)), file.path(out_dir, "counts_cells.txt"))
  has_corrected <- TRUE
} else {
  has_corrected <- FALSE
}

meta <- list(
  seed = seed,
  n_cells = n_cells,
  n_genes = n_genes,
  vst_flavor = vst_flavor,
  residual_type = residual_type,
  do_correct_umi = do_correct_umi,
  clip_range = as.numeric(clip_range),
  has_corrected_counts = has_corrected,
  n_variable_features = length(residual.features),
  vst_method = if (!is.null(vst.out$arguments$method)) vst.out$arguments$method else NA_character_,
  model_str = vst.out$model_str
)
write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "meta.json"))
