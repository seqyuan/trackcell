#!/usr/bin/env Rscript
# Export SCT VST artifacts for Trackcell's offline parity benchmark.
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
out_dir <- file.path(benchmark_root, "steps", "sct_benchmark", "r", paste0("seed_", seed), sample)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
counts <- Read10X(file.path(data_dir, sample), gene.column = 2)
obj <- CreateSeuratObject(counts, project = sample, min.cells = 100, min.features = 50)
obj <- RenameCells(obj, new.names = paste0(sample, "_", colnames(obj)))
umi <- GetAssayData(obj, assay = "RNA", slot = "counts")
writeLines(rownames(umi), file.path(ref_dir, paste0(sample, "_genes_min_cells.txt")))
cell_attr <- data.frame(log_umi = log10(obj$nCount_RNA), row.names = colnames(obj))
set.seed(seed)
vst.out <- vst(umi = umi, cell_attr = cell_attr, n_cells = ncells, n_genes = ngenes,
               vst.flavor = "v2", return_corrected_umi = FALSE, verbosity = 0)
write.csv(vst.out$model_pars, file.path(out_dir, "model_pars_step1.csv"))
write.csv(vst.out$model_pars_fit, file.path(out_dir, "model_pars_fit.csv"))
write.csv(vst.out$gene_attr, file.path(out_dir, "gene_attr.csv"))
writeLines(rownames(vst.out$gene_attr), file.path(out_dir, "genes_vst.txt"))
write.csv(data.frame(gmean = vst.out$genes_log_gmean), file.path(out_dir, "genes_log_gmean.csv"))
write.csv(data.frame(gmean = vst.out$genes_log_gmean_step1), file.path(out_dir, "genes_log_gmean_step1.csv"))
writeLines(vst.out$cells_step1, file.path(out_dir, "cells_step1.txt"))
writeLines(names(vst.out$genes_log_gmean_step1), file.path(out_dir, "genes_step1.txt"))
write.csv(vst.out$y, file.path(out_dir, "residuals.csv"))
hvg <- rownames(vst.out$gene_attr)[order(vst.out$gene_attr$residual_variance, decreasing = TRUE)]
writeLines(head(hvg, n_hvg), file.path(out_dir, "hvg_top3000.txt"))
min_var <- vst.out$arguments$min_variance
if (is.null(min_var) || !is.numeric(min_var)) min_var <- (median(umi@x) / 5)^2
write_json(list(sample = sample, seed = seed, n_cells = ncol(umi), n_genes = nrow(umi),
                min_variance = min_var, bw_adjust = 3, gmean_eps = 1),
           file.path(out_dir, "meta.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Exported", sample, "seed", seed, "to", out_dir, "\n")
