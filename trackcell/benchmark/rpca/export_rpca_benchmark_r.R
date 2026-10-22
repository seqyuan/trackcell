#!/usr/bin/env Rscript
# Offline Seurat v4 SCT + RPCA parity export for Trackcell.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(jsonlite); library(sctransform)
})
root <- Sys.getenv("RPCA_BENCHMARK_ROOT", unset = "/Volumes/process/tmp/tcl_test")
seed <- as.integer(Sys.getenv("RPCA_SEED", unset = "1448145"))
samples <- c("GSM8779707", "GSM8779708", "GSM8779709")
out <- file.path(root, "steps", "rpca_benchmark", "r", paste0("seed_", seed))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
write_json_file <- function(x, path) write(toJSON(x, auto_unbox = TRUE, pretty = TRUE), path)
objs <- list()
for (sample in samples) {
  counts <- Read10X(file.path(root, "data", sample), gene.column = 2)
  obj <- CreateSeuratObject(counts, project = sample, min.cells = 100, min.features = 50)
  obj[["orig.ident"]] <- sample
  obj <- RenameCells(obj, new.names = paste0(sample, "_", colnames(obj)))
  objs[[sample]] <- obj
  dir.create(file.path(out, "01_filter"), recursive = TRUE, showWarnings = FALSE)
  writeLines(colnames(obj), file.path(out, "01_filter", paste0(sample, "_cells.txt")))
}
step_sct <- file.path(out, "02_sct"); dir.create(step_sct, recursive = TRUE, showWarnings = FALSE)
objs <- lapply(objs, function(x) SCTransform(x, vst.flavor = "v2", variable.features.n = 3000,
                                             verbose = FALSE, seed.use = seed))
for (sample in samples) {
  obj <- objs[[sample]]
  writeLines(VariableFeatures(obj), file.path(step_sct, paste0(sample, "_hvg.txt")))
  model <- slot(obj@assays$SCT@SCTModel.list$model1, "feature.attributes")
  model <- model[, c("theta", "(Intercept)", "log_umi"), drop = FALSE]
  colnames(model) <- c("theta", "Intercept", "log_umi")
  write.csv(model, file.path(step_sct, paste0(sample, "_model_pars.csv")))
  write.csv(as.matrix(GetAssayData(obj, assay = "SCT", slot = "scale.data")),
            file.path(step_sct, paste0(sample, "_residuals_hvg.csv")))

  counts <- Read10X(file.path(root, "data", sample), gene.column = 2)
  raw <- CreateSeuratObject(counts, project = sample, min.cells = 100, min.features = 50)
  raw <- RenameCells(raw, new.names = paste0(sample, "_", colnames(raw)))
  umi <- GetAssayData(raw, assay = "RNA", slot = "counts")
  set.seed(seed)
  fit <- sctransform::vst(umi, vst.flavor = "v2", n_genes = 2000, n_cells = 2000,
                          verbosity = 0)
  writeLines(as.character(fit$cells_step1), file.path(step_sct, paste0(sample, "_cells_step1.txt")))
  writeLines(rownames(fit$model_pars), file.path(step_sct, paste0(sample, "_genes_step1.txt")))
  row_gmean <- function(mat, eps = 1) exp(rowMeans(log(mat + eps))) - eps
  all_gmean <- log10(row_gmean(umi))
  write.csv(data.frame(gene = names(all_gmean), log10_gmean = as.numeric(all_gmean)),
            file.path(step_sct, paste0(sample, "_genes_log_gmean_all.csv")), row.names = FALSE)
  write.csv(data.frame(gene = names(fit$genes_log_gmean_step1),
                       log10_gmean = as.numeric(fit$genes_log_gmean_step1)),
            file.path(step_sct, paste0(sample, "_genes_log_gmean_step1.csv")), row.names = FALSE)
  write.csv(as.data.frame(fit$model_pars),
            file.path(step_sct, paste0(sample, "_model_pars_step1.csv")))
}
features <- SelectIntegrationFeatures(objs, nfeatures = 2000)
dir.create(file.path(out, "03_features"), recursive = TRUE, showWarnings = FALSE)
writeLines(features, file.path(out, "03_features", "integration_features.txt"))
objs <- PrepSCTIntegration(objs, anchor.features = features, verbose = FALSE)
prep <- do.call(rbind, lapply(objs, function(obj) {
  mat <- t(as.matrix(GetAssayData(obj, assay = "SCT", slot = "scale.data")[features, , drop = FALSE]))
  rownames(mat) <- colnames(obj); mat
}))
dir.create(file.path(out, "04_prep"), recursive = TRUE, showWarnings = FALSE)
write.csv(prep, file.path(out, "04_prep", "prep_residuals.csv"))
objs <- lapply(objs, function(x) RunPCA(x, features = features, npcs = 30, verbose = FALSE))
anchors <- FindIntegrationAnchors(objs, anchor.features = features, normalization.method = "SCT",
                                  reduction = "rpca", dims = 1:30, k.anchor = 5,
                                  k.score = 30, n.trees = 50, verbose = FALSE)
dir.create(file.path(out, "05_anchors"), recursive = TRUE, showWarnings = FALSE)
write.csv(slot(anchors, "anchors"), file.path(out, "05_anchors", "anchors.csv"), row.names = FALSE)
integrated <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:30,
                            k.weight = 100, verbose = FALSE)
int_mat <- t(as.matrix(GetAssayData(integrated, assay = "integrated", slot = "scale.data")))
rownames(int_mat) <- colnames(integrated)
dir.create(file.path(out, "06_integrated"), recursive = TRUE, showWarnings = FALSE)
write.csv(int_mat, file.path(out, "06_integrated", "integrated_residuals.csv"))
write_json_file(list(seurat_version = as.character(packageVersion("Seurat")),
                     sctransform_version = as.character(packageVersion("sctransform")),
                     seed = seed, samples = samples, n_integration_features = length(features),
                     n_anchors = nrow(slot(anchors, "anchors"))),
                     n_integrated_cells = nrow(int_mat)), file.path(out, "pipeline_summary.json"))
