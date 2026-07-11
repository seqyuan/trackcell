#!/usr/bin/env Rscript
# Seurat v4 SCT + RPCA integration benchmark (GSE288946, 3 samples).
#
# Filtering matches SCT benchmark:
#   min.cells = 100, min.features = 50
#
# Env: RPCA_BENCHMARK_ROOT, RPCA_SEED (default 1448145)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(jsonlite)
  library(sctransform)
})

benchmark_root <- Sys.getenv("RPCA_BENCHMARK_ROOT", unset = "/Volumes/process/tmp/tcl_test")
seed <- as.integer(Sys.getenv("RPCA_SEED", unset = "1448145"))
options(future.globals.maxSize = 8000 * 1024^2)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
samples <- c("GSM8779707", "GSM8779708", "GSM8779709")
data_dir <- file.path(benchmark_root, "data")

min_cells_gene <- 100L
min_features <- 50L
N_INTEGRATION_FEATURES <- 2000L
N_SCT_HVG <- 3000L
DIMS <- 1:30
K_ANCHOR <- 5L
K_SCORE <- 30L
K_WEIGHT <- 100L
N_TREES <- 50L

step_root <- file.path(benchmark_root, "steps", "rpca_benchmark", "r", paste0("seed_", seed))
dir.create(step_root, recursive = TRUE, showWarnings = FALSE)

cell_id <- function(sample, barcode) paste0(sample, "_", barcode)

save_json <- function(obj, path) {
  write(toJSON(obj, auto_unbox = TRUE, pretty = TRUE, null = "null"), path)
}

cat(sprintf("=== RPCA benchmark seed=%d root=%s ===\n", seed, step_root))

cat("Step 01: load + QC filter\n")
obj.list <- list()
filter_summary <- list(
  min_cells_gene = min_cells_gene,
  min_features = min_features,
  samples = list()
)

for (s in samples) {
  counts <- Read10X(file.path(data_dir, s), gene.column = 2)
  obj <- CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = min_cells_gene,
    min.features = min_features
  )
  obj[["orig.ident"]] <- s
  n_before <- ncol(obj)
  obj <- RenameCells(
    obj,
    new.names = vapply(colnames(obj), function(bc) cell_id(s, bc), character(1))
  )
  n_after <- ncol(obj)
  filter_summary$samples[[s]] <- list(n_before = n_before, n_after = n_after)

  step_dir <- file.path(step_root, "01_filter")
  dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(colnames(obj), file.path(step_dir, paste0(s, "_cells.txt")))
  obj.list[[s]] <- obj
  cat(sprintf("  %s: %d cells\n", s, n_after))
}
save_json(filter_summary, file.path(step_root, "01_filter", "summary.json"))

cat("Step 02: SCTransform (vst.flavor='v2') per sample\n")
step_sct <- file.path(step_root, "02_sct")
dir.create(step_sct, recursive = TRUE, showWarnings = FALSE)

obj.list <- lapply(obj.list, function(x) {
  SCTransform(
    x,
    vst.flavor = "v2",
    variable.features.n = N_SCT_HVG,
    verbose = FALSE,
    seed.use = seed
  )
})

for (s in samples) {
  obj <- obj.list[[s]]
  hvg <- VariableFeatures(obj)
  writeLines(hvg, file.path(step_sct, paste0(s, "_hvg.txt")))

  model <- obj@assays$SCT@SCTModel.list$model1
  pars <- slot(model, "feature.attributes")
  pars <- pars[, c("theta", "(Intercept)", "log_umi"), drop = FALSE]
  colnames(pars) <- c("theta", "Intercept", "log_umi")
  write.csv(pars, file.path(step_sct, paste0(s, "_model_pars.csv")))

  scale_data <- GetAssayData(obj, assay = "SCT", slot = "scale.data")
  write.csv(
    as.matrix(scale_data),
    file.path(step_sct, paste0(s, "_residuals_hvg.csv"))
  )
  cat(sprintf("  %s: %d HVG, %d cells\n", s, length(hvg), ncol(obj)))
}

# Export step-1 details for P0 parity (cells_step1, genes_step1, genes_log_gmean, model_pars_step1)
cat("Step 02b: export step-1 details for SCT parity\n")
for (s in samples) {
  counts <- Read10X(file.path(data_dir, s), gene.column = 2)
  obj_raw <- CreateSeuratObject(counts = counts, project = s,
                                min.cells = min_cells_gene, min.features = min_features)
  obj_raw <- RenameCells(obj_raw,
    new.names = vapply(colnames(obj_raw), function(bc) cell_id(s, bc), character(1)))
  umi_mat <- GetAssayData(obj_raw, assay = "RNA", slot = "counts")

  set.seed(seed)
  vst_out <- sctransform::vst(
    umi = umi_mat,
    vst.flavor = "v2",
    n_genes = 2000L,
    n_cells = 2000L,
    verbosity = 0L
  )

  cells_step1 <- vst_out$cells_step1
  genes_step1 <- rownames(vst_out$model_pars)
  writeLines(as.character(cells_step1), file.path(step_sct, paste0(s, "_cells_step1.txt")))
  writeLines(genes_step1, file.path(step_sct, paste0(s, "_genes_step1.txt")))

  row_gmean <- function(mat, eps = 1) exp(rowMeans(log(mat + eps))) - eps
  all_gmean <- log10(row_gmean(umi_mat))
  write.csv(data.frame(gene = names(all_gmean), log10_gmean = as.numeric(all_gmean)),
            file.path(step_sct, paste0(s, "_genes_log_gmean_all.csv")), row.names = FALSE)

  write.csv(data.frame(gene = names(vst_out$genes_log_gmean_step1),
                       log10_gmean = as.numeric(vst_out$genes_log_gmean_step1)),
            file.path(step_sct, paste0(s, "_genes_log_gmean_step1.csv")), row.names = FALSE)

  write.csv(as.data.frame(vst_out$model_pars),
            file.path(step_sct, paste0(s, "_model_pars_step1.csv")))

  cat(sprintf("  %s: cells_step1=%d, genes_step1=%d\n",
              s, length(cells_step1), length(genes_step1)))
}

cat("Step 03: SelectIntegrationFeatures\n")
step_feat <- file.path(step_root, "03_features")
dir.create(step_feat, recursive = TRUE, showWarnings = FALSE)
int_features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = N_INTEGRATION_FEATURES
)
writeLines(int_features, file.path(step_feat, "integration_features.txt"))
save_json(
  list(n_features = length(int_features), seed = seed),
  file.path(step_feat, "summary.json")
)
cat(sprintf("  %d integration features\n", length(int_features)))

cat("Step 04: PrepSCTIntegration\n")
step_prep <- file.path(step_root, "04_prep")
dir.create(step_prep, recursive = TRUE, showWarnings = FALSE)
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = int_features,
  verbose = FALSE
)

prep_blocks <- list()
for (s in samples) {
  obj <- obj.list[[s]]
  mat <- t(as.matrix(GetAssayData(obj, assay = "SCT", slot = "scale.data")[int_features, , drop = FALSE]))
  rownames(mat) <- colnames(obj)
  prep_blocks[[s]] <- mat
}
prep_all <- do.call(rbind, prep_blocks)
write.csv(prep_all, file.path(step_prep, "prep_residuals.csv"), row.names = TRUE)
save_json(
  list(n_features = length(int_features), n_cells = nrow(prep_all)),
  file.path(step_prep, "summary.json")
)
cat(sprintf("  prep matrix: %d cells x %d genes\n", nrow(prep_all), length(int_features)))

cat("Step 04b: RunPCA per sample\n")
obj.list <- lapply(obj.list, function(x) {
  RunPCA(x, features = int_features, npcs = max(DIMS), verbose = FALSE)
})

cat("Step 05: FindIntegrationAnchors (RPCA, SCT)\n")
step_anchors <- file.path(step_root, "05_anchors")
dir.create(step_anchors, recursive = TRUE, showWarnings = FALSE)
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = int_features,
  normalization.method = "SCT",
  reduction = "rpca",
  dims = DIMS,
  k.anchor = K_ANCHOR,
  k.score = K_SCORE,
  n.trees = N_TREES,
  verbose = FALSE
)
anchor_df <- slot(anchors, "anchors")
write.csv(anchor_df, file.path(step_anchors, "anchors.csv"), row.names = FALSE)
save_json(
  list(
    n_anchors = nrow(anchor_df),
    params = list(
      dims = max(DIMS),
      k_anchor = K_ANCHOR,
      k_score = K_SCORE,
      n_trees = N_TREES
    )
  ),
  file.path(step_anchors, "summary.json")
)
cat(sprintf("  %d anchors\n", nrow(anchor_df)))

cat("Step 06: IntegrateData (SCT)\n")
step_int <- file.path(step_root, "06_integrated")
dir.create(step_int, recursive = TRUE, showWarnings = FALSE)
integrated <- IntegrateData(
  anchors,
  normalization.method = "SCT",
  dims = DIMS,
  k.weight = K_WEIGHT,
  verbose = FALSE
)
int_mat <- t(as.matrix(GetAssayData(integrated, assay = "integrated", slot = "scale.data")))
rownames(int_mat) <- colnames(integrated)
write.csv(int_mat, file.path(step_int, "integrated_residuals.csv"), row.names = TRUE)
save_json(
  list(n_cells = nrow(int_mat), n_features = ncol(int_mat)),
  file.path(step_int, "summary.json")
)
cat(sprintf("  integrated: %d cells x %d genes\n", nrow(int_mat), ncol(int_mat)))

save_json(
  list(
    seurat_version = as.character(packageVersion("Seurat")),
    sctransform_version = as.character(packageVersion("sctransform")),
    seed = seed,
    samples = samples,
    filter = filter_summary,
    n_integration_features = length(int_features),
    n_anchors = nrow(anchor_df),
    n_integrated_cells = nrow(int_mat)
  ),
  file.path(step_root, "pipeline_summary.json")
)
cat("RPCA benchmark R export done.\n")
