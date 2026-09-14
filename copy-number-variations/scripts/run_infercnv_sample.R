#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  BiocManager::install("rhdf5")
}
if (!requireNamespace("infercnv", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  BiocManager::install("infercnv")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  BiocManager::install("jsonlite")
}

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(infercnv)
  library(jsonlite)
})

cat("Reading input arguments.")

option_list <- list(
  make_option("--counts_h5", type = "character"),
  make_option("--annotations_tsv", type = "character"),
  make_option("--gene_order_file", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--ref_groups", type = "character", default = ""),
  make_option("--cutoff", type = "double", default = 0.1),
  make_option("--cluster_by_groups", type = "logical", default = TRUE),
  make_option("--denoise", type = "logical", default = TRUE),
  make_option("--HMM", type = "logical", default = FALSE),
  make_option("--num_threads", type = "integer", default = 4)
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading HDF5 counts: ", opt$counts_h5)
counts <- h5read(opt$counts_h5, "counts")
genes <- as.character(h5read(opt$counts_h5, "genes"))
barcodes <- as.character(h5read(opt$counts_h5, "barcodes"))

cat("counts dim: ", paste(dim(counts), collapse = " x "))
cat("n genes: ", length(genes))
cat("genes[1:5]: ", paste(head(genes, 5), collapse = ", "))
cat("n barcodes: ", length(barcodes))

if (length(dim(counts)) != 2) {
  stop("H5 dataset 'counts' is not 2-dimensional.")
}

# handle either genes x spots or spots x genes
if (nrow(counts) == length(genes) && ncol(counts) == length(barcodes)) {
  cat("Counts matrix orientation looks correct: genes x spots")
} else if (nrow(counts) == length(barcodes) && ncol(counts) == length(genes)) {
  cat("Counts matrix appears transposed: spots x genes. Transposing now.")
  counts <- t(counts)
} else {
  stop(
    "Counts dimensions do not match genes/barcodes.\n",
    "counts dim = ", paste(dim(counts), collapse = " x "), "\n",
    "length(genes) = ", length(genes), "\n",
    "length(barcodes) = ", length(barcodes)
  )
}

rownames(counts) <- genes
colnames(counts) <- barcodes

cat("Reading annotations: ", opt$annotations_tsv)
ann <- read.table(
  opt$annotations_tsv,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(ann) < 2) {
  stop("Annotation file must have at least 2 columns: barcode and group.")
}

colnames(ann)[1:2] <- c("barcode", "group")

# keep only shared barcodes, preserve annotation order
shared <- intersect(ann$barcode, colnames(counts))
if (length(shared) == 0) {
  stop("No overlap between annotation barcodes and H5 barcodes.")
}

ann <- ann[ann$barcode %in% shared, , drop = FALSE]
counts <- counts[, ann$barcode, drop = FALSE]

# write a reordered annotation file for reproducibility
ann_reordered_path <- file.path(opt$out_dir, "annotations_reordered.tsv")
write.table(
  ann[, c("barcode", "group")],
  file = ann_reordered_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

ref_groups_string <- opt$ref_groups
ref_groups_string <- gsub("__SPACE__", " ", ref_groups_string, fixed = TRUE)
ref_groups <- character(0)

if (nzchar(ref_groups_string)) {
  ref_groups <- trimws(strsplit(ref_groups_string, ",", fixed = TRUE)[[1]])
}

cat("Raw ref_groups string: ", ref_groups_string, "\n", sep = "")
cat("Parsed ref_groups: ", paste(ref_groups, collapse = " | "), "\n", sep = "")

if (length(ref_groups) > 0) {
  missing_refs <- setdiff(ref_groups, unique(ann$group))
  if (length(missing_refs) > 0) {
    stop(
      "These ref_groups were requested but are not present in annotations: ",
      paste(missing_refs, collapse = ", ")
    )
  }
}

export_infercnv_for_custom_plot <- function(infercnv_obj, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr <- infercnv_obj@expr.data
  genes <- rownames(expr)
  cells <- colnames(expr)

  # cell metadata from grouped indices
  cell_meta <- data.frame(
    cell = cells,
    cell_index = seq_along(cells),
    role = NA_character_,
    group = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(infercnv_obj@reference_grouped_cell_indices) > 0) {
    for (grp in names(infercnv_obj@reference_grouped_cell_indices)) {
      idx <- infercnv_obj@reference_grouped_cell_indices[[grp]]
      cell_meta$role[idx] <- "reference"
      cell_meta$group[idx] <- grp
    }
  }

  if (length(infercnv_obj@observation_grouped_cell_indices) > 0) {
    for (grp in names(infercnv_obj@observation_grouped_cell_indices)) {
      idx <- infercnv_obj@observation_grouped_cell_indices[[grp]]
      cell_meta$role[idx] <- "observation"
      cell_meta$group[idx] <- grp
    }
  }

  # gene metadata from gene_order if available
  gene_meta <- data.frame(
    gene = genes,
    gene_index = seq_along(genes),
    stringsAsFactors = FALSE
  )

  # infercnv objects usually keep genomic positions in gene_order
  if (!is.null(infercnv_obj@gene_order)) {
    go <- infercnv_obj@gene_order
    go <- as.data.frame(go, stringsAsFactors = FALSE)
    # try to standardize common column names
    colnames(go)[1:min(4, ncol(go))] <- c("gene", "chr", "start", "end")[1:min(4, ncol(go))]
    gene_meta <- merge(
      gene_meta,
      go,
      by = "gene",
      all.x = TRUE,
      sort = FALSE
    )
    gene_meta <- gene_meta[match(genes, gene_meta$gene), , drop = FALSE]
  }

  # write exact plotted matrix
  matrix_path <- file.path(out_dir, "heatmap_matrix.tsv.gz")
  con <- gzfile(matrix_path, "wt")
  write.table(
    expr,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
  close(con)

  write.table(
    cell_meta,
    file = file.path(out_dir, "cell_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    gene_meta,
    file = file.path(out_dir, "gene_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("Creating infercnv object")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = ann_reordered_path,
  delim = "\t",
  gene_order_file = opt$gene_order_file,
  ref_group_names = if (length(ref_groups) > 0) ref_groups else NULL
)

cat("Running infercnv")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = opt$cutoff,
  out_dir = opt$out_dir,
  cluster_by_groups = opt$cluster_by_groups,
  denoise = opt$denoise,
  HMM = opt$HMM,
  num_threads = opt$num_threads
)

saveRDS(infercnv_obj, file = file.path(opt$out_dir, "infercnv_obj.rds"))

cat("Exporting infercnv results.")
export_infercnv_for_custom_plot(infercnv_obj, opt$out_dir)

cat("Done: ", opt$out_dir)
