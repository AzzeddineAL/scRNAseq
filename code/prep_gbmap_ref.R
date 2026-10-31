# Builds a light Seurat reference from the GBmap core h5ad file.
# Run it once from the code/ folder : Rscript prep_gbmap_ref.R
# It uses ~6 GB of RAM, close other applications first.

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(Seurat)
})

h5ad <- "../data/gbmap_core/861acfd8-25f0-418b-a445-aa96da232827.h5ad"
out <- "../data/gbmap_core/gbmap_core_ref.rds"
max_per_type <- 1000  # max number of cells per annotation_level_4 type
set.seed(42)

# categorical columns are stored as integer codes + categories
read_cat <- function(path) {
  codes <- as.integer(h5read(h5ad, paste0(path, "/codes")))
  cats <- as.character(h5read(h5ad, paste0(path, "/categories")))
  cats[codes + 1L]
}

meta <- data.frame(
  row.names = as.character(h5read(h5ad, "obs/_index")),
  level1 = read_cat("obs/annotation_level_1"),
  level2 = read_cat("obs/annotation_level_2"),
  level3 = read_cat("obs/annotation_level_3"),
  level4 = read_cat("obs/annotation_level_4"),
  assay = read_cat("obs/assay"),
  donor = read_cat("obs/donor_id"))

# we keep rare cell types and take 10x cells first
is_10x <- grepl("^10x", meta$assay)
keep <- unlist(lapply(split(seq_len(nrow(meta)), meta$level4), function(idx) {
  idx <- sample(idx)
  idx <- idx[order(!is_10x[idx])]
  head(idx, max_per_type)
}))
keep <- sort(keep)

# raw/X holds raw counts stored by cell (CSR), so we read only the kept cells
indptr <- as.numeric(h5read(h5ad, "raw/X/indptr"))
genes <- read_cat("raw/var/feature_name")
fid <- H5Fopen(h5ad, flags = "H5F_ACC_RDONLY")
d_data <- H5Dopen(fid, "raw/X/data")
d_idx <- H5Dopen(fid, "raw/X/indices")
sp_data <- H5Dget_space(d_data)
sp_idx <- H5Dget_space(d_idx)

read_row <- function(d, sp, start, n) {
  H5Sselect_hyperslab(sp, "H5S_SELECT_SET", start = start, count = n)
  m <- H5Screate_simple(n)
  on.exit(H5Sclose(m))
  as.vector(H5Dread(d, h5spaceFile = sp, h5spaceMem = m))
}

x <- vector("list", length(keep))
i <- x
for (k in seq_along(keep)) {
  r <- keep[k]
  s <- indptr[r] + 1
  n <- indptr[r + 1] - indptr[r]
  x[[k]] <- read_row(d_data, sp_data, s, n)
  i[[k]] <- read_row(d_idx, sp_idx, s, n) + 1L
}
h5closeAll()

counts <- sparseMatrix(i = unlist(i), j = rep(seq_along(x), lengths(x)), x = unlist(x),
                       dims = c(length(genes), length(keep)),
                       dimnames = list(make.unique(genes), rownames(meta)[keep]))
rm(x, i)
gc()

gbmap <- CreateSeuratObject(counts, meta.data = meta[keep, ], project = "GBmap_core")
rm(counts)
gc()

gbmap <- NormalizeData(gbmap, verbose = FALSE)
gbmap <- FindVariableFeatures(gbmap, nfeatures = 3000, verbose = FALSE)
gbmap <- ScaleData(gbmap, verbose = FALSE)
gbmap <- RunPCA(gbmap, npcs = 30, verbose = FALSE)

saveRDS(gbmap, out)
print(gbmap)
