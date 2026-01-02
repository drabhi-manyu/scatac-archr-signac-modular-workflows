# Optional Signac entry point from Cell Ranger ATAC outputs.
# This is meant as a simple portfolio example, and can be used to generate a CSV of labels for ArchR flavor 3.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(readr)
  library(hdf5r)
})

set.seed(1234)

# Provide a CSV with columns:
# sample_id,outs_dir
# where outs_dir contains: filtered_peak_bc_matrix.h5, singlecell.csv, fragments.tsv.gz
samples <- readr::read_csv("config/signac_samples.csv", show_col_types = FALSE)

objs <- list()

for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  outs <- samples$outs_dir[i]

  counts <- Read10X_h5(filename = file.path(outs, "filtered_peak_bc_matrix.h5"))
  metadata <- read.csv(file = file.path(outs, "singlecell.csv"), header = TRUE, row.names = 1)

  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = file.path(outs, "fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
  )

  obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
  obj$sample_id <- sid
  objs[[sid]] <- obj
}

saveRDS(objs, file = file.path("results", "signac_objects.rds"))
message("Saved Signac objects: results/signac_objects.rds")
