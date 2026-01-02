source("scripts/00_archr_setup.R")

suppressPackageStartupMessages({
  library(readr)
})

proj <- loadArchRProject(file.path(out_base, "02_save_harmony"))

# External annotation CSV schema example:
# cell_id,barcode,predicted_label
# Where barcode is the 10x barcode (ends with -1), and cell_id can be ignored.
anno_csv <- file.path("data", "external_annotations.csv")
if (!file.exists(anno_csv)) {
  stop("Provide data/external_annotations.csv with columns: barcode, predicted_label")
}

annotations <- readr::read_csv(anno_csv, show_col_types = FALSE)

# Extract barcode from ArchR cell names, assumes pattern: <anything>#<10x_barcode>
proj$barcode <- sub(".*#", "", rownames(proj@cellColData))

match_idx <- match(proj$barcode, annotations$barcode)
proj$external_label <- annotations$predicted_label[match_idx]
proj$external_label[is.na(proj$external_label)] <- "Unknown"

p <- plotEmbedding(proj, colorBy = "cellColData", name = "external_label", embedding = "UMAP_Harmony")
plot_dir <- file.path(out_base, "03_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plot_dir, "umap_external_label.png"), p, width = 6, height = 5, dpi = 300)

saveArchRProject(proj, outputDirectory = file.path(out_base, "03_save_annotated"), load = FALSE)
