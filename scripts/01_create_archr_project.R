source("scripts/00_archr_setup.R")

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
})

samples <- readr::read_csv("config/samples.csv", show_col_types = FALSE)

arrow_files <- samples$arrow_file
stopifnot(all(file.exists(arrow_files)))

sample_names <- samples$sample_id

proj_dir <- file.path(out_base, "01_proj_archr")
proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = proj_dir,
  copyArrows = TRUE
)

# Optional metadata attachment
if (!is.null(cfg$metadata_xlsx) && file.exists(cfg$metadata_xlsx)) {
  meta <- readxl::read_excel(cfg$metadata_xlsx)

  sid_col <- cfg$metadata_sample_id_col
  meta[[sid_col]] <- as.character(meta[[sid_col]])

  # Example: extract a short sample id from proj$Sample if needed
  proj$Sample_ID <- as.character(proj$Sample)

  meta <- meta[meta[[sid_col]] %in% unique(proj$Sample_ID), , drop = FALSE]

  for (field in cfg$metadata_fields) {
    if (field %in% colnames(meta)) {
      m <- setNames(meta[[field]], meta[[sid_col]])
      proj[[field]] <- m[proj$Sample_ID]
    }
  }
}

# QC plots
qc_dir <- file.path(out_base, "01_qc_plots")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

p1 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
ggsave(file.path(qc_dir, "tss_enrichment_violin.png"), p1, width = 6, height = 5, dpi = 300)

p2 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
ggsave(file.path(qc_dir, "nfrags_violin.png"), p2, width = 6, height = 5, dpi = 300)

saveArchRProject(proj, outputDirectory = file.path(out_base, "01_save_proj"), load = FALSE)
message("Saved: ", file.path(out_base, "01_save_proj"))
