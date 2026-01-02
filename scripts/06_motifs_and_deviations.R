source("scripts/00_archr_setup.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(cowplot)
})

proj <- loadArchRProject(file.path(out_base, "05_save_with_peaks"))

markersPeaks_path <- file.path(out_base, "markers_peaks_da.rds")
if (!file.exists(markersPeaks_path)) stop("Missing markers_peaks_da.rds from flavor 5")
markersPeaks <- readRDS(markersPeaks_path)

# Motif annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = cfg$motif_set, name = "Motif")

# Motif enrichment on DA peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = cfg$motif_cutoff
)

hm <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
plot_dir <- file.path(out_base, "06_motif_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

plotPDF(hm, name = "motif_enrichment_heatmap.pdf", ArchRProj = proj, width = 8, height = 6)

# Deviations (chromVAR style)
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif", force = TRUE)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "variable_motif_deviation_scores.pdf", ArchRProj = proj, width = 5, height = 5)

saveArchRProject(proj, outputDirectory = file.path(out_base, "06_save_with_motifs"), load = FALSE)
message("Saved: ", file.path(out_base, "06_save_with_motifs"))
