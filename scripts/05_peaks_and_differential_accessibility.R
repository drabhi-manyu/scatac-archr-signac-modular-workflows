source("scripts/00_archr_setup.R")

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(openxlsx)
})

proj <- loadArchRProject(file.path(out_base, "04_save_downstream"))

# Choose the grouping column for peak calling and DA
group_col <- cfg$group_by_for_peaks
if (!group_col %in% colnames(proj@cellColData)) {
  stop("Grouping column not found in cellColData: ", group_col)
}

# Add group coverages
proj <- addGroupCoverages(proj, groupBy = group_col)

# Peak calling
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = group_col,
  pathToMacs2 = cfg$path_to_macs2
)

proj <- addPeakMatrix(proj)

saveArchRProject(proj, outputDirectory = file.path(out_base, "05_save_with_peaks"), load = FALSE)

# Differential accessibility
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = group_col,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = cfg$da_test_method
)

saveRDS(markersPeaks, file = file.path(out_base, "markers_peaks_da.rds"))

# Export annotated peaks per group
peakSet <- getPeakSet(proj)
peakAnno <- as.data.frame(peakSet, row.names = NULL)
colnames(peakAnno)[1:3] <- c("seqnames", "start", "end")

out_dir <- file.path(out_base, "05_da_peaks_by_group")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()
groups <- unique(proj[[group_col]])

for (g in groups) {
  da_gr <- getMarkers(markersPeaks, cutOff = cfg$da_cutoff, returnGR = TRUE)[[as.character(g)]]
  if (length(da_gr) == 0) next

  da_df <- as.data.frame(da_gr)
  merged_df <- merge(
    da_df,
    peakAnno[, c("seqnames", "start", "end", "nearestGene", "distToTSS", "peakType")],
    by = c("seqnames", "start", "end"),
    all.x = TRUE
  )

  write.csv(merged_df, file = file.path(out_dir, paste0("DA_Peaks_", gsub("[^A-Za-z0-9]", "_", g), ".csv")), row.names = FALSE)

  sheet <- substr(gsub("[^A-Za-z0-9]", "_", g), 1, 31)
  addWorksheet(wb, sheetName = sheet)
  writeData(wb, sheet = sheet, x = merged_df)
}

saveWorkbook(wb, file = file.path(out_base, "DA_Peaks_Annotated_AllGroups.xlsx"), overwrite = TRUE)
message("Done. DA outputs in: ", out_dir)
