source("scripts/00_archr_setup.R")

proj <- loadArchRProject(file.path(out_base, "03_save_annotated"))

# Keep only labeled cells if you want a clean downstream view
proj2 <- proj[proj$external_label != "Unknown", ]

proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_2",
  iterations = 2,
  clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
  varFeatures = 15000,
  dimsToUse = 1:30,
  force = TRUE
)

proj2 <- addHarmony(proj2, reducedDims = "IterativeLSI_2", name = "Harmony_2", groupBy = "Sample")

proj2 <- addUMAP(
  ArchRProj = proj2,
  reducedDims = "Harmony_2",
  name = "UMAP_Harmony2",
  nNeighbors = 30,
  minDist = 0.5,
  dimsToUse = 1:30,
  force = TRUE
)

proj2 <- addClusters(
  input = proj2,
  reducedDims = "Harmony_2",
  method = "Seurat",
  name = "Clusters_Harmony2",
  resolution = cfg$cluster_resolution,
  force = TRUE
)

plot_dir <- file.path(out_base, "04_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

p1 <- plotEmbedding(proj2, colorBy = "cellColData", name = "external_label", embedding = "UMAP_Harmony2")
p2 <- plotEmbedding(proj2, colorBy = "cellColData", name = "Clusters_Harmony2", embedding = "UMAP_Harmony2")

ggsave(file.path(plot_dir, "umap_harmony2_by_label.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "umap_harmony2_by_clusters.png"), p2, width = 6, height = 5, dpi = 300)

saveArchRProject(proj2, outputDirectory = file.path(out_base, "04_save_downstream"), load = FALSE)
