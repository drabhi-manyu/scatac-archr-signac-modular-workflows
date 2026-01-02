source("scripts/00_archr_setup.R")

proj <- loadArchRProject(file.path(out_base, "01_save_proj"))

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = cfg$lsi_iterations,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = cfg$lsi_var_features,
  dimsToUse = seq(cfg$lsi_dims[1], cfg$lsi_dims[2]),
  force = TRUE
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  resolution = cfg$cluster_resolution,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

# Harmony
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = cfg$harmony_group_by
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_Harmony",
  resolution = cfg$cluster_resolution,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

plot_dir <- file.path(out_base, "02_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters_Harmony", embedding = "UMAP_Harmony")

ggsave(file.path(plot_dir, "umap_harmony_by_sample.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "umap_harmony_by_clusters.png"), p2, width = 6, height = 5, dpi = 300)

# Filter very small clusters if desired
tab <- table(proj$Clusters_Harmony)
keep <- names(tab[tab >= cfg$min_cluster_size])
cells_to_keep <- proj$cellNames[proj$Clusters_Harmony %in% keep]
proj <- proj[cells_to_keep, ]

saveArchRProject(proj, outputDirectory = file.path(out_base, "02_save_harmony"), load = FALSE)
message("Saved: ", file.path(out_base, "02_save_harmony"))
