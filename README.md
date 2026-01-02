# atac-archr-signac-modular-workflows

Modular ATAC-seq analysis portfolio repo built from real working scripts, refactored to be GitHub safe.

This repo is intentionally organized as 6 analysis flavors you can run independently depending on what you want to demonstrate.

## Inputs
Option A: ArchR starting point
- Arrow files for each sample (produced elsewhere)
- Optional: a metadata table keyed by sample id

Option B: Cell Ranger ATAC outputs for Signac
- `filtered_peak_bc_matrix.h5`
- `singlecell.csv`
- `fragments.tsv.gz`

## Outputs
- Versioned ArchR project directories per step
- QC plots (TSS enrichment, fragments, cell counts)
- Harmony corrected embeddings, clusters
- Differential accessibility peaks tables
- Motif enrichment and deviation outputs
- Optional: transfer labels from RNA to ATAC then carry them into ArchR

## Configuration
- `config/config.yaml`: analysis settings and placeholder paths
- `config/samples.csv`: sample level paths and sample ids

## Analysis flavors (5 to 6 to pin as submodules)
Flavor 1: Build project and attach metadata
- `scripts/01_create_archr_project.R`

Flavor 2: LSI, UMAP, Harmony, clustering, filter small clusters
- `scripts/02_lsi_harmony_clustering.R`

Flavor 3: Add external annotations (for example from Signac label transfer)
- `scripts/03_add_external_annotations.R`

Flavor 4: Re-embed after filtering to annotated cells only, make clean UMAP plots
- `scripts/04_downstream_reembed_and_plots.R`

Flavor 5: Peak calling and differential accessibility
- `scripts/05_peaks_and_differential_accessibility.R`

Flavor 6: Motif enrichment and chromVAR deviations
- `scripts/06_motifs_and_deviations.R`

## Optional Signac helper
If you want a small Signac based entry point from Cell Ranger ATAC outputs, use:
- `scripts/signac_00_build_batch_objects.R`

## Quick start
1. Edit `config/config.yaml` and `config/samples.csv`
2. Run flavors in order, or just run the subset you want to showcase.
