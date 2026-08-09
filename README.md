# scATAC-seq Analysis with ArchR and Signac

A modular and reproducible workflow for single-cell ATAC-seq analysis using **ArchR** and **Signac**, adapted from pipelines developed for research datasets.

The repository demonstrates major stages of a typical scATAC-seq analysis, including project construction, dimensionality reduction, batch integration, clustering, annotation, peak calling, differential accessibility analysis, transcription-factor motif analysis, chromVAR-based deviations, and optional RNA-to-ATAC label transfer.

The workflow is organized into independent modules so that individual components can be run, adapted, or extended depending on the starting data and analytical objective.

## Workflow Overview

The main ArchR workflow includes:

1. Project construction and metadata integration
2. Iterative LSI dimensionality reduction
3. Harmony batch correction and integration
4. UMAP visualization and clustering
5. Filtering and re-embedding of selected populations
6. Integration of externally generated cell annotations
7. Peak calling and differential accessibility analysis
8. Transcription-factor motif enrichment and chromVAR deviation analysis

An optional Signac workflow provides an entry point from Cell Ranger ATAC outputs and can be used for RNA-to-ATAC annotation transfer before importing annotations into ArchR.

## Inputs

### Option A: Existing ArchR Arrow Files

The main ArchR workflow starts from:

* ArchR Arrow files generated for individual samples
* Optional sample-level metadata keyed by sample identifier

### Option B: Cell Ranger ATAC Outputs

The Signac helper workflow can start from standard Cell Ranger ATAC outputs, including:

* `filtered_peak_bc_matrix.h5`
* `singlecell.csv`
* `fragments.tsv.gz`

Input paths and analysis parameters are defined through configuration files rather than hard-coded into the analysis scripts.

## Configuration

Analysis settings and sample information are stored in:

* `config/config.yaml` — analysis parameters and input/output paths
* `config/samples.csv` — sample identifiers and sample-specific file paths

These files are intended to be edited before running the workflow on a new dataset.

## Analysis Modules

### 1. Create ArchR Project and Attach Metadata

`scripts/01_create_archr_project.R`

Initializes the ArchR project from existing Arrow files and attaches sample-level metadata required for downstream analyses.

Major tasks include:

* Loading Arrow files
* Creating the ArchR project
* Adding sample metadata
* Evaluating basic project-level information
* Saving the project for downstream analysis

### 2. Iterative LSI, Harmony Integration, Clustering, and UMAP

`scripts/02_lsi_harmony_clustering.R`

Performs dimensionality reduction and clustering while accounting for technical or sample-level variation.

Major tasks include:

* Iterative latent semantic indexing (LSI)
* Harmony-based batch integration
* Graph-based clustering
* UMAP embedding
* Evaluation of cluster sizes
* Filtering of small clusters when appropriate
* Generation of dimensionality-reduction plots
* Saving a versioned ArchR project

### 3. Add External Cell Annotations

`scripts/03_add_external_annotations.R`

Adds externally generated annotations to an existing ArchR project.

This can be used, for example, to incorporate cell-type labels generated through RNA-to-ATAC transfer in Signac or another annotation framework.

### 4. Re-embed Annotated Cells and Generate UMAPs

`scripts/04_downstream_reembed_and_plots.R`

Performs downstream dimensionality reduction after restricting the dataset to selected or successfully annotated cells.

Major tasks include:

* Filtering to annotated cell populations
* Recalculating reduced dimensions
* Generating updated UMAP embeddings
* Visualizing cell-type or metadata distributions
* Exporting publication-ready plots

### 5. Peak Calling and Differential Accessibility

`scripts/05_peaks_and_differential_accessibility.R`

Identifies accessible chromatin regions and tests for differential accessibility across biological groups or cell populations.

Major tasks include:

* Generation of group coverages
* Construction of reproducible peak sets
* Peak matrix generation
* Differential accessibility testing
* Correction for relevant technical covariates such as TSS enrichment and fragment depth
* Export of marker and differential-accessibility results

### 6. Motif Enrichment and chromVAR Deviations

`scripts/06_motifs_and_deviations.R`

Extends differential accessibility analysis to transcription-factor regulatory inference.

Major tasks include:

* Motif annotation
* Transcription-factor motif enrichment
* chromVAR-based motif deviation analysis
* Comparison of regulatory activity across cell populations or experimental conditions

## Optional Signac Workflow

`scripts/signac_00_build_batch_objects.R`

Provides a Signac-based entry point from Cell Ranger ATAC outputs.

This workflow can be used to:

* Construct sample-level Signac objects
* Perform initial quality-control processing
* Prepare datasets for downstream integration
* Support RNA-to-ATAC label transfer
* Generate annotations that can subsequently be incorporated into ArchR

## Outputs

Depending on the modules used, the workflow generates:

* Versioned ArchR project directories
* Quality-control summaries and plots
* Iterative LSI and Harmony-corrected embeddings
* Cluster assignments
* UMAP visualizations
* Externally transferred cell annotations
* Reproducible peak sets
* Differential accessibility tables
* Transcription-factor motif enrichment results
* chromVAR deviation scores
* Publication-ready figures
* CSV and other structured outputs for downstream interpretation

## Computational Environment

The workflows are written primarily in **R** and are designed for use in reproducible research-computing environments, including high-performance computing systems.

Core frameworks include:

* ArchR
* Signac
* Harmony
* chromVAR
* Seurat-based integration where applicable

The modular organization allows individual analysis stages to be modified without requiring the complete workflow to be rerun from the beginning.

## Reproducibility

The repository separates:

* Analysis code
* Configuration parameters
* Sample metadata
* Intermediate project states
* Downstream outputs

This structure is intended to facilitate reproducibility, troubleshooting, and adaptation of the workflow to independent scATAC-seq datasets.

Intermediate ArchR projects can be saved after major analytical stages, allowing analyses to resume from defined checkpoints rather than rebuilding the complete project for every downstream analysis.

## Quick Start

1. Clone the repository.
2. Edit `config/config.yaml` with the appropriate analysis parameters and paths.
3. Populate `config/samples.csv` with sample identifiers and input files.
4. Run the required analysis modules in numerical order.
5. Review quality-control and intermediate outputs before proceeding to downstream analyses.

The complete workflow can be run sequentially, or individual modules can be adapted depending on the available starting data and research question.

## Notes

This repository contains generalized versions of workflows developed and used for biological research. Dataset-specific file paths, identifiers, metadata, and restricted research data are not included.

The code is intended to demonstrate the analytical structure, computational approaches, and reproducible workflow design used for scATAC-seq analysis rather than provide a single dataset-specific pipeline.
