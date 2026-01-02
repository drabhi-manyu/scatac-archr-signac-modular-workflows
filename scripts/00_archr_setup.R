suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(yaml)
  library(ggplot2)
})

cfg <- yaml::read_yaml("config/config.yaml")

setwd(cfg$work_dir)

addArchRGenome(cfg$genome)
addArchRThreads(threads = cfg$threads)

out_base <- cfg$out_base
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
