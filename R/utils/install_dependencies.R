#!/usr/bin/env Rscript

# Install required ProteoLab package dependencies from DESCRIPTION.
# By default this installs non-base Imports only.
# Use --with-suggests to also install Suggests entries.
dependencies<-c(
  "arrow",
  "assertthat",
  "ComplexHeatmap",
  "ComplexUpset",
  "circlize",
  "clusterProfiler",
  "cowplot",
  "decoupleR",
  "data.table",
  "dplyr",
  "furrr",
  "factoextra",
  "fgsea",
  "GGally",
  "ggplot2",
  "ggpubr",
  "ggsci",
  "grid",
  "here",
  "imputeLCMD",
  "janitor",
  "limma",
  "magrittr",
  "methods",
  "msigdbr",
  "OmnipathR",
  "plyr",
  "purrr",
  "patchwork",
  "reshape2",
  "Rtsne",
  "S4Vectors",
  "shiny",
  "stats",
  "stringr",
  "SummarizedExperiment",
  "sva",
  "tibble",
  "tidyr",
  "umap",
  "utils",
  "vsn",
  "EnhancedVolcano",
  "lubridate",
  "patchwork"
)

install_dependencies <- function(r_packages) {

  base_r_packages <- c("grid", "methods", "stats", "utils")
  required_pkgs <- setdiff(r_packages, base_r_packages)

  target_pkgs <- required_pkgs
  if (isTRUE(with_suggests)) {
    target_pkgs <- unique(c(target_pkgs, suggests))
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  options(repos = BiocManager::repositories())

  installed <- rownames(installed.packages())
  to_install <- setdiff(target_pkgs, installed)

  message(sprintf("Required (non-base Imports): %d", length(required_pkgs)))
  if (isTRUE(with_suggests)) {
    message(sprintf("Including Suggests, total requested: %d", length(target_pkgs)))
  }
  message(sprintf("Already installed: %d", length(setdiff(target_pkgs, to_install))))
  message(sprintf("To install: %d", length(to_install)))

  if (length(to_install) == 0L) {
    message("All requested dependencies are already installed.")
    return(invisible(TRUE))
  }

  bioc_available <- rownames(BiocManager::available())
  bioc_pkgs <- intersect(to_install, bioc_available)
  cran_pkgs <- setdiff(to_install, bioc_pkgs)

  if (length(cran_pkgs) > 0L) {
    message(sprintf("Installing CRAN packages (%d)...", length(cran_pkgs)))
    install.packages(cran_pkgs)
  }

  if (length(bioc_pkgs) > 0L) {
    message(sprintf("Installing Bioconductor packages (%d)...", length(bioc_pkgs)))
    BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
  }

  message("Dependency installation finished.")
  invisible(TRUE)
}

# Only execute when called as a script (Rscript R/install_dependencies.R ...).
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  with_suggests <- "--with-suggests" %in% args
  install_dependencies(r_packages = dependencies)
}
