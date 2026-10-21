# One-time installation for SARTools_edgeR_Cancer_vs_Healthy.R.
# Run only if required libraries are not already installed.
# Internet access is required for first-time installation.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran_packages <- c("ggplot2", "ggforce", "remotes", "knitr", "xtable", "rmarkdown")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

bioc_packages <- c("edgeR", "DESeq2", "genefilter")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("SARTools", quietly = TRUE))
  remotes::install_github("PF2-pasteur-fr/SARTools", build_opts = "--no-resave-data", upgrade = "never")

for (pkg in c("SARTools", "edgeR", "ggplot2", "ggforce")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required package was not installed successfully: ", pkg)
}
message("Required analysis packages are available. Run the main R script next.")
