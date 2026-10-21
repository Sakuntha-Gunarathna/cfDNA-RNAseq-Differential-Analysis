
###############################################################################
# SARTools + edgeR: differential counts (Cancer vs Healthy example)
# Put counts.csv and target.txt in your R working directory, then run this.
# counts.csv: feature ID in first column; raw integer counts in other columns.
# target.txt: tab-separated columns label and group (sample labels must match).
# Designed for UNPAIRED samples. Required packages must already be installed.
###############################################################################

###############################################################################
# BLOCK 1: USER SETTINGS AND COLORS
###############################################################################
COUNTS_FILE <- "counts.csv"
TARGET_FILE <- "target.txt"
PROJECT_NAME <- "Cancer_vs_Healthy"
AUTHOR <- "Sakuntha Gunarathna"
REFERENCE_GROUP <- "Healthy"      # Positive log2FC means higher in Cancer
COMPARISON_GROUP <- "Cancer"
CPM_CUTOFF <- 0.1
FDR_CUTOFF <- 0.05
FONT_FAMILY <- "sans"              # Change to "Arial" if installed
SHOW_ELLIPSES <- FALSE            # Optional; works best with 3+ samples/group
MAKE_HTML_REPORT <- TRUE          # Set FALSE to omit optional SARTools report
GROUP_COLORS <- c(Healthy = "#0067A5", Cancer = "#BE0032")

###############################################################################
# BLOCK 2: LOAD INSTALLED LIBRARIES
###############################################################################
library(SARTools)
library(edgeR)
library(limma)
library(ggplot2)
library(ggforce)

# Keep directory changes inside a function; restore directory on any error.
run_analysis <- function() {
  start_dir <- getwd()
  on.exit(setwd(start_dir), add = TRUE)

  #############################################################################
  # BLOCK 3: IMPORT AND CHECK INPUTS
  #############################################################################
  if (!file.exists(COUNTS_FILE) || !file.exists(TARGET_FILE))
    stop("Input files not found in ", start_dir, ". Use setwd() first.")
  x <- read.csv(COUNTS_FILE, check.names = FALSE)
  if (ncol(x) < 3 || nrow(x) < 2) stop("Invalid counts.csv dimensions.")
  ids <- as.character(x[[1]])
  x <- x[, -1, drop = FALSE]
  samples <- names(x)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Feature IDs must be unique and nonempty.")
  if (anyDuplicated(samples) || any(!nzchar(samples)))
    stop("Sample names must be unique and nonempty.")
  if (!all(vapply(x, is.numeric, logical(1))))
    stop("Sample count columns must be numeric.")
  m <- as.matrix(x)
  if (any(!is.finite(m)) || any(m < 0) || any(m != round(m)))
    stop("Input counts must be nonnegative integers, not normalized counts.")
  rownames(x) <- ids

  t <- read.delim(TARGET_FILE, stringsAsFactors = FALSE)
  if (!all(c("label", "group") %in% names(t)))
    stop("target.txt needs the columns 'label' and 'group'.")
  if (anyDuplicated(t$label) || !setequal(samples, t$label))
    stop("Sample labels in counts.csv and target.txt do not match.")
  t <- t[match(samples, t$label), , drop = FALSE]
  groups <- c(REFERENCE_GROUP, COMPARISON_GROUP)
  if (length(unique(groups)) != 2 || !setequal(unique(t$group), groups))
    stop("The groups in target.txt must match BLOCK 1 settings.")
  if (any(table(t$group) < 2)) stop("At least two replicates/group are needed.")
  if (!all(groups %in% names(GROUP_COLORS)))
    stop("GROUP_COLORS needs a named color for each group.")

  #############################################################################
  # BLOCK 4: CREATE OUTPUT FOLDERS AND SARTools COUNT FILES
  #############################################################################
  projectName <- gsub("[^A-Za-z0-9_-]", "_", PROJECT_NAME)
  workDir <- file.path(start_dir, paste0("Results_", projectName))
  for (d in c("", "raw_counts", "tables", "figures"))
    dir.create(file.path(workDir, d), recursive = TRUE, showWarnings = FALSE)
  t$files <- sprintf("Sample_%02d.txt", seq_along(samples))
  for (i in seq_along(samples)) {
    write.table(data.frame(ids, x[[i]]),
                file = file.path(workDir, "raw_counts", t$files[i]),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  write.table(t[, c("label", "files", "group")],
              file = file.path(workDir, "target.txt"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  setwd(workDir)

  targetFile <- "target.txt"; rawDir <- "raw_counts"
  varInt <- "group"; condRef <- REFERENCE_GROUP; batch <- NULL
  alpha <- FDR_CUTOFF; cpmCutoff <- CPM_CUTOFF
  pAdjustMethod <- "BH"; normalizationMethod <- "TMM"
  gene.selection <- "pairwise"; featuresToRemove <- character(0)
  colors <- unname(GROUP_COLORS[groups]); author <- AUTHOR

  #############################################################################
  # BLOCK 5: SARTools / edgeR DIFFERENTIAL ANALYSIS
  #############################################################################
  checkParameters.edgeR(projectName = projectName, author = author,
    targetFile = targetFile, rawDir = rawDir,
    featuresToRemove = featuresToRemove, varInt = varInt,
    condRef = condRef, batch = batch, alpha = alpha,
    pAdjustMethod = pAdjustMethod, cpmCutoff = cpmCutoff,
    gene.selection = gene.selection, normalizationMethod = normalizationMethod,
    colors = colors)

  target <- loadTargetFile(targetFile = targetFile, varInt = varInt,
                           condRef = condRef, batch = batch)
  counts <- loadCountData(target = target, rawDir = rawDir,
                          featuresToRemove = featuresToRemove)
  out.edgeR <- run.edgeR(counts = counts, target = target, varInt = varInt,
    condRef = condRef, batch = batch, cpmCutoff = cpmCutoff,
    normalizationMethod = normalizationMethod, pAdjustMethod = pAdjustMethod)

  saveRDS(out.edgeR, "tables/edgeR_model.rds")
  write.csv(out.edgeR$dge$samples, "tables/TMM_factors.csv")

  #############################################################################
  # BLOCK 6: COMPLETE AND SIGNIFICANT RESULTS TABLES
  #############################################################################
  all_results <- exportResults.edgeR(out.edgeR, group = target[, varInt],
                                      counts = counts, alpha = alpha, export = FALSE)
  if (!is.list(all_results) || length(all_results) != 1 ||
      !is.data.frame(all_results[[1]])) stop("Unexpected SARTools results format.")
  results <- all_results[[1]]
  if (!all(c("log2FoldChange", "pvalue", "padj") %in% names(results)))
    stop("Differential table is missing expected columns.")
  if (!"Feature_ID" %in% names(results)) {
    if ("Id" %in% names(results)) {
      names(results)[names(results) == "Id"] <- "Feature_ID"
    } else if ("ID" %in% names(results)) {
      names(results)[names(results) == "ID"] <- "Feature_ID"
    } else {
      results <- data.frame(Feature_ID = rownames(results), results,
                            row.names = NULL, check.names = FALSE)
    }
  }
  results <- results[order(results$padj, results$pvalue, na.last = TRUE), , drop = FALSE]
  sig <- !is.na(results$padj) & is.finite(results$padj) & results$padj < alpha
  fc <- results$log2FoldChange
  comparison_higher <- results[sig & is.finite(fc) & fc > 0, , drop = FALSE]
  reference_higher <- results[sig & is.finite(fc) & fc < 0, , drop = FALSE]
  write.csv(results, "tables/Complete_results.csv", row.names = FALSE)
  write.csv(comparison_higher, "tables/Cancer_higher_FDR.csv", row.names = FALSE)
  write.csv(reference_higher, "tables/Healthy_higher_FDR.csv", row.names = FALSE)
  message("Differential analysis and CSV exports completed.")

  # Optional plotting/reporting failures are logged without stopping analysis.
  optional_errors <- character(0)
  optional <- function(block, action) tryCatch(force(action), error = function(e) {
    optional_errors <<- c(optional_errors, paste(block, conditionMessage(e), sep = ": "))
    message("Optional output failed [", block, "]: ", conditionMessage(e))
    NULL
  })
  pdf_output <- function(path, action) {
    pdf(path, width = 7, height = 6)
    on.exit(dev.off(), add = TRUE)
    action()
  }

  #############################################################################
  # BLOCK 7: MDS PLOTS (WITH AND WITHOUT SAMPLE LABELS)
  #############################################################################
  optional("MDS", {
    # IMPORTANT: use limma::plotMDS, not edgeR::plotMDS.
    mds <- limma::plotMDS(out.edgeR$dge, plot = FALSE,
      top = min(500, nrow(out.edgeR$dge$counts)),
      gene.selection = gene.selection)
    sn <- colnames(out.edgeR$dge$counts)
    df <- data.frame(MDS1 = mds$x, MDS2 = mds$y, sample = sn,
      group = as.character(target[match(sn, target$label), varInt]))
    write.csv(df, "tables/MDS_coordinates.csv", row.names = FALSE)
    p <- ggplot(df, aes(MDS1, MDS2, color = group)) +
      geom_point(size = 3) + scale_color_manual(values = GROUP_COLORS) +
      labs(title = "MDS Plot", x = "Leading logFC dimension 1",
           y = "Leading logFC dimension 2") +
      theme_classic(base_size = 12, base_family = FONT_FAMILY)
    if (SHOW_ELLIPSES && all(table(df$group) >= 3))
      p <- p + ggforce::geom_mark_ellipse(aes(fill = group, group = group),
                alpha = 0.1, show.legend = FALSE) +
                scale_fill_manual(values = GROUP_COLORS)
    labeled <- p + geom_text(aes(label = sample), vjust = -1, size = 3)
    ggsave("figures/MDS_with_labels.pdf", labeled, width = 6, height = 6)
    ggsave("figures/MDS_with_labels.png", labeled, width = 6, height = 6, dpi = 600)
    ggsave("figures/MDS_without_labels.pdf", p, width = 6, height = 6)
    ggsave("figures/MDS_without_labels.png", p, width = 6, height = 6, dpi = 600)
  })

  #############################################################################
  # BLOCK 8: SARTools DESCRIPTION PLOTS
  #############################################################################
  majSequences <- optional("Description", pdf_output(
    "figures/SARTools_description.pdf",
    function() descriptionPlots(counts = counts, group = target[, varInt], col = colors)))

  #############################################################################
  # BLOCK 9: SARTools CLUSTERING
  #############################################################################
  optional("Clustering", pdf_output("figures/SARTools_clustering.pdf",
    function() exploreCounts(object = out.edgeR$dge, group = target[, varInt],
                             gene.selection = gene.selection, col = colors)))

  #############################################################################
  # BLOCK 10: SARTools SUMMARY FIGURES / TABLES
  #############################################################################
  summaryResults <- optional("Summary", pdf_output("figures/SARTools_summary.pdf",
    function() summarizeResults.edgeR(out.edgeR, group = target[, varInt],
                                      counts = counts, alpha = alpha, col = colors)))

  #############################################################################
  # BLOCK 11: OPTIONAL SARTools HTML REPORT
  #############################################################################
  if (MAKE_HTML_REPORT && !is.null(majSequences) && !is.null(summaryResults)) {
    optional("HTML report", writeReport.edgeR(
      target = target, counts = counts, out.edgeR = out.edgeR,
      summaryResults = summaryResults, majSequences = majSequences,
      workDir = workDir, projectName = projectName, author = author,
      targetFile = targetFile, rawDir = rawDir,
      featuresToRemove = featuresToRemove, varInt = varInt,
      condRef = condRef, batch = batch, alpha = alpha,
      pAdjustMethod = pAdjustMethod, cpmCutoff = cpmCutoff,
      colors = colors, gene.selection = gene.selection,
      normalizationMethod = normalizationMethod))
  } else if (MAKE_HTML_REPORT) {
    message("HTML report skipped because a SARTools diagnostic did not complete.")
  }

  #############################################################################
  # BLOCK 12: SAVE SUMMARY AND SESSION INFORMATION
  #############################################################################
  write.csv(data.frame(Comparison = paste(COMPARISON_GROUP, "vs", REFERENCE_GROUP),
    Features = nrow(results), Significant_FDR = sum(sig),
    Cancer_higher = nrow(comparison_higher), Healthy_higher = nrow(reference_higher)),
    "tables/Analysis_summary.csv", row.names = FALSE)
  writeLines(capture.output(sessionInfo()), "tables/R_sessionInfo.txt")
  if (length(optional_errors)) {
    writeLines(optional_errors, "tables/Optional_output_errors.txt")
    message("Some optional steps failed; see tables/Optional_output_errors.txt")
  } else {
    message("All steps completed.")
  }
  message("Results saved to: ", workDir)
  invisible(workDir)
}

run_analysis()

