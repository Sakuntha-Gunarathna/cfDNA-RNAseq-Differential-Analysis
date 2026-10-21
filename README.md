# SARTools + edgeR — Cancer vs Healthy example

This is a **general example pipeline** for differential analysis of unpaired Cancer and Healthy samples using **raw integer counts**. The R script is organized into numbered output blocks and includes editable colors. It uses SARTools to run edgeR with TMM normalization, CPM cutoff 0.1, and Benjamini–Hochberg FDR threshold 0.05.

## Files

- `SARTools_edgeR_Cancer_vs_Healthy.R` — analysis script. Edit settings and plot colors at the top.
- `counts.csv` — **400 SYNTHETIC regions × 6 samples** for demonstrating input format and test-running the pipeline. **Not biological data.**
- `target.txt` — sample labels/groups matched to synthetic counts.
- `install_dependencies.R` — optional one-time installation of required libraries.

## Running

1. Place the R script, `counts.csv`, and `target.txt` in the same folder.
2. In R, set the working directory to that folder.
3. Install dependencies once if the required packages are not already installed:

   ```r
   source("install_dependencies.R")
   ```

4. Run:

   ```r
   source("SARTools_edgeR_Cancer_vs_Healthy.R")
   ```

5. Review `Results_Cancer_vs_Healthy/` for tables, MDS PDF/PNG plots, SARTools diagnostics, and an optional HTML report (if generation succeeds). **This is an example package; results on synthetic data have no biological meaning.**

To use your own data, replace `counts.csv` with a raw integer-count CSV containing feature IDs in the first column and one uniquely named sample per following column. Supply `target.txt` as a TAB-separated file containing `label` and `group`; sample labels must exactly match the CSV column headers. The script generates the SARTools `files` column itself. For the unchanged example settings, the only group names are `Healthy` and `Cancer`, with at least two independent biological replicates per group.

## Results

- `tables/Complete_results.csv` and group-higher FDR tables; `tables/TMM_factors.csv` and `tables/R_sessionInfo.txt`.
- `figures/MDS_with_labels.{pdf,png}` and `figures/MDS_without_labels.{pdf,png}`.
- SARTools QC/summary figures and HTML report.

`Healthy` is the reference condition; positive log2 fold change means higher in Cancer. Adjust colors in `GROUP_COLORS` and font in `FONT_FAMILY` near the top of the R file.

## Scope and validation

**Study design:** independent, unpaired samples only (not a matched before/after model). Inputs must be raw, unnormalized nonnegative integer counts, not fractions, CPM, or normalized signal.

The intended SARTools workflow follows the developer-provided edgeR template for SARTools 1.8.2; versions should be recorded via `sessionInfo()`.

SARTools source and documentation: https://github.com/PF2-pasteur-fr/SARTools
