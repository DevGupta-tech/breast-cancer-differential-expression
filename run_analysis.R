# GSE7904 breast cancer capstone
# Reproducible probe-level differential-expression analysis.
# Run from the project root: source("scripts/run_analysis.R")

options(stringsAsFactors = FALSE)

input_file <- "data/GSE7904_series_matrix.txt.gz"
if (!file.exists(input_file)) {
  stop("Cannot find data/GSE7904_series_matrix.txt.gz. Run this script from the project root.")
}

# Read GEO series-matrix text and locate its expression table.
con <- gzfile(input_file, open = "rt")
geo_lines <- readLines(con, warn = FALSE)
close(con)

begin <- match("!series_matrix_table_begin", geo_lines)
end <- match("!series_matrix_table_end", geo_lines)
if (is.na(begin) || is.na(end) || end <= begin) {
  stop("GEO series matrix table markers were not found.")
}

# Extract sample characteristics in the same order as matrix columns.
char_line_index <- grep("^!Sample_characteristics_ch1\\t", geo_lines)
if (length(char_line_index) != 1) {
  stop("Expected one !Sample_characteristics_ch1 line; inspect GEO metadata.")
}
char_line <- sub("^!Sample_characteristics_ch1\\t", "", geo_lines[char_line_index])
sample_labels <- strsplit(char_line, "\t", fixed = TRUE)[[1]]
sample_labels <- gsub('^"|"$', "", sample_labels)

table_text <- geo_lines[(begin + 1):(end - 1)]
tab <- read.delim(
  text = paste(table_text, collapse = "\n"),
  sep = "\t", header = TRUE, quote = '"',
  check.names = FALSE, stringsAsFactors = FALSE
)

if (ncol(tab) - 1 != length(sample_labels)) {
  stop("Sample metadata count does not match expression matrix columns.")
}

probe_ids <- tab[[1]]
raw_mat <- as.matrix(tab[, -1, drop = FALSE])
storage.mode(raw_mat) <- "numeric"
rownames(raw_mat) <- probe_ids

# Group definitions from the original GEO characteristics:
# Basal, BRCA1 and non-BLC = tumour; NB = normal breast; NO = excluded.
sample_group <- ifelse(
  sample_labels %in% c("Basal", "BRCA1", "non-BLC"), "Tumour",
  ifelse(sample_labels == "NB", "Normal", "Excluded")
)

sample_info <- data.frame(
  GSM_ID = colnames(raw_mat),
  Original_characteristic = sample_labels,
  Analysis_group = sample_group
)

keep <- sample_group %in% c("Tumour", "Normal")
group <- factor(sample_group[keep], levels = c("Normal", "Tumour"))
mat <- raw_mat[, keep, drop = FALSE]

if (sum(group == "Tumour") != 43 || sum(group == "Normal") != 7) {
  stop("Unexpected sample counts. Expected 43 tumour and 7 normal breast samples.")
}

# GEO labels these as dChip signal values. Some are non-positive, so floor
# values below 1 before log2 transformation. This preprocessing choice is
# documented in REPORT.md and should be considered a limitation.
log_mat <- log2(pmax(mat, 1))

# Vectorized two-sided Welch t-tests for each probe.
normal_idx <- which(group == "Normal")
tumour_idx <- which(group == "Tumour")
n0 <- length(normal_idx)
n1 <- length(tumour_idx)

normal_mat <- log_mat[, normal_idx, drop = FALSE]
tumour_mat <- log_mat[, tumour_idx, drop = FALSE]

mean0 <- rowMeans(normal_mat)
mean1 <- rowMeans(tumour_mat)
var0 <- apply(normal_mat, 1, var)
var1 <- apply(tumour_mat, 1, var)

a <- var0 / n0
b <- var1 / n1
se2 <- a + b
t_stat <- rep(0, nrow(log_mat))
valid_se <- is.finite(se2) & se2 > 0
t_stat[valid_se] <- (mean1[valid_se] - mean0[valid_se]) / sqrt(se2[valid_se])

denom <- (a^2 / (n0 - 1)) + (b^2 / (n1 - 1))
welch_df <- rep(Inf, nrow(log_mat))
valid_df <- is.finite(denom) & denom > 0
welch_df[valid_df] <- se2[valid_df]^2 / denom[valid_df]

p_value <- rep(1, nrow(log_mat))
p_value[valid_se] <- 2 * pt(-abs(t_stat[valid_se]), df = welch_df[valid_se])
zero_se <- is.finite(se2) & se2 == 0
p_value[zero_se & (mean1 != mean0)] <- 0
p_value[!is.finite(p_value)] <- NA_real_

# BH correction is applied across all tested probes.
adj_p <- p.adjust(p_value, method = "BH")

results <- data.frame(
  ProbeID = rownames(log_mat),
  Mean_log2_Normal = mean0,
  Mean_log2_Tumour = mean1,
  log2FoldChange_Tumour_minus_Normal = mean1 - mean0,
  Welch_t_statistic = t_stat,
  Welch_df = welch_df,
  P_value = p_value,
  BH_adjusted_P_value = adj_p
)
results$Significant_FDR_0.05 <- !is.na(adj_p) & adj_p < 0.05
results$Direction <- ifelse(
  results$log2FoldChange_Tumour_minus_Normal > 0, "Higher in tumour",
  ifelse(results$log2FoldChange_Tumour_minus_Normal < 0, "Lower in tumour", "No difference")
)
results <- results[order(results$P_value, results$ProbeID, na.last = TRUE), ]

significant <- results[results$Significant_FDR_0.05, , drop = FALSE]
n_up <- sum(significant$log2FoldChange_Tumour_minus_Normal > 0)
n_down <- sum(significant$log2FoldChange_Tumour_minus_Normal < 0)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write.csv(sample_info, "results/sample_classification.csv", row.names = FALSE)
write.csv(results, "results/all_probe_results.csv", row.names = FALSE)
write.csv(significant, "results/significant_probes_FDR_0.05.csv", row.names = FALSE)
write.csv(head(results, 20), "results/top_20_probes.csv", row.names = FALSE)

# Volcano plot
png("figures/volcano_plot.png", width = 1200, height = 850, res = 140)
y <- -log10(pmax(results$BH_adjusted_P_value, .Machine$double.xmin))
cols <- ifelse(results$Significant_FDR_0.05, "firebrick", "grey65")
plot(
  results$log2FoldChange_Tumour_minus_Normal, y,
  pch = 20, cex = 0.55, col = cols,
  xlab = "Log2 fold change (Tumour - Normal)",
  ylab = "-Log10(BH-adjusted p-value)",
  main = "GSE7904 differential expression (probe-level)"
)
abline(h = -log10(0.05), lty = 2)
abline(v = 0, lty = 3)
legend("topright", legend = c("FDR < 0.05", "Not significant"),
       col = c("firebrick", "grey65"), pch = 20, bty = "n")
dev.off()

# PCA of the 500 most variable probes; exploratory QC only.
probe_var <- apply(log_mat, 1, var)
top_idx <- order(probe_var, decreasing = TRUE)[seq_len(min(500, length(probe_var)))]
pca <- prcomp(t(log_mat[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
pc_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

png("figures/pca_samples.png", width = 1100, height = 850, res = 140)
plot(
  pca$x[, 1], pca$x[, 2],
  col = ifelse(group == "Tumour", "firebrick", "steelblue"),
  pch = 19, xlab = paste0("PC1 (", round(pc_var[1], 2), "%)"),
  ylab = paste0("PC2 (", round(pc_var[2], 2), "%)"),
  main = "GSE7904 PCA: top 500 variable probes"
)
legend("topright", legend = c("Tumour", "Normal breast"),
       col = c("firebrick", "steelblue"), pch = 19, bty = "n")
dev.off()

summary_lines <- c(
  "# GSE7904 analysis summary", "",
  "## Sample groups",
  paste0("- Tumour: ", n1),
  paste0("- Normal breast: ", n0),
  paste0("- Excluded NO samples: ", sum(sample_group == "Excluded")),
  "",
  "## Method",
  "- Processed dChip signal values; values below 1 floored at 1 before log2 transformation.",
  "- Probe-level two-sided Welch t-tests.",
  "- Benjamini-Hochberg multiple-testing correction.",
  "- Significance threshold: BH-adjusted p-value < 0.05.",
  "",
  "## Results",
  paste0("- Probes tested: ", nrow(results)),
  paste0("- Significant probes: ", nrow(significant)),
  paste0("- Higher in tumour: ", n_up),
  paste0("- Lower in tumour: ", n_down),
  "",
  "See REPORT.md for interpretation, confounders and limitations."
)
writeLines(summary_lines, "results/analysis_summary.md")

cat("\nANALYSIS COMPLETE\n")
cat("Tumour samples:", n1, "\n")
cat("Normal breast samples:", n0, "\n")
cat("Probes tested:", nrow(results), "\n")
cat("Significant probes (BH FDR < 0.05):", nrow(significant), "\n")
cat("Higher in tumour:", n_up, "\n")
cat("Lower in tumour:", n_down, "\n")
cat("Outputs saved in results/ and figures/.\n")
