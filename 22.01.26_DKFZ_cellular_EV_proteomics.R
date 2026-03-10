################################################################################
############################## 1️⃣ FILTERING ###################################
################################################################################

library(clValid)
library(openxlsx)  # for saving Excel files

# --- Set working directory ---
setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics")

# --- Load and preprocess data ---
data <- read.table("DKFZ_HNSCC_massspec _Cellular_EV.txt", 
header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Check available columns ---
cat("Columns in dataset:\n")
print(colnames(data))

# --- Make sure target columns exist ---
required_cols <- c("Only identified by site", "Reverse", "Potential contaminant", 
"PG.Organisms", "PG.ProteinDescriptions")
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
warning(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# --- Apply filtering safely ---
data_filtered <- data

# Remove rows marked by MaxQuant flags
if ("Only identified by site" %in% colnames(data_filtered)) {
data_filtered <- data_filtered[data_filtered$`Only identified by site` != "+", ]
}
if ("Reverse" %in% colnames(data_filtered)) {
data_filtered <- data_filtered[data_filtered$Reverse != "+", ]
}
if ("Potential contaminant" %in% colnames(data_filtered)) {
data_filtered <- data_filtered[data_filtered$`Potential contaminant` != "+", ]
}

# Remove rows with "Potential contaminant" in PG.Organisms
if ("PG.Organisms" %in% colnames(data_filtered)) {
data_filtered <- data_filtered[!grepl("Potential contaminant", data_filtered$PG.Organisms, ignore.case = TRUE), ]
}

# Remove rows with "Bos taurus" in PG.ProteinDescriptions
if ("PG.ProteinDescriptions" %in% colnames(data_filtered)) {
data_filtered <- data_filtered[!grepl("Bos taurus", data_filtered$PG.ProteinDescriptions, ignore.case = TRUE, fixed = FALSE), ]
}

# --- Report summary ---
cat("\nSummary:\n")
cat("Original rows:", nrow(data), "\n")
cat("Filtered rows:", nrow(data_filtered), "\n")
cat("Removed rows:", nrow(data) - nrow(data_filtered), "\n")

# --- Reset row names ---
row.names(data_filtered) <- NULL

# --- Save filtered data (TXT) ---
write.table(data_filtered,
"DKFZ_HNSCC_massspec _Cellular_EV_filtered.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save filtered data (Excel) ---
write.xlsx(data_filtered,
"DKFZ_HNSCC_massspec _Cellular_EV_filtered.xlsx",
overwrite = TRUE)

cat("\nFiltering complete. Files saved.\n")


################################################################################
############################## 2️⃣ LOG10 TRANSFORM ##############################
################################################################################

# --- Identify LFQ intensity columns (original names) ---
lfq_cols <- grep("^LFQ_", colnames(data_filtered), value = TRUE)

# --- Copy before transforming ---
log10_data_filtered <- data_filtered

# --- Robust log10 transform of LFQ columns ---
log10_data_filtered[lfq_cols] <- lapply(data_filtered[lfq_cols], function(x) {
# 1) force to character, strip thousands separators, then to numeric
x_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
# 2) treat zeros/negatives as missing (will be imputed later)
x_num[x_num <= 0] <- NA_real_
# 3) log10 and clean non-finite
y <- log10(x_num)
y[!is.finite(y)] <- NA_real_
y
})

cat("Log10 transform complete on", length(lfq_cols), "LFQ columns.\n")

# --- Save (TXT) ---
write.table(log10_data_filtered,
"DKFZ_HNSCC_massspec _Cellular_EV_filtered_log10LFQ.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save (Excel) ---
write.xlsx(log10_data_filtered,
"DKFZ_HNSCC_massspec _Cellular_EV_filtered_log10LFQ.xlsx",
overwrite = TRUE)



################################################################################
############################## 4️⃣ QC VALIDATION ###############################
################################################################################

# --- Check: same number of rows & columns ---
cat("Rows (filtered):", nrow(log10_data_filtered), "\n")
cat("Rows (reordered):", nrow(log10_data_filtered), "\n")

cat("Columns (filtered):", ncol(log10_data_filtered), "\n")
cat("Columns (reordered):", ncol(log10_data_filtered), "\n")

# --- Verify identical content ignoring column order ---
same_values <- all.equal(
log10_data_filtered[, sort(colnames(log10_data_filtered))],
log10_data_filtered[, sort(colnames(log10_data_filtered))],
check.attributes = FALSE
)

if (isTRUE(same_values)) {
cat("\n✅ Data check passed: All numeric values match exactly (only column order changed).\n")
} else {
cat("\n⚠️ Data mismatch detected! Details:\n")
print(same_values)
}

# --- Optional: check that metadata columns match row-by-row ---
meta_cols <- c("PG.ProteinGroups", "PG.Genes", "PG.Organisms", "PG.ProteinDescriptions")
meta_same <- all.equal(log10_data_filtered[meta_cols], log10_data_filtered[meta_cols])
if (isTRUE(meta_same)) {
cat("✅ Metadata columns identical row-by-row.\n")
} else {
cat("⚠️ Metadata mismatch detected!\n")
print(meta_same)
}

################################################################################
################################################################################
################################################################################

################################################################################
############################## 5️⃣ 50% FILTERING ###############################
################################################################################

library(openxlsx)

# --- Prepare output folders ---
dir.create("SCC1", showWarnings = FALSE)
dir.create("SCC6", showWarnings = FALSE)
dir.create("SCC47", showWarnings = FALSE)
dir.create("SCC90", showWarnings = FALSE)
dir.create("Group_comparison", showWarnings = FALSE)

# --- Function: filter proteins with ≥50% valid values in either group ---
filter_50 <- function(df, group1_cols, group2_cols, meta_cols) {
keep <- apply(df, 1, function(row) {
g1 <- row[group1_cols]
g2 <- row[group2_cols]
valid_g1 <- sum(!is.na(g1))
valid_g2 <- sum(!is.na(g2))
# Keep if ≥50% non-NA in either group
valid_g1 >= length(g1) * 0.5 || valid_g2 >= length(g2) * 0.5
})
df[keep, c(meta_cols, group1_cols, group2_cols)]
}

# --- Metadata columns ---
meta_cols <- c("PG.ProteinGroups", "PG.Genes", "PG.Organisms", "PG.ProteinDescriptions")

# --- Detect LFQ columns ---
lfq_cols <- grep("^LFQ_", colnames(log10_data_filtered), value = TRUE)

# --- Helper function to extract columns by pattern ---
get_cols <- function(cell, dose) {
grep(paste0("LFQ_", cell, "_", dose, "_"), lfq_cols, value = TRUE)
}

################################################################################
# Individual cell line comparisons
################################################################################

cell_lines <- c("SCC1", "SCC6", "SCC47", "SCC90")

for (cell in cell_lines) {
cols_0Gy <- get_cols(cell, "0Gy")
cols_6Gy <- get_cols(cell, "6Gy")

if (length(cols_0Gy) == 0 | length(cols_6Gy) == 0) {
cat("⚠️ Skipping", cell, "- missing one of the dose groups.\n")
next
}

filtered_50percent <- filter_50(log10_data_filtered, cols_0Gy, cols_6Gy, meta_cols)

# --- Save outputs ---
write.table(filtered_50percent,
file.path(cell, paste0(cell, "_0Gy_vs_6Gy_filtered_50percent.txt")),
sep = "\t", quote = FALSE, row.names = FALSE)

write.xlsx(filtered_50percent,
file.path(cell, paste0(cell, "_0Gy_vs_6Gy_filtered_50percent.xlsx")),
overwrite = TRUE)

cat("✅", cell, ": filtered_50percent created (", nrow(filtered_50percent), "proteins kept)\n")
}

################################################################################
# Group comparison (all 0Gy vs 6Gy, excluding SCC6)
################################################################################

# remove SCC6 samples first
lfq_cols_noSCC6 <- lfq_cols[!grepl("SCC6", lfq_cols)]

# split into all 0Gy and 6Gy (excluding SCC6)
cols_0Gy_all <- grep("_0Gy_", lfq_cols_noSCC6, value = TRUE)
cols_6Gy_all <- grep("_6Gy_", lfq_cols_noSCC6, value = TRUE)

filtered_50percent <- filter_50(log10_data_filtered, cols_0Gy_all, cols_6Gy_all, meta_cols)

# --- Save outputs ---
write.table(filtered_50percent,
"Group_comparison/All_0Gy_vs_All_6Gy_noSCC6_filtered_50percent.txt",
sep = "\t", quote = FALSE, row.names = FALSE)

write.xlsx(filtered_50percent,
"Group_comparison/All_0Gy_vs_All_6Gy_noSCC6_filtered_50percent.xlsx",
overwrite = TRUE)

cat("✅ Group comparison (no SCC6): filtered_50percent created (", nrow(filtered_50percent), "proteins kept)\n")


################################################################################
################################################################################
################################################################################

################################################################################
############################## 6️⃣ CORRELATION QC ###############################
################################################################################

library(openxlsx)
library(pheatmap)
library(ggplot2)
library(dplyr)

################################################################################
# --- Function to analyse one results folder ---
################################################################################

analyze_correlation_folder <- function(folder) {
cat("\n=== Processing", folder, "===\n")

# Output subfolder (renamed to Step6)
out_dir <- file.path(folder, "Step6_Correlation_QC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Locate 50% filtered matrix (Step 5 output)
filtered_files <- list.files(folder, pattern = "_filtered_50percent\\.txt$", full.names = TRUE)
if (length(filtered_files) == 0) {
cat("⚠️ No _filtered_50percent.txt found in", folder, "\n")
return(NULL)
}
filtered_file <- filtered_files[1]
filtered_50 <- read.table(filtered_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Identify LFQ columns for this dataset
cols_group <- grep("^LFQ_", colnames(filtered_50), value = TRUE)
if (length(cols_group) < 2) {
cat("⚠️ Not enough samples for correlation in", folder, "\n")
return(NULL)
}

df_group <- filtered_50[cols_group]

# Convert numeric safely (important if LFQ are characters)
df_group <- as.data.frame(lapply(df_group, function(x) suppressWarnings(as.numeric(x))))

# --- Pearson correlation matrix ---
corr_mat <- cor(df_group, use = "pairwise.complete.obs", method = "pearson")

# --- Save numeric outputs ---
write.table(
corr_mat,
file.path(out_dir, paste0(folder, "_correlation_matrix.txt")),
sep = "\t", quote = FALSE, col.names = NA
)
write.xlsx(
corr_mat,
file.path(out_dir, paste0(folder, "_correlation_matrix.xlsx")),
overwrite = TRUE
)

# --- Plot heatmap (PNG + PDF) ---
heat_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
heat_name <- paste("Sample Correlation –", folder)
heat_fn_base <- file.path(out_dir, paste0(folder, "_correlation_heatmap"))

# PNG
pheatmap(
corr_mat,
main = heat_name, cluster_rows = TRUE, cluster_cols = TRUE,
display_numbers = TRUE, number_format = "%.2f",
color = heat_palette, fontsize_number = 8, border_color = NA,
filename = paste0(heat_fn_base, ".png"),
width = 8, height = 7, dpi = 300
)

# PDF (vector)
pdf(paste0(heat_fn_base, ".pdf"), width = 8, height = 7)
pheatmap(
corr_mat,
main = heat_name, cluster_rows = TRUE, cluster_cols = TRUE,
display_numbers = TRUE, number_format = "%.2f",
color = heat_palette, fontsize_number = 8, border_color = NA
)
dev.off()

# --- Mean correlation per sample (exclude self-corr ~1) ---
mean_corr <- apply(corr_mat, 1, function(x) mean(x[!is.na(x) & x < 0.999999], na.rm = TRUE))

# --- Classification ---
summary_df <- data.frame(
Sample = names(mean_corr),
Mean_Correlation = round(mean_corr, 3),
Status = case_when(
mean_corr >= 0.85 ~ "✅ OK (≥0.85)",
mean_corr >= 0.7  ~ "⚠ Borderline (0.7–0.85)",
TRUE              ~ "❌ Outlier (<0.7)"
),
stringsAsFactors = FALSE
)

# --- Save summary ---
write.table(
summary_df,
file.path(out_dir, paste0(folder, "_correlation_summary.txt")),
sep = "\t", quote = FALSE, row.names = FALSE
)
write.xlsx(
summary_df,
file.path(out_dir, paste0(folder, "_correlation_summary.xlsx")),
overwrite = TRUE
)

# --- Console report ---
n_ok      <- sum(summary_df$Mean_Correlation >= 0.85)
n_border  <- sum(summary_df$Mean_Correlation >= 0.7 & summary_df$Mean_Correlation < 0.85)
n_outlier <- sum(summary_df$Mean_Correlation < 0.7)

cat("Summary for", folder, ":\n")
cat("  ✅ OK:", n_ok, "\n")
cat("  ⚠ Borderline:", n_border, "\n")
cat("  ❌ Outliers:", n_outlier, "\n")

if (n_outlier > 0 || n_border > 0) print(summary_df)

return(summary_df)
}

################################################################################
# --- Run for each results folder ---
################################################################################

folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")
all_summaries <- lapply(folders, analyze_correlation_folder)

# Combine all summaries
summary_combined <- do.call(rbind, all_summaries)
if (!is.null(summary_combined)) {
write.xlsx(summary_combined, "All_Correlation_QC_Summary.xlsx", overwrite = TRUE)
}

cat("\n✅ Step 6 complete — correlation QC results (TXT, XLSX, PNG, PDF) saved in each folder’s Step6_Correlation_QC subfolder.\n")
cat("Combined summary: All_Correlation_QC_Summary.xlsx\n")



################################################################################
############################## 7️⃣ REMOVE LOW-CORR SAMPLES ######################
################################################################################

library(openxlsx)
library(dplyr)

# --- User setting: choose whether to filter based on correlation threshold ---
corr_filter    <- "yes"   # "yes" to filter, "no" to skip
corr_threshold <- 0.70    # correlation cutoff for exclusion

if (!corr_filter %in% c("yes", "no")) {
stop("⚠️ Please set corr_filter to either 'yes' or 'no'.")
}

apply_filter <- corr_filter == "yes"

if (apply_filter) {
cat("\n✅ Filtering enabled — samples with correlation < ", corr_threshold, " will be removed.\n", sep = "")
} else {
cat("\n⚠️ Filtering skipped — all samples will be retained (but still saved as new files).\n")
}

# --- Folders to process ---
folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

# --- Optional: track removed samples ---
removed_summary <- data.frame(
Folder = character(),
Removed_Sample = character(),
Correlation = numeric(),
stringsAsFactors = FALSE
)

for (folder in folders) {
cat("\n=== Processing", folder, "===\n")

# Path to correlation summary (from Step 6)
corr_summary_file <- file.path(folder, "Step6_Correlation_QC", paste0(folder, "_correlation_summary.txt"))
if (!file.exists(corr_summary_file)) {
cat("⚠️ No correlation summary found in", folder, "\n")
next
}
corr_summary <- read.table(corr_summary_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Identify keep/remove samples
if (apply_filter) {
keep_samples   <- corr_summary$Sample[corr_summary$Mean_Correlation >= corr_threshold]
remove_samples <- corr_summary$Sample[corr_summary$Mean_Correlation < corr_threshold]
} else {
keep_samples   <- corr_summary$Sample
remove_samples <- character(0)
}

# Report
if (apply_filter && length(remove_samples) > 0) {
cat("❌ Removing", length(remove_samples), "low-correlation samples (<", corr_threshold, "):\n")
print(remove_samples)

removed_summary <- rbind(
removed_summary,
data.frame(
Folder = folder,
Removed_Sample = remove_samples,
Correlation = corr_summary$Mean_Correlation[corr_summary$Sample %in% remove_samples],
stringsAsFactors = FALSE
)
)
} else if (apply_filter) {
cat("✅ All samples in", folder, "passed correlation threshold (≥", corr_threshold, ")\n")
} else {
cat("➡ Filtering skipped for", folder, "(but file will still be saved as _corrFiltered).\n")
}

# Locate the correct 50% filtered file (Step 5 output)
f50 <- list.files(folder, pattern = "_filtered_50percent\\.txt$", full.names = TRUE)
if (length(f50) == 0) {
cat("⚠️ No *_filtered_50percent.txt found in", folder, "\n")
next
}
f50 <- f50[1]

# Load the 50% filtered matrix
filtered_50percent <- read.table(f50, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Identify LFQ columns
lfq_cols <- grep("^LFQ_", colnames(filtered_50percent), value = TRUE)

# Keep only samples meeting criteria (or all if skipping)
keep_cols <- setdiff(colnames(filtered_50percent), lfq_cols)
keep_cols <- c(keep_cols, intersect(lfq_cols, keep_samples))

filtered_clean <- filtered_50percent[, keep_cols, drop = FALSE]

# --- Always save with a single suffix (_corrFiltered.[txt/xlsx]) ---
base_name <- gsub("\\.txt$", "", basename(f50))
suffix <- "_corrFiltered"

out_txt  <- file.path(folder, paste0(base_name, suffix, ".txt"))
out_xlsx <- file.path(folder, paste0(base_name, suffix, ".xlsx"))

write.table(filtered_clean, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(filtered_clean, out_xlsx, overwrite = TRUE)

cat("✅ Output saved to:", out_txt, "\n")
}

# --- Save summary of removed samples if any ---
if (apply_filter && nrow(removed_summary) > 0) {
write.xlsx(removed_summary, "Step7_Removed_Samples_Summary.xlsx", overwrite = TRUE)
cat("\n📄 Summary of removed samples saved to: Step7_Removed_Samples_Summary.xlsx\n")
}

cat("\n✅ Step 7 complete. Filtering applied:", ifelse(apply_filter, "YES", "NO"), "\n")
cat("All outputs saved with simplified names (_corrFiltered) for Step 8–9 plots.\n")



################################################################################
############################## 8️⃣ PROTEIN IDS PLOT #############################
################################################################################

library(dplyr)
library(ggplot2)

# Prefer Step 7 matrix (_corrFiltered). If missing, fallback to _filtered_50percent.
get_input_matrix_file <- function(folder) {
f_corr <- list.files(folder, pattern = "_corrFiltered\\.txt$", full.names = TRUE)
if (length(f_corr) > 0) return(f_corr[1])

f50 <- list.files(folder, pattern = "_filtered_50percent\\.txt$", full.names = TRUE)
if (length(f50) > 0) return(f50[1])

return(NA_character_)
}

plot_step8_for_folder <- function(folder) {
if (!dir.exists(folder)) {
cat("⚠️ Folder does not exist:", folder, "\n")
return(invisible(NULL))
}

# Locate matrix
f_in <- get_input_matrix_file(folder)
if (is.na(f_in) || !file.exists(f_in)) {
cat("⚠️ No input matrix found (_corrFiltered or _filtered_50percent) in", folder, "\n")
return(invisible(NULL))
}

# load matrix
dat <- read.table(f_in, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# identify LFQ columns
sample_cols <- grep("^LFQ_", colnames(dat), value = TRUE)
if (length(sample_cols) == 0) {
cat("⚠️ No LFQ_ columns in", folder, "\n")
return(invisible(NULL))
}

# Protein IDs per sample (with cumulative/shared)
protein_matrix <- !is.na(as.matrix(dat[sample_cols]))
protein_counts_per_sample <- colSums(protein_matrix)
sample_names <- sample_cols

cumulative_protein <- sapply(seq_along(sample_names), function(i) {
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) > 0)
})
shared_protein <- sapply(seq_along(sample_names), function(i) {
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) == i)
})

df8 <- data.frame(
Sample       = factor(sample_names, levels = sample_names),
ProteinCount = protein_counts_per_sample,
Cumulative   = cumulative_protein,
Shared       = shared_protein
)

y_max <- max(df8$ProteinCount, df8$Cumulative, df8$Shared, na.rm = TRUE) * 1.10
n_samples  <- length(sample_names)
plot_width <- max(12, n_samples * 0.25)

p8 <- ggplot(df8, aes(x = Sample)) +
geom_bar(aes(y = ProteinCount), stat = "identity", fill = "#0065bd", width = 0.85) +
geom_line(aes(y = Cumulative, group = 1, color = "Cumulative"), linewidth = 1.2) +
geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2.2) +
geom_line(aes(y = Shared, group = 1, color = "Shared"), linewidth = 1.2) +
geom_point(aes(y = Shared, color = "Shared"), size = 2.2) +
scale_color_manual(values = c("Cumulative" = "green4", "Shared" = "orange2")) +
labs(
y = "Number of proteins (ID)",
x = "Samples",
title = "Protein IDs per Sample with Shared and Cumulative Trends"
) +
scale_x_discrete(labels = sample_names) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
theme_classic(base_size = 22) +
theme(
legend.title = element_blank(),
legend.position.inside = c(0.15, 0.15),
axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = 14),
axis.text.y = element_text(size = 22, color = "black", face = "bold"),
axis.title.x = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
axis.title.y = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
plot.title   = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22, t = 35)),
plot.background  = element_rect(fill = "white", color = NA),
panel.background = element_rect(fill = "white", color = NA)
)

# outputs (include source label)
source_tag <- if (grepl("_corrFiltered\\.txt$", f_in)) "corrFiltered" else "filtered_50percent"

ggsave(
file.path(folder, paste0("step8_protein_IDs_per_sample_trends_", source_tag, ".png")),
p8, width = 10, height = 10, dpi = 300, bg = "white"
)
ggsave(
file.path(folder, paste0("step8_protein_IDs_per_sample_trends_", source_tag, ".pdf")),
p8, width = 10, height = 10, dpi = 300, bg = "white", device = cairo_pdf
)

cat("✅ Step 8 completed for:", folder, " | input:", basename(f_in), "\n")
}

folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")
invisible(lapply(folders, plot_step8_for_folder))

################################################################################
############################## 9️⃣ PROTEIN TOTALS PLOT ##########################
################################################################################

library(dplyr)
library(ggplot2)
library(openxlsx)

plot_step9_for_folder <- function(folder) {
  if (!dir.exists(folder)) {
    cat("⚠️ Folder does not exist:", folder, "\n")
    return(invisible(NULL))
  }
  
  # Prefer Step 7 matrix (_corrFiltered). If missing, fallback to _filtered_50percent.
  f_in <- get_input_matrix_file(folder)
  if (is.na(f_in) || !file.exists(f_in)) {
    cat("⚠️ No input matrix found (_corrFiltered or _filtered_50percent) in", folder, "\n")
    return(invisible(NULL))
  }
  
  # load data
  dat <- read.table(f_in, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
  
  # identify LFQ columns
  sample_cols <- grep("^LFQ_", colnames(dat), value = TRUE)
  if (length(sample_cols) == 0) {
    cat("⚠️ No LFQ_ columns in", folder, "\n")
    return(invisible(NULL))
  }
  
  # convert numeric safely
  dat[sample_cols] <- lapply(dat[sample_cols], function(x) suppressWarnings(as.numeric(x)))
  
  # count identified proteins (LFQ > 0)
  protein_counts_mat <- sapply(sample_cols, function(col) sum(dat[[col]] > 0, na.rm = TRUE))
  
  # extract sample and group from column names (e.g., LFQ_SCC1_0Gy_1 → SCC1, 0Gy)
  sample_info <- do.call(rbind, lapply(sample_cols, function(nm) {
    m <- regexec("^LFQ_([^_]+)_([^_]+)_.*$", nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 3) c(Cell = parts[2], Dose = parts[3]) else c(Cell = NA, Dose = NA)
  }))
  sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  rownames(sample_info) <- sample_cols
  
  df_protein <- data.frame(
    Cell = sample_info$Cell,
    Dose = sample_info$Dose,
    ProteinCount = as.numeric(protein_counts_mat),
    stringsAsFactors = FALSE
  )
  
  # --- enforce stable order if present ---
  dose_order <- c("0Gy", "6Gy")
  if (all(dose_order %in% unique(df_protein$Dose))) {
    df_protein$Dose <- factor(df_protein$Dose, levels = dose_order)
  } else {
    df_protein$Dose <- factor(df_protein$Dose, levels = sort(unique(df_protein$Dose)))
  }
  
  # calculate mean ± SD per dose group
  summary_protein <- df_protein %>%
    group_by(Dose) %>%
    summarise(
      n    = sum(!is.na(ProteinCount)),
      mean = mean(ProteinCount, na.rm = TRUE),
      sd   = sd(ProteinCount, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      `Mean±SD` = paste0(round(mean, 1), " ± ", round(sd, 1))
    )
  
  y_max <- max(summary_protein$mean + summary_protein$sd,
               df_protein$ProteinCount, na.rm = TRUE) * 1.10
  
  # plot: bar (mean ± SD) + points for samples
  p9 <- ggplot(summary_protein, aes(x = Dose, y = mean, fill = Dose)) +
    geom_bar(stat = "identity", width = 0.75, color = "black") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  width = 0.3, linewidth = 1.2, color = "black") +
    geom_point(
      data = df_protein, aes(x = Dose, y = ProteinCount),
      position = position_jitter(width = 0.13, height = 0),
      size = 3.5, color = "black", fill = "white", shape = 21, stroke = 1.2, alpha = 1
    ) +
    labs(
      x = NULL,
      y = "Number of identified proteins",
      title = NULL
    ) +
    scale_fill_manual(values = c("0Gy" = "#4575b4", "6Gy" = "#d73027")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
    theme_classic(base_size = 20) +
    theme(
      legend.position = "none",
      axis.line        = element_line(color = "black", linewidth = 1.2),
      axis.ticks       = element_line(color = "black", linewidth = 1),
      axis.text.x      = element_text(size = 32, color = "black", face = "bold"),
      axis.text.y      = element_text(size = 32, color = "black", face = "bold"),
      axis.title.x     = element_text(size = 32, color = "black", face = "bold", margin = margin(t = 18)),
      axis.title.y     = element_text(size = 32, color = "black", face = "bold", margin = margin(r = 18)),
      plot.title       = element_text(size = 32, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22, t = 35)),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  source_tag <- if (grepl("_corrFiltered\\.txt$", f_in)) "corrFiltered" else "filtered_50percent"
  
  # save both PNG + PDF (use pdf device for best portability)
  ggsave(
    file.path(folder, paste0("step9_protein_totals_by_group_opencircle_", source_tag, ".png")),
    p9, width = 10, height = 10, dpi = 300, bg = "white"
  )
  ggsave(
    file.path(folder, paste0("step9_protein_totals_by_group_opencircle_", source_tag, ".pdf")),
    p9, width = 10, height = 10, dpi = 300, bg = "white", device = "pdf"
  )
  
  # ==========================
  # ✅ EXPORT SUMMARY TO EXCEL
  # ==========================
  # one Excel per folder (contains mean, SD, n, and Mean±SD)
  out_xlsx <- file.path(folder, paste0("step9_protein_totals_summary_meanSD_", source_tag, ".xlsx"))
  
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_protein %>% mutate(Dose = as.character(Dose)))
  
  addWorksheet(wb, "Per_Sample")
  writeData(wb, "Per_Sample", df_protein %>% mutate(Dose = as.character(Dose)))
  
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  
  cat("✅ Step 9 completed for:", folder,
      " | input:", basename(f_in),
      " | summary:", basename(out_xlsx), "\n")
  
  invisible(list(
    folder = folder,
    input_file = f_in,
    source_tag = source_tag,
    summary = summary_protein,
    per_sample = df_protein,
    summary_xlsx = out_xlsx
  ))
}

folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

# run + collect results
res_list <- lapply(folders, plot_step9_for_folder)
res_list <- Filter(Negate(is.null), res_list)

# OPTIONAL (bonus): also write ONE combined Excel across all folders
if (length(res_list) > 0) {
  combined_summary <- bind_rows(lapply(res_list, function(x) {
    x$summary %>%
      mutate(
        Folder = x$folder,
        InputFile = basename(x$input_file),
        SourceTag = x$source_tag
      )
  })) %>%
    dplyr::select(Folder, SourceTag, InputFile, dplyr::everything())
  
  combined_per_sample <- bind_rows(lapply(res_list, function(x) {
    x$per_sample %>%
      mutate(
        Folder = x$folder,
        InputFile = basename(x$input_file),
        SourceTag = x$source_tag
      )
  })) %>%
    dplyr::select(Folder, SourceTag, InputFile, dplyr::everything())
  
  out_all <- file.path(getwd(), "step9_ALL_folders_proteinTotals_summary_meanSD.xlsx")
  
  wb_all <- createWorkbook()
  addWorksheet(wb_all, "Summary_AllFolders")
  writeData(wb_all, "Summary_AllFolders", combined_summary)
  
  addWorksheet(wb_all, "PerSample_AllFolders")
  writeData(wb_all, "PerSample_AllFolders", combined_per_sample)
  
  saveWorkbook(wb_all, out_all, overwrite = TRUE)
  
  cat("\n✅ Combined Step 9 summary Excel saved:\n", out_all, "\n")
} else {
  cat("\n⚠️ No folders produced results, so combined Excel not written.\n")
}




################################################################################
############################## 🔟 IMPUTATION (PERSEUS STYLE) ###################
################################################################################

library(openxlsx)

# --- Fixed seed for reproducibility ---
set.seed(1)

# --- Folders to process (same as before) ---
folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

for (folder in folders) {
cat("\n=== Step 10: Imputation for", folder, "===\n")

# --- Find Step 9 output file (_corrFiltered.txt) ---
step9_files <- list.files(
folder,
pattern = "_filtered_50percent_corrFiltered\\.txt$",
full.names = TRUE
)

if (length(step9_files) == 0) {
cat("⚠️ No Step 9 (_corrFiltered.txt) file found in", folder, "— skipping.\n")
next
}

infile <- step9_files[1]
base_name <- sub("\\.txt$", "", basename(infile))
cat("📄 Using:", basename(infile), "\n")

# --- Load the data ---
df <- read.table(infile, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Identify LFQ columns ---
expr_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (!length(expr_cols)) {
cat("⚠️ No LFQ_ columns found in", folder, "— skipping.\n")
next
}

# --- Convert to numeric and clean ---
X <- as.matrix(df[expr_cols])
storage.mode(X) <- "double"
X[!is.finite(X)] <- NA_real_

# --- Compute μ and σ from all valid values (Perseus style) ---
vals <- X[!is.na(X)]
if (length(vals) < 2) {
cat("⚠️ Not enough valid values to impute in", folder, "— skipping.\n")
next
}
mu  <- mean(vals)
sig <- sd(vals)
if (!is.finite(sig) || sig == 0) {
cat("⚠️ σ invalid in", folder, "— skipping.\n")
next
}

# --- Perseus-style parameters ---
mu_imp <- mu - 1.8 * sig
sd_imp <- 0.3 * sig

# --- Impute ---
miss <- which(is.na(X))
if (length(miss)) {
X[miss] <- rnorm(length(miss), mean = mu_imp, sd = sd_imp)
cat("🧩 Imputed", length(miss), "missing values in", folder, "\n")
} else {
cat("✅ No missing values detected in", folder, "\n")
}

# --- Replace data and save ---
imputed_data <- df
imputed_data[expr_cols] <- X

out_txt  <- file.path(folder, paste0(base_name, "_imputed_fixedseed.txt"))
out_xlsx <- file.path(folder, paste0(base_name, "_imputed_fixedseed.xlsx"))

write.table(imputed_data, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(imputed_data, out_xlsx, overwrite = TRUE)

cat("💾 Saved:\n  -", out_txt, "\n  -", out_xlsx, "\n")
}

################################################################################
cat("\n✅ Step 10 complete: Perseus-style imputation done for all Step 9 outputs.\n")
cat("Each folder now contains its imputed files.\n")
################################################################################


################################################################################
############################## 1️⃣1️⃣ QUANTILE NORMALIZATION ####################
################################################################################

library(openxlsx)

# Folders to process
folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

quantile_normalize <- function(X) {
# Remove rows with all NA
all_na <- apply(X, 1, function(x) all(is.na(x)))
X2 <- X[!all_na, , drop = FALSE]

# Rank each column, ignoring NA
ranks <- apply(X2, 2, function(col) rank(col, ties.method = "min", na.last = "keep"))

# Sort each column, keeping NA at end
sorted_X <- apply(X2, 2, function(col) sort(col, na.last = TRUE))

# For each rank, get mean (excluding NAs at that rank)
mean_sorted <- numeric(nrow(X2))
for (i in seq_len(nrow(X2))) {
mean_sorted[i] <- mean(sorted_X[i, ], na.rm = TRUE)
}

# Now assign means to each value according to its rank
X_norm <- X2
for (j in seq_len(ncol(X2))) {
for (i in seq_len(nrow(X2))) {
if (!is.na(ranks[i, j])) {
X_norm[i, j] <- mean_sorted[ranks[i, j]]
} else {
X_norm[i, j] <- NA_real_
}
}
}
# Re-insert all-NA rows
if (any(all_na)) {
res <- matrix(NA, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
res[!all_na, ] <- X_norm
colnames(res) <- colnames(X)
rownames(res) <- rownames(X)
X_norm <- res
}
return(X_norm)
}

for (folder in folders) {
cat("\n=== Step 11: Quantile normalization for", folder, "===\n")

# Find the imputed file from Step 10
step10_file <- list.files(
folder,
pattern = "_imputed_fixedseed\\.txt$",
full.names = TRUE
)
if (length(step10_file) == 0) {
cat("⚠️  No imputed file found in", folder, "— skipping.\n")
next
}

infile <- step10_file[1]
imputed_data <- read.table(infile, header = TRUE, sep = "\t",
quote = "", check.names = FALSE)

# Identify LFQ columns
all_lfq_cols <- grep("^LFQ_", colnames(imputed_data), value = TRUE)
cell_name <- folder
if (cell_name == "Group_comparison") {
expr_cols <- all_lfq_cols
} else {
expr_cols <- grep(paste0("^LFQ_", cell_name, "_"), all_lfq_cols, value = TRUE)
}

cat("🔎 Using", length(expr_cols), "LFQ columns in", folder, "\n")
if (length(expr_cols) == 0) {
cat("❌  No matching LFQ columns — skipping.\n")
next
}

# Coerce LFQ columns to numeric safely
for (cn in expr_cols) {
imputed_data[[cn]] <- suppressWarnings(as.numeric(imputed_data[[cn]]))
}

# Build numeric matrix
X <- as.matrix(imputed_data[expr_cols])
storage.mode(X) <- "double"

if (!all(is.finite(X) | is.na(X))) {
bad <- sum(!is.finite(X) & !is.na(X))
stop(paste("Non-finite values (", bad, ") found in LFQ columns for", folder,
"- check Step 10 output.", sep = " "))
}

# ---- Apply **proper quantile normalization** across samples ----
cat("🔄 Performing quantile normalization on", ncol(X), "LFQ columns...\n")
X_norm <- quantile_normalize(X)

quantile_norm_data <- imputed_data
quantile_norm_data[expr_cols] <- X_norm

# ---- Output ----
base_in <- sub("\\.txt$", "", basename(infile))
base_out <- sub("_filtered_50percent_corrFiltered_imputed_fixedseed$", "", base_in, perl = TRUE)
base_out <- sub("_corrFiltered_imputed_fixedseed$", "", base_out, perl = TRUE)
base_out <- sub("_imputed_fixedseed$", "", base_out, perl = TRUE)

out_txt  <- file.path(folder, paste0(base_out, "_quantile_normalized.txt"))
out_xlsx <- file.path(folder, paste0(base_out, "_quantile_normalized.xlsx"))

write.table(quantile_norm_data, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(quantile_norm_data, out_xlsx, overwrite = TRUE)

cat("💾 Saved:\n  -", out_txt, "\n  -", out_xlsx, "\n")
}

cat("\n✅ Step 11 complete — quantile normalization finished for all datasets.\n")



################################################################################
############################## 1️⃣2️⃣ sPLS-DA ANALYSIS ##########################
################################################################################

################################################################################
########################### sPLS-DA with ellipse ###############################
################################################################################

# ---- Load required packages ----
library(ggplot2)
library(ggrepel)
library(mixOmics)
library(ellipse)
library(dplyr)
library(openxlsx)

# ---- Main function ----
run_splsda_proteomics <- function(
cell_prefix,
ncomp = 2,
keepX = c(100, 100),
colors = c("0Gy" = "#1f78b4", "6Gy" = "#d73027"),
ellipse_level = 0.9,
seed = 123
) {
set.seed(seed)
out_dir <- paste0("sPLSDA_", cell_prefix)
dir.create(out_dir, showWarnings = FALSE)

# ---- Load normalized data ----
norm_pattern <- "_quantile_normalized\\.txt$"
folder <- cell_prefix
norm_file <- list.files(folder, pattern = norm_pattern, full.names = TRUE)[1]
if (is.na(norm_file)) {
warning("⚠️ No normalized file found for: ", cell_prefix)
return(invisible(NULL))
}
data <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# ---- Identify sample columns ----
regex <- if (cell_prefix == "Group_comparison") "^LFQ_.*_(0Gy|6Gy)_" else
paste0("^LFQ_", cell_prefix, "_(0Gy|6Gy)_")
sample_cols <- grep(regex, colnames(data), value = TRUE)
if (length(sample_cols) < 3) {
warning("⚠️ Too few samples for: ", cell_prefix)
return(invisible(NULL))
}

# ---- Prepare data ----
X <- as.data.frame(t(data[, sample_cols, drop = FALSE]))
Y <- factor(ifelse(grepl("_0Gy_", rownames(X)), "0Gy", "6Gy"), levels = c("0Gy", "6Gy"))
cat("\n🧩", cell_prefix, " – samples per group:\n"); print(table(Y))

# ---- Remove zero-variance features ----
X <- X[, apply(X, 2, sd) > 0, drop = FALSE]

# ---- Fit model ----
splsda_model <- suppressWarnings(
splsda(X, Y, ncomp = ncomp, keepX = keepX, scale = TRUE)
)
expl_var <- apply(splsda_model$variates$X^2, 2, sum) / sum(splsda_model$X^2)

# ---- Save keepX plot (y-axis limit = 250) ----
keep_df <- data.frame(Component = paste0("Comp", seq_along(keepX)), keepX = keepX)
keep_plot <- ggplot(keep_df, aes(x = Component, y = keepX)) +
geom_bar(stat = "identity", fill = "#0065bd", width = 0.6) +
geom_text(aes(label = keepX), vjust = -0.5, size = 8) +
labs(title = "sPLS-DA", x = "Component", y = "Variables Selected") +
scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0.05))) +
theme_minimal(base_size = 28) +
theme(
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
axis.title = element_text(size = 28, face = "bold"),
axis.text = element_text(size = 24)
)
ggsave(file.path(out_dir, paste0(cell_prefix, "_keepX_plot.png")),
keep_plot, dpi = 300, width = 7, height = 6, bg = "white")
ggsave(file.path(out_dir, paste0(cell_prefix, "_keepX_plot.pdf")),
keep_plot, dpi = 300, width = 7, height = 6, bg = "white")

# ---- Scores and centroids ----
scores <- data.frame(
Comp1 = splsda_model$variates$X[, 1],
Comp2 = splsda_model$variates$X[, 2],
Group = Y
)
centroids <- scores |>
group_by(Group) |>
summarise(Comp1 = mean(Comp1), Comp2 = mean(Comp2), n = n(), .groups = "drop")

# ---- Build ellipse coordinates ----
build_ellipses <- function(df_by_group, level = 0.9) {
out <- lapply(split(df_by_group, df_by_group$Group), function(df) {
if (nrow(df) < 3) return(NULL)
S <- cov(df[, c("Comp1","Comp2")])
if (any(!is.finite(S))) return(NULL)
if (det(S) <= 1e-12) S <- S + diag(2) * 1e-6
el <- ellipse::ellipse(S,
centre = colMeans(df[, c("Comp1","Comp2")]),
level = level, npoints = 200)
el <- as.data.frame(el)
colnames(el) <- c("x","y")
el$Group <- df$Group[1]
el
})
dplyr::bind_rows(Filter(Negate(is.null), out))
}
ellipse_data <- build_ellipses(scores, level = ellipse_level)

# ---- Plot sPLS-DA ----
xlab <- paste0("Component 1 (", round(expl_var[1] * 100, 1), "%)")
ylab <- paste0("Component 2 (", round(expl_var[2] * 100, 1), "%)")

p <- ggplot(scores, aes(x = Comp1, y = Comp2, color = Group, fill = Group))

if (nrow(ellipse_data) > 0) {
p <- p +
geom_polygon(data = ellipse_data, aes(x = x, y = y, fill = Group),
alpha = 0.25, color = NA) +
geom_path(data = ellipse_data, aes(x = x, y = y, color = Group),
linewidth = 1.2)
} else {
p <- p + stat_ellipse(type = "t", level = ellipse_level, alpha = 0.25, linewidth = 1.2)
}

p <- p +
geom_point(size = 7, alpha = 0.9) +
geom_text_repel(
data = centroids,
aes(label = paste0(Group, "\n(n=", n, ")")),
color = "black", size = 9, fontface = "bold"
) +
scale_color_manual(values = colors) +
scale_fill_manual(values = colors) +
labs(title = NULL, x = xlab, y = ylab, color = "Group", fill = "Group") +
coord_equal() +
theme_minimal(base_size = 30) +
theme(
plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
legend.position = "noon",
legend.title = element_text(size = 28, face = "bold"),
legend.text = element_text(size = 26),
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
axis.title = element_text(size = 30, face = "bold"),
axis.text = element_text(size = 26)
)

# ---- Save sPLS-DA plot ----
ggsave(file.path(out_dir, paste0(cell_prefix, "_sPLSDA_plot.png")),
p, dpi = 300, width = 10, height = 10, bg = "white")
ggsave(file.path(out_dir, paste0(cell_prefix, "_sPLSDA_plot.pdf")),
p, dpi = 300, width = 10, height = 10, bg = "white")

cat("✅ Completed sPLS-DA for:", cell_prefix, "\n")
}

# ---- Run all comparisons ----
comparisons <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")
for (grp in comparisons) run_splsda_proteomics(grp)

cat("\n✅ All groups processed: sPLS-DA analyses and plots completed.\n")
################################################################################



################################################################################
############################## PCA ANALYSIS (Final) ############################
################################################################################

# ---- Load only necessary packages ----
library(mixOmics)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(ellipse)

run_pca_proteomics <- function(
cell_prefix,
colors = c("0Gy" = "#1f78b4", "6Gy" = "#d73027"),
ellipse_level = 0.9,
variance_target = 0.80,
max_ncomp_cap = 10,
seed = 1
) {
set.seed(seed)
out_dir <- paste0("PCA_", cell_prefix)
dir.create(out_dir, showWarnings = FALSE)

# ---- Locate and load normalized file ----
norm_pattern <- "_quantile_normalized\\.txt$"
folder <- cell_prefix
norm_file <- list.files(folder, pattern = norm_pattern, full.names = TRUE)[1]
if (is.na(norm_file)) {
warning("⚠️ No normalized file found for: ", cell_prefix)
return(invisible(NULL))
}
df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
cat("📄 Using input file:", norm_file, "\n")

# ---- Detect LFQ columns ----
regex <- if (cell_prefix == "Group_comparison") "^LFQ_.*_(0Gy|6Gy)_" else
paste0("^LFQ_", cell_prefix, "_(0Gy|6Gy)_")
sample_cols <- grep(regex, colnames(df), value = TRUE)
if (length(sample_cols) < 3) {
warning("⚠️ Too few samples for:", cell_prefix)
return(invisible(NULL))
}

# ---- Prepare data ----
X <- as.data.frame(t(df[, sample_cols, drop = FALSE]))
Y <- factor(ifelse(grepl("_0Gy_", rownames(X)), "0Gy", "6Gy"), levels = c("0Gy", "6Gy"))
cat("\n🧩", cell_prefix, "– samples per group:\n"); print(table(Y))

# ---- Remove zero-variance proteins ----
X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

# ---- Determine number of components ----
total_samples <- nrow(X)
max_possible <- min(total_samples - 1, ncol(X), max_ncomp_cap)

pca_tmp <- suppressWarnings(
pca(X, ncomp = max_possible, center = TRUE, scale = TRUE)
)
prop <- as.numeric(pca_tmp$prop_expl_var$X)
cumprop <- cumsum(prop)
ncomp_auto <- max(2, min(which(cumprop >= variance_target)[1], max_possible))
cat("ℹ️ PCA components chosen:", ncomp_auto,
"(cumulative variance =", round(cumprop[ncomp_auto] * 100, 2), "%)\n")

# ---- Final PCA ----
pca_fit <- suppressWarnings(pca(X, ncomp = ncomp_auto, center = TRUE, scale = TRUE))

# ---- Scree plot ----
scree_df <- data.frame(
PC = factor(paste0("PC", seq_along(prop)), levels = paste0("PC", seq_along(prop))),
Variance = prop * 100,
Cumulative = cumprop * 100
)

gg_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
geom_bar(stat = "identity", width = 0.7, fill = "#0065bd") +
geom_text(aes(label = sprintf("%.1f%%", Variance)), vjust = -0.4, size = 8) +
scale_y_continuous(limits = c(0, 100)) +
labs(title = "PCA", x = "Principal Component", y = "Explained Variance (%)") +
theme_minimal(base_size = 28) +
theme(
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
axis.title = element_text(size = 28, face = "bold"),
axis.text = element_text(size = 24)
)

ggsave(file.path(out_dir, paste0(cell_prefix, "_PCA_scree.png")),
gg_scree, dpi = 300, width = 9, height = 7, bg = "white")
ggsave(file.path(out_dir, paste0(cell_prefix, "_PCA_scree.pdf")),
gg_scree, dpi = 300, width = 9, height = 7, bg = "white")

# ---- PC1 vs PC2 scatter ----
scores <- data.frame(
PC1 = pca_fit$variates$X[, 1],
PC2 = pca_fit$variates$X[, 2],
Group = Y
)
centroids <- scores |>
group_by(Group) |>
summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = dplyr::n(), .groups = "drop")

# ---- Ellipse builder ----
build_ellipses <- function(df_by_group, level = 0.9) {
out <- lapply(split(df_by_group, df_by_group$Group), function(df) {
if (nrow(df) < 3) return(NULL)
S <- cov(df[, c("PC1","PC2")])
if (any(!is.finite(S))) return(NULL)
if (det(S) <= 1e-12) S <- S + diag(2) * 1e-6
el <- ellipse::ellipse(S,
centre = colMeans(df[, c("PC1","PC2")]),
level = level, npoints = 200)
el <- as.data.frame(el)
colnames(el) <- c("x","y")
el$Group <- df$Group[1]
el
})
dplyr::bind_rows(Filter(Negate(is.null), out))
}
ellipse_data <- build_ellipses(scores, level = ellipse_level)

# ---- Plot PCA ----
xlab <- paste0("PC1 (", round(prop[1] * 100, 1), "%)")
ylab <- paste0("PC2 (", round(prop[2] * 100, 1), "%)")

p <- ggplot(scores, aes(x = PC1, y = PC2, color = Group, fill = Group))
if (nrow(ellipse_data) > 0) {
p <- p +
geom_polygon(data = ellipse_data, aes(x = x, y = y, fill = Group),
alpha = 0.25, color = NA) +
geom_path(data = ellipse_data, aes(x = x, y = y, color = Group),
linewidth = 1.2)
} else {
p <- p + stat_ellipse(type = "t", level = ellipse_level, alpha = 0.25, linewidth = 1.2)
}

p <- p +
geom_point(size = 7, alpha = 0.9) +
geom_text_repel(
data = centroids,
aes(label = paste0(Group, "\n(n=", n, ")")),
color = "black", size = 9, fontface = "bold"
) +
scale_color_manual(values = colors) +
scale_fill_manual(values = colors) +
labs(title = "PCA", x = xlab, y = ylab, color = "Group", fill = "Group") +
coord_equal() +
theme_minimal(base_size = 30) +
theme(
plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
legend.position = "right",
legend.title = element_text(size = 28, face = "bold"),
legend.text = element_text(size = 26),
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
axis.title = element_text(size = 30, face = "bold"),
axis.text = element_text(size = 26)
)

# ---- Save PCA plot ----
ggsave(file.path(out_dir, paste0(cell_prefix, "_PCA_PC1_vs_PC2.png")),
p, dpi = 300, width = 10, height = 9, bg = "white")
ggsave(file.path(out_dir, paste0(cell_prefix, "_PCA_PC1_vs_PC2.pdf")),
p, dpi = 300, width = 10, height = 9, bg = "white")

# ---- Save summary ----
write.xlsx(list(
Variance = data.frame(
PC = paste0("PC", seq_along(prop)),
Explained_percent = round(prop * 100, 2),
Cumulative_percent = round(cumprop * 100, 2)
)
), file.path(out_dir, paste0(cell_prefix, "_PCA_summary.xlsx")), overwrite = TRUE)

cat("✅ PCA completed for", cell_prefix, "→", out_dir, "\n")
}

# ---- Run PCA for all groups ----
groups <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")
for (g in groups) run_pca_proteomics(g)

cat("\n✅ Step 13 complete — PCA analysis finished for all normalized datasets.\n")
################################################################################

################################################################################
############################## 1️⃣4️⃣ VOLCANO PLOTS (ROBUST) ###################
################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(grid)          # for unit()
library(Cairo)         # for cairo_pdf device

# ========================== SETTINGS ==========================================
pval_thr  <- 0.05            # p-value threshold (e.g., 0.05)
fc_thr    <- log10(1.5)      # log10 fold-change threshold (e.g., 1.5-fold)
test_type <- "welch"         # "" = auto (F-test); "welch" = Welch; "student" = Student; "wilcoxon" = Wilcoxon
n_label   <- 10              # Max number of gene labels per direction
# ==============================================================================

folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

for (folder in folders) {
cat("\n🚀 Processing:", folder, "\n")

# Find quantile-normalized file in this folder
qn_file <- list.files(folder, pattern = "_quantile_normalized\\.txt$",
full.names = TRUE, recursive = TRUE)
if (length(qn_file) == 0) {
cat("⚠️ No quantile-normalized file found in", folder, "\n")
next
}

quantile_norm_data <- read.table(qn_file[1], header = TRUE, sep = "\t",
   quote = "", check.names = FALSE)
df <- quantile_norm_data

# Detect LFQ columns
sample_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(sample_cols) < 4) {
cat("⚠️ Too few LFQ columns — skipping", folder, "\n")
next
}

# Always assign group columns by dose, not group order!
dose_labels <- sapply(sample_cols, function(x) {
if (grepl("_0Gy_", x)) return("0Gy")
if (grepl("_6Gy_", x)) return("6Gy")
return(NA_character_)
})

g1_cols <- sample_cols[dose_labels == "0Gy"]
g2_cols <- sample_cols[dose_labels == "6Gy"]
g1 <- "0Gy"
g2 <- "6Gy"

if (length(g1_cols) < 2 || length(g2_cols) < 2) {
cat("⚠️ Not enough samples for comparison — skipping.\n")
next
}

cat("🧬 Groups detected:", g1, "vs", g2, "\n")
cat("📊 n =", length(g1_cols), "/", length(g2_cols), "samples\n")

# Prepare result data
gene_names <- if ("PG.Genes" %in% colnames(df)) df$`PG.Genes` else df$`PG.ProteinGroups`
results <- data.frame(
Gene = gene_names,
log10FC = NA_real_,
p_value = NA_real_,
mean_g1 = NA_real_,
mean_g2 = NA_real_,
stringsAsFactors = FALSE
)

# Statistical test per protein
for (i in seq_len(nrow(df))) {
v1 <- as.numeric(df[i, g1_cols])
v2 <- as.numeric(df[i, g2_cols])

if (sum(!is.na(v1)) >= 2 && sum(!is.na(v2)) >= 2) {
mytype <- tolower(test_type)

test <- tryCatch({
if (mytype == "wilcoxon") {
wilcox.test(v1, v2)
} else if (mytype == "student") {
t.test(v1, v2, paired = FALSE, var.equal = TRUE)
} else if (mytype == "welch") {
t.test(v1, v2, paired = FALSE, var.equal = FALSE)
} else if (mytype == "") {
# Auto: F-test for equal variance
ft <- tryCatch(var.test(v1, v2), error = function(e) NULL)
vequal <- if (!is.null(ft) && !is.na(ft$p.value) && ft$p.value > 0.05) TRUE else FALSE
t.test(v1, v2, paired = FALSE, var.equal = vequal)
} else {
stop("Invalid test type.")
}
}, error = function(e) NULL)

if (!is.null(test)) {
results$log10FC[i] <- mean(v2, na.rm = TRUE) - mean(v1, na.rm = TRUE) # 6Gy - 0Gy
results$p_value[i] <- test$p.value
results$mean_g1[i] <- mean(v1, na.rm = TRUE)
results$mean_g2[i] <- mean(v2, na.rm = TRUE)
}
}
}

results <- na.omit(results)
results$negLog10P <- -log10(results$p_value)

# Significance classification
results$group <- "Non-significant"
results$group[results$log10FC >  fc_thr & results$p_value < pval_thr] <- "Higher in 6Gy"
results$group[results$log10FC < -fc_thr & results$p_value < pval_thr] <- "Higher in 0Gy"

# Make group a factor so legend/order are stable
results$group <- factor(results$group, levels = c("Higher in 6Gy", "Higher in 0Gy", "Non-significant"))

# Only label significant top genes
top_up <- results %>%
filter(group == "Higher in 6Gy" & negLog10P > -log10(pval_thr) & abs(log10FC) > fc_thr) %>%
arrange(p_value) %>% head(n_label)

top_down <- results %>%
filter(group == "Higher in 0Gy" & negLog10P > -log10(pval_thr) & abs(log10FC) > fc_thr) %>%
arrange(p_value) %>% head(n_label)

# IMPORTANT FIX: explicitly map Non-significant to true black
color_map <- c(
"Higher in 6Gy"      = "#D73027",
"Higher in 0Gy"      = "#0065BD",
"Non-significant"    = "#000000"
)

base_name  <- sub("_quantile_normalized\\.txt$", "", basename(qn_file[1]))
plot_title <- paste0("Volcano Plot: ", folder)

volcano_plot <- ggplot(results, aes(x = log10FC, y = negLog10P)) +
geom_point(aes(color = group), size = 1, alpha = 1) +   # alpha=1 ensures perfect black
scale_color_manual(
values = color_map,
breaks = c("Higher in 6Gy", "Higher in 0Gy", "Non-significant")
) +
geom_hline(yintercept = -log10(pval_thr), linetype = "dashed", color = "black", linewidth = 1.2) +
geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", color = "black", linewidth = 1.2) +
geom_text_repel(data = top_up, aes(label = Gene), color = "#D73027", size = 10, max.overlaps = 20) +
geom_text_repel(data = top_down, aes(label = Gene), color = "#0065BD", size = 10, max.overlaps = 20) +
labs(
title = NULL,
x = expression(Log[10]~Fold~Change~"(6Gy - 0Gy)"),
y = expression(-Log[10]~p~value),
color = "Significance"
) +
theme_minimal(base_size = 28) +
theme(
plot.title   = element_text(size = 36, face = "bold", hjust = 0.5),
axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 20)),
axis.title.y = element_text(size = 32, face = "bold", margin = margin(r = 20)),
axis.text.x  = element_text(size = 24, face = "bold"),
axis.text.y  = element_text(size = 24, face = "bold"),
legend.title = element_text(size = 28, face = "bold"),
legend.text  = element_text(size = 24, face = "bold"),
legend.key.size = unit(2, "lines"),
legend.position = "bottom",
plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
)

# Save outputs
out_png  <- file.path(folder, paste0("Volcano_", base_name, ".png"))
out_pdf  <- file.path(folder, paste0("Volcano_", base_name, ".pdf"))
out_txt  <- file.path(folder, paste0("Volcano_data_", base_name, ".txt"))
out_xlsx <- file.path(folder, paste0("Volcano_data_", base_name, ".xlsx"))

ggsave(out_png, volcano_plot, dpi = 300, width = 14, height = 10, bg = "white")
ggsave(out_pdf, volcano_plot, dpi = 300, width = 14, height = 10, bg = "white", device = cairo_pdf)

write.table(results, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(results, out_xlsx, overwrite = TRUE)

cat("✅ Volcano PNG/PDF and tables saved in:", folder, "\n")
}

################################################################################
cat("\n✅ Step 14 complete — volcano plots and statistics saved in their folders.\n")
################################################################################


################################################################################
############################## 1️⃣5️⃣ ANNOTATION ################################
################################################################################
# Annotates significant proteins (Higher in 6 Gy / Higher in 0 Gy)
# Keeps only significant proteins with LFQ values + GO, KEGG, Reactome annotations
# All functions are fully namespaced (package::function) to avoid conflicts.
################################################################################

# ============================ LIBRARIES =======================================
library(dplyr)
library(tidyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(ReactomePA)
library(clusterProfiler)
library(biomaRt)
library(openxlsx)

# ============================ FOLDERS =========================================
folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")

# ============================ HELPERS =========================================
split_genes <- function(gene_list) {
df <- data.frame(Original = gene_list, stringsAsFactors = FALSE)
df <- df %>%
dplyr::mutate(Gene_Single = strsplit(Original, ";")) %>%
tidyr::unnest(Gene_Single) %>%
dplyr::mutate(Gene_Single = base::trimws(Gene_Single))
df
}

map_to_entrez <- function(symbols) {
mapping <- AnnotationDbi::select(org.Hs.eg.db,
keys = symbols,
columns = "ENTREZID",
keytype = "SYMBOL")
mapping <- mapping[!duplicated(mapping$SYMBOL), ]

# Fallback via biomaRt
if (any(is.na(mapping$ENTREZID))) {
missing_syms <- mapping$SYMBOL[is.na(mapping$ENTREZID)]
mart <- tryCatch({
biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
}, error = function(e) NULL)

if (!is.null(mart)) {
bm <- tryCatch({
biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
filters = "hgnc_symbol",
values = missing_syms,
mart = mart)
}, error = function(e) NULL)
if (!is.null(bm) && nrow(bm) > 0) {
bm <- bm[!is.na(bm$entrezgene_id) & bm$entrezgene_id != "", ]
if (nrow(bm) > 0) {
cat("🧭 Added", nrow(bm), "Entrez IDs via biomaRt fallback.\n")
for (i in seq_len(nrow(bm))) {
mapping$ENTREZID[mapping$SYMBOL == bm$hgnc_symbol[i]] <-
as.character(bm$entrezgene_id[i])
}
}
}
}
}

mapping <- mapping[!is.na(mapping$ENTREZID) & mapping$ENTREZID != "", ]
mapping
}

# ============================ MAIN LOOP =======================================
for (folder in folders) {
cat("\n==========================================================\n")
cat("Processing:", folder, "\n")
cat("==========================================================\n")

volcano_files <- list.files(folder, pattern = "^Volcano_data_.*\\.txt$", full.names = TRUE)
qn_file       <- list.files(folder, pattern = "_quantile_normalized\\.txt$", full.names = TRUE)

if (!length(volcano_files) || !length(qn_file)) {
cat("⚠️ Missing volcano or normalized data in", folder, "\n")
next
}

cat("✅ Volcano file found:", basename(volcano_files[1]), "\n")
cat("✅ Quantile-normalized file found:", basename(qn_file[1]), "\n")

volcano_df <- utils::read.table(volcano_files[1], header = TRUE, sep = "\t", quote = "",
check.names = FALSE)
norm_df    <- utils::read.table(qn_file[1], header = TRUE, sep = "\t", quote = "",
check.names = FALSE)

up_genes   <- stats::na.omit(unique(volcano_df$Gene[volcano_df$group == "Higher in 6Gy"]))
down_genes <- stats::na.omit(unique(volcano_df$Gene[volcano_df$group == "Higher in 0Gy"]))

cat("🔹 Significant 6Gy-high genes:", length(up_genes), "\n")
cat("🔹 Significant 0Gy-high genes:", length(down_genes), "\n")

if (length(up_genes) == 0 && length(down_genes) == 0) {
cat("⚠️ No significant genes in", folder, "\n")
next
}

# -------------------- ANNOTATION FUNCTION -----------------------------------
annotate_genes <- function(gene_vec) {
if (length(gene_vec) == 0) return(NULL)

gene_map <- split_genes(gene_vec)
gene_list_unique <- unique(stats::na.omit(base::trimws(gene_map$Gene_Single)))

gene_entrez <- map_to_entrez(gene_list_unique)
combined_entrez <- unique(gene_entrez$ENTREZID)
if (length(combined_entrez) < 2) {
cat("⚠️ Too few Entrez IDs — skipping enrichment.\n")
return(NULL)
}

# ---------- GO ----------
go_results <- tryCatch({
go_anno <- AnnotationDbi::select(org.Hs.eg.db, keys = combined_entrez,
columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
go_terms <- AnnotationDbi::select(GO.db, keys = unique(go_anno$GO),
columns = "TERM", keytype = "GOID")
go_merged <- base::merge(go_anno, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
go_merged <- dplyr::left_join(go_merged, gene_entrez, by = "ENTREZID",
relationship = "many-to-many")
go_merged$GO_Pathway <- paste(go_merged$GO, go_merged$ONTOLOGY, go_merged$TERM, sep = " | ")
dplyr::select(go_merged, SYMBOL, GO_Pathway)
}, error = function(e) NULL)

go_agg <- if (!is.null(go_results) && nrow(go_results) > 0) {
go_results %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(GO_Pathway = paste(unique(GO_Pathway), collapse = "; "),
.groups = "drop")
} else data.frame(SYMBOL = character(), GO_Pathway = character())

# ---------- KEGG ----------
kegg_agg <- tryCatch({
# Enrichment via KEGG
kegg_enrich <- clusterProfiler::enrichKEGG(
gene = combined_entrez,
organism = "hsa",
keyType = "kegg",
pvalueCutoff = 1,
qvalueCutoff = 1
)

if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
kegg_df <- as.data.frame(kegg_enrich)
kegg_agg <- kegg_df %>%
dplyr::select(KEGG_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep = "/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(
KEGG_Pathway = paste(unique(KEGG_Pathway), collapse = "; "),
.groups = "drop"
)
} else {
# Fallback: direct mapping if enrichment returns empty
kegg_link <- tryCatch(
clusterProfiler::bitr_kegg(combined_entrez, fromType = "kegg", toType = "path", organism = "hsa"),
error = function(e) NULL
)
if (!is.null(kegg_link) && nrow(kegg_link) > 0) {
kegg_desc <- tryCatch({
desc <- clusterProfiler::keggList("pathway", "hsa")
data.frame(KEGG_Pathway = sub("path:", "", names(desc)),
KEGG_Description = unname(desc))
}, error = function(e) NULL)

kegg_results <- kegg_link %>%
dplyr::rename(ENTREZID = from, KEGG_Pathway = to) %>%
dplyr::mutate(KEGG_Pathway = sub("path:", "", KEGG_Pathway)) %>%
dplyr::left_join(kegg_desc, by = "KEGG_Pathway") %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::mutate(KEGG_Annotation = paste(KEGG_Description, KEGG_Pathway, sep = " | ")) %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(
KEGG_Pathway = paste(unique(KEGG_Annotation), collapse = "; "),
.groups = "drop"
)
} else {
data.frame(SYMBOL = character(), KEGG_Pathway = character())
}
}
}, error = function(e) {
cat("⚠️ KEGG enrichment skipped (", length(combined_entrez), " genes)\n")
data.frame(SYMBOL = character(), KEGG_Pathway = character())
})

# ---------- REACTOME ----------
reactome_agg <- tryCatch({
enr <- ReactomePA::enrichPathway(gene = combined_entrez,
organism = "human",
readable = FALSE,
pvalueCutoff = 1,
qvalueCutoff = 1)
if (is.null(enr) || nrow(as.data.frame(enr)) == 0)
return(data.frame(SYMBOL = character(), Reactome_Pathway = character()))
df <- as.data.frame(enr)
df %>%
dplyr::select(Reactome_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep = "/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(Reactome_Pathway =
paste(unique(Reactome_Pathway), collapse = "; "),
.groups = "drop")
}, error = function(e) {
cat("⚠️ Reactome skipped (", length(combined_entrez), " genes)\n")
data.frame(SYMBOL = character(), Reactome_Pathway = character())
})

# ---------- MERGE ----------
combined_annots <- gene_map %>%
dplyr::left_join(go_agg, by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(reactome_agg, by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(kegg_agg, by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::group_by(Original) %>%
dplyr::summarise(
GO_Pathway       = paste(stats::na.omit(unique(GO_Pathway)), collapse = "; "),
Reactome_Pathway = paste(stats::na.omit(unique(Reactome_Pathway)), collapse = "; "),
KEGG_Pathway     = paste(stats::na.omit(unique(KEGG_Pathway)), collapse = "; "),
.groups = "drop"
)
combined_annots
}

# -------------------- SAVE SIGNIFICANT ONLY ---------------------------------
save_sig_annot <- function(annot_df, gene_vec, label) {
if (is.null(annot_df) || nrow(annot_df) == 0) return(NULL)
sig_norm <- norm_df[norm_df$PG.Genes %in% gene_vec, ]
annotated <- sig_norm %>% dplyr::left_join(annot_df, by = c("PG.Genes" = "Original"))
outfile_txt  <- file.path(folder, paste0(folder, "_", label, "_annotated.txt"))
outfile_xlsx <- file.path(folder, paste0(folder, "_", label, "_annotated.xlsx"))
utils::write.table(annotated, outfile_txt, sep = "\t", row.names = FALSE, quote = FALSE)
openxlsx::write.xlsx(annotated, outfile_xlsx, overwrite = TRUE)
cat("💾 Saved:", basename(outfile_txt), "\n")
}

up_annot   <- annotate_genes(up_genes)
down_annot <- annotate_genes(down_genes)

save_sig_annot(up_annot, up_genes, "6Gy_high")
save_sig_annot(down_annot, down_genes, "0Gy_high")

cat("✅ Annotation done for", folder, "\n")
}

cat("\n🎉 Step 15 complete — annotation files saved (GO/KEGG/Reactome ready).\n")
################################################################################


################################################################################
############################## 1️⃣6️⃣ FISHER TEST ################################
################################################################################
# Fisher Exact Test Enrichment for Step 15 annotated results
# - Per folder background from its quantile-normalized file (PG.Genes/PG.ProteinGroups)
# - Runs on <folder>_0Gy_high_annotated.txt and <folder>_6Gy_high_annotated.txt
# - Saves individual Fisher outputs (GO / KEGG / Reactome) as before
# - Additionally saves ONE merged file per condition that combines GO+KEGG+Reactome:
#     Merged_Fisher_<folder>_0Gy.xlsx
#     Merged_Fisher_<folder>_6Gy.xlsx
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# ========================== SETTINGS ==========================================
folders <- c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison")
# ==============================================================================

# --------------------- Fisher helper (one annotation column) ------------------
fisher_on_annotation <- function(ann_df, annot_col, background_genes, out_xlsx) {
if (!(annot_col %in% colnames(ann_df))) return(invisible(NULL))
if (!("PG.Genes" %in% colnames(ann_df))) {
cat("❌ Missing PG.Genes column; skipping.\n")
return(invisible(NULL))
}

# Expand multi-annotation cells and clean
df <- ann_df %>%
filter(!is.na(.data[[annot_col]]) & .data[[annot_col]] != "") %>%
transmute(PG.Genes, Pathway = .data[[annot_col]]) %>%
mutate(Pathway = strsplit(as.character(Pathway), ";")) %>%
unnest(Pathway) %>%
mutate(Pathway = str_trim(Pathway)) %>%
filter(!is.na(Pathway) & Pathway != "" & Pathway != "NA")
df <- df[!grepl("^NA\\s*\\|\\s*NA\\s*\\|\\s*NA$", df$Pathway, ignore.case = TRUE), ]
df$Pathway_Name <- ifelse(grepl(" - Homo sapiens", df$Pathway),
sub(" - Homo sapiens.*", "", df$Pathway),
df$Pathway)

# Significant set = all PG.Genes in the annotated file
sig_genes <- unique(trimws(as.character(ann_df$PG.Genes)))
sig_genes <- sig_genes[sig_genes != "" & !is.na(sig_genes)]

terms <- unique(df$Pathway_Name)
if (length(terms) == 0 || length(sig_genes) == 0 || length(background_genes) == 0) {
return(invisible(NULL))
}

out_list <- vector("list", length(terms))
idx <- 0

for (term in terms) {
term_genes <- unique(df$PG.Genes[df$Pathway_Name == term & df$PG.Genes %in% background_genes])

a <- length(intersect(term_genes, sig_genes))                         # sig in term
b <- length(setdiff(sig_genes, term_genes))                           # sig not in term
c <- length(setdiff(term_genes, intersect(term_genes, sig_genes)))    # non-sig in term
d <- length(setdiff(background_genes, union(sig_genes, term_genes)))  # non-sig not in term
if (any(c(a,b,c,d) < 0)) next

ft <- tryCatch(fisher.test(matrix(c(a,b,c,d), nrow = 2)), error = function(e) NULL)
if (is.null(ft)) next

gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
bg_ratio   <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
ef         <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)

idx <- idx + 1
out_list[[idx]] <- data.frame(
Pathway               = term,
Sig_protein_volcano   = a,
Sig_NotIn_Pathway     = b,
Nonsig_In_Pathway     = c,
Nonsig_NotIn_Pathway  = d,
Genes                 = paste(intersect(term_genes, sig_genes), collapse = ";"),
p_value               = ft$p.value,
enrichment_factor     = ef,
stringsAsFactors = FALSE
)
}

res <- dplyr::bind_rows(out_list)
if (is.null(res) || nrow(res) == 0) return(invisible(NULL))

res$p_adj <- p.adjust(res$p_value, method = "BH")
res <- res %>% arrange(p_adj)

# Save individual file
write.xlsx(res, out_xlsx, overwrite = TRUE)
cat("💾 Saved:", basename(out_xlsx), "(", nrow(res), "terms)\n")

invisible(res)
}

# --------------------------- MAIN LOOP ---------------------------------------
for (folder in folders) {
cat("\n==================================================================\n")
cat("📂 Processing:", folder, "\n")
cat("==================================================================\n")

# Per-folder background from quantile-normalized file
qn_file <- list.files(folder, pattern = "_quantile_normalized\\.txt$", full.names = TRUE)
if (!length(qn_file)) {
cat("⚠️ No quantile-normalized file → skipping.\n")
next
}
norm_df <- read.table(qn_file[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE)
bg_vec <- if ("PG.Genes" %in% colnames(norm_df)) norm_df$`PG.Genes` else norm_df$`PG.ProteinGroups`
background_genes <- unique(trimws(as.character(bg_vec)))
background_genes <- background_genes[background_genes != "" & !is.na(background_genes)]
cat("✅ Background size:", length(background_genes), "\n")

# Step 15 outputs
annot_files <- list.files(
folder,
pattern = paste0("^", folder, "_(6Gy|0Gy)_high_annotated\\.txt$"),
full.names = TRUE
)
if (!length(annot_files)) {
cat("⚠️ No annotated files found in", folder, "→ skipping.\n")
next
}

out_dir <- file.path(folder, "Fisher")
dir.create(out_dir, showWarnings = FALSE)

# Containers to build merged-per-condition (GO+KEGG+Reactome)
merged_condition <- list(
`0Gy` = list(),
`6Gy` = list()
)

for (annot_file in annot_files) {
base <- tools::file_path_sans_ext(basename(annot_file))
condition <- if (grepl("0Gy", annot_file)) "0Gy" else "6Gy"
cat("🔬 File:", basename(annot_file), "\n")

ann_df <- read.table(annot_file, header = TRUE, sep = "\t", quote = "",
stringsAsFactors = FALSE, check.names = FALSE)

# Individual outputs (kept)
res_go <- res_kegg <- res_rea <- NULL

if ("GO_Pathway" %in% colnames(ann_df)) {
res_go <- fisher_on_annotation(
ann_df, "GO_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_GO_", base, ".xlsx"))
)
if (!is.null(res_go)) res_go <- res_go %>% mutate(Pathway_Type = "GO", Source_File = base)
}
if ("KEGG_Pathway" %in% colnames(ann_df)) {
res_kegg <- fisher_on_annotation(
ann_df, "KEGG_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_KEGG_", base, ".xlsx"))
)
if (!is.null(res_kegg)) res_kegg <- res_kegg %>% mutate(Pathway_Type = "KEGG", Source_File = base)
}
if ("Reactome_Pathway" %in% colnames(ann_df)) {
res_rea <- fisher_on_annotation(
ann_df, "Reactome_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_Reactome_", base, ".xlsx"))
)
if (!is.null(res_rea)) res_rea <- res_rea %>% mutate(Pathway_Type = "Reactome", Source_File = base)
}

# Append any available results into the merged-per-condition list
merged_piece <- dplyr::bind_rows(res_go, res_kegg, res_rea)
if (!is.null(merged_piece) && nrow(merged_piece) > 0) {
merged_condition[[condition]][[base]] <- merged_piece
}
}

# ---- Write ONE merged file per condition (GO+KEGG+Reactome together) ----
for (cond in c("0Gy", "6Gy")) {
merged_df <- dplyr::bind_rows(merged_condition[[cond]], .id = "Annotated_Set")
if (!is.null(merged_df) && nrow(merged_df) > 0) {
# Order by adjusted p-value, then by Pathway_Type for readability
merged_df <- merged_df %>%
arrange(p_adj, Pathway_Type, Pathway)

out_merged <- file.path(out_dir, paste0("Merged_Fisher_", folder, "_", cond, ".xlsx"))
openxlsx::write.xlsx(merged_df, out_merged, overwrite = TRUE)
cat("📘 Merged (GO+KEGG+Reactome) for", folder, cond, "→", basename(out_merged),
" (", nrow(merged_df), "rows)\n")
}
}

cat("✅ Fisher enrichment done for", folder, "\n")
}

cat("\n🎉 Step 16 complete — individual (GO/KEGG/Reactome) + merged-per-condition saved.\n")
################################################################################


################################################################################
###################### 1️⃣7️⃣ VOLCANO ENRICHMENT HEATMAPS #######################
################################################################################
# Fully automatic width & height scaling:
# - Each gene & pathway has at least 30 px visual space (~0.76 inches)
# - Dynamic figure width and height for full label visibility
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)
library(scales)
})

# ---- USER SETTINGS ----
volcano_root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics"
P_ADJ_CUTOFF <- 0.01
MIN_PROTEINS <- 2
MAX_PATHWAYS_PER_PAGE <- 15
PATHWAY_WRAP_WIDTH <- 40
BASE_FONTSIZE <- 14
LEGEND_TITLE_SIZE <- 28
LEGEND_TEXT_SIZE <- 22
LEGEND_KEY_HEIGHT_PT <- 60
LEGEND_KEY_WIDTH_PT <- 24

# ---- KEYWORDS ----
global_keywords <- c(
# --- Cell movement ---
"migration","cell migration","chemotaxis","cell motility","invasion",
"metastasis","locomotion","cell movement","wound healing",

# --- Cell death / survival ---
"survival","cell survival","apoptosis","cell death","autophagy",
"cell cycle","mitosis","mitotic","senescence","cell differentiation",
"DNA replication",

# --- ECM / Matrix ---
"matrix","extracellular matrix","ECM","adhesion","cell adhesion",
"focal adhesion","angiogenesis","vasculature","vascular development",
"collagen","fibronectin","basement membrane","extracellular structure",
"integrin",

# --- Cancer / HPV ---
"cancer","carcinoma","tumor","tumour","neoplasm","oncogene",
"tumor suppressor","p53","radiation","HPV","human papillomavirus",
"malignancy","metastatic","tumorigenesis","carcinogenesis",

# --- Stress / DNA damage ---
"stress","cellular stress","genotoxic","oxidative stress",
"redox","oxidation","reduction","ROS","chaperone","heat shock",
"protein folding","UPR",

# --- Mitochondrial / Metabolic ---
"mitochondria","mitochondrial","metabolic","metabolism","lipid",
"fatty acid","glycolysis","gluconeogenesis","ATP",
"oxidative phosphorylation","TCA cycle",

# --- Vesicle / Exosome ---
"exosome","extracellular vesicle","vesicle","multivesicular body",
"endosome","endocytosis","exocytosis","secretion","microvesicle",
"vesicular transport","platelet",

# --- Platelet ---
"platelet","megakaryocyte","thrombocyte","thrombopoiesis",
"platelet activation","platelet aggregation","platelet degranulation",
"platelet derived","platelet factor","thrombus","thrombosis",
"coagulation","clot",

# --- Signaling ---
"MAPK","ERK","PI3K","AKT","mTOR","JAK","STAT","WNT","Notch","TGF"
)

# ---- HELPERS ----
first_col <- function(df, cand) {
c <- cand[cand %in% names(df)]
if (length(c)) c[1] else NA_character_
}
count_genes <- function(x) {
if (is.null(x) || is.na(x)) return(0L)
toks <- unlist(stringr::str_split(as.character(x), "[;,/\\s]+"))
toks <- stringr::str_trim(toks)
toks <- toks[toks != ""]
length(unique(toks))
}
clean_pathway <- function(pw) {
pw <- stringr::str_remove(pw, "^GO:\\d+\\s*\\|\\s*")
pw <- stringr::str_replace(pw, "(^|\\|\\s*)(BP|CC|MF)\\s*\\|\\s*", "")
pw <- stringr::str_replace(pw, "^KEGG\\s*\\|\\s*", "")
pw <- stringr::str_replace(pw, "^Reactome\\s*\\|\\s*", "")
stringr::str_trim(pw)
}
add_source_tag <- function(name, src) {
src_tag <- ifelse(grepl("GO", src, ignore.case = TRUE), "[GO]",
ifelse(grepl("KEGG", src, ignore.case = TRUE), "[KEGG]",
ifelse(grepl("Reactome", src, ignore.case = TRUE), "[Reactome]", "")))
paste0(name, " ", src_tag)
}

# ---- FIND ALL FISHER DIRECTORIES ----
all_dirs <- list.dirs(volcano_root, recursive = TRUE, full.names = TRUE)
fisher_dirs <- all_dirs[grepl("[/\\\\]Fisher$", all_dirs)]
if (!length(fisher_dirs)) stop("❌ No 'Fisher' folders found under: ", normalizePath(volcano_root))

# ---- MAIN LOOP ----
for (fdir in fisher_dirs) {
cat("\n📂 Processing:", fdir, "\n")
files <- list.files(fdir, pattern = "^Fisher_(GO|KEGG|Reactome)_.+\\.xlsx$", full.names = TRUE)
if (!length(files)) {
cat("⚠️ No Fisher enrichment files in", fdir, "\n")
next
}

info <- stringr::str_match(basename(files), "^Fisher_(GO|KEGG|Reactome)_(.*)\\.xlsx$")
colnames(info) <- c("full", "Source", "Direction")
directions <- unique(info[, "Direction"])

for (direction in directions) {
cat("   🔎 Direction:", direction, "\n")
f_direction <- files[info[, "Direction"] == direction]
src <- info[info[, "Direction"] == direction, "Source"]
pieces <- list()

for (k in seq_along(f_direction)) {
fpath <- f_direction[k]
source_type <- src[k]
df0 <- tryCatch(openxlsx::read.xlsx(fpath), error = function(e) NULL)
if (is.null(df0) || !"Pathway" %in% names(df0)) next

pcol <- first_col(df0, c("p_adj", "p.adj", "padj", "adj_p_value", "adj.P.Val", "p_adjust"))
gcol <- first_col(df0, c("Genes", "Gene", "Gene_List", "GeneID", "Gene_Names", "Proteins", "Protein_List"))
if (is.na(pcol) || is.na(gcol)) next

df <- df0 %>%
dplyr::mutate(
.p = suppressWarnings(as.numeric(.data[[pcol]])),
.genes_raw = as.character(.data[[gcol]]),
.gene_n = vapply(.genes_raw, count_genes, integer(1))
) %>%
dplyr::filter(is.finite(.p), .p < P_ADJ_CUTOFF, .gene_n >= MIN_PROTEINS) %>%
dplyr::mutate(
Pathway = stringr::str_trim(Pathway),
Source = source_type,
log10_p = -log10(.p)
) %>%
dplyr::filter(
stringr::str_detect(
stringr::str_to_lower(Pathway),
paste(stringr::str_to_lower(global_keywords), collapse = "|")
)
) %>%
tidyr::separate_rows(dplyr::all_of(gcol), sep = "[;,/\\s]+") %>%
dplyr::mutate(Gene = stringr::str_trim(.data[[gcol]])) %>%
dplyr::filter(Gene != "") %>%
dplyr::mutate(
Clean_Pathway = clean_pathway(Pathway),
Pathway_Display = add_source_tag(
stringr::str_wrap(Clean_Pathway, width = PATHWAY_WRAP_WIDTH),
Source
)
)
pieces[[length(pieces) + 1]] <- df
}

side_df <- dplyr::bind_rows(pieces)
if (!nrow(side_df)) next

gene_order_by_freq <- side_df %>% dplyr::count(Gene, sort = TRUE) %>% dplyr::pull(Gene)
keep_genes <- unique(gene_order_by_freq)
side_df <- side_df %>% dplyr::filter(Gene %in% keep_genes)

pathway_n_proteins <- side_df %>%
dplyr::group_by(Pathway_Display) %>%
dplyr::summarise(n_protein = dplyr::n_distinct(Gene), .groups = "drop")

mat <- side_df %>%
dplyr::select(Pathway_Display, Gene, log10_p) %>%
tidyr::complete(Pathway_Display, Gene = keep_genes)

out_dir <- file.path(fdir, "Combined_Enrichment_Plots", direction)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pths_all <- unique(mat$Pathway_Display)
chunks <- split(pths_all, ceiling(seq_along(pths_all) / MAX_PATHWAYS_PER_PAGE))

for (i in seq_along(chunks)) {
sel <- chunks[[i]]
page_pathway_order <- pathway_n_proteins %>%
dplyr::filter(Pathway_Display %in% sel) %>%
dplyr::arrange(desc(n_protein)) %>%
dplyr::pull(Pathway_Display)

d <- mat %>%
dplyr::filter(Pathway_Display %in% sel) %>%
dplyr::mutate(Pathway_Display = factor(Pathway_Display, levels = rev(page_pathway_order)))

# --- AUTO WIDTH + HEIGHT WITH MINIMUM VISUAL SIZE ---
n_genes <- length(unique(d$Gene))
n_paths <- length(unique(d$Pathway_Display))

# Each gene/pathway gets at least 0.76 inches (~30 px)
W_auto <- max(n_genes * 0.8, 10, n_genes * 0.76)
H_auto <- max(n_paths * 0.45, 10, n_paths * 0.76)

# Adaptive fonts
font_x <- max(10, 28 - 0.03 * n_genes)
font_y <- max(12, min(36, 18 + 0.25 * n_paths))

cat("📊 Genes:", n_genes, " Pathways:", n_paths,
" => Width:", round(W_auto, 1), " Height:", round(H_auto, 1), "\n")

title_txt <- paste0(direction, " (", basename(dirname(fdir)), ") — Page ", i)
caption_txt <- "Filters: adj p < 0.01; proteins ≥ 2; all significant genes shown (super-wide mode)."

p <- ggplot(d, aes(x = Gene, y = Pathway_Display, fill = log10_p)) +
geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE, width = 0.95, height = 0.5) +
scale_fill_gradientn(colours = c("green", "yellow", "red"),
na.value = "black", name = "-log10(p.adj)") +
scale_x_discrete(position = "top", guide = guide_axis(check.overlap = TRUE)) +
labs(title = title_txt, x = NULL, y = "Enriched Pathway", caption = caption_txt) +
theme_minimal(base_size = BASE_FONTSIZE) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size = font_x, angle = 90, hjust = 1, vjust = 0.5),
axis.text.y = element_text(size = font_y),
axis.title.y = element_text(size = font_y + 2, face = "bold"),
plot.title = element_text(size = font_y + 6, face = "bold", hjust = 0.5),
legend.position = "right",
legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
legend.text = element_text(size = LEGEND_TEXT_SIZE),
legend.key.height = unit(LEGEND_KEY_HEIGHT_PT, "pt"),
legend.key.width = unit(LEGEND_KEY_WIDTH_PT, "pt"),
plot.caption = element_text(size = 10, hjust = 1)
) +
coord_cartesian(clip = "off", expand = FALSE)

out_base <- file.path(out_dir, paste0("Enrichment_", direction, "_",
basename(dirname(fdir)), "_Page", i))
ggsave(paste0(out_base, ".pdf"), p, width = W_auto, height = H_auto,
units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(out_base, ".png"), p, width = W_auto, height = H_auto,
units = "in", dpi = 600, bg = "white", limitsize = FALSE)

cat("   ✅ Saved:", out_base, "\n")
}
}
}

cat("\n🎉 All enrichment heatmaps created (super-wide + min-size mode complete).\n")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# ==============================================================================
# 📌 Fisher enrichment + Clustered Sankey (ggalluvial version, multi-color)
# ==============================================================================

# ------------------------------------------------------------------------------
# 📦 Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(ggalluvial)
library(ggplot2)
library(stringr)
library(dplyr)

# ------------------------------------------------------------------------------
# 1️⃣ Working directory
# ------------------------------------------------------------------------------
setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics")
cat("📁 Working directory:", getwd(), "\n")

# ------------------------------------------------------------------------------
# 2️⃣ Settings
# ------------------------------------------------------------------------------
folders      <- c("SCC1","SCC6","SCC47","SCC90","Group_comparison")

colors <- c(
"0Gy" = "#1f78b4",   # blue
"6Gy" = "#d73027"    # red
)

group_order <- c("0Gy", "6Gy")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ------------------------------------------------------------------------------
# GLOBAL KEYWORDS
# ------------------------------------------------------------------------------
global_keywords <- c(
# --- Cell movement ---
"migration","cell migration","chemotaxis","cell motility","invasion",
"metastasis","locomotion","cell movement","wound healing",

# --- Cell death / survival ---
"survival","cell survival","apoptosis","cell death","autophagy",
"cell cycle","mitosis","mitotic","senescence","cell differentiation",
"DNA replication",

# --- ECM / Matrix ---
"matrix","extracellular matrix","ECM","adhesion","cell adhesion",
"focal adhesion","angiogenesis","vasculature","vascular development",
"collagen","fibronectin","basement membrane","extracellular structure",
"integrin",

# --- Cancer / HPV ---
"cancer","carcinoma","tumor","tumour","neoplasm","oncogene",
"tumor suppressor","p53","radiation","HPV","human papillomavirus",
"malignancy","metastatic","tumorigenesis","carcinogenesis",

# --- Stress / DNA damage ---
"stress","cellular stress","genotoxic","oxidative stress",
"redox","oxidation","reduction","ROS","chaperone","heat shock",
"protein folding","UPR",

# --- Mitochondrial / Metabolic ---
"mitochondria","mitochondrial","metabolic","metabolism","lipid",
"fatty acid","glycolysis","gluconeogenesis","ATP",
"oxidative phosphorylation","TCA cycle",

# --- Vesicle / Exosome ---
"exosome","extracellular vesicle","vesicle","multivesicular body",
"endosome","endocytosis","exocytosis","secretion","microvesicle",
"vesicular transport","platelet",

# --- Platelet ---
"platelet","megakaryocyte","thrombocyte","thrombopoiesis",
"platelet activation","platelet aggregation","platelet degranulation",
"platelet derived","platelet factor","thrombus","thrombosis",
"coagulation","clot",

# --- Signaling ---
"MAPK","ERK","PI3K","AKT","mTOR","JAK","STAT","WNT","Notch","TGF"
)

# ------------------------------------------------------------------------------
# CLUSTER DEFINITIONS + COLORS
# ------------------------------------------------------------------------------
clusters <- list(
"Cell movement" = c(
"migration","cell migration","chemotaxis","cell motility",
"invasion","metastasis","locomotion","cell movement",
"wound healing"
),

"Cell survival / death" = c(
"survival","cell survival","apoptosis","cell death",
"autophagy","cell cycle","mitosis","mitotic","senescence",
"cell differentiation","DNA replication"
),

"ECM / Matrix" = c(
"matrix","extracellular matrix","ECM","adhesion","cell adhesion",
"focal adhesion","angiogenesis","vasculature",
"vascular development","collagen","fibronectin",
"basement membrane","extracellular structure","integrin"
),

"Cancer / HPV" = c(
"cancer","carcinoma","tumor","tumour","neoplasm","oncogene",
"tumor suppressor","p53","radiation","HPV","human papillomavirus",
"malignancy","metastatic","tumorigenesis","carcinogenesis"
),

"Stress / DNA Damage" = c(
"stress","cellular stress","genotoxic","oxidative stress",
"redox","oxidation","reduction","ROS","chaperone",
"heat shock","protein folding","UPR"
),

"Mitochondrial / Metabolic" = c(
"mitochondria","mitochondrial","metabolic","metabolism",
"lipid","fatty acid","glycolysis","gluconeogenesis",
"ATP","oxidative phosphorylation","TCA cycle"
),

"Vesicle / Exosome" = c(
"exosome","extracellular vesicle","vesicle","multivesicular body",
"endosome","endocytosis","exocytosis","secretion",
"microvesicle","vesicular transport"
),

"Platelet" = c(
"platelet","megakaryocyte","thrombocyte","thrombopoiesis",
"platelet activation","platelet aggregation",
"platelet degranulation","platelet derived",
"platelet factor","thrombus","thrombosis",
"coagulation","clot"
),

"Signaling pathways" = c(
"MAPK","ERK","PI3K","AKT","mTOR","JAK","STAT",
"WNT","Notch","TGF"
),

"Other" = character(0)
)


cluster_colors <- c(
"Cell movement"             = "#1b9e77",
"Cell survival / death"     = "#d95f02",
"ECM / Matrix"              = "#7570b3",
"Cancer / HPV"              = "#e7298a",
"Stress / DNA Damage"       = "#66a61e",
"Mitochondrial / Metabolic" = "#e6ab02",
"Vesicle / Exosome"         = "#a6761d",
"Platelet"                  = "#8c564b",
"Signaling pathways"        = "#666666",
"Other"                     = "#999999"
)

# ------------------------------------------------------------------------------
# CLEAN PATHWAY TEXT
# ------------------------------------------------------------------------------
clean_pathway <- function(p) {
p <- str_replace_all(p, "GO:\\d+", "")
p <- str_replace_all(p, "\\b(BP|CC|MF)\\b", "")
p <- str_replace_all(p, "KEGG\\s*\\|?", "")
p <- str_replace_all(p, "Reactome\\s*\\|?", "")
p <- str_replace_all(p, "[\\|]", " ")
str_squish(p)
}

# ------------------------------------------------------------------------------
# READ ONE FISHER RESULT
# ------------------------------------------------------------------------------
read_fisher_one <- function(file, group_label) {

if (!file.exists(file)) return(NULL)
df <- tryCatch(read.xlsx(file), error = function(e) NULL)
if (is.null(df) || !"Pathway" %in% names(df)) return(NULL)

pcol <- intersect(c("p_adj","p.adj","padj","adj_p_value","adj.P.Val"), names(df))[1]
gcol <- intersect(c("Genes","Gene","Gene_List","Proteins","Protein_List"), names(df))[1]
if (is.na(pcol) || is.na(gcol)) return(NULL)

df$p_adj <- suppressWarnings(as.numeric(df[[pcol]]))

df <- df %>%
filter(!is.na(p_adj), p_adj < FDR_CUTOFF,
!is.na(.data[[gcol]]), .data[[gcol]] != "") %>%
mutate(Pathway = clean_pathway(Pathway)) %>%
filter(str_detect(tolower(Pathway), paste(tolower(global_keywords), collapse="|"))) %>%
separate_rows(all_of(gcol), sep="[;,/\\s]+") %>%
mutate(Gene = str_trim(.data[[gcol]])) %>%
filter(Gene != "") %>%
group_by(Pathway) %>%
mutate(n_prot = n_distinct(Gene)) %>%
ungroup() %>%
filter(n_prot >= MIN_PROTEINS) %>%
mutate(Group = group_label) %>%
distinct(Pathway, Gene, p_adj, Group)

return(df)
}

# ------------------------------------------------------------------------------
# MAIN LOOP (ggalluvial)
# ------------------------------------------------------------------------------
cat("🔎 Generating Sankey (alluvial) diagrams...\n")

for (folder in folders) {

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("📂 Folder:", folder, "\n")

fisher_dir <- file.path(folder, "Fisher")
if (!dir.exists(fisher_dir)) next

f0 <- list.files(fisher_dir, pattern="Merged_Fisher_.*_0Gy\\.xlsx$", full.names=TRUE)
f6 <- list.files(fisher_dir, pattern="Merged_Fisher_.*_6Gy\\.xlsx$", full.names=TRUE)
if (!length(f0) || !length(f6)) next

d0 <- read_fisher_one(f0[1], "0Gy")
d6 <- read_fisher_one(f6[1], "6Gy")

dat <- bind_rows(d0, d6)
if (!nrow(dat)) next

keep_paths <- dat %>%
count(Pathway) %>% arrange(desc(n)) %>% head(TOP_N) %>% pull(Pathway)
dat <- dat %>% filter(Pathway %in% keep_paths)

# ==============================================================================
# 6. ASSIGN CLUSTERS
# ==============================================================================
# Assign each pathway to the appropriate biological cluster based on keyword matching.

dat$Cluster <- sapply(dat$Pathway, function(p) {
p.l <- tolower(p)
for (cl in names(clusters)) {
if (any(str_detect(p.l, tolower(clusters[[cl]])))) return(cl)
}
return("Other")
})

# ==============================================================================
# 7. PREPARE DATA FOR ALLUVIAL
# ==============================================================================
# Build frequency table Cluster → Pathway → Group + set factor levels.

sankey_df <- dat %>%
count(Cluster, Pathway, Group, name = "Freq") %>%
ungroup()

sankey_df$Cluster <- factor(sankey_df$Cluster, levels = names(cluster_colors))
sankey_df$Group   <- factor(sankey_df$Group, levels = group_order)

# Order pathways inside clusters
pathway_order <- sankey_df %>%
group_by(Cluster, Pathway) %>%
summarise(total = sum(Freq), .groups = "drop") %>%
arrange(
factor(Cluster, levels = names(cluster_colors)),
Pathway
) %>%
pull(Pathway)

sankey_df$Pathway <- factor(sankey_df$Pathway, levels = unique(pathway_order))

# ==============================================================================
# 8. PLOT (as in cell line script)
# ==============================================================================
# Construct alluvial plot with cluster-colored base layer and group-colored overlay.
# Add stratum labels and custom top/bottom borders only.

p <- ggplot(
sankey_df,
aes(axis1 = Cluster, axis3 = Pathway, axis5 = Group, y = Freq)
) +
geom_alluvium(aes(fill = Cluster),
width = 1/12, alpha = 0.3,
color = "black", size = 1) +
geom_alluvium(aes(fill = Group),
width = 1/12, alpha = 0.3,
linetype = 0, stat = "alluvium", na.rm = TRUE) +
geom_stratum(width = 0.09, fill = NA, color = NA) +
geom_text(stat = "stratum", aes(label = after_stat(stratum)),
size = 15, color = "black") +
scale_x_discrete(expand = c(.1, .1),
labels = c("Cluster", "Pathway", "Group")) +
scale_fill_manual(values = c(cluster_colors, colors)) +
theme_minimal(base_size = 24) +
theme(
legend.position = "none",
axis.text.y = element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
panel.grid = element_blank(),
plot.title = element_text(hjust = 0.5, size= 50)
) +
ggtitle(paste0("Clustered Sankey Diagram — ", folder, " (0 Gy vs 6 Gy)"))

pb <- ggplot_build(p)
idx <- which(
sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum"))
)
stratum_data <- pb$data[[idx]]

p <- p +
geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
inherit.aes = FALSE, color = "black", size = 7
) +
geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
inherit.aes = FALSE, color = "black", size = 7
)

# SAVE ---------------------------------------------------------------------
out_dir <- file.path(folder, "Plots_Sankey")
if (!dir.exists(out_dir)) dir.create(out_dir)

write.xlsx(sankey_df,
file.path(out_dir, paste0("Sankey_", folder, "_Clustered.xlsx")),
rownames = FALSE)

ggsave(
file.path(out_dir, paste0("Sankey_", folder, "_Clustered.png")),
p, width = 32, height = 25, dpi = 300, bg = "white"
)

ggsave(
file.path(out_dir, paste0("Sankey_", folder, "_Clustered.pdf")),
p, width = 32, height = 25, bg = "white"
)

cat("✅ Saved Sankey for:", folder, "\n")
}

cat("\n🎉 ALL DONE — Sankey diagrams with TOP/BOTTOM stratum borders only!\n")


################################################################################
################################################################################
################################################################################

################################################################################
# 📌 CLEAN Sankey: Pathway → Group (0Gy vs 6Gy) — ONE PLOT PER FOLDER
# ✅ Updated to MATCH STYLE of the “All-Groups 2-Way Sankey” script:
#    - width/alpha/size settings
#    - stratum style (no fill), labels (no wrap)
#    - pathway order A -> Z (after TOP_N filtering)
#    - margins + coord_cartesian(clip="off")
#    - black top/bottom borders on strata (geom_segment trick)
#
# 🔧 FIX (requested):
#    - Pathway labels are placed OUTSIDE the left stratum (left of the graph)
#    - Group labels are placed OUTSIDE the right stratum (right of the graph)
#    - No overlap with the strata/flows
#    - Nothing else changed
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggalluvial)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

# ==============================================================================
# 1. WORKING DIRECTORY
# ==============================================================================
setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics")
cat("📁 Working directory:", getwd(), "\n")

# ==============================================================================
# 2. SETTINGS
# ==============================================================================

folders <- c("SCC1","SCC6","SCC47","SCC90","Group_comparison")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ---- 0Gy vs 6Gy ----
group_order <- c("0Gy", "6Gy")
group_colors <- c(
  "0Gy" = "#1f78b4",
  "6Gy" = "#d73027"
)

# ==============================================================================
# 3. GLOBAL KEYWORDS
# ==============================================================================

global_keywords <- c(
  # --- Cell movement ---
  "migration","cell migration","chemotaxis","cell motility","invasion",
  "metastasis","locomotion","cell movement","wound healing",
  
  # --- Cell death / survival ---
  "survival","cell survival","apoptosis","cell death","autophagy",
  "cell cycle","mitosis","mitotic","senescence","cell differentiation",
  "DNA replication",
  
  # --- ECM / Matrix ---
  "matrix","extracellular matrix","ECM","adhesion","cell adhesion",
  "focal adhesion","angiogenesis","vasculature","vascular development",
  "collagen","fibronectin","basement membrane","extracellular structure",
  "integrin",
  
  # --- Cancer / HPV ---
  "cancer","carcinoma","tumor","tumour","neoplasm","oncogene",
  "tumor suppressor","p53","radiation","HPV","human papillomavirus",
  "malignancy","metastatic","tumorigenesis","carcinogenesis",
  
  # --- Stress / DNA damage ---
  "stress","cellular stress","genotoxic","oxidative stress",
  "redox","oxidation","reduction","ROS","chaperone","heat shock",
  "protein folding","UPR",
  
  # --- Mitochondrial / Metabolic ---
  "mitochondria","mitochondrial","metabolic","metabolism","lipid",
  "fatty acid","glycolysis","gluconeogenesis","ATP",
  "oxidative phosphorylation","TCA cycle",
  
  # --- Vesicle / Exosome ---
  "exosome","extracellular vesicle","vesicle","multivesicular body",
  "endosome","endocytosis","exocytosis","secretion","microvesicle",
  "vesicular transport","platelet",
  
  # --- Platelet ---
  "platelet","megakaryocyte","thrombocyte","thrombopoiesis",
  "platelet activation","platelet aggregation","platelet degranulation",
  "platelet derived","platelet factor","thrombus","thrombosis",
  "coagulation","clot",
  
  # --- Signaling ---
  "MAPK","ERK","PI3K","AKT","mTOR","JAK","STAT","WNT","Notch","TGF"
)

# ==============================================================================
# 4. HELPERS
# ==============================================================================

clean_pathway <- function(p) {
  p <- str_replace_all(p, "GO:\\d+", "")
  p <- str_replace_all(p, "\\b(BP|CC|MF)\\b", "")
  p <- str_replace_all(p, "KEGG\\s*\\|?", "")
  p <- str_replace_all(p, "Reactome\\s*\\|?", "")
  p <- str_replace_all(p, "[|]", " ")
  str_squish(p)
}

first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[1] else NA_character_
}

# ==============================================================================
# 5. READ FISHER FILE (robust column detection like reference script)
# ==============================================================================

read_fisher <- function(file, group) {
  
  df <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  
  # pathway column robust
  if (!"Pathway" %in% names(df)) {
    alt <- first_col(df, c("Pathway_Name","Clean_Path","Term","Description","NAME","name"))
    if (is.na(alt)) return(NULL)
    df$Pathway <- df[[alt]]
  }
  
  # p_adj + genes/proteins list robust
  pcol <- first_col(df, c("p_adj","p.adj","padj","adj_p_value","adj.P.Val","p.value","p_value"))
  gcol <- first_col(df, c("IDs","Genes","Gene","Gene_List","Proteins","Protein_List"))
  
  if (is.na(pcol) || is.na(gcol)) return(NULL)
  
  df$p_adj <- suppressWarnings(as.numeric(df[[pcol]]))
  df[[gcol]] <- as.character(df[[gcol]])
  
  kw_pat <- paste(stringr::str_to_lower(global_keywords), collapse = "|")
  
  out <- df %>%
    dplyr::filter(
      !is.na(p_adj),
      is.finite(p_adj),
      p_adj < FDR_CUTOFF,
      !is.na(.data[[gcol]]),
      .data[[gcol]] != ""
    ) %>%
    dplyr::mutate(Pathway = clean_pathway(Pathway)) %>%
    dplyr::filter(stringr::str_detect(stringr::str_to_lower(Pathway), kw_pat)) %>%
    tidyr::separate_rows(dplyr::all_of(gcol), sep = "[;,/\\s]+") %>%
    dplyr::mutate(Gene = stringr::str_trim(.data[[gcol]])) %>%
    dplyr::filter(Gene != "" & !is.na(Gene) & Gene != "NA") %>%
    dplyr::group_by(Pathway) %>%
    dplyr::mutate(n_prot = dplyr::n_distinct(Gene)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_prot >= MIN_PROTEINS) %>%
    dplyr::mutate(Group = group) %>%
    dplyr::distinct(Pathway, Gene, p_adj, Group)
  
  if (!nrow(out)) return(NULL)
  out
}

# ==============================================================================
# 6. MAIN LOOP — ONE SANKEY PER FOLDER (STYLE-MATCHED)
# ==============================================================================

for (folder in folders) {
  
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("📂 Processing:", folder, "\n")
  
  fisher_dir <- file.path(folder, "Fisher")
  if (!dir.exists(fisher_dir)) {
    cat("⚠️ No Fisher folder found. Skipping.\n")
    next
  }
  
  # Detect Fisher files (same patterns you used)
  f0 <- list.files(fisher_dir, pattern = "Merged_Fisher_.*_0Gy\\.xlsx$", full.names = TRUE)
  f6 <- list.files(fisher_dir, pattern = "Merged_Fisher_.*_6Gy\\.xlsx$", full.names = TRUE)
  
  if (!length(f0) || !length(f6)) {
    cat("⚠️ Missing 0Gy or 6Gy Fisher file. Skipping.\n")
    next
  }
  
  # Load data
  d0 <- read_fisher(f0[1], "0Gy")
  d6 <- read_fisher(f6[1], "6Gy")
  
  dat <- dplyr::bind_rows(d0, d6)
  if (!nrow(dat)) {
    cat("⚠️ No pathways detected in folder (after filters). Skipping.\n")
    next
  }
  
  # Keep only the two groups (and enforce order)
  dat <- dat %>% dplyr::filter(Group %in% group_order)
  
  # TOP N pathways (by overall frequency)
  top_paths <- dat %>%
    dplyr::count(Pathway, sort = TRUE) %>%
    dplyr::slice_head(n = TOP_N) %>%
    dplyr::pull(Pathway)
  
  dat <- dat %>% dplyr::filter(Pathway %in% top_paths)
  
  # Prep Sankey data (2-WAY: Pathway <-> Group)
  sankey_df <- dat %>%
    dplyr::count(Pathway, Group, name = "Freq") %>%
    dplyr::ungroup()
  
  # Enforce group order
  sankey_df$Group <- factor(sankey_df$Group, levels = group_order)
  
  # Pathway order A -> Z (match reference script)
  pathway_order <- sankey_df %>%
    dplyr::distinct(Pathway) %>%
    dplyr::pull(Pathway) %>%
    sort()
  
  sankey_df$Pathway <- factor(sankey_df$Pathway, levels = pathway_order)
  
  # ==============================================================================
  # 7. PLOT (MATCHED STYLE)
  # ==============================================================================
  
  # --- build base plot WITHOUT labels (we will place labels outside using stratum geometry) ---
  p <- ggplot(
    sankey_df,
    aes(axis1 = Pathway, axis2 = Group, y = Freq)
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = Group),
      width = 0.5,
      alpha = 0.3,
      color = "black",
      size = 1
    ) +
    ggalluvial::geom_stratum(
      width = 0.2,
      fill = NA,
      color = NA
    ) +
    ggplot2::scale_x_discrete(
      expand = c(.25, .25),
      labels = c("Pathway", "Group")
    ) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::theme_minimal(base_size = 40) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 70),
      plot.margin = ggplot2::margin(t = 20, r = 5, b = 20, l = 1500, unit = "pt")
    ) +
    coord_cartesian(clip = "off") +
    ggplot2::ggtitle("")
  
  # --- extract stratum rectangles from built plot ---
  pb <- ggplot_build(p)
  idx <- which(sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum")))
  stratum_data <- pb$data[[idx]]
  
  # --- add black top/bottom borders to each stratum (same trick) ---
  p <- p +
    ggplot2::geom_segment(
      data = stratum_data,
      aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
      inherit.aes = FALSE,
      color = "black",
      size = 3
    ) +
    ggplot2::geom_segment(
      data = stratum_data,
      aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
      inherit.aes = FALSE,
      color = "black",
      size = 3
    )
  
  # --- place labels OUTSIDE strata: Pathway left, Group right (no overlap) ---
  left_lab <- stratum_data %>%
    dplyr::filter(x == 1) %>%
    dplyr::mutate(
      x_lab = xmin - 0.10,
      y_lab = (ymin + ymax) / 2
    )
  
  right_lab <- stratum_data %>%
    dplyr::filter(x == 2) %>%
    dplyr::mutate(
      x_lab = xmax + 0.10,
      y_lab = (ymin + ymax) / 2
    )
  
  p <- p +
    ggplot2::geom_text(
      data = left_lab,
      aes(x = x_lab, y = y_lab, label = stratum),
      inherit.aes = FALSE,
      hjust = 1,
      size = 25,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = right_lab,
      aes(x = x_lab, y = y_lab, label = stratum),
      inherit.aes = FALSE,
      hjust = 0,
      size = 30,
      color = "black"
    )
  
  # ==============================================================================
  # 8. SAVE INSIDE EACH FOLDER
  # ==============================================================================
  
  out_dir <- file.path(folder, "Plots_Sankey_Pathway")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  openxlsx::write.xlsx(
    sankey_df,
    file.path(out_dir, paste0("Sankey_Data_", folder, "_0Gy_vs_6Gy.xlsx")),
    rowNames = FALSE,
    overwrite = TRUE
  )
  
  ggplot2::ggsave(
    file.path(out_dir, paste0("Sankey_", folder, "_0Gy_vs_6Gy.png")),
    p, width = 45, height = 45, dpi = 300, bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    file.path(out_dir, paste0("Sankey_", folder, "_0Gy_vs_6Gy.pdf")),
    p, width = 45, height = 45, bg = "white",
    limitsize = FALSE
  )
  
  cat("✅ Saved Sankey for:", folder, "\n")
}

cat("\n🎉 ALL DONE — one keyword-filtered 2-WAY Sankey per folder (STYLE-MATCHED) created!\n")
################################################################################

################################################################################
# Step 8 (UPDATED): TARGET HEATMAPS of MEAN GROUPS — FULL (conflicted-safe + ht_opt fix)
# - X-axis shows ONLY: 6Gy and 0Gy (LEFT = 6Gy, RIGHT = 0Gy)
# - Legend moved farther RIGHT WITHOUT moving heatmap (HEATMAP_LEGEND_PADDING)
# - Legend title spaced further from the colorbar (newline in title)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================
# CONFIG (edit only this)
# =========================
cfg <- list(
  folders = c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison"),
  
  target_genes = c("CD9", "CD63", "CD81", "PDCD6IP", "SDCBP", "TSG101"),
  gene_labels = c(
    CD9     = "CD9",
    CD63    = "CD63",
    CD81    = "CD81",
    PDCD6IP = "Alix",
    SDCBP   = "Syntenin-1",
    TSG101  = "TSG101"
  ),
  
  # output + file picking
  out_subdir = "Step8_TargetProtein_Heatmaps_MeanDose",
  pattern_step7 = "_corrFiltered\\.txt$",
  pattern_step5 = "_filtered_50percent\\.txt$",
  
  # heatmap look
  heatmap_name = "mean_log10_LFQ",
  na_col = "black",
  colors = c("white", "#4F81BD", "#C00000"),   # low, mid, high
  
  # text + layout
  font_row = 20,
  font_col = 22,
  col_angle = 0,
  title_fontsize = 20,
  
  # legend text
  legend_title_fontsize  = 18,
  legend_labels_fontsize = 18,
  legend_title_position  = "topcenter",
  
  # (kept for config compatibility; spacing is implemented via newline in title)
  legend_title_gap_mm    = 30,
  
  # heatmap -> legend distance (does NOT move heatmap)
  legend_gap_right_mm    = 16,
  
  # column order (LEFT -> RIGHT)
  dose_order = c("6Gy", "0Gy"),
  
  # image export
  width_in = 8,
  height_in = 8,
  dpi = 300,
  
  # plot padding (outer whitespace). padding order: top, right, bottom, left
  margins = c(t = 14, r = 14, b = 14, l = 30),
  margin_unit = "mm",
  
  legend_side = "right"
)

# =========================
# Helpers
# =========================
split_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  g <- unlist(strsplit(as.character(x), ";"))
  g <- stringr::str_trim(g)
  g[g != ""]
}

row_has_target <- function(pg_genes, targets) {
  any(split_genes(pg_genes) %in% targets)
}

pick_input_file <- function(folder, pattern7, pattern5) {
  f7 <- list.files(folder, pattern = pattern7, full.names = TRUE)
  if (length(f7) > 0) return(f7[1])
  f5 <- list.files(folder, pattern = pattern5, full.names = TRUE)
  if (length(f5) > 0) return(f5[1])
  NA_character_
}

# Parse LFQ_<CellLine>_<Dose>_<Rep>  -> return Dose
parse_dose_from_lfq <- function(x) {
  m <- stringr::str_match(x, "^LFQ_[^_]+_([^_]+)_\\d+$")
  m[, 2]
}

# Build matrix with 2 columns named EXACTLY: 6Gy and 0Gy
build_target_mean_matrix <- function(df, lfq_cols, target_genes, gene_labels, dose_order) {
  
  long <- df %>%
    dplyr::select(PG.Genes, dplyr::all_of(lfq_cols)) %>%
    dplyr::mutate(GeneHit = sapply(PG.Genes, function(x) {
      hits <- base::intersect(split_genes(x), target_genes)  # ✅ conflicted fix
      if (length(hits) == 0) NA_character_ else hits[1]
    })) %>%
    dplyr::filter(!is.na(GeneHit)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(lfq_cols),
                        names_to = "Sample", values_to = "Value") %>%
    dplyr::mutate(
      Value = suppressWarnings(as.numeric(Value)),
      Dose  = parse_dose_from_lfq(Sample),
      GeneLabel = unname(gene_labels[GeneHit])
    ) %>%
    dplyr::filter(!is.na(Dose))
  
  # If multiple rows map to same GeneLabel, collapse per sample using MAX
  collapsed_sample <- long %>%
    dplyr::group_by(GeneLabel, Dose, Sample) %>%
    dplyr::summarise(
      Value = if (all(is.na(Value))) NA_real_ else base::max(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # MEAN per dose (across replicates/samples)
  collapsed_dose <- collapsed_sample %>%
    dplyr::group_by(GeneLabel, Dose) %>%
    dplyr::summarise(
      Value = if (all(is.na(Value))) NA_real_ else base::mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  wide <- collapsed_dose %>%
    dplyr::mutate(Dose = factor(Dose, levels = dose_order)) %>%
    tidyr::pivot_wider(names_from = Dose, values_from = Value)
  
  ordered_labels <- unname(gene_labels[target_genes])
  all_rows <- tibble::tibble(GeneLabel = ordered_labels)
  wide <- all_rows %>% dplyr::left_join(wide, by = "GeneLabel")
  
  mat <- as.data.frame(wide)
  rownames(mat) <- mat$GeneLabel
  mat$GeneLabel <- NULL
  
  # Ensure both columns exist and in correct order
  for (d in dose_order) if (!d %in% colnames(mat)) mat[[d]] <- NA_real_
  mat <- mat[, dose_order, drop = FALSE]
  
  as.matrix(mat)
}

save_target_heatmap <- function(mat, base_fn, cfg) {
  rng <- range(mat, na.rm = TRUE)
  if (!all(is.finite(rng))) rng <- c(0, 1)
  
  col_fun <- circlize::colorRamp2(
    c(rng[1], base::mean(rng), rng[2]),
    cfg$colors
  )
  
  # newline in title for guaranteed spacing
  legend_param <- list(
    title          = paste0(cfg$heatmap_name, "\n"), # change to "\n\n" if you want more space
    title_gp       = grid::gpar(fontsize = cfg$legend_title_fontsize, fontface = "bold"),
    labels_gp      = grid::gpar(fontsize = cfg$legend_labels_fontsize),
    title_position = cfg$legend_title_position
  )
  
  ht <- Heatmap(
    mat,
    name = cfg$heatmap_name,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = cfg$na_col,
    row_names_gp = grid::gpar(fontsize = cfg$font_row),
    column_names_gp = grid::gpar(fontsize = cfg$font_col),
    column_names_rot = cfg$col_angle,
    column_title = NULL,
    heatmap_legend_param = legend_param
  )
  
  pad <- grid::unit(
    c(cfg$margins["t"], cfg$margins["r"], cfg$margins["b"], cfg$margins["l"]),
    cfg$margin_unit
  )
  
  # ✅ FIX: ht_opt is NOT namespaced
  old_leg_pad <- ht_opt$HEATMAP_LEGEND_PADDING
  ht_opt$HEATMAP_LEGEND_PADDING <- grid::unit(cfg$legend_gap_right_mm, "mm")
  on.exit({ ht_opt$HEATMAP_LEGEND_PADDING <- old_leg_pad }, add = TRUE)
  
  png(paste0(base_fn, ".png"),
      width = cfg$width_in, height = cfg$height_in,
      units = "in", res = cfg$dpi)
  draw(ht, padding = pad, heatmap_legend_side = cfg$legend_side)
  dev.off()
  
  pdf(paste0(base_fn, ".pdf"),
      width = cfg$width_in, height = cfg$height_in)
  draw(ht, padding = pad, heatmap_legend_side = cfg$legend_side)
  dev.off()
}

# =========================
# Run
# =========================
for (folder in cfg$folders) {
  cat("\n=== Step 8 (MEAN, 2-COL): Target heatmap for", folder, "===\n")
  
  in_file <- pick_input_file(folder, cfg$pattern_step7, cfg$pattern_step5)
  if (is.na(in_file)) {
    cat("⚠️ No input file found in", folder, "— skipping.\n")
    next
  }
  
  cat("Reading:", in_file, "\n")
  df <- read.table(in_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
  
  if (!"PG.Genes" %in% colnames(df)) {
    cat("⚠️ PG.Genes not found in", folder, "— skipping.\n")
    next
  }
  
  lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
  if (length(lfq_cols) < 2) {
    cat("⚠️ Not enough LFQ columns in", folder, "— skipping.\n")
    next
  }
  
  df_t <- df[sapply(df$PG.Genes, row_has_target, targets = cfg$target_genes), , drop = FALSE]
  cat("Matched rows:", nrow(df_t), "\n")
  if (nrow(df_t) == 0) {
    cat("⚠️ No target genes found in", folder, "— skipping.\n")
    next
  }
  
  mat <- build_target_mean_matrix(
    df_t, lfq_cols,
    cfg$target_genes, cfg$gene_labels,
    cfg$dose_order
  )
  
  out_dir <- file.path(folder, cfg$out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  mat_df <- as.data.frame(mat) %>% tibble::rownames_to_column("Protein")
  write.table(
    mat_df,
    file.path(out_dir, paste0(folder, "_TargetProteins_MeanDose_log10LFQ_matrix.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  openxlsx::write.xlsx(
    mat_df,
    file.path(out_dir, paste0(folder, "_TargetProteins_MeanDose_log10LFQ_matrix.xlsx")),
    overwrite = TRUE
  )
  
  base_fn <- file.path(out_dir, paste0(folder, "_TargetProteins_MeanDose_log10LFQ_heatmap"))
  save_target_heatmap(mat, base_fn, cfg)
  
  cat("✅ Saved outputs to:", out_dir, "\n")
}

cat("\n✅ Step complete — per-folder MEAN-dose target heatmaps saved.\n")

################################################################################

################################################################################
# Step 8: TARGET HEATMAPS (config-driven) — FULL explicit :: (conflicted-safe)
# - Uses base::intersect / base::max / base::mean to avoid conflicted errors
# - Uses ComplexHeatmap::Heatmap / ComplexHeatmap::draw
# - IMPORTANT: ht_opt is NOT namespaced (must stay ht_opt)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

cfg <- list(
  folders = c("SCC1", "SCC6", "SCC47", "SCC90", "Group_comparison"),
  
  target_genes = c("CD9", "CD63", "CD81", "PDCD6IP", "SDCBP", "TSG101"),
  gene_labels = c(
    CD9     = "CD9",
    CD63    = "CD63",
    CD81    = "CD81",
    PDCD6IP = "Alix",
    SDCBP   = "Syntenin-1",
    TSG101  = "TSG101"
  ),
  
  out_subdir = "Step8_TargetProtein_Heatmaps",
  pattern_step7 = "_corrFiltered\\.txt$",
  pattern_step5 = "_filtered_50percent\\.txt$",
  
  heatmap_name = "log10 LFQ",
  na_col = "black",
  colors = c("white", "#4F81BD", "#C00000"),
  
  font_row = 20,
  font_col = 20,
  col_angle = 90,
  title_fontsize = 20,
  
  legend_title_fontsize  = 18,
  legend_labels_fontsize = 18,
  legend_title_position  = "topcenter",
  
  legend_gap_right_mm = 16,
  legend_side = "right",
  
  dose_order = c("6Gy", "0Gy"),
  cellline_order = c("SCC1", "SCC47", "SCC90"),
  
  width_in = 9,
  height_in = 9,
  dpi = 300,
  
  margins = c(t = 12, r = 8, b = 14, l = 10),
  margin_unit = "mm"
)

split_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  g <- unlist(strsplit(as.character(x), ";"))
  g <- stringr::str_trim(g)
  g[g != ""]
}

row_has_target <- function(pg_genes, targets) {
  any(split_genes(pg_genes) %in% targets)
}

pick_input_file <- function(folder, pattern7, pattern5) {
  f7 <- list.files(folder, pattern = pattern7, full.names = TRUE)
  if (length(f7) > 0) return(f7[1])
  f5 <- list.files(folder, pattern = pattern5, full.names = TRUE)
  if (length(f5) > 0) return(f5[1])
  NA_character_
}

parse_lfq_name <- function(x) {
  m <- stringr::str_match(x, "^LFQ_([^_]+)_([^_]+)_(\\d+)$")
  data.frame(
    Sample   = x,
    CellLine = m[, 2],
    Dose     = m[, 3],
    Rep      = suppressWarnings(as.integer(m[, 4])),
    stringsAsFactors = FALSE
  )
}

order_lfq_cols <- function(cols,
                           dose_order = c("6Gy", "0Gy"),
                           cellline_order = c("SCC1", "SCC47", "SCC90")) {
  info <- parse_lfq_name(cols)
  ok <- !is.na(info$CellLine) & !is.na(info$Dose) & !is.na(info$Rep)
  
  info$DoseRank <- match(info$Dose, dose_order)
  info$CellRank <- match(info$CellLine, cellline_order)
  info$DoseRank[is.na(info$DoseRank)] <- 999L
  info$CellRank[is.na(info$CellRank)] <- 999L
  
  idx_ok <- which(ok)
  idx_sorted <- idx_ok[order(info$DoseRank[idx_ok],
                             info$CellRank[idx_ok],
                             info$Rep[idx_ok],
                             info$Sample[idx_ok])]
  
  cols_ok  <- info$Sample[idx_sorted]
  cols_bad <- info$Sample[!ok]
  c(cols_ok, cols_bad)
}

build_target_matrix <- function(df, lfq_cols, target_genes, gene_labels) {
  long <- df %>%
    dplyr::select(PG.Genes, dplyr::all_of(lfq_cols)) %>%
    dplyr::mutate(GeneHit = sapply(PG.Genes, function(x) {
      hits <- base::intersect(split_genes(x), target_genes)  # ✅ conflicted-safe
      if (length(hits) == 0) NA_character_ else hits[1]
    })) %>%
    dplyr::filter(!is.na(GeneHit)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(lfq_cols),
                        names_to = "Sample", values_to = "Value") %>%
    dplyr::mutate(
      Value     = suppressWarnings(as.numeric(Value)),
      GeneLabel = unname(gene_labels[GeneHit])
    )
  
  collapsed <- long %>%
    dplyr::group_by(GeneLabel, Sample) %>%
    dplyr::summarise(
      Value = if (all(is.na(Value))) NA_real_ else base::max(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  wide <- collapsed %>%
    tidyr::pivot_wider(names_from = Sample, values_from = Value)
  
  ordered_labels <- unname(gene_labels[target_genes])
  wide$GeneLabel <- factor(wide$GeneLabel, levels = ordered_labels)
  wide <- wide %>% dplyr::arrange(GeneLabel)
  
  mat <- as.data.frame(wide)
  rownames(mat) <- as.character(mat$GeneLabel)
  mat$GeneLabel <- NULL
  as.matrix(mat)
}

save_target_heatmap <- function(mat, title, base_fn, cfg) {
  rng <- range(mat, na.rm = TRUE)
  if (!all(is.finite(rng))) rng <- c(0, 1)
  
  col_fun <- circlize::colorRamp2(
    c(rng[1], base::mean(rng), rng[2]),
    cfg$colors
  )
  
  legend_param <- list(
    title          = paste0(cfg$heatmap_name, "\n"),
    title_gp       = grid::gpar(fontsize = cfg$legend_title_fontsize, fontface = "bold"),
    labels_gp      = grid::gpar(fontsize = cfg$legend_labels_fontsize),
    title_position = cfg$legend_title_position
  )
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = cfg$heatmap_name,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = cfg$na_col,
    row_names_gp = grid::gpar(fontsize = cfg$font_row),
    column_names_gp = grid::gpar(fontsize = cfg$font_col),
    column_names_rot = cfg$col_angle,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = cfg$title_fontsize, fontface = "bold"),
    heatmap_legend_param = legend_param
  )
  
  pad <- grid::unit(
    c(cfg$margins["t"], cfg$margins["r"], cfg$margins["b"], cfg$margins["l"]),
    cfg$margin_unit
  )
  
  # IMPORTANT: ht_opt is NOT namespaced
  old_leg_pad <- ht_opt$HEATMAP_LEGEND_PADDING
  ht_opt$HEATMAP_LEGEND_PADDING <- grid::unit(cfg$legend_gap_right_mm, "mm")
  on.exit({ ht_opt$HEATMAP_LEGEND_PADDING <- old_leg_pad }, add = TRUE)
  
  png(paste0(base_fn, ".png"),
      width = cfg$width_in, height = cfg$height_in,
      units = "in", res = cfg$dpi)
  ComplexHeatmap::draw(ht, padding = pad, heatmap_legend_side = cfg$legend_side)
  dev.off()
  
  pdf(paste0(base_fn, ".pdf"),
      width = cfg$width_in, height = cfg$height_in)
  ComplexHeatmap::draw(ht, padding = pad, heatmap_legend_side = cfg$legend_side)
  dev.off()
}

for (folder in cfg$folders) {
  cat("\n=== Step 8: Target heatmap for", folder, "===\n")
  
  in_file <- pick_input_file(folder, cfg$pattern_step7, cfg$pattern_step5)
  if (is.na(in_file)) {
    cat("⚠️ No input file found in", folder, "— skipping.\n")
    next
  }
  
  cat("Reading:", in_file, "\n")
  df <- read.table(in_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
  
  if (!"PG.Genes" %in% colnames(df)) {
    cat("⚠️ PG.Genes not found in", folder, "— skipping.\n")
    next
  }
  
  lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
  if (length(lfq_cols) < 2) {
    cat("⚠️ Not enough LFQ columns in", folder, "— skipping.\n")
    next
  }
  
  df_t <- df[sapply(df$PG.Genes, row_has_target, targets = cfg$target_genes), , drop = FALSE]
  cat("Matched rows:", nrow(df_t), "\n")
  if (nrow(df_t) == 0) {
    cat("⚠️ No target genes found in", folder, "— skipping.\n")
    next
  }
  
  mat <- build_target_matrix(df_t, lfq_cols, cfg$target_genes, cfg$gene_labels)
  mat <- mat[, order_lfq_cols(colnames(mat), cfg$dose_order, cfg$cellline_order), drop = FALSE]
  
  out_dir <- file.path(folder, cfg$out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  mat_df <- as.data.frame(mat) %>% tibble::rownames_to_column("Protein")
  write.table(
    mat_df,
    file.path(out_dir, paste0(folder, "_TargetProteins_log10LFQ_matrix.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  openxlsx::write.xlsx(
    mat_df,
    file.path(out_dir, paste0(folder, "_TargetProteins_log10LFQ_matrix.xlsx")),
    overwrite = TRUE
  )
  
  heat_title <- paste0("EV markers (log10 LFQ) – ", folder)
  base_fn <- file.path(out_dir, paste0(folder, "_TargetProteins_log10LFQ_heatmap"))
  
  save_target_heatmap(mat, heat_title, base_fn, cfg)
  
  cat("✅ Saved outputs to:", out_dir, "\n")
}

cat("\n✅ Step complete — per-folder target heatmaps saved.\n")
