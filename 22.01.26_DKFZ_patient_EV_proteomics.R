###############################################################################
############################## 1️⃣ FILTERING ###################################
################################################################################

library(clValid)
library(openxlsx)  # for saving Excel files

# --- Set working directory ---
setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics")

# --- Load and preprocess data ---
data <- read.table(
"DKFZ_HNSCC_patient_EV.txt",
header = TRUE, sep = "\t", quote = "", check.names = FALSE
)

# --- Check available columns ---
cat("Columns in dataset:\n")
print(colnames(data))

# --- Make sure target columns exist ---
required_cols <- c(
"Only identified by site", "Reverse", "Potential contaminant",
"PG.Organisms", "PG.ProteinDescriptions"
)
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
warning(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# --- Apply filtering safely ---
data_filtered <- data

# Remove rows marked by MaxQuant flags
if ("Only identified by site" %in% colnames(data_filtered)) {
data_filtered <- subset(data_filtered, `Only identified by site` != "+")
}
if ("Reverse" %in% colnames(data_filtered)) {
data_filtered <- subset(data_filtered, Reverse != "+")
}
if ("Potential contaminant" %in% colnames(data_filtered)) {
data_filtered <- subset(data_filtered, `Potential contaminant` != "+")
}

# Remove rows with "Potential contaminant" in PG.Organisms
if ("PG.Organisms" %in% colnames(data_filtered)) {
data_filtered <- subset(
data_filtered,
!grepl("Potential contaminant", PG.Organisms, ignore.case = TRUE)
)
}

# Remove rows with "Bos taurus" in PG.ProteinDescriptions
if ("PG.ProteinDescriptions" %in% colnames(data_filtered)) {
data_filtered <- subset(
data_filtered,
!grepl("Bos taurus", PG.ProteinDescriptions, ignore.case = TRUE)
)
}

# --- Report summary ---
cat("\nSummary:\n")
cat("Original rows:", nrow(data), "\n")
cat("Filtered rows:", nrow(data_filtered), "\n")
cat("Removed rows:", nrow(data) - nrow(data_filtered), "\n")

# --- Reset row names ---
row.names(data_filtered) <- NULL

# --- Save filtered data (TXT & Excel) ---
write.table(
data_filtered,
"DKFZ_HNSCC_patient_EV_filtered.txt",
sep = "\t", quote = FALSE, row.names = FALSE
)

write.xlsx(
data_filtered,
"DKFZ_HNSCC_patient_EV_filtered.xlsx",
overwrite = TRUE
)

cat("\nFiltering complete for DKFZ_HNSCC_patient_EV. Files saved.\n")


################################################################################
############################## 2️⃣ LOG10 TRANSFORM ##############################
################################################################################

# --- Identify LFQ intensity columns (based on naming pattern) ---
lfq_cols <- grep("^LFQ_", colnames(data_filtered), value = TRUE)

if (length(lfq_cols) == 0) {
warning("No LFQ intensity columns detected. Check column names!")
} else {
cat("LFQ columns detected:\n")
print(lfq_cols)
}

# --- Copy before transforming ---
log10_data_filtered <- data_filtered

# --- Robust log10 transform of LFQ columns ---
log10_data_filtered[lfq_cols] <- lapply(data_filtered[lfq_cols], function(x) {
# 1️⃣ Force to character, remove commas (thousands separators), convert to numeric
x_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))

# 2️⃣ Replace zeros or negatives with NA (to avoid -Inf)
x_num[x_num <= 0] <- NA_real_

# 3️⃣ Apply log10, replace non-finite values with NA
y <- log10(x_num)
y[!is.finite(y)] <- NA_real_
y
})

cat("\nLog10 transform complete on", length(lfq_cols), "LFQ columns.\n")

# --- Save transformed data (TXT) ---
write.table(
log10_data_filtered,
"DKFZ_HNSCC_patient_EV_filtered_log10LFQ.txt",
sep = "\t", quote = FALSE, row.names = FALSE
)

# --- Save transformed data (Excel) ---
write.xlsx(
log10_data_filtered,
"DKFZ_HNSCC_patient_EV_filtered_log10LFQ.xlsx",
overwrite = TRUE
)

cat("\nLog10-transformed files saved successfully.\n")

################################################################################

################################################################################
############################## 3️⃣ QC VALIDATION ###############################
################################################################################

# QC on log10 transformed dataset (old reorder step removed)

df <- log10_data_filtered

# --- Identify LFQ columns ---
lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols) == 0) stop("❌ QC Step 3: No LFQ_ columns detected!")

cat("✅ QC Step 3: LFQ columns detected:", length(lfq_cols), "\n")

# --- Check metadata columns (optional but informative) ---
meta_cols_requested <- c("PG.ProteinGroups", "PG.Genes", "PG.Organisms", "PG.ProteinDescriptions")
meta_cols <- intersect(meta_cols_requested, colnames(df))

if (length(meta_cols) < length(meta_cols_requested)) {
warning(paste(
"QC Step 3: Missing metadata columns (continuing):",
paste(setdiff(meta_cols_requested, meta_cols), collapse = ", ")
))
} else {
cat("✅ QC Step 3: All requested metadata columns present.\n")
}

# --- Check row/column counts (basic sanity) ---
cat("\nDataset dimensions:\n")
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n\n")

# --- Ensure LFQ columns are numeric (robust check) ---
lfq_is_numeric <- sapply(df[, lfq_cols, drop = FALSE], is.numeric)
if (!all(lfq_is_numeric)) {
non_num <- names(lfq_is_numeric)[!lfq_is_numeric]
warning("⚠️ QC Step 3: These LFQ columns are NOT numeric:\n",
paste(non_num, collapse = ", "),
"\n(They may still behave fine if coerced later, but best to fix upstream.)")
} else {
cat("✅ QC Step 3: All LFQ columns are numeric.\n")
}

# --- Check for non-finite values (Inf / -Inf) ---
mat <- as.matrix(df[, lfq_cols, drop = FALSE])
non_finite_count <- sum(!is.finite(mat), na.rm = TRUE)

# Note: is.finite(NA) is NA, so na.rm=TRUE ignores NAs.
# We specifically check Inf/-Inf by testing finite on non-NA values:
inf_count <- sum(is.infinite(mat), na.rm = TRUE)

if (inf_count > 0) {
warning("❌ QC Step 3: Found Inf/-Inf values in LFQ matrix! Count = ", inf_count,
"\nThis usually means zeros/negatives were not converted to NA before log10.")
} else {
cat("✅ QC Step 3: No Inf/-Inf values in LFQ columns.\n")
}

# --- Report NA rates per group + overall ---
na_total <- sum(is.na(mat))
total_vals <- length(mat)
cat("\nNA summary (LFQ only):\n")
cat("Total LFQ values:", total_vals, "\n")
cat("NA values:", na_total, "\n")
cat("NA fraction:", round(na_total / total_vals, 4), "\n")

# --- Optional: show worst columns by NA fraction ---
na_frac_by_col <- colMeans(is.na(mat))
worst <- sort(na_frac_by_col, decreasing = TRUE)[1:min(10, length(na_frac_by_col))]

cat("\nTop NA-heavy LFQ columns (up to 10):\n")
print(round(worst, 3))

cat("\n✅ Step 3 QC complete.\n")


################################################################################
#################### 4️⃣ GROUP-wise 50% FILTER (robust) ########################
################################################################################

library(openxlsx)

# Use log10 transformed data
df <- log10_data_filtered

main_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"

# ---- Metadata columns (keep those that exist) ----
meta_cols_requested <- c("PG.ProteinGroups","PG.Genes","PG.Organisms","PG.ProteinDescriptions")
meta_cols <- intersect(meta_cols_requested, colnames(df))

if (length(meta_cols) < length(meta_cols_requested)) {
warning(paste(
"Step 4: Missing metadata columns (will continue):",
paste(setdiff(meta_cols_requested, meta_cols), collapse = ", ")
))
}

# ---- Identify LFQ columns ----
lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols) == 0) stop("❌ No LFQ columns detected!")

# ---- Hard check: enforce LFQ_<GROUP>_<repNumber> naming ----
bad <- lfq_cols[!grepl("^LFQ_.+_[0-9]+$", lfq_cols)]
if (length(bad) > 0) {
stop("❌ These LFQ columns don't match 'LFQ_<GROUP>_<repNumber>':\n",
paste(bad, collapse = "\n"))
}

# ---- Extract group names ----
get_group <- function(x) sub("^LFQ_(.+)_([0-9]+)$", "\\1", x)
all_groups <- unique(get_group(lfq_cols))

# Optional: enforce a stable order (Healthy -> Pre_RT -> Post_RT -> anything else)
stages <- c("Healthy", "Pre_RT", "Post_RT")
groups <- c(intersect(stages, all_groups), setdiff(all_groups, stages))

cat("Detected groups:\n")
print(groups)

# ---- LFQ columns per group ----
group_cols <- lapply(groups, function(g) {
grep(paste0("^LFQ_", g, "_[0-9]+$"), lfq_cols, value = TRUE)
})
names(group_cols) <- groups

# ---- Force expression matrix numeric ----
expr <- df[, lfq_cols, drop = FALSE]
expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(as.character(x))))

# ---- Group-wise 50% filter (OR logic): keep if ANY group has >=50% non-NA ----
keep <- rep(FALSE, nrow(expr))

for (g in groups) {
cols <- group_cols[[g]]
if (length(cols) == 0) next

n_valid <- rowSums(!is.na(expr[, cols, drop = FALSE]))
thr <- ceiling(length(cols) * 0.50)

keep <- keep | (n_valid >= thr)

cat(sprintf("Group %-10s : %d samples, threshold=%d, pass=%d\n",
g, length(cols), thr, sum(n_valid >= thr)))
}

cat("\nTotal proteins:", nrow(df), "\n")
cat("Proteins kept (group-wise ≥50% in ANY group):", sum(keep), "\n")

# ---- Filter dataset ----
df_50 <- df[keep, , drop = FALSE]

# ---- Save ----
out_txt  <- file.path(main_dir, "DKFZ_HNSCC_patient_EV_filtered_50percent.txt")
out_xlsx <- file.path(main_dir, "DKFZ_HNSCC_patient_EV_filtered_50percent.xlsx")

write.table(df_50, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(df_50, out_xlsx, overwrite = TRUE)

cat("\n✅ Step 4 complete — GROUP-wise 50% filtered dataset saved.\n")
cat("   →", out_txt, "\n")
cat("   →", out_xlsx, "\n")


################################################################################
############################## 5️⃣ PROTEIN IDS PLOT ############################
################################################################################

library(dplyr)
library(ggplot2)

wd <- getwd()

file50 <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent.txt")
if (!file.exists(file50)) stop("❌ ERROR: filtered_50percent file not found in working directory!")

df <- read.table(file50, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
cat("Loaded file:", file50, "with", nrow(df), "proteins.\n")

# --- Extract LFQ columns ---
sample_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(sample_cols) == 0) stop("❌ No LFQ columns found!")

# --- Order LFQ columns by stage + replicate (important for cumulative/shared) ---
get_group <- function(x) sub("^LFQ_(.+)_([0-9]+)$", "\\1", x)
get_rep   <- function(x) as.numeric(sub("^.*_([0-9]+)$", "\\1", x))

stages <- c("Healthy", "Pre_RT", "Post_RT")

grp <- get_group(sample_cols)
rep <- get_rep(sample_cols)

grp_order <- match(grp, stages)
grp_order[is.na(grp_order)] <- length(stages) + 1  # unknown groups go last

sample_cols <- sample_cols[order(grp_order, grp, rep)]

cat("Detected", length(sample_cols), "sample columns.\n")
cat("LFQ order used for cumulative/shared:\n")
print(sample_cols)

# --- Compute protein ID information ---
protein_matrix <- !is.na(as.matrix(df[, sample_cols, drop = FALSE]))
protein_counts <- colSums(protein_matrix)

cumulative <- sapply(seq_along(sample_cols), function(i)
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) > 0))

shared <- sapply(seq_along(sample_cols), function(i)
sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) == i))

df_plot <- data.frame(
Sample = factor(sample_cols, levels = sample_cols),
ProteinCount = protein_counts,
Cumulative = cumulative,
Shared = shared
)

y_max <- max(df_plot$ProteinCount, df_plot$Cumulative, df_plot$Shared) * 1.15
plot_width <- max(12, length(sample_cols) * 0.35)

# --- Plot ---
p <- ggplot(df_plot, aes(x = Sample)) +
geom_bar(aes(y = ProteinCount), stat = "identity", fill = "#0072B2", width = 0.8) +
geom_line(aes(y = Cumulative, color = "Cumulative"), group = 1, linewidth = 1.2) +
geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2) +
geom_line(aes(y = Shared, color = "Shared"), group = 1, linewidth = 1.2) +
geom_point(aes(y = Shared, color = "Shared"), size = 2) +
scale_color_manual(values = c("Cumulative" = "forestgreen", "Shared" = "orange")) +
labs(
y = "Number of Proteins",
x = "Samples",
title = "Protein IDs per Sample (50% Filtered)"
) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
theme_classic(base_size = 18) +
theme(
legend.title = element_blank(),
legend.position = c(0.15, 0.15),
axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
axis.text.y = element_text(size = 16),
axis.title.x = element_text(size = 18, face = "bold"),
axis.title.y = element_text(size = 18, face = "bold"),
plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
)

base_name <- "DKFZ_HNSCC_patient_EV_filtered_50percent"

ggsave(paste0(base_name, "_proteinIDs.png"),
p, width = plot_width, height = 6, dpi = 300, bg = "white")

# safer PDF (works without Cairo)
ggsave(paste0(base_name, "_proteinIDs.pdf"),
p, width = plot_width, height = 6, dpi = 300, bg = "white", device = "pdf")

cat("\n🎉 Step 5 done — saved as:\n",
paste0(base_name, "_proteinIDs.png"), "\n",
paste0(base_name, "_proteinIDs.pdf"), "\n")


################################################################################
############################## 6️⃣ PROTEIN TOTALS PLOT #########################
################################################################################

library(dplyr)
library(ggplot2)
library(openxlsx)

wd <- getwd()

file50 <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent.txt")
if (!file.exists(file50)) stop("❌ ERROR: filtered_50percent file not found in working directory!")

df <- read.table(file50, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
cat("Loaded:", file50, "with", nrow(df), "proteins\n")

# --- Identify LFQ columns ---
sample_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(sample_cols) == 0) stop("❌ No LFQ columns found!")

# --- Convert LFQ to numeric (safe) ---
df[sample_cols] <- lapply(df[sample_cols], function(x) suppressWarnings(as.numeric(as.character(x))))

# --- Count proteins per sample: non-NA (log10 data => NA means missing) ---
protein_counts <- sapply(sample_cols, function(col) sum(!is.na(df[[col]])))

# --- Extract stages from LFQ_<Stage>_<rep> ---
sample_info <- data.frame(
Sample = sample_cols,
Stage  = gsub("^LFQ_|_[0-9]+$", "", sample_cols),
stringsAsFactors = FALSE
)

df_protein <- data.frame(
Stage = sample_info$Stage,
ProteinCount = as.numeric(protein_counts)
)

# --- Force stage order ---
stage_order <- c("Healthy", "Pre_RT", "Post_RT")
df_protein$Stage <- factor(df_protein$Stage, levels = stage_order)

# --- Colors (requested) ---
stage_colors <- c(
"Healthy" = "#66bd63",   # choose any distinct Healthy color you like
"Pre_RT"  = "#4575b4",
"Post_RT" = "#d73027"
)

# --- Summary mean ± SD (ALL groups) ---
summary_protein <- df_protein %>%
group_by(Stage) %>%
summarise(
mean = mean(ProteinCount, na.rm = TRUE),
sd   = sd(ProteinCount, na.rm = TRUE),
.groups = "drop"
)

y_max <- max(summary_protein$mean + summary_protein$sd,
df_protein$ProteinCount, na.rm = TRUE) * 1.10

# --- Plot 1 (ORIGINAL: includes Healthy) ---
p6 <- ggplot(summary_protein, aes(x = Stage, y = mean, fill = Stage)) +
geom_bar(stat = "identity", width = 0.7, color = "black") +
geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
width = 0.3, linewidth = 1.1, color = "black") +
geom_point(data = df_protein,
aes(x = Stage, y = ProteinCount),
position = position_jitter(width = 0.12),
size = 3.5, color = "black", fill = "white",
shape = 21, stroke = 1.1) +
scale_fill_manual(values = stage_colors, drop = FALSE) +
labs(
x = NULL,
y = "Number of Identified Proteins",
title = NULL
) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
theme_classic(base_size = 20) +
theme(
axis.text.x = element_text(size = 20, color = "black", face = "bold"),
axis.text.y = element_text(size = 20, color = "black", face = "bold"),
axis.title.x = element_text(size = 22, face = "bold"),
axis.title.y = element_text(size = 22, face = "bold"),
plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
legend.position = "none"
)

base_name <- "DKFZ_HNSCC_patient_EV_filtered_50percent"

ggsave(paste0(base_name, "_proteinTotals.png"),
p6, width = 10, height = 10, dpi = 300, bg = "white")

ggsave(paste0(base_name, "_proteinTotals.pdf"),
p6, width = 10, height = 10, dpi = 300, bg = "white", device = "pdf")

cat("\n🎉 Step 6 done — saved:\n",
paste0(base_name, "_proteinTotals.png"), "\n",
paste0(base_name, "_proteinTotals.pdf"), "\n")


# ==============================================================================
# ✅ ADDITIONAL PLOT: EXCLUDING HEALTHY
# ==============================================================================

df_protein_noHealthy <- df_protein %>%
filter(Stage %in% c("Pre_RT", "Post_RT")) %>%
droplevels()

summary_protein_noHealthy <- df_protein_noHealthy %>%
group_by(Stage) %>%
summarise(
mean = mean(ProteinCount, na.rm = TRUE),
sd   = sd(ProteinCount, na.rm = TRUE),
.groups = "drop"
)

y_max_noHealthy <- max(summary_protein_noHealthy$mean + summary_protein_noHealthy$sd,
df_protein_noHealthy$ProteinCount, na.rm = TRUE) * 1.10

p6_noHealthy <- ggplot(summary_protein_noHealthy, aes(x = Stage, y = mean, fill = Stage)) +
geom_bar(stat = "identity", width = 0.7, color = "black") +
geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
width = 0.3, linewidth = 1.1, color = "black") +
geom_point(data = df_protein_noHealthy,
aes(x = Stage, y = ProteinCount),
position = position_jitter(width = 0.12),
size = 3.5, color = "black", fill = "white",
shape = 21, stroke = 1.1) +
scale_fill_manual(values = stage_colors, drop = FALSE) +
labs(
x = NULL,
y = "Number of Identified Proteins",
title = NULL
) +
scale_y_continuous(expand = c(0, 0), limits = c(0, y_max_noHealthy)) +
theme_classic(base_size = 20) +
theme(
axis.text.x = element_text(size = 32, color = "black", face = "bold"),
axis.text.y = element_text(size = 32, color = "black", face = "bold"),
axis.title.x = element_text(size = 32, face = "bold"),
axis.title.y = element_text(size = 32, face = "bold"),
plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
legend.position = "none"
)

ggsave(paste0(base_name, "_proteinTotals_noHealthy.png"),
p6_noHealthy, width = 10, height = 10, dpi = 300, bg = "white")

ggsave(paste0(base_name, "_proteinTotals_noHealthy.pdf"),
p6_noHealthy, width = 10, height = 10, dpi = 300, bg = "white", device = "pdf")

cat("\n🎉 Additional plot (No Healthy) saved:\n",
paste0(base_name, "_proteinTotals_noHealthy.png"), "\n",
paste0(base_name, "_proteinTotals_noHealthy.pdf"), "\n")


# ==============================================================================
# ✅ EXPORT SUMMARY TABLES (mean ± SD) TO EXCEL (with & without Healthy)
# ==============================================================================

# Make sure Stage is plain text in the exported tables
summary_protein_export <- summary_protein %>%
mutate(Stage = as.character(Stage)) %>%
rename(Mean_ProteinCount = mean, SD_ProteinCount = sd) %>%
mutate(`Mean±SD` = paste0(round(Mean_ProteinCount, 1), " ± ", round(SD_ProteinCount, 1)))

summary_protein_noHealthy_export <- summary_protein_noHealthy %>%
mutate(Stage = as.character(Stage)) %>%
rename(Mean_ProteinCount = mean, SD_ProteinCount = sd) %>%
mutate(`Mean±SD` = paste0(round(Mean_ProteinCount, 1), " ± ", round(SD_ProteinCount, 1)))

out_summary_xlsx <- file.path(wd, paste0(base_name, "_proteinTotals_summary_meanSD.xlsx"))

wb <- createWorkbook()

addWorksheet(wb, "With_Healthy")
writeData(wb, "With_Healthy", summary_protein_export)

addWorksheet(wb, "No_Healthy")
writeData(wb, "No_Healthy", summary_protein_noHealthy_export)

saveWorkbook(wb, out_summary_xlsx, overwrite = TRUE)

cat("\n✅ Summary Excel saved:\n", out_summary_xlsx, "\n")



################################################################################
############################## 7️⃣ GLOBAL CORRELATION QC ########################
################################################################################

library(openxlsx)
library(pheatmap)

wd <- getwd()
file50 <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent.txt")

if (!file.exists(file50)) {
stop("❌ ERROR: DKFZ_HNSCC_patient_EV_filtered_50percent.txt not found in working directory!")
}

cat("\nLoading file:", file50, "\n")
df <- read.table(file50, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Extract LFQ columns
lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols) == 0) stop("❌ ERROR: No LFQ_ columns found in the dataset!")

cat("\nTotal LFQ samples included in correlation:", length(lfq_cols), "\n")

# Ensure numeric matrix
expr_matrix <- df[, lfq_cols, drop = FALSE]
expr_matrix[] <- lapply(expr_matrix, function(x) suppressWarnings(as.numeric(as.character(x))))
expr_matrix <- as.matrix(expr_matrix)

# Compute Pearson correlation matrix
corr <- cor(expr_matrix, use = "pairwise.complete.obs", method = "pearson")

base_name <- "DKFZ_HNSCC_patient_EV_filtered_50percent"

# Save correlation matrix (TXT + XLSX)
write.table(
corr,
file.path(wd, paste0(base_name, "_correlation_matrix.txt")),
sep = "\t", quote = FALSE, col.names = NA
)

write.xlsx(
corr,
file.path(wd, paste0(base_name, "_correlation_matrix.xlsx")),
overwrite = TRUE
)

# Compute mean correlation per sample (exclude diagonal properly)
mean_corr <- sapply(seq_len(ncol(corr)), function(i) mean(corr[i, -i], na.rm = TRUE))
names(mean_corr) <- colnames(corr)

summary_df <- data.frame(
Sample = names(mean_corr),
Mean_Correlation = round(mean_corr, 3),
row.names = NULL
)

write.table(
summary_df,
file.path(wd, paste0(base_name, "_correlation_summary.txt")),
sep = "\t", quote = FALSE, row.names = FALSE
)

write.xlsx(
summary_df,
file.path(wd, paste0(base_name, "_correlation_summary.xlsx")),
overwrite = TRUE
)

# Heatmap (PNG)
pheatmap(
corr,
main = "Global LFQ Sample Correlation (50% Filtered)",
cluster_rows = TRUE,
cluster_cols = TRUE,
display_numbers = TRUE,
number_format = "%.2f",
color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
filename = file.path(wd, paste0(base_name, "_correlation_heatmap.png")),
width = 9, height = 9, dpi = 300
)

# Heatmap (PDF)
pdf(file.path(wd, paste0(base_name, "_correlation_heatmap.pdf")), width = 9, height = 9)
pheatmap(
corr,
main = "Global LFQ Sample Correlation (50% Filtered)",
cluster_rows = TRUE,
cluster_cols = TRUE,
display_numbers = TRUE,
number_format = "%.2f",
color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
dev.off()

cat("\n✅ Step 7 complete — Global correlation files saved in working directory.\n")

################################################################################
############################## 8️⃣ GLOBAL REMOVE LOW-CORR SAMPLES ##############
################################################################################

library(openxlsx)
library(dplyr)

corr_threshold <- 0.70
cat("\nSamples with mean correlation <", corr_threshold, "will be removed.\n")

wd <- getwd()

file50 <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent.txt")
if (!file.exists(file50)) stop("❌ ERROR: 50% filtered file not found!")

df <- read.table(file50, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

summary_file <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent_correlation_summary.txt")
if (!file.exists(summary_file)) stop("❌ ERROR: Global correlation summary file not found! Run Step 7 first.")

corr_summary <- read.table(summary_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# LFQ columns present in the dataset
lfq_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols) == 0) stop("❌ ERROR: No LFQ_ columns found in dataset!")

# Decide remove/keep from summary
remove_samples <- corr_summary$Sample[corr_summary$Mean_Correlation < corr_threshold]
keep_samples   <- corr_summary$Sample[corr_summary$Mean_Correlation >= corr_threshold]

# Ensure samples exist in df (prevents subsetting errors)
remove_samples <- intersect(remove_samples, lfq_cols)
keep_samples   <- intersect(keep_samples, lfq_cols)

cat("\nSamples to REMOVE (present in df):\n"); print(remove_samples)
cat("\nSamples to KEEP (present in df):\n");   print(keep_samples)

# Keep metadata + kept LFQs
meta_cols <- setdiff(colnames(df), lfq_cols)
keep_cols <- c(meta_cols, keep_samples)

df_clean <- df[, keep_cols, drop = FALSE]

cat("\nFinal number of samples kept:", length(keep_samples), "\n")

base_out <- "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered"

write.table(
df_clean,
file.path(wd, paste0(base_out, ".txt")),
sep = "\t", quote = FALSE, row.names = FALSE
)

write.xlsx(
df_clean,
file.path(wd, paste0(base_out, ".xlsx")),
overwrite = TRUE
)

# Save removed sample summary if anything removed
if (length(remove_samples) > 0) {
removed_df <- corr_summary[corr_summary$Sample %in% remove_samples, , drop = FALSE]

write.xlsx(
removed_df,
file.path(wd, "Step8_Removed_Samples_Summary.xlsx"),
overwrite = TRUE
)

cat("\nRemoved samples summary saved.\n")
} else {
cat("\nNo samples removed at this threshold.\n")
}

cat("\n🎉 Step 8 complete — global low-correlation samples removed.\n")


################################################################################
############################## 9️⃣ IMPUTATION (PER-COLUMN, Perseus-like) ########
################################################################################

library(openxlsx)
set.seed(1)

wd <- getwd()
infile <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered.txt")
if (!file.exists(infile)) stop("❌ ERROR: corrFiltered file not found. Run Step 8 first.")

cat("\n📄 Loading:", basename(infile), "\n")
df <- read.table(infile, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

expr_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(expr_cols) == 0) stop("❌ ERROR: No LFQ columns found for imputation.")
cat("Detected", length(expr_cols), "LFQ columns.\n")

# Force LFQ columns to numeric safely
for (cn in expr_cols) {
df[[cn]] <- suppressWarnings(as.numeric(gsub(",", "", as.character(df[[cn]]))))
df[[cn]][!is.finite(df[[cn]])] <- NA_real_
}

# Matrix with stable dimnames
X <- as.matrix(df[, expr_cols, drop = FALSE])
storage.mode(X) <- "double"

shift <- 1.8
width <- 0.3

imputed_n <- 0L
param_table <- data.frame(
Sample = expr_cols,
Mean = NA_real_,
SD = NA_real_,
Imp_Mean = NA_real_,
Imp_SD = NA_real_,
row.names = NULL
)

for (j in seq_along(expr_cols)) {
col <- X[, j]
obs <- col[!is.na(col)]

if (length(obs) < 2) {
warning("Sample ", expr_cols[j], ": <2 observed values. Skipping imputation for this column.")
next
}

mu  <- mean(obs)
sig <- sd(obs)
if (!is.finite(sig) || sig == 0) {
warning("Sample ", expr_cols[j], ": invalid SD. Skipping imputation for this column.")
next
}

mu_imp <- mu - shift * sig
sd_imp <- width * sig

miss_idx <- which(is.na(col))
if (length(miss_idx) > 0) {
col[miss_idx] <- rnorm(length(miss_idx), mean = mu_imp, sd = sd_imp)
X[, j] <- col
imputed_n <- imputed_n + length(miss_idx)
}

param_table$Mean[j]     <- mu
param_table$SD[j]       <- sig
param_table$Imp_Mean[j] <- mu_imp
param_table$Imp_SD[j]   <- sd_imp
}

cat("\n🧩 Total imputed values:", imputed_n, "\n")

df_imputed <- df
df_imputed[, expr_cols] <- X

base_name <- "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_imputed_fixedseed"
out_txt    <- file.path(wd, paste0(base_name, ".txt"))
out_xlsx   <- file.path(wd, paste0(base_name, ".xlsx"))
out_params <- file.path(wd, "Step9_Imputation_Parameters.xlsx")

write.table(df_imputed, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(df_imputed, out_xlsx, overwrite = TRUE)
write.xlsx(param_table, out_params, overwrite = TRUE)

cat("\n💾 Saved imputed files:\n")
cat("  -", out_txt, "\n")
cat("  -", out_xlsx, "\n")
cat("  -", out_params, " (imputation params per sample)\n")
cat("\n🎉 Step 9 complete — Perseus-like per-column imputation finished.\n")

################################################################################
############################## 1️⃣1️⃣ QUANTILE NORMALIZATION ####################
################################################################################

################################################################################
############################## 🔟 QUANTILE NORMALIZATION (SAFER) ################
################################################################################

library(openxlsx)
library(ggplot2)
library(reshape2)

set.seed(1)

quantile_normalize <- function(X) {
all_na <- apply(X, 1, function(x) all(is.na(x)))
X2 <- X[!all_na, , drop = FALSE]
if (nrow(X2) == 0) stop("ERROR: all rows are NA. Cannot quantile-normalize.")

ranks <- apply(X2, 2, function(col) rank(col, ties.method = "average", na.last = "keep"))
sorted_X <- apply(X2, 2, function(col) sort(col, na.last = TRUE))
mean_sorted <- rowMeans(sorted_X, na.rm = TRUE)

X_norm <- X2
for (j in seq_len(ncol(X2))) {
rj <- ranks[, j]
ok <- !is.na(rj)
X_norm[ok, j] <- mean_sorted[rj[ok]]
X_norm[!ok, j] <- NA_real_
}

if (any(all_na)) {
X_full <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
X_full[!all_na, ] <- X_norm
X_norm <- X_full
}

X_norm
}

wd <- getwd()

infile <- file.path(wd, "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_imputed_fixedseed.txt")
if (!file.exists(infile)) stop("❌ ERROR: Imputed file from Step 9 not found.")

cat("\n=== Step 10: Quantile normalization (global) ===\n")
cat("Using file:", infile, "\n")

imputed_data <- read.table(infile, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

expr_cols <- grep("^LFQ_", colnames(imputed_data), value = TRUE)
if (length(expr_cols) == 0) stop("❌ ERROR: No LFQ_ columns found in the imputed dataset.")
cat("Found", length(expr_cols), "LFQ columns.\n")

for (cn in expr_cols) {
imputed_data[[cn]] <- suppressWarnings(as.numeric(gsub(",", "", as.character(imputed_data[[cn]]))))
imputed_data[[cn]][!is.finite(imputed_data[[cn]])] <- NA_real_
}

X <- as.matrix(imputed_data[, expr_cols, drop = FALSE])
storage.mode(X) <- "double"

if (!all(is.finite(X) | is.na(X))) {
bad <- sum(!is.finite(X) & !is.na(X))
stop(paste("ERROR:", bad, "invalid numeric values detected in expression matrix."))
}

# QC BEFORE
df_before <- melt(imputed_data[, expr_cols, drop = FALSE])
colnames(df_before) <- c("Sample", "Expression")

g_before <- ggplot(df_before, aes(x = Sample, y = Expression)) +
geom_boxplot(outlier.size = 0.5) +
theme_bw(base_size = 12) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
ggtitle("Before Quantile Normalization – Imputed Data")

ggsave(file.path(wd, "Step10_Boxplot_BEFORE.png"), g_before, width = 10, height = 6, dpi = 300)
cat("Saved: Step10_Boxplot_BEFORE.png\n")

# Normalize
cat("Performing quantile normalization across", ncol(X), "samples...\n")
X_norm <- quantile_normalize(X)

quantile_norm_data <- imputed_data
quantile_norm_data[, expr_cols] <- X_norm

# Flat proteins check
row_identical <- apply(X_norm, 1, function(r) {
v <- r[!is.na(r)]
if (length(v) <= 1) return(TRUE)
all(abs(v - v[1]) < 1e-12)
})
num_identical <- sum(row_identical)

if (num_identical > 0) {
cat("\nWARNING:", num_identical,
"proteins have IDENTICAL values across all samples after normalization.\n")

flat_file <- file.path(wd, "Step10_Flat_Proteins.txt")
protein_col <- if ("PG.ProteinGroups" %in% colnames(imputed_data)) {
imputed_data$PG.ProteinGroups
} else if ("Protein.IDs" %in% colnames(imputed_data)) {
imputed_data$Protein.IDs
} else {
seq_len(nrow(imputed_data))
}

flat_df <- data.frame(
RowIndex = which(row_identical),
Protein  = protein_col[row_identical],
stringsAsFactors = FALSE
)

write.table(flat_df, flat_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Flat protein list saved to:", flat_file, "\n")
} else {
cat("\nNo flat proteins detected (all proteins show variability across samples).\n")
}

# QC AFTER
df_after <- melt(quantile_norm_data[, expr_cols, drop = FALSE])
colnames(df_after) <- c("Sample", "Expression")

g_after <- ggplot(df_after, aes(x = Sample, y = Expression)) +
geom_boxplot(outlier.size = 0.5) +
theme_bw(base_size = 12) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
ggtitle("After Quantile Normalization – Imputed Data")

ggsave(file.path(wd, "Step10_Boxplot_AFTER.png"), g_after, width = 10, height = 6, dpi = 300)
cat("Saved: Step10_Boxplot_AFTER.png\n")

# Save normalized output
base_in  <- sub("\\.txt$", "", basename(infile))
base_out <- sub("_imputed_fixedseed$", "", base_in)

out_txt  <- file.path(wd, paste0(base_out, "_quantile_normalized.txt"))
out_xlsx <- file.path(wd, paste0(base_out, "_quantile_normalized.xlsx"))

write.table(quantile_norm_data, out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(quantile_norm_data, out_xlsx, overwrite = TRUE)

cat("\nSaved normalized files:\n  -", out_txt, "\n  -", out_xlsx, "\n")
cat("\n✅ Step 10 complete — quantile normalization + QC finished.\n")

################################################################################
############################## 1️⃣1️⃣ GLOBAL sPLS-DA ###########################
##############################      (NO LEGEND)     ###########################
################################################################################

library(ggplot2)
library(ggrepel)
library(mixOmics)
library(ellipse)
library(dplyr)

set.seed(1)

wd <- getwd()

# ✅ FIXED INPUT FILE (from Step 10)
norm_file <- file.path(
  wd,
  "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(norm_file))
  stop("❌ Normalized file from Step 10 not found! Expected:\n", norm_file)

cat("📄 Using normalized file:\n", norm_file, "\n\n")

df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Extract LFQ columns
expr_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(expr_cols) < 4) stop("❌ Not enough LFQ samples for sPLS-DA!")

cat("✔ Detected", length(expr_cols), "samples.\n")

# Group labels
group_labels <- gsub("^LFQ_|_[0-9]+$", "", expr_cols)
Y <- factor(group_labels)

cat("\n🧩 Samples per group:\n")
print(table(Y))

if (length(unique(Y)) < 2)
  stop("❌ Only one group present! sPLS-DA requires ≥2 groups.")

# Ensure numeric LFQs
df[, expr_cols] <- lapply(df[, expr_cols, drop = FALSE], function(x)
  suppressWarnings(as.numeric(as.character(x)))
)

# Data matrix: samples x proteins
X <- as.data.frame(t(df[, expr_cols, drop = FALSE]))

# Remove zero-variance proteins
X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
if (ncol(X) < 2) stop("❌ Too few variable proteins after filtering zero-variance.")

# Run sPLS-DA
ncomp <- 2
keepX <- c(100, 100)

cat("\n⚙ Running sPLS-DA with keepX =", paste(keepX, collapse = ", "), "\n")
splsda_model <- splsda(X, Y, ncomp = ncomp, keepX = keepX, scale = TRUE)

# Explained variance (safe)
expl_var <- NULL
if (!is.null(splsda_model$prop_expl_var$X)) {
  expl_var <- splsda_model$prop_expl_var$X[1:ncomp]
}

# Output directory
out_dir <- file.path(wd, "sPLSDA_Global")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comp_name <- "Global_sPLSDA"

# keepX bar chart
keep_df <- data.frame(
  Component = paste0("Comp", seq_along(keepX)),
  keepX = keepX
)

max_y <- max(150, max(keepX) + 10)

p_keep <- ggplot(keep_df, aes(Component, keepX)) +
  geom_bar(stat = "identity", fill = "#0065bd", width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 7) +
  labs(
    title = "sPLS-DA Variable Selection (keepX)",
    x = "Component", y = "Variables Selected"
  ) +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal(base_size = 22) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

ggsave(file.path(out_dir, paste0(comp_name, "_keepX_plot.png")),
       p_keep, dpi = 300, width = 7, height = 6, bg = "white")
ggsave(file.path(out_dir, paste0(comp_name, "_keepX_plot.pdf")),
       p_keep, dpi = 300, width = 7, height = 6, bg = "white")

# Scores + centroids
scores <- data.frame(
  Comp1 = splsda_model$variates$X[, 1],
  Comp2 = splsda_model$variates$X[, 2],
  Group = Y
)

centroids <- scores %>%
  group_by(Group) %>%
  summarise(
    Comp1 = mean(Comp1),
    Comp2 = mean(Comp2),
    n = n(),
    .groups = "drop"
  )

# Ellipse builder
build_ellipses <- function(df, level = 0.9) {
  out <- lapply(split(df, df$Group), function(d) {
    if (nrow(d) < 3) return(NULL)
    
    S <- cov(d[, c("Comp1", "Comp2")])
    if (det(S) <= 1e-12) S <- S + diag(2) * 1e-6
    
    e <- ellipse(
      S,
      centre = colMeans(d[, c("Comp1", "Comp2")]),
      level = level,
      npoints = 200
    )
    
    e <- as.data.frame(e)
    colnames(e) <- c("x", "y")
    e$Group <- d$Group[1]
    e
  })
  
  bind_rows(Filter(Negate(is.null), out))
}

ellipse_data <- build_ellipses(scores)

# Colors
palette <- c("#4575b4", "#fdae61", "#d73027", "#66c2a5", "#abdda4")
colors <- setNames(palette[seq_along(unique(Y))], unique(Y))

# ✅ FIXED if/else blocks (prevents "unexpected else")
xlab <- if (!is.null(expl_var)) {
  paste0("Component 1 (", round(expl_var[1] * 100, 1), "%)")
} else {
  "Component 1"
}

ylab <- if (!is.null(expl_var)) {
  paste0("Component 2 (", round(expl_var[2] * 100, 1), "%)")
} else {
  "Component 2"
}

# ===========================
# FINAL sPLS-DA PLOT (NO LEGEND)
# ===========================
p <- ggplot(scores, aes(Comp1, Comp2, color = Group, fill = Group))

if (nrow(ellipse_data) > 0) {
  p <- p +
    geom_polygon(
      data = ellipse_data,
      aes(x = x, y = y, fill = Group),
      alpha = 0.25, color = NA
    ) +
    geom_path(
      data = ellipse_data,
      aes(x = x, y = y, color = Group),
      linewidth = 1.1
    )
} else {
  p <- p + stat_ellipse(type = "t", level = 0.9)
}

p <- p +
  geom_point(size = 6, alpha = 0.9) +
  geom_text_repel(
    data = centroids,
    aes(label = paste0(Group, "\n(n=", n, ")")),
    color = "black", size = 7, fontface = "bold"
  ) +
  labs(title = NULL, x = xlab, y = ylab) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(color = "none", fill = "none") +
  theme_minimal(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "Global_sPLSDA_plot.png"),
       p, dpi = 300, width = 10, height = 10)
ggsave(file.path(out_dir, "Global_sPLSDA_plot.pdf"),
       p, dpi = 300, width = 10, height = 10)

cat("\n🎉 GLOBAL sPLS-DA completed (legend removed).\n🔥 Output saved to:\n", out_dir, "\n")
################################################################################

################################################################################
########################### 1️⃣1️⃣ sPLS-DA (Pre_RT vs Post_RT only) ############
###########################         (NO LEGEND + NO ELSE ERROR) ################
################################################################################
# ✅ This script:
#   - Loads the Step10 quantile-normalized file
#   - Keeps ONLY Pre_RT and Post_RT samples (Healthy removed)
#   - Runs sPLS-DA (2 components by default)
#   - Saves ALL outputs into a NEW folder
# ✅ Fixes:
#   - Removes legend properly
#   - Fixes "unexpected else" by using { } in if/else blocks
#   - Uses legend.position = "none" (NOT NULL)

library(ggplot2)
library(ggrepel)
library(mixOmics)
library(ellipse)
library(dplyr)
library(openxlsx)

set.seed(1)

wd <- getwd()

# --- INPUT FILE from Step10 ---
norm_file <- file.path(
  wd,
  "DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(norm_file)) {
  stop("❌ Normalized file from Step 10 not found! Expected:\n", norm_file)
}

cat("📄 Using normalized file:\n", norm_file, "\n\n")

df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Identify LFQ columns ---
all_expr_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(all_expr_cols) < 4) stop("❌ Not enough LFQ samples for sPLS-DA!")

# --- Extract group labels from column names ---
all_group_labels <- gsub("^LFQ_|_[0-9]+$", "", all_expr_cols)

# --- Keep ONLY Pre_RT and Post_RT (remove Healthy and anything else) ---
keep_groups <- c("Pre_RT", "Post_RT")

keep_mask <- all_group_labels %in% keep_groups
expr_cols <- all_expr_cols[keep_mask]
group_labels <- all_group_labels[keep_mask]

cat("✔ Total LFQ samples detected:", length(all_expr_cols), "\n")
cat("✔ Samples kept for sPLS-DA:", length(expr_cols), "\n")

if (length(expr_cols) < 4) {
  stop("❌ After filtering to Pre_RT/Post_RT there are <4 samples. Check your LFQ column names.")
}

# --- Define Y with fixed order ---
Y <- factor(group_labels, levels = keep_groups)

cat("\n🧩 Samples per group (after filtering):\n")
print(table(Y))

if (length(unique(Y)) < 2) stop("❌ Only one group present after filtering! Need Pre_RT and Post_RT.")

# --- Ensure numeric LFQs (robust coercion) ---
df[, expr_cols] <- lapply(df[, expr_cols, drop = FALSE], function(x) {
  x_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  x_num[!is.finite(x_num)] <- NA_real_
  x_num
})

# --- Build X matrix: samples x proteins ---
# df is proteins x samples -> transpose to samples x proteins
X <- as.data.frame(t(df[, expr_cols, drop = FALSE]))

# --- Remove proteins with zero variance (important for sPLS-DA) ---
sd_vec <- apply(X, 2, sd, na.rm = TRUE)
X <- X[, is.finite(sd_vec) & sd_vec > 0, drop = FALSE]

if (ncol(X) < 2) stop("❌ Too few variable proteins after removing zero-variance features.")

cat("\n✔ Matrix used for sPLS-DA:\n")
cat("   Samples:", nrow(X), "\n")
cat("   Proteins:", ncol(X), "\n\n")

################################################################################
# --- OUTPUT FOLDER (NEW) ---
################################################################################
out_dir <- file.path(wd, "sPLSDA_PreRT_PostRT_only")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("📁 Output folder:\n", out_dir, "\n\n")

################################################################################
# --- (Optional) keepX tuning grid ---
################################################################################
do_tune <- FALSE

ncomp <- 2
keepX <- c(100, 100)

if (do_tune) {
  cat("🔎 Tuning keepX (this can take some time)...\n")
  
  test_keepX <- c(10, 25, 50, 75, 100, 150, 200)
  tune <- tune.splsda(
    X, Y,
    ncomp = ncomp,
    test.keepX = replicate(ncomp, test_keepX, simplify = FALSE),
    validation = "Mfold",
    folds = 5,
    nrepeat = 50,
    dist = "max.dist",
    measure = "BER",
    progressBar = TRUE
  )
  
  keepX <- tune$choice.keepX
  cat("✅ Best keepX chosen:\n")
  print(keepX)
}

################################################################################
# --- Run sPLS-DA ---
################################################################################
cat("\n⚙ Running sPLS-DA with ncomp =", ncomp, "and keepX =", paste(keepX, collapse = ", "), "\n")

splsda_model <- splsda(X, Y, ncomp = ncomp, keepX = keepX, scale = TRUE)

# Explained variance (safe)
expl_var <- NULL
if (!is.null(splsda_model$prop_expl_var$X)) {
  expl_var <- splsda_model$prop_expl_var$X[1:ncomp]
}

################################################################################
# --- Save model + key tables ---
################################################################################
saveRDS(splsda_model, file.path(out_dir, "splsda_model_PreRT_vs_PostRT.rds"))

# Selected variables per component
sel_vars <- lapply(1:ncomp, function(k) selectVar(splsda_model, comp = k)$name)
names(sel_vars) <- paste0("Comp", 1:ncomp)

# Write selected features to Excel
sel_df <- data.frame(
  Component = rep(names(sel_vars), times = sapply(sel_vars, length)),
  Protein   = unlist(sel_vars, use.names = FALSE),
  row.names = NULL
)
write.xlsx(sel_df, file.path(out_dir, "Selected_Proteins_by_Component.xlsx"), overwrite = TRUE)

################################################################################
# --- keepX bar plot (NO LEGEND) ---
################################################################################
keep_df <- data.frame(
  Component = paste0("Comp", seq_along(keepX)),
  keepX = keepX
)

max_y <- max(150, max(keepX) + 10)

p_keep <- ggplot(keep_df, aes(Component, keepX)) +
  geom_bar(stat = "identity", fill = "#0065bd", width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 7) +
  labs(
    title = "sPLS-DA Variable Selection (keepX)\nPre_RT vs Post_RT only",
    x = "Component", y = "Variables Selected"
  ) +
  scale_y_continuous(limits = c(0, max_y)) +
  theme_minimal(base_size = 22) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "PrePost_sPLSDA_keepX_plot.png"),
       p_keep, dpi = 300, width = 8, height = 6, bg = "white")
ggsave(file.path(out_dir, "PrePost_sPLSDA_keepX_plot.pdf"),
       p_keep, dpi = 300, width = 8, height = 6, bg = "white")

################################################################################
# --- Scores plot + centroids + ellipses (NO LEGEND) ---
################################################################################
scores <- data.frame(
  Comp1 = splsda_model$variates$X[, 1],
  Comp2 = splsda_model$variates$X[, 2],
  Group = Y
)

centroids <- scores %>%
  group_by(Group) %>%
  summarise(
    Comp1 = mean(Comp1),
    Comp2 = mean(Comp2),
    n = n(),
    .groups = "drop"
  )

build_ellipses <- function(df, level = 0.9) {
  out <- lapply(split(df, df$Group), function(d) {
    if (nrow(d) < 3) return(NULL)
    
    S <- cov(d[, c("Comp1", "Comp2")])
    if (det(S) <= 1e-12) S <- S + diag(2) * 1e-6
    
    e <- ellipse(
      S,
      centre = colMeans(d[, c("Comp1", "Comp2")]),
      level = level, npoints = 200
    )
    
    e <- as.data.frame(e)
    colnames(e) <- c("x", "y")
    e$Group <- d$Group[1]
    e
  })
  
  bind_rows(Filter(Negate(is.null), out))
}

ellipse_data <- build_ellipses(scores)

# Colors (fixed for 2 groups)
colors <- c("Pre_RT" = "#4575b4", "Post_RT" = "#d73027")

# ✅ FIXED if/else blocks to prevent "unexpected else"
xlab <- if (!is.null(expl_var)) {
  paste0("Component 1 (", round(expl_var[1] * 100, 1), "%)")
} else {
  "Component 1"
}

ylab <- if (!is.null(expl_var)) {
  paste0("Component 2 (", round(expl_var[2] * 100, 1), "%)")
} else {
  "Component 2"
}

p <- ggplot(scores, aes(Comp1, Comp2, color = Group, fill = Group))

if (nrow(ellipse_data) > 0) {
  p <- p +
    geom_polygon(
      data = ellipse_data,
      aes(x = x, y = y, fill = Group),
      alpha = 0.25, color = NA
    ) +
    geom_path(
      data = ellipse_data,
      aes(x = x, y = y, color = Group),
      linewidth = 1.1
    )
} else {
  p <- p + stat_ellipse(type = "t", level = 0.9)
}

p <- p +
  geom_point(size = 6, alpha = 0.9) +
  geom_text_repel(
    data = centroids,
    aes(label = paste0(Group, "\n(n=", n, ")")),
    color = "black", size = 7, fontface = "bold"
  ) +
  labs(
    title = NULL,
    x = xlab, y = ylab
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(color = "none", fill = "none") +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "PrePost_sPLSDA_plot.png"),
       p, dpi = 300, width = 10, height = 10, bg = "white")
ggsave(file.path(out_dir, "PrePost_sPLSDA_plot.pdf"),
       p, dpi = 300, width = 10, height = 10, bg = "white")

################################################################################
# --- Performance (CV) summary (recommended QC) ---
################################################################################
perf_res <- perf(
  splsda_model,
  validation = "Mfold",
  folds = 5,
  nrepeat = 50,
  progressBar = TRUE,
  auc = TRUE
)

saveRDS(perf_res, file.path(out_dir, "perf_PrePost_sPLSDA.rds"))

# --- Save error rates table (ROBUST across mixOmics versions) ---
err <- perf_res$error.rate

# pick distance name safely
dist_name <- if ("max.dist" %in% names(err)) "max.dist" else names(err)[1]
err_mat <- err[[dist_name]]

cat("\nUsing distance:", dist_name, "\n")
cat("Error-rate matrix columns:\n")
print(colnames(err_mat))

# Convert to data.frame safely
err_df <- as.data.frame(err_mat)
err_df$ncomp <- seq_len(nrow(err_mat))
err_df <- err_df[, c("ncomp", setdiff(colnames(err_df), "ncomp")), drop = FALSE]

# Standardize a couple common column names if present (optional)
if ("Overall" %in% colnames(err_df) && !"overall" %in% colnames(err_df)) {
  colnames(err_df)[colnames(err_df) == "Overall"] <- "overall"
}

# Keep only the most useful columns if they exist
want_cols <- c("ncomp", "overall", "BER")
have_cols <- intersect(want_cols, colnames(err_df))

if (length(have_cols) >= 2) {
  err_df_out <- err_df[, have_cols, drop = FALSE]
} else {
  # fallback: save everything available (often includes per-class error rates)
  err_df_out <- err_df
}

write.xlsx(err_df_out, file.path(out_dir, "CV_ErrorRates.xlsx"), overwrite = TRUE)
cat("✅ Saved CV error rates to CV_ErrorRates.xlsx\n")

cat("\n🎉 sPLS-DA Pre_RT vs Post_RT completed (legend removed, no else error).\n🔥 Output saved to:\n", out_dir, "\n")
################################################################################



################################################################################
############################## 1️⃣2️⃣ GLOBAL PCA ##############################
################################################################################

library(mixOmics)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ellipse)

set.seed(1)

wd <- getwd()

# ✅ FIXED INPUT FILE (from Step 10)
norm_file <- file.path(
wd,
"DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(norm_file))
stop("❌ Normalized file from Step 10 not found! Expected:\n", norm_file)

cat("📄 Using normalized file:\n", norm_file, "\n\n")

df <- read.table(norm_file, header = TRUE, sep = "\t",
quote = "", check.names = FALSE)

# Identify LFQ columns
sample_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(sample_cols) < 3)
stop("❌ Not enough samples for PCA!")

# Safe numeric coercion
df[, sample_cols] <- lapply(df[, sample_cols, drop = FALSE], function(x) {
x <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
x[!is.finite(x)] <- NA_real_
x
})

# Groups
Y <- factor(gsub("^LFQ_|_[0-9]+$", "", sample_cols))

cat("\n🧩 PCA — samples per group:\n")
print(table(Y))

# Matrix: samples × proteins
X <- t(as.matrix(df[, sample_cols, drop = FALSE]))
storage.mode(X) <- "double"

# Remove zero-variance / all-NA proteins
sds <- apply(X, 2, sd, na.rm = TRUE)
keep <- is.finite(sds) & sds > 0
X <- X[, keep, drop = FALSE]

if (ncol(X) < 2)
stop("❌ Not enough variable proteins for PCA!")

# PCA
pca_res <- mixOmics::pca(
X,
ncomp = 10,
center = TRUE,
scale = TRUE
)

expl_var <- pca_res$prop_expl_var$X
cum_var <- cumsum(expl_var)

# Scree plot
scree_df <- data.frame(
PC = factor(paste0("PC", seq_along(expl_var)),
levels = paste0("PC", seq_along(expl_var))),
Variance = expl_var * 100
)

p_scree <- ggplot(scree_df, aes(PC, Variance)) +
geom_bar(stat = "identity", fill = "#1565c0") +
geom_text(aes(label = sprintf("%.1f%%", Variance)),
vjust = -0.5, size = 6) +
labs(title = "Global PCA Scree Plot",
x = "Principal Component",
y = "Explained Variance (%)") +
theme_minimal(base_size = 24) +
theme(panel.border = element_rect(color = "black", fill = NA))

# PC1 vs PC2
scores <- data.frame(
PC1 = pca_res$variates$X[, 1],
PC2 = pca_res$variates$X[, 2],
Group = Y
)

centroids <- scores |>
group_by(Group) |>
summarise(
PC1 = mean(PC1),
PC2 = mean(PC2),
n = n(),
.groups = "drop"
)

build_ellipses <- function(df, level = 0.9) {
out <- lapply(split(df, df$Group), function(d) {
if (nrow(d) < 3) return(NULL)
S <- cov(d[, c("PC1", "PC2")])
if (!is.finite(det(S)) || det(S) <= 1e-12)
S <- S + diag(2) * 1e-6
e <- ellipse(S, centre = colMeans(d[, c("PC1", "PC2")]),
level = level, npoints = 200)
e <- as.data.frame(e)
colnames(e) <- c("x", "y")
e$Group <- d$Group[1]
e
})
out <- Filter(Negate(is.null), out)
if (length(out) == 0) return(NULL)
bind_rows(out)
}

ellipse_data <- build_ellipses(scores)

default_colors <- c(
"Healthy" = "#f9d739",
"Pre_RT"  = "#1565c0",
"Post_RT" = "#d73027"
)
colors <- default_colors[names(default_colors) %in% levels(Y)]

xlab <- paste0("PC1 (", round(expl_var[1] * 100, 1), "%)")
ylab <- paste0("PC2 (", round(expl_var[2] * 100, 1), "%)")

p_pc12 <- ggplot(scores, aes(PC1, PC2, color = Group, fill = Group)) +
geom_point(size = 5) +
geom_text_repel(
data = centroids,
aes(label = paste0(Group, "\n(n=", n, ")")),
size = 6, fontface = "bold"
) +
scale_color_manual(values = colors, drop = FALSE) +
scale_fill_manual(values = colors, drop = FALSE) +
labs(title = "Global PCA", x = xlab, y = ylab) +
coord_equal() +
theme_minimal(base_size = 26) +
theme(panel.border = element_rect(color = "black", fill = NA))

if (!is.null(ellipse_data)) {
p_pc12 <- p_pc12 +
geom_polygon(data = ellipse_data,
aes(x = x, y = y, fill = Group),
alpha = 0.25, color = NA) +
geom_path(data = ellipse_data,
aes(x = x, y = y, color = Group),
linewidth = 1.2)
}

out_dir <- file.path(wd, "PCA_results")
dir.create(out_dir, showWarnings = FALSE)

ggsave(file.path(out_dir, "Global_PCA_scree.png"),
p_scree, dpi = 300, width = 10, height = 7, bg = "white")
ggsave(file.path(out_dir, "Global_PCA_PC1_PC2.png"),
p_pc12, dpi = 300, width = 10, height = 9, bg = "white")

ggsave(file.path(out_dir, "Global_PCA_scree.pdf"),
p_scree, dpi = 300, width = 10, height = 7, bg = "white", device = "pdf")
ggsave(file.path(out_dir, "Global_PCA_PC1_PC2.pdf"),
p_pc12, dpi = 300, width = 10, height = 9, bg = "white", device = "pdf")

cat("\n🎉 Step 12 complete — Global PCA finished.\n")

################################################################################
############################## 🔵 STEP 13 — VOLCANO PLOTS #######################
################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(openxlsx)

# ========================== SETTINGS ==========================================
pval_thr  <- 0.05        # p-value threshold (raw p)
use_fdr   <- FALSE       # TRUE = use BH-adjusted p-values for threshold + y-axis
fc_thr    <- log10(1.5)  # threshold on log10 ratio (1.5x)
test_type <- "welch"     # welch / student / wilcoxon
n_label   <- 12          # number of top gene labels per side
# ==============================================================================

wd <- getwd()

# ✅ FIXED INPUT FILE (from Step 10)
norm_file <- file.path(
wd,
"DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(norm_file))
stop("❌ Normalized file not found in working directory! Expected:\n", norm_file)

cat("📄 Using normalized file:\n", norm_file, "\n\n")

df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

###############################################################################
# 1 — DETECT & SORT LFQ COLUMNS BY STAGE + REPLICATE NUMBER
###############################################################################

sample_cols <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(sample_cols) == 0) stop("❌ No LFQ_ columns found in normalized file!")

# Robust numeric LFQ conversion
df[, sample_cols] <- lapply(df[, sample_cols, drop = FALSE], function(x) {
x <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
x[!is.finite(x)] <- NA_real_
x
})

parse_col <- function(x) {
stage <- sub("^LFQ_(.*)_([0-9]+)$", "\\1", x)
num   <- as.numeric(sub(".*_([0-9]+)$", "\\1", x))
data.frame(col = x, stage = stage, num = num)
}

meta <- do.call(rbind, lapply(sample_cols, parse_col))

stage_order <- c("Healthy", "Pre_RT", "Post_RT")
meta$stage <- factor(meta$stage, levels = stage_order)

# Unknown stages go last (won't break)
meta$stage2 <- as.character(meta$stage)
meta$stage2[is.na(meta$stage)] <- meta$stage[is.na(meta$stage)] <- "ZZZ_OTHER"

meta <- meta[order(meta$stage, meta$num), ]
sample_cols <- meta$col

###############################################################################
# 2 — STAGE LABELS
###############################################################################

clean_stage <- function(x) {
x <- gsub("^LFQ_", "", x)
x <- sub("_[0-9]+$", "", x)

if (grepl("^Pre_RT$", x, ignore.case = TRUE)) return("Pre_RT")
if (grepl("^Post_RT$", x, ignore.case = TRUE)) return("Post_RT")
if (grepl("^Healthy$", x, ignore.case = TRUE)) return("Healthy")
return(x)
}

stage_labels <- sapply(sample_cols, clean_stage)

cat("✔ Detected stages:\n")
print(table(stage_labels))

comparisons <- list(
Pre_RT_vs_Healthy    = c("Pre_RT", "Healthy"),
Healthy_vs_Post_RT   = c("Healthy", "Post_RT"),
Pre_RT_vs_Post_RT    = c("Pre_RT", "Post_RT")
)

###############################################################################
# MAIN LOOP
###############################################################################

for (comp_name in names(comparisons)) {

groups <- comparisons[[comp_name]]
g1 <- groups[1]
g2 <- groups[2]

cat("\n\n🚀 Processing:", comp_name, "|", g1, "vs", g2, "\n")

out_dir <- file.path(wd, comp_name)
dir.create(out_dir, showWarnings = FALSE)

g1_cols <- sample_cols[stage_labels == g1]
g2_cols <- sample_cols[stage_labels == g2]

if (length(g1_cols) < 2 || length(g2_cols) < 2) {
cat("⚠️ Not enough samples for:", comp_name, "\n")
next
}

gene_names <- if ("PG.Genes" %in% colnames(df)) df$PG.Genes else df$PG.ProteinGroups

results <- data.frame(
Gene = gene_names,
log10FC = NA_real_,
p_value = NA_real_,
stringsAsFactors = FALSE
)

# Per-protein stats
for (i in seq_len(nrow(df))) {

v1 <- as.numeric(df[i, g1_cols])
v2 <- as.numeric(df[i, g2_cols])

if (sum(!is.na(v1)) >= 2 && sum(!is.na(v2)) >= 2) {

pval <- tryCatch({
if (tolower(test_type) == "welch") {
t.test(v1, v2, var.equal = FALSE)$p.value
} else if (tolower(test_type) == "student") {
t.test(v1, v2, var.equal = TRUE)$p.value
} else if (tolower(test_type) == "wilcoxon") {
wilcox.test(v1, v2)$p.value
} else {
stop("Unknown test_type: ", test_type)
}
}, error = function(e) NA_real_)

if (is.finite(pval)) {
results$log10FC[i] <- mean(v2, na.rm = TRUE) - mean(v1, na.rm = TRUE)
results$p_value[i] <- pval
}
}
}

results <- results[is.finite(results$p_value) & is.finite(results$log10FC), , drop = FALSE]

if (nrow(results) == 0) {
cat("⚠️ No valid proteins for stats in:", comp_name, "\n")
next
}

# Multiple testing (BH)
results$p_adj <- p.adjust(results$p_value, method = "BH")

# Choose p column for threshold + y-axis
p_used <- if (use_fdr) results$p_adj else results$p_value

# Significance labels
results$negLog10P <- -log10(p_used)
results$group <- "Non-significant"

results$group[results$log10FC >  fc_thr & p_used < pval_thr] <- paste("Higher_in", g2, sep = "_")
results$group[results$log10FC < -fc_thr & p_used < pval_thr] <- paste("Higher_in", g1, sep = "_")

# Save tables
write.xlsx(results, file.path(out_dir, paste0("Volcano_data_", comp_name, ".xlsx")), overwrite = TRUE)
write.table(results, file.path(out_dir, paste0("Volcano_data_", comp_name, ".txt")),
sep = "\t", quote = FALSE, row.names = FALSE)

higher_g2 <- results %>% filter(group == paste("Higher_in", g2, sep = "_"))
higher_g1 <- results %>% filter(group == paste("Higher_in", g1, sep = "_"))

write.table(higher_g2, file.path(out_dir, paste0("Higher_in_", g2, ".txt")),
sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(higher_g2, file.path(out_dir, paste0("Higher_in_", g2, ".xlsx")), overwrite = TRUE)

write.table(higher_g1, file.path(out_dir, paste0("Higher_in_", g1, ".txt")),
sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(higher_g1, file.path(out_dir, paste0("Higher_in_", g1, ".xlsx")), overwrite = TRUE)

# Plot labels
top_up   <- higher_g2 %>% arrange(p_value) %>% head(n_label)
top_down <- higher_g1 %>% arrange(p_value) %>% head(n_label)

colors <- setNames(
c("#D73027", "#0065BD", "black"),
c(paste("Higher_in", g2, sep="_"), paste("Higher_in", g1, sep="_"), "Non-significant")
)

y_thr <- -log10(pval_thr)

volcano_plot <- ggplot(results, aes(x = log10FC, y = negLog10P)) +
geom_point(aes(color = group), size = 1, alpha = 1) +
scale_color_manual(values = colors) +
geom_hline(yintercept = y_thr, linetype = "dashed") +
geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed") +
geom_text_repel(data = top_up, aes(label = Gene),
color = "#D73027", size = 10) +
geom_text_repel(data = top_down, aes(label = Gene),
color = "#0065BD", size = 10) +
labs(
title = NULL,
x = paste0("Log10 Fold Change (", g2, " – ", g1, ")"),
y = if (use_fdr)
expression(-Log[10]~FDR)
else
expression(-Log[10]~p~value)
) +
theme_minimal(base_size = 28) +
theme(
plot.title   = element_text(size = 36, face = "plain", hjust = 0.5),
axis.title.x = element_text(size = 32, face = "plain", margin = margin(t = 20)),
axis.title.y = element_text(size = 32, face = "plain", margin = margin(r = 20)),
axis.text.x  = element_text(size = 24, face = "plain"),
axis.text.y  = element_text(size = 24, face = "plain"),
legend.title = element_text(size = 28, face = "plain"),
legend.text  = element_text(size = 24, face = "plain"),
legend.key.size = unit(2, "lines"),
legend.position = "bottom",
plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
)


ggsave(file.path(out_dir, paste0("Volcano_", comp_name, ".png")),
volcano_plot, dpi = 300, width = 14, height = 10, bg = "white")

ggsave(file.path(out_dir, paste0("Volcano_", comp_name, ".pdf")),
volcano_plot, dpi = 300, width = 14, height = 10, device = "pdf")

cat("✅ Volcano plot saved for:", comp_name, "\n")
}

cat("\n🎉 Step 13 complete — Volcano + Higher_in_* files saved.\n")


################################################################################
############################## 🔵 STEP 14 — ANNOTATION (FIXED) ##################
################################################################################
# ✅ Uses global normalized file produced by Step 10:
#   DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt
#
# ✅ Reads Higher_in_*.txt inside:
#   Pre_RT_vs_Healthy
#   Healthy_vs_Post_RT
#   Pre_RT_vs_Post_RT
#
# ✅ Writes annotated outputs back into the same comparison folder.
################################################################################

suppressPackageStartupMessages({
library(dplyr); library(tidyr); library(stringr)
library(AnnotationDbi); library(org.Hs.eg.db)
library(GO.db); library(ReactomePA); library(clusterProfiler)
library(KEGGREST); library(biomaRt)
library(openxlsx)
})

wd <- getwd()
cat("🏠 Working directory:\n", wd, "\n")

# ✅ FIXED normalized file name (matches Step 10 / Steps 11–13)
norm_file <- file.path(
wd,
"DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(norm_file)) {
stop("❌ Cannot find normalized file in working directory:\n", norm_file)
}

cat("✔ Using global normalized file:\n", norm_file, "\n")

norm_df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# ✅ Robust: allow either PG.Genes or PG.ProteinGroups / Protein.IDs as key
key_col <- NULL
if ("PG.Genes" %in% colnames(norm_df)) {
key_col <- "PG.Genes"
} else if ("PG.ProteinGroups" %in% colnames(norm_df)) {
key_col <- "PG.ProteinGroups"
} else if ("Protein.IDs" %in% colnames(norm_df)) {
key_col <- "Protein.IDs"
} else {
stop("❌ None of these ID columns exist in normalized file: PG.Genes, PG.ProteinGroups, Protein.IDs")
}

cat("🔑 Using ID column for merge:", key_col, "\n")

comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")
folders <- file.path(wd, comparisons)
folders <- folders[dir.exists(folders)]

if (length(folders) == 0) {
stop("❌ No comparison folders found in working directory.")
}

cat("📂 Found comparison folders:\n")
print(folders)

# ---------------- Helper functions ----------------

split_genes <- function(gene_list) {
data.frame(Original = gene_list, stringsAsFactors = FALSE) %>%
mutate(Gene_Single = strsplit(Original, ";")) %>%
tidyr::unnest(Gene_Single) %>%
mutate(Gene_Single = trimws(Gene_Single)) %>%
filter(!is.na(Gene_Single) & Gene_Single != "")
}

map_to_entrez <- function(symbols) {
symbols <- unique(na.omit(symbols))
if (length(symbols) == 0) return(data.frame())

mapping <- tryCatch(
AnnotationDbi::select(
org.Hs.eg.db,
keys = symbols,
columns = c("SYMBOL", "ENTREZID", "UNIPROT"),
keytype = "SYMBOL"
),
error = function(e) data.frame()
)

mapping <- mapping[!duplicated(mapping$SYMBOL), ]

# Fallback: if "symbols" are actually UniProt IDs
if (nrow(mapping) == 0 || all(is.na(mapping$ENTREZID))) {
uni_map <- tryCatch(
clusterProfiler::bitr(symbols, fromType = "UNIPROT",
toType = c("SYMBOL", "ENTREZID"),
OrgDb = "org.Hs.eg.db"),
error = function(e) NULL
)
if (!is.null(uni_map) && nrow(uni_map) > 0) mapping <- uni_map
}

# Fallback: biomaRt for missing Entrez (may fail offline)
if (nrow(mapping) > 0 && any(is.na(mapping$ENTREZID))) {
miss <- mapping$SYMBOL[is.na(mapping$ENTREZID)]
mart <- tryCatch(useEnsembl("genes", "hsapiens_gene_ensembl"), error = function(e) NULL)

if (!is.null(mart) && length(miss) > 0) {
bm <- tryCatch(
getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
filters = "hgnc_symbol", values = miss, mart = mart),
error = function(e) NULL
)
if (!is.null(bm) && nrow(bm) > 0) {
for (i in seq_len(nrow(bm))) {
mapping$ENTREZID[mapping$SYMBOL == bm$hgnc_symbol[i]] <- as.character(bm$entrezgene_id[i])
}
}
}
}

mapping <- mapping[!is.na(mapping$ENTREZID) & mapping$ENTREZID != "", ]
unique(mapping)
}

safe_kegg <- function(entrez_ids, gene_entrez) {
entrez_ids <- as.character(entrez_ids)

k <- tryCatch(
enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "ncbi-geneid",
pvalueCutoff = 1, qvalueCutoff = 1),
error = function(e) NULL
)

if (!is.null(k) && nrow(as.data.frame(k)) > 0) {
dfk <- as.data.frame(k)
return(
dfk %>%
dplyr::select(KEGG_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep = "/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(KEGG_Pathway = paste(unique(KEGG_Pathway), collapse = "; "),
.groups = "drop")
)
}

# Fallback: KEGGREST
cat("📦 KEGGREST fallback\n")
path_list <- KEGGREST::keggLink("pathway", "hsa")
path_list <- path_list[gsub("hsa:", "", names(path_list)) %in% entrez_ids]
if (length(path_list) == 0) return(data.frame())

desc <- KEGGREST::keggList("pathway", "hsa")

df <- data.frame(
ENTREZID = sub("hsa:", "", names(path_list)),
PATHWAY  = sub("path:", "", path_list),
stringsAsFactors = FALSE
)

df <- merge(df,
data.frame(PATHWAY = sub("path:", "", names(desc)),
KEGG_Pathway = desc,
stringsAsFactors = FALSE),
by = "PATHWAY")

df %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(KEGG_Pathway = paste(unique(KEGG_Pathway), collapse = "; "),
.groups = "drop")
}

# ---------------- MAIN LOOP ----------------

for (folder in folders) {

cat("\n============================================================\n")
cat("📂 Processing comparison folder:", folder, "\n")
cat("============================================================\n")

higher_files <- list.files(folder, pattern = "^Higher_in_.*\\.txt$", full.names = TRUE)

if (!length(higher_files)) {
cat("⚠️ No Higher_in_*.txt in", folder, "→ skipping.\n")
next
}

cat("✔ Found", length(higher_files), "Higher_in files in", basename(folder), ":\n")
print(basename(higher_files))

for (hf in higher_files) {

grp <- gsub("^Higher_in_|\\.txt$", "", basename(hf))
cat("\n🔹 Annotating:", basename(folder), "— Higher_in_", grp, "\n", sep = "")

higher_df <- read.table(hf, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

if (!("Gene" %in% colnames(higher_df))) {
cat("⚠️ No 'Gene' column in", basename(hf), "→ skip.\n")
next
}

genes_original <- unique(na.omit(higher_df$Gene))
if (length(genes_original) == 0) {
cat("⚠️ No genes in", basename(hf), "→ skip.\n")
next
}

# Split "A;B;C" only for mapping
map <- split_genes(genes_original)
symbols <- unique(map$Gene_Single)

gene_entrez <- map_to_entrez(symbols)
entrez <- unique(gene_entrez$ENTREZID)

if (length(entrez) < 2) {
cat("⚠️ Too few Entrez IDs for", grp, "→ skip.\n")
next
}

# GO annotation
go <- tryCatch({
go_anno <- AnnotationDbi::select(
org.Hs.eg.db,
keys = entrez,
columns = c("GO", "ONTOLOGY"),
keytype = "ENTREZID"
)

go_terms <- AnnotationDbi::select(
GO.db,
keys = unique(go_anno$GO),
columns = "TERM",
keytype = "GOID"
)

gm <- merge(go_anno, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
gm <- dplyr::left_join(gm, gene_entrez, by = "ENTREZID")
gm$GO_Pathway <- paste(gm$GO, gm$ONTOLOGY, gm$TERM, sep = " | ")

gm %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(GO_Pathway = paste(unique(GO_Pathway), collapse = "; "),
.groups = "drop")
}, error = function(e) {
cat("⚠️ GO failed:", e$message, "\n")
data.frame()
})

# KEGG
kegg <- safe_kegg(entrez, gene_entrez)

# Reactome
rea <- tryCatch({
enr <- ReactomePA::enrichPathway(
gene = as.numeric(entrez),
organism = "human",
pvalueCutoff = 1, qvalueCutoff = 1
)
dfrea <- as.data.frame(enr)

dfrea %>%
dplyr::select(Reactome_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep = "/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by = "ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(Reactome_Pathway = paste(unique(Reactome_Pathway), collapse = "; "),
.groups = "drop")
}, error = function(e) {
cat("⚠️ Reactome failed:", e$message, "\n")
data.frame()
})

# Merge annotations back to ORIGINAL gene strings
ann <- map %>%
dplyr::left_join(go,   by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(kegg, by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::left_join(rea,  by = c("Gene_Single" = "SYMBOL")) %>%
dplyr::group_by(Original) %>%
dplyr::summarise(
GO_Pathway       = paste(na.omit(unique(GO_Pathway)), collapse = "; "),
KEGG_Pathway     = paste(na.omit(unique(KEGG_Pathway)), collapse = "; "),
Reactome_Pathway = paste(na.omit(unique(Reactome_Pathway)), collapse = "; "),
.groups = "drop"
)

# ✅ Subset normalized table using the SAME IDs stored in Higher_in_*.txt (Gene column)
sig_norm <- norm_df[norm_df[[key_col]] %in% genes_original, , drop = FALSE]

if (nrow(sig_norm) == 0) {
cat("⚠️ No matching rows found in normalized table using key column '", key_col,
"' for Higher_in_", grp, " → skip saving.\n", sep = "")
next
}

annotated <- dplyr::left_join(sig_norm, ann, by = setNames("Original", key_col))

out_pref <- file.path(folder, paste0(basename(folder), "_Higher_in_", grp, "_annotated"))

write.table(annotated, paste0(out_pref, ".txt"),
sep = "\t", quote = FALSE, row.names = FALSE)
openxlsx::write.xlsx(annotated, paste0(out_pref, ".xlsx"), overwrite = TRUE)

cat("💾 Saved:\n  ", out_pref, ".txt\n  ", out_pref, ".xlsx\n", sep = "")
}
}

cat("\n🎉 STEP 14 FINISHED — annotations written into each comparison folder.\n")
################################################################################

################################################################################
############################## 🔵 STEP 15 — FISHER TEST #########################
################################################################################
# Fisher Exact Test Enrichment for annotated results (Step 14/15)
#
# - GLOBAL background file:
#     DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt
# - Annotation files inside comparison folders:
#     <comparison>/*_annotated.txt
# - Fisher enrichment on:
#     GO_Pathway, KEGG_Pathway, Reactome_Pathway
# - Outputs:
#     Fisher_GO_<base>.xlsx
#     Fisher_KEGG_<base>.xlsx
#     Fisher_Reactome_<base>.xlsx
#     Merged_Fisher_<comparison>_<condition>.xlsx
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
})

# ==============================================================================
# 1) Define comparison folders (YOUR PATH)
# ==============================================================================

main_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")
folders <- file.path(main_dir, comparisons)

# ==============================================================================
# 2) Load GLOBAL background file (matches Step 10 fixed output)
# ==============================================================================

global_background_file <- file.path(
main_dir,
"DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(global_background_file)) {
stop("❌ Global normalized file not found: ", global_background_file)
}

global_norm_df <- read.table(global_background_file, header = TRUE, sep = "\t",
quote = "", check.names = FALSE)

# Pick best available ID column for background
bg_col <- NULL
if ("PG.Genes" %in% colnames(global_norm_df)) {
bg_col <- "PG.Genes"
} else if ("PG.ProteinGroups" %in% colnames(global_norm_df)) {
bg_col <- "PG.ProteinGroups"
} else if ("Protein.IDs" %in% colnames(global_norm_df)) {
bg_col <- "Protein.IDs"
} else {
stop("❌ Background file missing ID columns: PG.Genes / PG.ProteinGroups / Protein.IDs")
}

background_genes <- unique(trimws(as.character(global_norm_df[[bg_col]])))
background_genes <- background_genes[background_genes != "" & !is.na(background_genes)]

cat("✅ Using GLOBAL background:", length(background_genes), "IDs from", bg_col, "\n")

# ==============================================================================
# 3) Helper: pick ID column from annotated file
# ==============================================================================

pick_id_col <- function(df) {
if ("PG.Genes" %in% colnames(df)) return("PG.Genes")
if ("PG.ProteinGroups" %in% colnames(df)) return("PG.ProteinGroups")
if ("Protein.IDs" %in% colnames(df)) return("Protein.IDs")
return(NULL)
}

# ==============================================================================
# 4) Fisher test function (robust to different ID columns)
# ==============================================================================

fisher_on_annotation <- function(ann_df, id_col, annot_col, background_genes, out_xlsx) {

if (is.null(id_col) || !(id_col %in% colnames(ann_df))) return(invisible(NULL))
if (!(annot_col %in% colnames(ann_df))) return(invisible(NULL))

df <- ann_df %>%
filter(!is.na(.data[[annot_col]]) & .data[[annot_col]] != "") %>%
transmute(ID = trimws(as.character(.data[[id_col]])),
Pathway = .data[[annot_col]]) %>%
mutate(Pathway = strsplit(as.character(Pathway), ";")) %>%
unnest(Pathway) %>%
mutate(Pathway = str_trim(Pathway)) %>%
filter(!is.na(ID) & ID != "" & !is.na(Pathway) & Pathway != "" & Pathway != "NA")

# remove fully empty GO rows "NA | NA | NA"
df <- df[!grepl("^NA\\s*\\|\\s*NA\\s*\\|\\s*NA$", df$Pathway, ignore.case = TRUE), ]

# clean " - Homo sapiens" suffix (KEGG)
df$Pathway_Name <- ifelse(
grepl(" - Homo sapiens", df$Pathway),
sub(" - Homo sapiens.*", "", df$Pathway),
df$Pathway
)

sig_ids <- unique(trimws(as.character(ann_df[[id_col]])))
sig_ids <- sig_ids[sig_ids != "" & !is.na(sig_ids)]

terms <- unique(df$Pathway_Name)
if (length(terms) == 0 || length(sig_ids) == 0 || length(background_genes) == 0)
return(invisible(NULL))

out_list <- list()

for (term in terms) {

term_ids <- unique(df$ID[df$Pathway_Name == term & df$ID %in% background_genes])

a <- length(intersect(term_ids, sig_ids))
b <- length(setdiff(sig_ids, term_ids))
c <- length(setdiff(term_ids, intersect(term_ids, sig_ids)))
d <- length(setdiff(background_genes, union(sig_ids, term_ids)))

if (any(c(a, b, c, d) < 0)) next

ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)), error = function(e) NULL)
if (is.null(ft)) next

gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
bg_ratio   <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)

out_list[[term]] <- data.frame(
Pathway               = term,
Sig_protein_volcano   = a,
Sig_NotIn_Pathway     = b,
Nonsig_In_Pathway     = c,
Nonsig_NotIn_Pathway  = d,
IDs                   = paste(intersect(term_ids, sig_ids), collapse = ";"),
p_value               = ft$p.value,
enrichment_factor     = enrichment_factor,
stringsAsFactors      = FALSE
)
}

res <- dplyr::bind_rows(out_list)
if (is.null(res) || nrow(res) == 0) return(invisible(NULL))

res$p_adj <- p.adjust(res$p_value, method = "BH")
res <- res %>% arrange(p_adj)

openxlsx::write.xlsx(res, out_xlsx, overwrite = TRUE)
cat("💾 Saved:", basename(out_xlsx), "(", nrow(res), "terms)\n")

invisible(res)
}

# ==============================================================================
# 5) Main loop over comparison folders
# ==============================================================================

for (folder in folders) {
cat("\n==================================================================\n")
cat("📂 Processing:", folder, "\n")
cat("==================================================================\n")

if (!dir.exists(folder)) {
cat("⚠️ Folder not found:", folder, "→ skipping.\n")
next
}

annot_files <- list.files(folder, pattern = "_annotated\\.txt$", full.names = TRUE)

if (!length(annot_files)) {
cat("⚠️ No annotated files found in", folder, "→ skipping.\n")
next
}

out_dir <- file.path(folder, "Fisher")
dir.create(out_dir, showWarnings = FALSE)

merged_condition <- list()

for (annot_file in annot_files) {

base <- tools::file_path_sans_ext(basename(annot_file))

# ✅ robust condition extraction: expects "..._Higher_in_<COND>_annotated.txt"
condition <- sub("^.*_Higher_in_", "", basename(annot_file))
condition <- sub("_annotated\\.txt$", "", condition)

cat("🔬 File:", basename(annot_file), "→ condition:", condition, "\n")

ann_df <- read.table(
annot_file, header = TRUE, sep = "\t",
quote = "", stringsAsFactors = FALSE, check.names = FALSE
)

id_col <- pick_id_col(ann_df)
if (is.null(id_col)) {
cat("❌ No ID column found (PG.Genes / PG.ProteinGroups / Protein.IDs) → skip\n")
next
}

res_go <- res_kegg <- res_rea <- NULL

if ("GO_Pathway" %in% colnames(ann_df))
res_go <- fisher_on_annotation(
ann_df, id_col, "GO_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_GO_", base, ".xlsx"))
)

if ("KEGG_Pathway" %in% colnames(ann_df))
res_kegg <- fisher_on_annotation(
ann_df, id_col, "KEGG_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_KEGG_", base, ".xlsx"))
)

if ("Reactome_Pathway" %in% colnames(ann_df))
res_rea <- fisher_on_annotation(
ann_df, id_col, "Reactome_Pathway", background_genes,
file.path(out_dir, paste0("Fisher_Reactome_", base, ".xlsx"))
)

merged_piece <- dplyr::bind_rows(
if (!is.null(res_go))   dplyr::mutate(res_go,   Pathway_Type = "GO"),
if (!is.null(res_kegg)) dplyr::mutate(res_kegg, Pathway_Type = "KEGG"),
if (!is.null(res_rea))  dplyr::mutate(res_rea,  Pathway_Type = "Reactome")
)

if (!is.null(merged_piece) && nrow(merged_piece) > 0) {
if (is.null(merged_condition[[condition]])) merged_condition[[condition]] <- list()
merged_condition[[condition]][[base]] <- merged_piece
}
}

# Write merged per condition
for (cond in names(merged_condition)) {
merged_df <- dplyr::bind_rows(merged_condition[[cond]], .id = "Annotated_Set")

if (!is.null(merged_df) && nrow(merged_df) > 0) {
merged_df <- merged_df %>% arrange(p_adj, Pathway_Type, Pathway)

out_merged <- file.path(
out_dir,
paste0("Merged_Fisher_", basename(folder), "_", cond, ".xlsx")
)

openxlsx::write.xlsx(merged_df, out_merged, overwrite = TRUE)

cat("📘 Merged Fisher (GO+KEGG+Reactome):", cond, "→",
basename(out_merged), " (", nrow(merged_df), "rows)\n")
}
}

cat("✅ Fisher enrichment complete for", folder, "\n")
}

cat("\n🎉 Step 15 complete — individual (GO/KEGG/Reactome) + merged-per-condition files saved.\n")
################################################################################


################################################################################
######################## STEP 16 — ENRICHMENT HEATMAPS #########################
################################################################################
# - Uses Fisher outputs from Step 15 (gene-based recommended) or old (IDs)
# - Auto-detects Fisher directories
# - Reads Fisher_(GO|KEGG|Reactome)_*.xlsx
# - Generates publication-grade heatmaps with large genes & pathways
# - HARDENED: all dplyr verbs explicitly namespaced with dplyr:: (no masking)
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)
library(grid)
})

# ======================== SETTINGS ============================================

P_ADJ_CUTOFF <- 0.01
MIN_GENES <- 2
MAX_PATHWAYS_PER_PAGE <- 15
PATHWAY_WRAP_WIDTH <- 45

# BIG TEXT SIZES
BASE_FONTSIZE <- 20
GENE_TEXT_SIZE <- 28
PATHWAY_TEXT_SIZE <- 30
TITLE_TEXT_SIZE <- 40
YLAB_TEXT_SIZE <- 34
LEGEND_TITLE_SIZE <- 30
LEGEND_TEXT_SIZE  <- 26

# Output sizes
FIG_WIDTH  <- 20   # inches
FIG_HEIGHT <- 16   # inches

# IMPORTANT: set project root so list.dirs(".") finds your comparison folders
main_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
setwd(main_dir)
cat("📌 Working dir:", getwd(), "\n")

# Keywords of interest
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

# ======================= HELPERS ==============================================

first_col <- function(df, cand) {
c <- intersect(cand, names(df))
if (length(c)) c[1] else NA_character_
}

count_genes <- function(x) {
if (is.null(x) || is.na(x)) return(0L)
toks <- unlist(stringr::str_split(as.character(x), "[;,/\\s]+"))
toks <- stringr::str_trim(toks)
toks <- toks[toks != "" & !is.na(toks) & toks != "NA"]
length(unique(toks))
}

clean_pathway <- function(pw) {
pw <- stringr::str_remove(pw, "^GO:\\d+\\s*\\|\\s*")
pw <- stringr::str_replace(pw, "(BP|CC|MF)\\s*\\|\\s*", "")
pw <- stringr::str_replace(pw, "^KEGG\\s*\\|\\s*", "")
pw <- stringr::str_replace(pw, "^Reactome\\s*\\|\\s*", "")
stringr::str_trim(pw)
}

add_source_tag <- function(name, src) {
tag <- ifelse(grepl("GO", src, ignore.case = TRUE), "[GO]",
ifelse(grepl("KEGG", src, ignore.case = TRUE), "[KEGG]",
ifelse(grepl("Reactome", src, ignore.case = TRUE), "[Reactome]", "")))
paste0(name, " ", tag)
}

# ====================== FIND FISHER DIRECTORIES ===============================

all_dirs <- list.dirs(".", recursive = TRUE, full.names = TRUE)
fisher_dirs <- all_dirs[grepl("[/\\\\]Fisher$", all_dirs)]

cat("🔎 Fisher dirs found:", length(fisher_dirs), "\n")
if (length(fisher_dirs)) print(fisher_dirs)

if (!length(fisher_dirs)) stop("❌ No Fisher/ folders found under: ", getwd())

# ============================= MAIN LOOP ======================================

for (fdir in fisher_dirs) {

cat("\n📂 Processing Fisher directory:", fdir, "\n")

comparison <- basename(dirname(fdir))

files <- list.files(
fdir,
pattern = "^Fisher_(GO|KEGG|Reactome)_.+\\.xlsx$",
full.names = TRUE
)

if (!length(files)) {
cat("⚠️ No Fisher files in:", fdir, "\n")
next
}

meta <- stringr::str_match(basename(files),
"^Fisher_(GO|KEGG|Reactome)_(.*)\\.xlsx$")
colnames(meta) <- c("full", "Source", "Direction")

directions <- unique(meta[, "Direction"])

# ----------------------------------------------------------------------
# Loop per annotated "Higher_in_*" direction
# ----------------------------------------------------------------------
for (direction in directions) {

cat("   🔎 Direction:", direction, "\n")

f_direction <- files[meta[, "Direction"] == direction]
src_vec     <- meta[meta[, "Direction"] == direction, "Source"]

pieces <- list()

# ---- Load Fisher files ----
for (k in seq_along(f_direction)) {

fp <- f_direction[k]
src <- src_vec[k]

df0 <- tryCatch(openxlsx::read.xlsx(fp), error = function(e) NULL)
if (is.null(df0)) next
if (!"Pathway" %in% names(df0)) next

# p_adj column
pcol <- first_col(df0, c("p_adj", "padj", "adj_p_value", "p.adj"))

# gene list column (new gene-based Fisher writes "Genes"; old writes "IDs")
gcol <- first_col(df0, c("Genes", "IDs", "Gene", "GeneID"))

if (is.na(pcol) || is.na(gcol)) next

df <- df0 %>%
dplyr::mutate(
p_adj_raw = suppressWarnings(as.numeric(.data[[pcol]])),
gene_raw  = as.character(.data[[gcol]]),
gene_n    = vapply(gene_raw, count_genes, integer(1))
) %>%
dplyr::filter(
is.finite(p_adj_raw),
p_adj_raw < P_ADJ_CUTOFF,
gene_n >= MIN_GENES
) %>%
dplyr::mutate(
Source     = src,
Pathway    = stringr::str_trim(Pathway),
log10_p    = -log10(p_adj_raw),
Clean_Path = clean_pathway(Pathway)
) %>%
dplyr::filter(
stringr::str_detect(
stringr::str_to_lower(Clean_Path),
paste(stringr::str_to_lower(global_keywords), collapse = "|")
)
) %>%
tidyr::separate_rows(gene_raw, sep = "[;,/\\s]+") %>%
dplyr::mutate(Gene = stringr::str_trim(gene_raw)) %>%
dplyr::filter(Gene != "" & !is.na(Gene) & Gene != "NA") %>%
dplyr::mutate(
Pathway_Display = add_source_tag(
stringr::str_wrap(Clean_Path, width = PATHWAY_WRAP_WIDTH),
Source
)
)

pieces[[length(pieces) + 1]] <- df
}

side_df <- dplyr::bind_rows(pieces)

if (!nrow(side_df)) {
cat("   ⚠️ No significant pathways for:", direction, "\n")
next
}

# ---- Gene ordering ----
gene_order <- side_df %>%
dplyr::count(Gene, sort = TRUE) %>%
dplyr::pull(Gene)

# ---- Pathway ordering ----
pw_counts <- side_df %>%
dplyr::group_by(Pathway_Display) %>%
dplyr::summarise(n = dplyr::n_distinct(Gene), .groups = "drop")

# ---- Build long tile matrix ----
mat <- side_df %>%
dplyr::select(Pathway_Display, Gene, log10_p) %>%
tidyr::complete(Pathway_Display, Gene = gene_order)

# ---- Output directory ----
out_dir <- file.path(fdir, "Combined_Enrichment_Plots", direction)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_pathways <- unique(mat$Pathway_Display)
pages <- split(all_pathways, ceiling(seq_along(all_pathways) / MAX_PATHWAYS_PER_PAGE))

# ========================== PLOT EACH PAGE ===================================
for (i in seq_along(pages)) {

page_pw <- pages[[i]]

pw_order <- pw_counts %>%
dplyr::filter(Pathway_Display %in% page_pw) %>%
dplyr::arrange(desc(n)) %>%
dplyr::pull(Pathway_Display)

dplot <- mat %>%
dplyr::filter(Pathway_Display %in% page_pw) %>%
dplyr::mutate(
Pathway_Display = factor(Pathway_Display, levels = rev(pw_order)),
Gene = factor(Gene, levels = gene_order)
)

title_txt <- paste0(comparison)

p <- ggplot2::ggplot(dplot, ggplot2::aes(x = Gene, y = Pathway_Display, fill = log10_p)) +
ggplot2::geom_tile(
color = "white",
linewidth = 0.3,
width = 0.95,
height = 0.65
) +
ggplot2::scale_fill_gradientn(
colours = c("green", "yellow", "red"),
na.value = "black",
name = "-log10(p.adj)"
) +
ggplot2::scale_x_discrete(position = "top") +
ggplot2::labs(
title = NULL,
x = NULL,
y = "Enriched Pathway"
) +
ggplot2::theme_minimal(base_size = BASE_FONTSIZE) +
ggplot2::theme(
panel.grid.major = ggplot2::element_blank(),
panel.grid.minor = ggplot2::element_blank(),

axis.text.x = ggplot2::element_text(
size = GENE_TEXT_SIZE,
angle = 90,
vjust = 0.5,
hjust = 1,
face = "bold"
),
axis.text.y = ggplot2::element_text(
size = PATHWAY_TEXT_SIZE,
face = "bold"
),
axis.title.y = ggplot2::element_text(
size = YLAB_TEXT_SIZE,
face = "bold"
),
plot.title = ggplot2::element_text(
size = TITLE_TEXT_SIZE,
face = "bold",
hjust = 0.5
),

legend.position = "right",
legend.title = ggplot2::element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
legend.text  = ggplot2::element_text(size = LEGEND_TEXT_SIZE),
legend.key.height = grid::unit(40, "pt"),
legend.key.width  = grid::unit(18, "pt")
)

out_base <- file.path(out_dir, paste0("Enrichment_", comparison, "_", direction, "_Page", i))

ggplot2::ggsave(
paste0(out_base, ".pdf"),
p,
width = FIG_WIDTH,
height = FIG_HEIGHT,
units = "in",
device = cairo_pdf,
bg = "white"
)

ggplot2::ggsave(
paste0(out_base, ".png"),
p,
width = FIG_WIDTH,
height = FIG_HEIGHT,
units = "in",
dpi = 500,
bg = "white"
)

cat("   ✅ Saved:", out_base, "\n")
}
}
}

cat("\n🎉 Step 16 complete — large-text enrichment heatmaps generated.\n")


################################################################################
################################################################################
################################################################################
################################################################################
### All-Groups Clustered Sankey (GGALLUVIAL) — ALL CLUSTERS, Healthy First    ###
### Reads Step 15 Fisher_(GO|KEGG|Reactome)_..._annotated.xlsx                ###
### Robust to Fisher columns: Pathway/Term/Description and IDs/Genes          ###
################################################################################
################################################################################
################################################################################
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
# 1. SETTINGS
# ==============================================================================

root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ---- Healthy, Pre_RT, Post_RT (Healthy always on top) ----
group_order <- c("Healthy", "Pre_RT", "Post_RT")
group_colors <- c(
Healthy = "#33a02c",
Pre_RT  = "#1f78b4",
Post_RT = "#d73027"
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

# ------------------------------------------------------------------------------
# GLOBAL KEYWORDS (ALL CLUSTERS)
# ------------------------------------------------------------------------------
global_keywords <- unique(unlist(clusters[names(clusters) != "Other"]))

# ==============================================================================
# 2. HELPERS
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

read_one <- function(file, grp) {

df <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
if (is.null(df) || !nrow(df)) return(NULL)

# ---- pathway column robust ----
if (!"Pathway" %in% names(df)) {
alt <- first_col(df, c("Pathway_Name","Clean_Path","Term","Description","NAME","name"))
if (is.na(alt)) return(NULL)
df$Pathway <- df[[alt]]
}

# ---- p_adj and gene/protein list column robust ----
pcol <- first_col(df, c("p_adj","p.adj","padj","adj_p_value","adj.P.Val","p.value","p_value"))
# Step 15 writes "IDs" (important!)
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
dplyr::mutate(Group = grp) %>%
dplyr::distinct(Pathway, Gene, p_adj, Group)

if (!nrow(out)) return(NULL)
out
}

# ==============================================================================
# 3. READ ALL DATA
# ==============================================================================

all_dat <- list()

for (cmp in comparisons) {

f_dir <- file.path(root, cmp, "Fisher")
if (!dir.exists(f_dir)) next

# Step 15 example:
# Fisher_Reactome_Healthy_vs_Post_RT_Higher_in_Healthy_annotated.xlsx
files <- list.files(
f_dir,
pattern = "^Fisher_(GO|KEGG|Reactome)_.+_annotated\\.xlsx$",
full.names = TRUE,
ignore.case = TRUE
)

for (f in files) {

bn <- basename(f)

grp <- dplyr::case_when(
stringr::str_detect(bn, stringr::regex("Higher_in_Healthy", ignore_case = TRUE))  ~ "Healthy",
stringr::str_detect(bn, stringr::regex("Higher_in_Pre_RT",  ignore_case = TRUE))  ~ "Pre_RT",
stringr::str_detect(bn, stringr::regex("Higher_in_Post_RT", ignore_case = TRUE))  ~ "Post_RT",
TRUE ~ NA_character_
)

if (!is.na(grp)) {
df <- read_one(f, grp)
if (!is.null(df)) all_dat[[length(all_dat) + 1]] <- df
}
}
}

dat <- dplyr::bind_rows(all_dat)
if (!nrow(dat)) stop("No pathways found. (Check FDR_CUTOFF / MIN_PROTEINS / keyword filter)")

# ==============================================================================
# 4. KEEP ALL GROUPS (order: Healthy, Pre_RT, Post_RT)
# ==============================================================================

dat <- dat %>% dplyr::filter(Group %in% group_order)

# ==============================================================================
# 5. FILTER TO TOP PATHWAYS
# ==============================================================================

top_paths <- dat %>%
dplyr::count(Pathway, sort = TRUE) %>%
dplyr::slice_head(n = TOP_N) %>%
dplyr::pull(Pathway)

dat <- dat %>% dplyr::filter(Pathway %in% top_paths)

# ==============================================================================
# 6. ASSIGN CLUSTERS
# ==============================================================================

dat$Cluster <- vapply(dat$Pathway, function(p) {
pl <- stringr::str_to_lower(p)
for (cl in names(clusters)) {
if (any(stringr::str_detect(pl, stringr::str_to_lower(clusters[[cl]])))) return(cl)
}
"Other"
}, character(1))

# ==============================================================================
# 7. PREP ALLUVIAL DATA
# ==============================================================================

sankey_df <- dat %>%
dplyr::count(Cluster, Pathway, Group, name = "Freq") %>%
dplyr::ungroup()

# Set levels to enforce group order everywhere
sankey_df$Cluster <- factor(sankey_df$Cluster, levels = names(cluster_colors))
sankey_df$Group   <- factor(sankey_df$Group,   levels = group_order)

pathway_order <- sankey_df %>%
dplyr::group_by(Cluster, Pathway) %>%
dplyr::summarise(total = sum(Freq), .groups = "drop") %>%
dplyr::arrange(factor(Cluster, levels = names(cluster_colors)), Pathway) %>%
dplyr::pull(Pathway)

sankey_df$Pathway <- factor(sankey_df$Pathway, levels = unique(pathway_order))

# ==============================================================================
# 8. PLOT
# ==============================================================================

p <- ggplot(
sankey_df,
aes(axis1 = Cluster, axis3 = Pathway, axis5 = Group, y = Freq)
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
fill = NA, color = NA
) +
ggplot2::geom_text(
stat = "stratum",
aes(label = after_stat(stratum)),
size = 20,
color = "black"
) +
ggplot2::scale_x_discrete(
expand = c(.1, .1),
labels = c("Cluster", "Pathway", "Group")
) +
ggplot2::scale_fill_manual(values = group_colors) +
ggplot2::theme_minimal(base_size = 40) +
ggplot2::theme(
legend.position = "none",
axis.text.y = ggplot2::element_blank(),
axis.ticks = ggplot2::element_blank(),
axis.title = ggplot2::element_blank(),
panel.grid = ggplot2::element_blank(),
plot.title = ggplot2::element_text(hjust = 0.5, size = 70)
) +
ggplot2::ggtitle("Clustered Sankey Diagram — Healthy → Pre_RT → Post_RT")

# --- add black top/bottom borders to each stratum ---
pb <- ggplot_build(p)
idx <- which(sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum")))
stratum_data <- pb$data[[idx]]

p <- p +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
inherit.aes = FALSE, color = "black", size = 3
) +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
inherit.aes = FALSE, color = "black", size = 3
)

# ==============================================================================
# 9. SAVE
# ==============================================================================

out_dir <- file.path(root, "AllGroups_Sankey_GGALLUVIAL_WithHealthy")
if (!dir.exists(out_dir)) dir.create(out_dir)

openxlsx::write.xlsx(
sankey_df,
file.path(out_dir, "Sankey_AllGroups_Clustered_WithHealthy.xlsx"),
rowNames = FALSE,
overwrite = TRUE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_Clustered_WithHealthy.png"),
p, width = 32, height = 45, dpi = 300, bg = "white"
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_Clustered_WithHealthy.pdf"),
p, width = 32, height = 45, bg = "white"
)

cat("\n🎉 ALL DONE — ALL CLUSTERS, Healthy first!\n")


################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
### All-Groups Clustered Sankey (GGALLUVIAL) — NO HEALTHY IN PLOT            ###
### Reads Step 15 Fisher_(GO|KEGG|Reactome)_..._annotated.xlsx               ###
### Groups shown: Pre_RT → Post_RT (Healthy removed from Sankey)             ###
################################################################################
################################################################################
################################################################################
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
# 1. SETTINGS
# ==============================================================================

root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ---- ONLY Pre_RT, Post_RT (Healthy removed) ----
group_order <- c("Pre_RT", "Post_RT")
group_colors <- c(
Pre_RT  = "#1f78b4",
Post_RT = "#d73027"
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

# ------------------------------------------------------------------------------
# GLOBAL KEYWORDS (ALL CLUSTERS)
# ------------------------------------------------------------------------------
global_keywords <- unique(unlist(clusters[names(clusters) != "Other"]))

# ==============================================================================
# 2. HELPERS
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

read_one <- function(file, grp) {

df <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
if (is.null(df) || !nrow(df)) return(NULL)

# ---- pathway column robust ----
if (!"Pathway" %in% names(df)) {
alt <- first_col(df, c("Pathway_Name","Clean_Path","Term","Description","NAME","name"))
if (is.na(alt)) return(NULL)
df$Pathway <- df[[alt]]
}

# ---- p_adj and gene/protein list column robust ----
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
dplyr::mutate(Group = grp) %>%
dplyr::distinct(Pathway, Gene, p_adj, Group)

if (!nrow(out)) return(NULL)
out
}

# ==============================================================================
# 3. READ ALL DATA
# ==============================================================================

all_dat <- list()

for (cmp in comparisons) {

f_dir <- file.path(root, cmp, "Fisher")
if (!dir.exists(f_dir)) next

# Step 15 example:
# Fisher_Reactome_Healthy_vs_Post_RT_Higher_in_Healthy_annotated.xlsx
files <- list.files(
f_dir,
pattern = "^Fisher_(GO|KEGG|Reactome)_.+_annotated\\.xlsx$",
full.names = TRUE,
ignore.case = TRUE
)

for (f in files) {

bn <- basename(f)

grp <- dplyr::case_when(
stringr::str_detect(bn, stringr::regex("Higher_in_Healthy", ignore_case = TRUE))  ~ "Healthy",
stringr::str_detect(bn, stringr::regex("Higher_in_Pre_RT",  ignore_case = TRUE))  ~ "Pre_RT",
stringr::str_detect(bn, stringr::regex("Higher_in_Post_RT", ignore_case = TRUE))  ~ "Post_RT",
TRUE ~ NA_character_
)

if (!is.na(grp)) {
df <- read_one(f, grp)
if (!is.null(df)) all_dat[[length(all_dat) + 1]] <- df
}
}
}

dat <- dplyr::bind_rows(all_dat)
if (!nrow(dat)) stop("No pathways found. (Check FDR_CUTOFF / MIN_PROTEINS / keyword filter)")

# ==============================================================================
# 4. KEEP ONLY Pre_RT and Post_RT (Healthy removed)
# ==============================================================================

dat <- dat %>% dplyr::filter(Group %in% group_order)

# ==============================================================================
# 5. FILTER TO TOP PATHWAYS
# ==============================================================================

top_paths <- dat %>%
dplyr::count(Pathway, sort = TRUE) %>%
dplyr::slice_head(n = TOP_N) %>%
dplyr::pull(Pathway)

dat <- dat %>% dplyr::filter(Pathway %in% top_paths)

# ==============================================================================
# 6. ASSIGN CLUSTERS
# ==============================================================================

dat$Cluster <- vapply(dat$Pathway, function(p) {
pl <- stringr::str_to_lower(p)
for (cl in names(clusters)) {
if (any(stringr::str_detect(pl, stringr::str_to_lower(clusters[[cl]])))) return(cl)
}
"Other"
}, character(1))

# ==============================================================================
# 7. PREP ALLUVIAL DATA
# ==============================================================================

sankey_df <- dat %>%
dplyr::count(Cluster, Pathway, Group, name = "Freq") %>%
dplyr::ungroup()

sankey_df$Cluster <- factor(sankey_df$Cluster, levels = names(cluster_colors))
sankey_df$Group   <- factor(sankey_df$Group,   levels = group_order)

pathway_order <- sankey_df %>%
dplyr::group_by(Cluster, Pathway) %>%
dplyr::summarise(total = sum(Freq), .groups = "drop") %>%
dplyr::arrange(factor(Cluster, levels = names(cluster_colors)), Pathway) %>%
dplyr::pull(Pathway)

sankey_df$Pathway <- factor(sankey_df$Pathway, levels = unique(pathway_order))

# ==============================================================================
# 8. PLOT (2 AXES NOW: Cluster -> Pathway -> Group (Pre_RT/Post_RT))
# ==============================================================================

p <- ggplot(
sankey_df,
aes(axis1 = Cluster, axis3 = Pathway, axis5 = Group, y = Freq)
) +
ggalluvial::geom_alluvium(
aes(fill = Group),
width = 0.5,
alpha = 0.3,
color = "black",
size  = 1
) +
ggalluvial::geom_stratum(
width = 0.2,
fill = NA, color = NA
) +
ggplot2::geom_text(
stat = "stratum",
aes(label = after_stat(stratum)),
size = 20,
color = "black"
) +
ggplot2::scale_x_discrete(
expand = c(.1, .1),
labels = c("Cluster", "Pathway", "Group")
) +
ggplot2::scale_fill_manual(values = group_colors) +
ggplot2::theme_minimal(base_size = 40) +
ggplot2::theme(
legend.position = "none",
axis.text.y = ggplot2::element_blank(),
axis.ticks = ggplot2::element_blank(),
axis.title = ggplot2::element_blank(),
panel.grid = ggplot2::element_blank(),
plot.title = ggplot2::element_text(hjust = 0.5, size = 70)
) +
ggplot2::ggtitle("Clustered Sankey Diagram — Pre_RT → Post_RT (Healthy removed)")

# --- add black top/bottom borders to each stratum ---
pb <- ggplot_build(p)
idx <- which(sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum")))
stratum_data <- pb$data[[idx]]

p <- p +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
inherit.aes = FALSE, color = "black", size = 3
) +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
inherit.aes = FALSE, color = "black", size = 3
)

# ==============================================================================
# 9. SAVE
# ==============================================================================

out_dir <- file.path(root, "AllGroups_Sankey_GGALLUVIAL_NoHealthy")
if (!dir.exists(out_dir)) dir.create(out_dir)

openxlsx::write.xlsx(
sankey_df,
file.path(out_dir, "Sankey_AllGroups_Clustered_NoHealthy.xlsx"),
rowNames = FALSE,
overwrite = TRUE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_Clustered_NoHealthy.png"),
p, width = 32, height = 45, dpi = 300, bg = "white"
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_Clustered_NoHealthy.pdf"),
p, width = 32, height = 45, bg = "white"
)

cat("\n🎉 ALL DONE — ALL CLUSTERS, Healthy removed from Sankey!\n")
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
### All-Groups 2-Way Sankey (GGALLUVIAL) — NO CLUSTERS, Healthy First         ###
### Reads Step 15 Fisher_(GO|KEGG|Reactome)_..._annotated.xlsx                ###
### Robust to Fisher columns: Pathway/Term/Description and IDs/Genes          ###
################################################################################
################################################################################
################################################################################
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
# 1. SETTINGS
# ==============================================================================

root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ---- Healthy, Pre_RT, Post_RT (Healthy always on top) ----
group_order <- c("Healthy", "Pre_RT", "Post_RT")
group_colors <- c(
Healthy = "#33a02c",
Pre_RT  = "#1f78b4",
Post_RT = "#d73027"
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

# ------------------------------------------------------------------------------
# GLOBAL KEYWORDS (ALL CLUSTERS)
# ------------------------------------------------------------------------------
global_keywords <- unique(unlist(clusters[names(clusters) != "Other"]))

# ==============================================================================
# 2. HELPERS
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

read_one <- function(file, grp) {

df <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
if (is.null(df) || !nrow(df)) return(NULL)

# ---- pathway column robust ----
if (!"Pathway" %in% names(df)) {
alt <- first_col(df, c("Pathway_Name","Clean_Path","Term","Description","NAME","name"))
if (is.na(alt)) return(NULL)
df$Pathway <- df[[alt]]
}

# ---- p_adj and gene/protein list column robust ----
pcol <- first_col(df, c("p_adj","p.adj","padj","adj_p_value","adj.P.Val","p.value","p_value"))
# Step 15 writes "IDs" (important!)
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
dplyr::mutate(Group = grp) %>%
dplyr::distinct(Pathway, Gene, p_adj, Group)

if (!nrow(out)) return(NULL)
out
}

# ==============================================================================
# 3. READ ALL DATA
# ==============================================================================

all_dat <- list()

for (cmp in comparisons) {

f_dir <- file.path(root, cmp, "Fisher")
if (!dir.exists(f_dir)) next

# Step 15 example:
# Fisher_Reactome_Healthy_vs_Post_RT_Higher_in_Healthy_annotated.xlsx
files <- list.files(
f_dir,
pattern = "^Fisher_(GO|KEGG|Reactome)_.+_annotated\\.xlsx$",
full.names = TRUE,
ignore.case = TRUE
)

for (f in files) {

bn <- basename(f)

grp <- dplyr::case_when(
stringr::str_detect(bn, stringr::regex("Higher_in_Healthy", ignore_case = TRUE))  ~ "Healthy",
stringr::str_detect(bn, stringr::regex("Higher_in_Pre_RT",  ignore_case = TRUE))  ~ "Pre_RT",
stringr::str_detect(bn, stringr::regex("Higher_in_Post_RT", ignore_case = TRUE))  ~ "Post_RT",
TRUE ~ NA_character_
)

if (!is.na(grp)) {
df <- read_one(f, grp)
if (!is.null(df)) all_dat[[length(all_dat) + 1]] <- df
}
}
}

dat <- dplyr::bind_rows(all_dat)
if (!nrow(dat)) stop("No pathways found. (Check FDR_CUTOFF / MIN_PROTEINS / keyword filter)")

# ==============================================================================
# 4. KEEP ALL GROUPS (order: Healthy, Pre_RT, Post_RT)
# ==============================================================================

dat <- dat %>% dplyr::filter(Group %in% group_order)

# ==============================================================================
# 5. FILTER TO TOP PATHWAYS
# ==============================================================================

top_paths <- dat %>%
dplyr::count(Pathway, sort = TRUE) %>%
dplyr::slice_head(n = TOP_N) %>%
dplyr::pull(Pathway)

dat <- dat %>% dplyr::filter(Pathway %in% top_paths)

# ==============================================================================
# 6. ASSIGN CLUSTERS
#    (kept for full script parity; NOT USED in the 2-way Sankey)
# ==============================================================================

dat$Cluster <- vapply(dat$Pathway, function(p) {
pl <- stringr::str_to_lower(p)
for (cl in names(clusters)) {
if (any(stringr::str_detect(pl, stringr::str_to_lower(clusters[[cl]])))) return(cl)
}
"Other"
}, character(1))

# ==============================================================================
# 7. PREP ALLUVIAL DATA (2-WAY: Pathway <-> Group)
# ==============================================================================

sankey_df <- dat %>%
dplyr::count(Pathway, Group, name = "Freq") %>%
dplyr::ungroup()

# Set levels to enforce group order everywhere
sankey_df$Group <- factor(sankey_df$Group, levels = group_order)

# Keep pathway order stable (by total frequency, then name)
pathway_order <- sankey_df %>%
dplyr::group_by(Pathway) %>%
dplyr::summarise(total = sum(Freq), .groups = "drop") %>%
dplyr::arrange(dplyr::desc(total), Pathway) %>%
dplyr::pull(Pathway)

sankey_df$Pathway <- factor(sankey_df$Pathway, levels = unique(pathway_order))

# ==============================================================================
# 8. PLOT (2 AXES ONLY: Pathway, Group)
# ==============================================================================

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
fill = NA, color = NA
) +
ggplot2::geom_text(
stat = "stratum",
aes(label = after_stat(stratum)),
size = 20,
color = "black"
) +
ggplot2::scale_x_discrete(
# >>> FIX: give the left side more breathing room too <<<
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

# >>> FIX: increase margins (left is the key) <<<
# units are "pt" by default; large left margin prevents clipping
plot.margin = ggplot2::margin(t = 20, r = 40, b = 20, l = 100, unit = "pt")
) +

# >>> FIX: allow drawing outside panel (so labels can extend into the margin) <<<
coord_cartesian(clip = "off") +

ggplot2::ggtitle("Sankey Diagram (Healthy → Pre_RT → Post_RT)")

# --- add black top/bottom borders to each stratum ---
pb <- ggplot_build(p)
idx <- which(sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum")))
stratum_data <- pb$data[[idx]]

p <- p +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
inherit.aes = FALSE, color = "black", size = 3
) +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
inherit.aes = FALSE, color = "black", size = 3
)

# ==============================================================================
# 9. SAVE
# ==============================================================================

out_dir <- file.path(root, "2Way_AllGroups_Sankey_GGALLUVIAL_WithHealthy")
if (!dir.exists(out_dir)) dir.create(out_dir)

openxlsx::write.xlsx(
sankey_df,
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_WithHealthy.xlsx"),
rowNames = FALSE,
overwrite = TRUE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_WithHealthy.png"),
p, width = 35, height = 45, dpi = 300, bg = "white",
limitsize = FALSE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_WithHealthy.pdf"),
p, width = 35, height = 45, bg = "white",
limitsize = FALSE
)

cat("\n🎉 ALL DONE — 2-WAY Sankey (Pathway ↔ Group), NO CLUSTERS, Healthy first!\n")
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
### All-Groups 2-Way Sankey (GGALLUVIAL) — NO CLUSTERS, NO HEALTHY            ###
### Reads Step 15 Fisher_(GO|KEGG|Reactome)_..._annotated.xlsx                ###
### Robust to Fisher columns: Pathway/Term/Description and IDs/Genes          ###
### Pathway order: A -> Z                                                    ###
################################################################################
################################################################################
################################################################################
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
# 1. SETTINGS
# ==============================================================================

root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
comparisons <- c("Pre_RT_vs_Healthy", "Healthy_vs_Post_RT", "Pre_RT_vs_Post_RT")

FDR_CUTOFF   <- 0.01
MIN_PROTEINS <- 2
TOP_N        <- 20

# ---- Pre_RT, Post_RT (NO Healthy) ----
group_order <- c("Pre_RT", "Post_RT")
group_colors <- c(
Pre_RT  = "#1f78b4",
Post_RT = "#d73027"
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

# ------------------------------------------------------------------------------
# GLOBAL KEYWORDS (ALL CLUSTERS)
# ------------------------------------------------------------------------------
global_keywords <- unique(unlist(clusters[names(clusters) != "Other"]))

# ==============================================================================
# 2. HELPERS
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

read_one <- function(file, grp) {

df <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
if (is.null(df) || !nrow(df)) return(NULL)

# ---- pathway column robust ----
if (!"Pathway" %in% names(df)) {
alt <- first_col(df, c("Pathway_Name","Clean_Path","Term","Description","NAME","name"))
if (is.na(alt)) return(NULL)
df$Pathway <- df[[alt]]
}

# ---- p_adj and gene/protein list column robust ----
pcol <- first_col(df, c("p_adj","p.adj","padj","adj_p_value","adj.P.Val","p.value","p_value"))
# Step 15 writes "IDs" (important!)
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
dplyr::mutate(Group = grp) %>%
dplyr::distinct(Pathway, Gene, p_adj, Group)

if (!nrow(out)) return(NULL)
out
}

# ==============================================================================
# 3. READ ALL DATA
# ==============================================================================

all_dat <- list()

for (cmp in comparisons) {

f_dir <- file.path(root, cmp, "Fisher")
if (!dir.exists(f_dir)) next

# Step 15 example:
# Fisher_Reactome_Healthy_vs_Post_RT_Higher_in_Healthy_annotated.xlsx
files <- list.files(
f_dir,
pattern = "^Fisher_(GO|KEGG|Reactome)_.+_annotated\\.xlsx$",
full.names = TRUE,
ignore.case = TRUE
)

for (f in files) {

bn <- basename(f)

grp <- dplyr::case_when(
stringr::str_detect(bn, stringr::regex("Higher_in_Pre_RT",  ignore_case = TRUE))  ~ "Pre_RT",
stringr::str_detect(bn, stringr::regex("Higher_in_Post_RT", ignore_case = TRUE))  ~ "Post_RT",
TRUE ~ NA_character_
)

if (!is.na(grp)) {
df <- read_one(f, grp)
if (!is.null(df)) all_dat[[length(all_dat) + 1]] <- df
}
}
}

dat <- dplyr::bind_rows(all_dat)
if (!nrow(dat)) stop("No pathways found. (Check FDR_CUTOFF / MIN_PROTEINS / keyword filter)")

# ==============================================================================
# 4. KEEP ONLY Pre_RT & Post_RT (order: Pre_RT, Post_RT)
# ==============================================================================

dat <- dat %>% dplyr::filter(Group %in% group_order)

# ==============================================================================
# 5. FILTER TO TOP PATHWAYS (by overall frequency), THEN ORDER PATHWAYS A->Z
# ==============================================================================

top_paths <- dat %>%
dplyr::count(Pathway, sort = TRUE) %>%
dplyr::slice_head(n = TOP_N) %>%
dplyr::pull(Pathway)

dat <- dat %>% dplyr::filter(Pathway %in% top_paths)

# ==============================================================================
# 6. ASSIGN CLUSTERS
#    (kept for full script parity; NOT USED in the 2-way Sankey)
# ==============================================================================

dat$Cluster <- vapply(dat$Pathway, function(p) {
pl <- stringr::str_to_lower(p)
for (cl in names(clusters)) {
if (any(stringr::str_detect(pl, stringr::str_to_lower(clusters[[cl]])))) return(cl)
}
"Other"
}, character(1))

# ==============================================================================
# 7. PREP ALLUVIAL DATA (2-WAY: Pathway <-> Group)
# ==============================================================================

sankey_df <- dat %>%
dplyr::count(Pathway, Group, name = "Freq") %>%
dplyr::ungroup()

# Set levels to enforce group order everywhere
sankey_df$Group <- factor(sankey_df$Group, levels = group_order)

# --- Pathway order A -> Z ---
pathway_order <- sankey_df %>%
dplyr::distinct(Pathway) %>%
dplyr::pull(Pathway) %>%
sort()

sankey_df$Pathway <- factor(sankey_df$Pathway, levels = pathway_order)

# ==============================================================================
# 8. PLOT (2 AXES ONLY: Pathway, Group)
# ==============================================================================

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
fill = NA, color = NA
) +
ggplot2::geom_text(
stat = "stratum",
aes(label = after_stat(stratum)),
size = 20,
color = "black"
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
plot.margin = ggplot2::margin(t = 20, r = 40, b = 20, l = 100, unit = "pt")
) +
coord_cartesian(clip = "off") +
ggplot2::ggtitle("Sankey Diagram (Pre_RT → Post_RT)")

# --- add black top/bottom borders to each stratum ---
pb <- ggplot_build(p)
idx <- which(sapply(pb$plot$layers, function(l) inherits(l$geom, "GeomStratum")))
stratum_data <- pb$data[[idx]]

p <- p +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
inherit.aes = FALSE, color = "black", size = 3
) +
ggplot2::geom_segment(
data = stratum_data,
aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
inherit.aes = FALSE, color = "black", size = 3
)

# ==============================================================================
# 9. SAVE
# ==============================================================================

out_dir <- file.path(root, "2Way_AllGroups_Sankey_GGALLUVIAL_NoHealthy")
if (!dir.exists(out_dir)) dir.create(out_dir)

openxlsx::write.xlsx(
sankey_df,
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_NoHealthy.xlsx"),
rowNames = FALSE,
overwrite = TRUE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_NoHealthy.png"),
p, width = 35, height = 45, dpi = 300, bg = "white",
limitsize = FALSE
)

ggplot2::ggsave(
file.path(out_dir, "Sankey_AllGroups_2Way_Pathway_Group_NoHealthy.pdf"),
p, width = 35, height = 45, bg = "white",
limitsize = FALSE
)

cat("\n🎉 ALL DONE — 2-WAY Sankey (Pathway ↔ Group), NO CLUSTERS, NO Healthy! (Pathways A→Z)\n")



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
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(broom)

# --- Set seed for reproducibility ---
set.seed(1)

# --- Set working directory ---
setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics")

# --- Load and preprocess data ---
data <- read.table("DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt", 
header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Reshape to long format ---
data_long <- data %>%
pivot_longer(
cols = grep("LFQ_", colnames(data)),   # selects all LFQ columns
names_to = "Condition",
values_to = "Expression"
) %>%
mutate(Group = case_when(
grepl("Healthy", Condition) ~ "Healthy",
grepl("Pre_RT", Condition) ~ "Pre_RT",
grepl("Post_RT", Condition) ~ "Ende_RT",
TRUE ~ NA_character_
)) %>%
mutate(Group = factor(Group, levels = c("Healthy", "Pre_RT", "Ende_RT")),
Group_num = as.numeric(Group))

# --- Calculate mean expression for plotting and trend logic ---
mean_expression <- data_long %>%
group_by(PG.Genes, Group) %>%
summarize(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# --- Pearson correlation per gene ---
trend_results <- data_long %>%
group_by(PG.Genes) %>%
summarize(
pearson_r = cor(Expression, Group_num, method = "pearson", use = "complete.obs"),
pearson_p = cor.test(Expression, Group_num, method = "pearson")$p.value,
.groups = "drop"
)

# --- Determine U-shaped or inverted-U-shaped trends ---
trend_pattern <- mean_expression %>%
pivot_wider(names_from = Group, values_from = mean_expr) %>%
mutate(
trend = case_when(
Healthy > Pre_RT & Ende_RT > Pre_RT ~ "Tumor suppressor protein",
Healthy < Pre_RT & Ende_RT < Pre_RT ~ "Oncoprotein",
TRUE ~ "No Trend"
)
) %>%
dplyr::select("PG.Genes", "trend")   # FIXED HERE

# --- Merge all results ---
final_results <- left_join(trend_pattern, trend_results, by = "PG.Genes")

# --- Save results ---
write_xlsx(final_results, 
path = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/Pearson_dose_dependent_results.xlsx")

write.table(final_results, 
file = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/Pearson_dose_dependent_results.txt",
sep = "\t", row.names = FALSE, quote = FALSE)

# --- Plotting ---
plot <- ggplot(mean_expression, aes(x = Group, y = mean_expr, color = PG.Genes, group = PG.Genes)) +
geom_point() +
geom_line() +
theme_minimal() +
labs(title = "Mean Expression Trends by Group (Pearson)", y = "Mean Expression")

print(plot)

# --- Merge trends back to original data ---
data_with_trends <- left_join(data, final_results, by = "PG.Genes")

# Save full matrix
write_xlsx(data_with_trends, 
path = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/Pearson_dose_dependent_results_with_trends.xlsx")

write.table(data_with_trends, 
file = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/Pearson_dose_dependent_results_with_trends.txt",
sep = "\t", row.names = FALSE, quote = FALSE)

# View preview
head(data_with_trends)

################################################################################

###############################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
###########################         SOTA        ################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
# FIXED SOTA SCRIPT — FULL VERSION WITH CORRECT MERGING
################################################################################

library(openxlsx)
library(clValid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# ============================================
# GLOBAL SEED FOR REPRODUCIBILITY
# ============================================
set.seed(1)

setwd("C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics")

# Load annotated Pearson file (Step 1 result)
sota_data <- read.table(
"Pearson_dose_dependent_results_with_trends.txt",
header = TRUE, sep = "\t", quote = "", check.names = FALSE
)

trend_groups <- c("Oncoprotein", "Tumor suppressor protein")

# Output root folder
root_dir <- "SOTA_Cluster"
dir.create(root_dir, showWarnings = FALSE)

################################################################################
# Antibody filter
################################################################################

antibody_keywords <- c(
"IGHV","IGLV","IGKV","IGHA","IGHD","IGHE","IGHG","IGLC","IGK","IGH","IGL","IGKV1"
)
antibody_pattern <- paste(antibody_keywords, collapse = "|")

contains_antibody_gene <- function(gene_string) {
genes <- unlist(strsplit(gene_string, ";"))
any(grepl(antibody_pattern, genes, ignore.case = TRUE))
}

################################################################################
# Expand multi-symbol PG.Genes into single gene rows
################################################################################

expand_gene_rows <- function(df) {
df$PG.Genes <- gsub('"|\'', "", df$PG.Genes)
df$PG.Genes <- trimws(df$PG.Genes)

split_rows <- strsplit(df$PG.Genes, ";")

expanded_df <- do.call(rbind, lapply(seq_along(split_rows), function(i) {
genes <- trimws(split_rows[[i]])
rows  <- df[rep(i, length(genes)), , drop = FALSE]
rows$PG.Genes <- genes
rows
}))

rownames(expanded_df) <- make.unique(expanded_df$PG.Genes)
expanded_df
}

################################################################################
# START PROCESSING
################################################################################

for (trend_type in trend_groups) {

message("Processing trend group: ", trend_type)

# ---------------- FILTER ----------------
filtered_sota_data <- sota_data[sota_data$trend == trend_type, ]
filtered_sota_data <- filtered_sota_data[
!is.na(filtered_sota_data$PG.Genes) & filtered_sota_data$PG.Genes != "", ]

# Remove antibody proteins
filtered_sota_data <- filtered_sota_data[
!sapply(filtered_sota_data$PG.Genes, contains_antibody_gene), ]

# Expand rows (IMPORTANT!)
expanded_data <- expand_gene_rows(filtered_sota_data)

if (nrow(expanded_data) < 3) {
message("Too few genes for ", trend_type)
next
}

# ---------------- CREATE TREND FOLDER ----------------
trend_folder <- file.path(root_dir, gsub(" ", "_", trend_type))
dir.create(trend_folder, recursive = TRUE, showWarnings = FALSE)

# Save expanded table
write.table(
expanded_data,
file = file.path(trend_folder,
paste0(gsub(" ", "_", trend_type), "_Filtered_Full_Table.txt")),
sep = "\t", quote = FALSE, row.names = FALSE
)

# ---------------- EXTRACT LFQ MATRIX ----------------
treatment_cols <- grep(
"LFQ_Healthy_|LFQ_Pre_RT_|LFQ_Post_RT_",
colnames(expanded_data),
value = TRUE
)

express <- expanded_data[, treatment_cols]
gene_names <- make.unique(expanded_data$PG.Genes)
rownames(expanded_data) <- gene_names
rownames(express) <- gene_names

express_clean <- express[complete.cases(express), ]
expanded_clean <- expanded_data[rownames(express_clean), ]
express_matrix <- as.matrix(express_clean)

if (nrow(express_matrix) < 3) {
message("Not enough genes after NA filtering for ", trend_type)
next
}

# ---------------- ELBOW PLOT ----------------
# ============================================
# GLOBAL SEED FOR REPRODUCIBILITY
# ============================================
set.seed(1)
max_k <- min(20, nrow(express_matrix) - 1)
wss <- sapply(2:max_k, function(k)
kmeans(express_matrix, centers = k, nstart = 10)$tot.withinss)

png(file.path(trend_folder, "elbow_plot.png"), 800, 600)
plot(2:max_k, wss, type = "b", main = paste("Elbow -", trend_type),
xlab = "k", ylab = "WSS")
dev.off()

# ---------------- SOTA CLUSTERING ----------------
set.seed(1)
optimal_clusters <- 4
sotaCl <- sota(
data = express_matrix,
maxCycles = optimal_clusters,
maxEpochs = 2000,
distance = "euclidean",
wcell = 0.01,
pcell = 0.005,
scell = 0.001,
delta = 1e-4,
neighb.level = 0,
maxDiversity = 0.5,
unrest.growth = TRUE
)

cluster_labels <- sotaCl$clust
names(cluster_labels) <- rownames(express_matrix)

# ---------------- PER-CLUSTER EXPORT ----------------
cluster_list <- split(names(cluster_labels), cluster_labels)
wb <- createWorkbook()

for (cl_id in names(cluster_list)) {

genes <- cluster_list[[cl_id]]

cluster_dir <- file.path(trend_folder, paste0("Cluster_", cl_id))
dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

# Extract expression
cluster_matrix <- express_matrix[genes, , drop = FALSE]
cluster_mean <- colMeans(cluster_matrix, na.rm = TRUE)

# ---------------- FIXED MERGE ----------------
# Merge MUST USE expanded_clean (single-gene rows, correct LFQs)
gene_df <- data.frame(PG.Genes = genes, stringsAsFactors = FALSE)

cluster_full <- merge(
expanded_clean,         # FIX: not sota_data!
gene_df,
by = "PG.Genes",
all.y = TRUE,
sort = FALSE
)

# reorder
cluster_full$PG.Genes <- factor(cluster_full$PG.Genes, levels = genes)
cluster_full <- cluster_full[order(cluster_full$PG.Genes), ]
cluster_full$PG.Genes <- as.character(cluster_full$PG.Genes)

# save txt
write.table(cluster_full,
file = file.path(cluster_dir, paste0("Cluster_", cl_id, "_genes.txt")),
sep = "\t", quote = FALSE, row.names = FALSE)

# save excel
write.xlsx(cluster_full,
file = file.path(cluster_dir, paste0("Cluster_", cl_id, "_genes.xlsx")),
overwrite = TRUE)

# write gene-only list
addWorksheet(wb, paste0("Cluster_", cl_id))
writeData(wb, sheet = paste0("Cluster_", cl_id), genes, colNames = FALSE)
}

saveWorkbook(wb,
file = file.path(trend_folder, "Cluster_Protein_Names_List.xlsx"),
overwrite = TRUE)

message("Completed ", trend_type)
}

cat("All trend groups processed.\n")

################################################################################

################################################################################
# SUMMARY MULTIPANEL PLOT (READS Cluster_*_genes.txt FROM FOLDERS)
# Folder layout:
# SOTA_Cluster/<Trend>/Cluster_<id>/Cluster_<id>_genes.txt
################################################################################

root_dir <- "SOTA_Cluster"

# extract cluster number from "Cluster_12"
get_cl_id <- function(cluster_dir_name) {
as.integer(sub("^Cluster_([0-9]+)$", "\\1", cluster_dir_name))
}

plot_trend_folder <- function(trend_folder) {

trend_name <- basename(trend_folder)

# list Cluster_* subfolders (direct children)
cl_dirs <- list.dirs(trend_folder, recursive = FALSE, full.names = TRUE)
cl_dirs <- cl_dirs[grepl("/Cluster_[0-9]+$", cl_dirs) | grepl("\\\\Cluster_[0-9]+$", cl_dirs)]

if (length(cl_dirs) == 0) {
message("[", trend_name, "] No Cluster_* folders found. Skipping.")
return(invisible(NULL))
}

# order cluster folders numerically
cl_base <- basename(cl_dirs)
cl_ids <- get_cl_id(cl_base)
ord <- order(cl_ids)
cl_dirs <- cl_dirs[ord]
cl_ids <- cl_ids[ord]

mats <- list()
all_vals <- c()
x_labels <- NULL

for (i in seq_along(cl_dirs)) {

cl_id <- cl_ids[i]
f_txt <- file.path(cl_dirs[i], paste0("Cluster_", cl_id, "_genes.txt"))

if (!file.exists(f_txt)) {
message("[", trend_name, "] Missing file: ", f_txt, " (skipping cluster ", cl_id, ")")
next
}

df <- read.table(f_txt, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

lfq_cols <- grep("LFQ_Healthy_|LFQ_Pre_RT_|LFQ_Post_RT_", colnames(df), value = TRUE)
if (length(lfq_cols) == 0) {
message("[", trend_name, "] No LFQ cols in cluster ", cl_id, " (skipping)")
next
}

if (is.null(x_labels)) x_labels <- lfq_cols

mat <- as.matrix(df[, lfq_cols, drop = FALSE])
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- df$PG.Genes

# keep only complete trajectories (same as your original plotting)
mat <- mat[complete.cases(mat), , drop = FALSE]
if (nrow(mat) < 1) {
message("[", trend_name, "] Cluster ", cl_id, " has 0 complete rows after NA filter (skipping)")
next
}

mats[[as.character(cl_id)]] <- mat
all_vals <- c(all_vals, as.numeric(mat))
}

if (length(mats) == 0) {
message("[", trend_name, "] No usable clusters. Skipping plot.")
return(invisible(NULL))
}

# global y-range for comparability across panels
global_y <- range(all_vals, na.rm = TRUE)

# panel layout
cluster_ids <- sort(as.integer(names(mats)))
k <- length(cluster_ids)
ncol_panels <- ceiling(sqrt(k))
nrow_panels <- ceiling(k / ncol_panels)

out_png <- file.path(trend_folder, paste0("SOTA_", trend_name, "_AllClusters_Multipanel.png"))

png(out_png, width = 1800, height = 1000, res = 150)
par(mfrow = c(nrow_panels, ncol_panels),
mar = c(7, 4, 3, 1),
oma = c(0, 0, 2, 0))

x <- seq_along(x_labels)

for (cl_id in cluster_ids) {

mat <- mats[[as.character(cl_id)]]
mu  <- colMeans(mat, na.rm = TRUE)
sdv <- apply(mat, 2, sd, na.rm = TRUE)

plot(x, rep(NA, length(x)),
type = "n",
ylim = global_y,
xaxt = "n",
xlab = "",
ylab = "LFQ",
main = paste0("Cluster ", cl_id, " (n=", nrow(mat), ")"))

# mean ± 1 SD ribbon
polygon(c(x, rev(x)), c(mu - sdv, rev(mu + sdv)),
border = NA, col = adjustcolor("red", alpha.f = 0.12))

# gene lines
apply(mat, 1, function(v)
lines(x, v, col = adjustcolor("grey30", alpha.f = 0.18), lwd = 1)
)

# mean line
lines(x, mu, col = "red", lwd = 3)

# x labels rotated
axis(1, at = x, labels = FALSE)
text(x,
par("usr")[3] - 0.03 * diff(par("usr")[3:4]),
labels = x_labels,
srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
}

mtext(paste0("SOTA trend: ", trend_name, " (from Cluster_* files)"),
outer = TRUE, cex = 1.1, line = 0.5)

dev.off()
message("[", trend_name, "] Saved: ", out_png)
}

# Run for BOTH trend folders (direct children of SOTA_Cluster)
trend_folders <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
for (tf in trend_folders) {
plot_trend_folder(tf)
}



################################################################################
################################################################################
################################################################################

################################################################################
# Welch t-test based cluster plots — MEAN + SEM VERSION
################################################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(openxlsx)
library(ggsignif)
library(patchwork)
library(scales)
library(tidyr)

# ------------- Load Data ----------------------
data_full <- read.table(
"C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt",
header = TRUE, sep = "\t", quote = "", check.names = FALSE
)

expr_cols <- grep(
"LFQ_Healthy_|LFQ_Pre_RT_|LFQ_Post_RT_",
colnames(data_full), value = TRUE
)

root_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/SOTA_Cluster"

# ------------- Main Function ------------------
process_trend_clusters <- function(trend_label) {

cat("\n==============================\n")
cat(" Processing TREND:", trend_label, "\n")
cat("==============================\n")

trend_folder <- file.path(root_dir, gsub(" ", "_", trend_label))

cluster_dirs <- list.dirs(trend_folder, recursive = FALSE, full.names = TRUE)
cluster_dirs <- cluster_dirs[grepl("Cluster_", basename(cluster_dirs))]

if (length(cluster_dirs) == 0) {
message("⚠️ No cluster folders found: ", trend_label)
return(NULL)
}

all_group_summaries <- list()
all_pvals <- list()

for (cl_dir in cluster_dirs) {

cluster_id <- basename(cl_dir)
cat("\n--- Cluster:", cluster_id, "\n")

gene_file <- file.path(cl_dir, paste0(cluster_id, "_genes.txt"))
if (!file.exists(gene_file)) {
message("⚠️ Missing gene list: ", gene_file)
next
}

genes <- read.table(gene_file, header = TRUE, sep = "\t")$PG.Genes
sub <- data_full[data_full$PG.Genes %in% genes, expr_cols, drop = FALSE]

if (nrow(sub) < 3) {
message("⚠️ Skipping ", cluster_id, " — not enough genes")
next
}

# ------------------------------------------------------------------
# MEAN (instead of MEDIAN)
# ------------------------------------------------------------------
sample_means <- apply(sub, 2, mean, na.rm = TRUE)

df_plot <- data.frame(
Sample = names(sample_means),
Value = as.numeric(sample_means)
) %>%
mutate(
Group = case_when(
grepl("Healthy", Sample) ~ "Healthy",
grepl("Pre_RT", Sample) ~ "Pre_RT",
grepl("Post_RT", Sample) ~ "Ende_RT",
TRUE ~ NA_character_
),
Group = factor(Group, levels = c("Healthy", "Pre_RT", "Ende_RT"))
)

# --- Welch t-tests unchanged ---
combs <- combn(levels(df_plot$Group), 2, simplify = FALSE)
test_results <- lapply(combs, function(pair) {
g1 <- df_plot$Value[df_plot$Group == pair[1]]
g2 <- df_plot$Value[df_plot$Group == pair[2]]
p <- tryCatch(t.test(g1, g2, var.equal = FALSE)$p.value, error = function(e) NA)
data.frame(Comparison = paste(pair, collapse = "-"), p_value = p)
})

t_df <- bind_rows(test_results) %>%
mutate(
p_adj = p.adjust(p_value, method = "BH"),
stars = case_when(
p_adj < 0.0001 ~ "****",
p_adj < 0.001 ~ "***",
p_adj < 0.01 ~ "**",
p_adj < 0.05 ~ "*",
TRUE ~ ""
),
Cluster = cluster_id
)

signif_pairs <- strsplit(t_df$Comparison, "-")
signif_labels <- t_df$stars
signif_pairs <- signif_pairs[signif_labels != ""]
signif_labels <- signif_labels[signif_labels != ""]

# ------------------------------------------------------------------
# GROUP SUMMARY — MEAN + SEM
# ------------------------------------------------------------------
group_summary <- df_plot %>%
group_by(Group) %>%
summarise(
Mean = mean(Value, na.rm = TRUE),
SEM = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
.groups = "drop"
) %>%
mutate(Cluster = cluster_id)

all_group_summaries[[cluster_id]] <- group_summary
all_pvals[[cluster_id]] <- t_df

# =================== VISUALIZATION ========================

base_y <- 0.2
step   <- 0.5
y_positions <- base_y + step * seq_along(signif_labels)

blank_df <- data.frame(
Group = factor(levels(df_plot$Group), levels = levels(df_plot$Group)),
dummy = 1
)

p_top <- ggplot(blank_df, aes(x = Group, y = dummy)) +
geom_blank() +
theme_void() +
coord_cartesian(clip = "off") +
theme(plot.margin = margin(1, 5, 0, 5))

if (length(signif_labels) > 0) {
p_top <- p_top +
geom_signif(
comparisons = signif_pairs,
annotations = signif_labels,
y_position = y_positions,
tip_length = 0.01,
textsize = 7,
vjust = 0.3
)
}

# ---- MAIN PLOT (MEAN + SEM) ----
p_main <- ggplot(df_plot, aes(Group, Value)) +
geom_jitter(width = 0.15, size = 10, alpha = 0.6, color = "#222222") +

geom_errorbar(
data = group_summary,
aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),
width = 0.18, size = 2.5, color = "#0072B2"
) +

geom_point(
data = group_summary,
aes(y = Mean),
size = 10, shape = 18, color = "#D55E00"
) +

geom_line(
data = group_summary,
aes(x = as.numeric(Group), y = Mean, group = 1),
size = 1.2, color = "#0072B2"
) +

theme_minimal(base_size = 36) +
labs(y = "LFQ Intensity (mean ± SEM)", x = NULL) +
theme(
plot.margin = margin(12, 12, 28, 12),
axis.text.x = element_text(size = 32),
axis.title.y = element_text(size = 32),
axis.line = element_line(color = "grey30", linewidth = 0.8),
panel.grid.minor = element_blank(),
panel.grid.major.x = element_blank()
)

full_plot <- (p_top / p_main + plot_layout(heights = c(0.05, 1))) +
plot_annotation(
title = paste0(trend_label, " – ", cluster_id),
theme = theme(
plot.title = element_text(
size = 32, face = "bold", hjust = 0.5,
margin = margin(b = 20)
),
plot.margin = margin(10, 12, 10, 12)
)
)

png_file <- file.path(cl_dir, paste0(cluster_id, "_Welch_TwoPanel_mean.png"))
pdf_file <- file.path(cl_dir, paste0(cluster_id, "_Welch_TwoPanel_mean.pdf"))
ggsave(png_file, full_plot, width = 14, height = 11, dpi = 300)
ggsave(pdf_file, full_plot, width = 14, height = 11, device = cairo_pdf)

# ---- EXPORT ----
write.table(group_summary,
file = file.path(cl_dir, paste0(cluster_id, "_mean_SEM_per_group.txt")),
sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(t_df,
file = file.path(cl_dir, paste0(cluster_id, "_Welch_pairwise_pvalues.txt")),
sep = "\t", row.names = FALSE, quote = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "Mean_SEM")
writeData(wb, "Mean_SEM", group_summary)
addWorksheet(wb, "Pairwise_pvalues")
writeData(wb, "Pairwise_pvalues", t_df)
saveWorkbook(wb, file = file.path(cl_dir, paste0(cluster_id, "_SummaryStats.xlsx")), overwrite = TRUE)
}

# ---- SUMMARY EXPORTS ----
all_group_df <- bind_rows(all_group_summaries)
all_pvals_df <- bind_rows(all_pvals)

write.xlsx(all_group_df, file.path(trend_folder, "AllClusters_Mean_SEM.xlsx"), overwrite = TRUE)
write.xlsx(all_pvals_df, file.path(trend_folder, "AllClusters_Pvalues.xlsx"), overwrite = TRUE)

cat("✅ Exported global summaries for:", trend_label, "\n")

# ---- MEAN LFQ PER GENE PER GROUP ----
lfq_long <- data_full[, c("PG.Genes", expr_cols)] %>%
pivot_longer(cols = -PG.Genes, names_to = "Sample", values_to = "LFQ") %>%
mutate(
Group = case_when(
grepl("Healthy", Sample) ~ "Healthy",
grepl("Pre_RT", Sample) ~ "Pre_RT",
grepl("Post_RT", Sample) ~ "Ende_RT",
TRUE ~ NA_character_
)
) %>%
filter(!is.na(Group))

lfq_mean <- lfq_long %>%
group_by(PG.Genes, Group) %>%
summarise(mean_LFQ = mean(as.numeric(LFQ), na.rm = TRUE), .groups = "drop") %>%
pivot_wider(names_from = Group, values_from = mean_LFQ)

write.xlsx(lfq_mean, file.path(trend_folder, "LFQ_mean_per_gene_per_group.xlsx"), overwrite = TRUE)
write.table(lfq_mean, file.path(trend_folder, "LFQ_mean_per_gene_per_group.txt"),
sep = "\t", row.names = FALSE)

cat("🎯 Exported MEAN LFQ per gene per group for trend:", trend_label, "\n")
}

# ------------- RUN ------------------
process_trend_clusters("Oncoprotein")
process_trend_clusters("Tumor suppressor protein")

cat("\n✔️ ALL Welch t-test cluster plots (mean version) completed.\n")

################################################################################

################################################################################
# Welch t-test based cluster plots — MEAN + SEM VERSION (NO HEALTHY)
################################################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(openxlsx)
library(ggsignif)
library(patchwork)
library(scales)
library(tidyr)

# ------------- Load Data ----------------------
data_full <- read.table(
"C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt",
header = TRUE, sep = "\t", quote = "", check.names = FALSE
)

# ONLY Pre and Post
expr_cols <- grep(
"LFQ_Pre_RT_|LFQ_Post_RT_",
colnames(data_full), value = TRUE
)

root_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/SOTA_Cluster"

# ------------- Main Function ------------------
process_trend_clusters <- function(trend_label) {

cat("\n==============================\n")
cat(" Processing TREND:", trend_label, "\n")
cat("==============================\n")

trend_folder <- file.path(root_dir, gsub(" ", "_", trend_label))

cluster_dirs <- list.dirs(trend_folder, recursive = FALSE, full.names = TRUE)
cluster_dirs <- cluster_dirs[grepl("Cluster_", basename(cluster_dirs))]

if (length(cluster_dirs) == 0) {
message("⚠️ No cluster folders found: ", trend_label)
return(NULL)
}

all_group_summaries <- list()
all_pvals <- list()

for (cl_dir in cluster_dirs) {

cluster_id <- basename(cl_dir)
cat("\n--- Cluster:", cluster_id, "\n")

gene_file <- file.path(cl_dir, paste0(cluster_id, "_genes.txt"))
if (!file.exists(gene_file)) {
message("⚠️ Missing gene list: ", gene_file)
next
}

genes <- read.table(gene_file, header = TRUE, sep = "\t")$PG.Genes
sub <- data_full[data_full$PG.Genes %in% genes, expr_cols, drop = FALSE]

if (nrow(sub) < 3) {
message("⚠️ Skipping ", cluster_id, " — not enough genes")
next
}

# ------------------------------------------------------------------
# MEAN (instead of MEDIAN)
# ------------------------------------------------------------------
sample_means <- apply(sub, 2, mean, na.rm = TRUE)

df_plot <- data.frame(
Sample = names(sample_means),
Value  = as.numeric(sample_means)
) %>%
mutate(
Group = case_when(
grepl("Pre_RT",  Sample) ~ "Pre_RT",
grepl("Post_RT", Sample) ~ "Ende_RT",
TRUE ~ NA_character_
),
Group = factor(Group, levels = c("Pre_RT", "Ende_RT"))
) %>%
filter(!is.na(Group))

# --- Welch t-tests (ONLY Pre vs Ende now) ---
combs <- combn(levels(df_plot$Group), 2, simplify = FALSE)
test_results <- lapply(combs, function(pair) {
g1 <- df_plot$Value[df_plot$Group == pair[1]]
g2 <- df_plot$Value[df_plot$Group == pair[2]]
p  <- tryCatch(t.test(g1, g2, var.equal = FALSE)$p.value, error = function(e) NA)
data.frame(Comparison = paste(pair, collapse = "-"), p_value = p)
})

t_df <- bind_rows(test_results) %>%
mutate(
p_adj = p.adjust(p_value, method = "BH"),
stars = case_when(
p_adj < 0.0001 ~ "****",
p_adj < 0.001  ~ "***",
p_adj < 0.01   ~ "**",
p_adj < 0.05   ~ "*",
TRUE ~ ""
),
Cluster = cluster_id
)

signif_pairs  <- strsplit(t_df$Comparison, "-")
signif_labels <- t_df$stars
signif_pairs  <- signif_pairs[signif_labels != ""]
signif_labels <- signif_labels[signif_labels != ""]

# ------------------------------------------------------------------
# GROUP SUMMARY — MEAN + SEM
# ------------------------------------------------------------------
group_summary <- df_plot %>%
group_by(Group) %>%
summarise(
Mean = mean(Value, na.rm = TRUE),
SEM  = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
.groups = "drop"
) %>%
mutate(Cluster = cluster_id)

all_group_summaries[[cluster_id]] <- group_summary
all_pvals[[cluster_id]] <- t_df

# =================== VISUALIZATION ========================

base_y <- 0.2
step   <- 0.5
y_positions <- base_y + step * seq_along(signif_labels)

blank_df <- data.frame(
Group = factor(levels(df_plot$Group), levels = levels(df_plot$Group)),
dummy = 1
)

p_top <- ggplot(blank_df, aes(x = Group, y = dummy)) +
geom_blank() +
theme_void() +
coord_cartesian(clip = "off") +
theme(plot.margin = margin(1, 5, 0, 5))

if (length(signif_labels) > 0) {
p_top <- p_top +
geom_signif(
comparisons = signif_pairs,
annotations = signif_labels,
y_position  = y_positions,
tip_length  = 0.01,
textsize    = 20,
vjust       = 0.3
)
}

# ---- MAIN PLOT (MEAN + SEM) ----
p_main <- ggplot(df_plot, aes(Group, Value)) +
geom_jitter(width = 0.15, size = 10, alpha = 0.6, color = "#222222") +

geom_errorbar(
data = group_summary,
aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),
width = 0.18, size = 5, color = "#0072B2"
) +

geom_point(
data = group_summary,
aes(y = Mean),
size = 10, shape = 18, color = "#D55E00"
) +

geom_line(
data = group_summary,
aes(x = as.numeric(Group), y = Mean, group = 1),
size = 1.2, color = "#0072B2"
) +

theme_minimal(base_size = 36) +
labs(y = "LFQ Intensity (mean ± SEM)", x = NULL) +
theme(
plot.margin = margin(12, 12, 28, 12),
axis.text.x = element_text(size = 36, color = "#222222"),
axis.title.y = element_text(size = 36, color = "#222222"),
axis.line = element_line(color = "#222222", linewidth = 0.8),
panel.grid.minor = element_blank(),
panel.grid.major.x = element_blank()
)

full_plot <- (p_top / p_main + plot_layout(heights = c(0.05, 1))) +
plot_annotation(
title = paste0(trend_label, " – ", cluster_id),
theme = theme(
plot.title = element_text(
size = 32, face = "bold", hjust = 0.5,
margin = margin(b = 20)
),
plot.margin = margin(10, 12, 10, 12)
)
)

# ---- SAVE with NEW NAMES ----
png_file <- file.path(cl_dir, paste0(cluster_id, "_Welch_TwoPanel_mean_NO_HEALTHY.png"))
pdf_file <- file.path(cl_dir, paste0(cluster_id, "_Welch_TwoPanel_mean_NO_HEALTHY.pdf"))
ggsave(png_file, full_plot, width = 14, height = 11, dpi = 300)
ggsave(pdf_file, full_plot, width = 14, height = 11, device = cairo_pdf)

# ---- EXPORT with NEW NAMES ----
write.table(
group_summary,
file = file.path(cl_dir, paste0(cluster_id, "_mean_SEM_per_group_NO_HEALTHY.txt")),
sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
t_df,
file = file.path(cl_dir, paste0(cluster_id, "_Welch_pairwise_pvalues_NO_HEALTHY.txt")),
sep = "\t", row.names = FALSE, quote = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "Mean_SEM")
writeData(wb, "Mean_SEM", group_summary)
addWorksheet(wb, "Pairwise_pvalues")
writeData(wb, "Pairwise_pvalues", t_df)
saveWorkbook(
wb,
file = file.path(cl_dir, paste0(cluster_id, "_SummaryStats_NO_HEALTHY.xlsx")),
overwrite = TRUE
)
}

# ---- SUMMARY EXPORTS (NEW NAMES) ----
all_group_df <- bind_rows(all_group_summaries)
all_pvals_df <- bind_rows(all_pvals)

write.xlsx(all_group_df, file.path(trend_folder, "AllClusters_Mean_SEM_NO_HEALTHY.xlsx"), overwrite = TRUE)
write.xlsx(all_pvals_df, file.path(trend_folder, "AllClusters_Pvalues_NO_HEALTHY.xlsx"), overwrite = TRUE)

cat("✅ Exported global summaries for:", trend_label, "(NO HEALTHY)\n")

# ---- MEAN LFQ PER GENE PER GROUP (NO HEALTHY) ----
lfq_long <- data_full[, c("PG.Genes", expr_cols)] %>%
pivot_longer(cols = -PG.Genes, names_to = "Sample", values_to = "LFQ") %>%
mutate(
Group = case_when(
grepl("Pre_RT",  Sample) ~ "Pre_RT",
grepl("Post_RT", Sample) ~ "Ende_RT",
TRUE ~ NA_character_
)
) %>%
filter(!is.na(Group))

lfq_mean <- lfq_long %>%
group_by(PG.Genes, Group) %>%
summarise(mean_LFQ = mean(as.numeric(LFQ), na.rm = TRUE), .groups = "drop") %>%
pivot_wider(names_from = Group, values_from = mean_LFQ)

write.xlsx(lfq_mean, file.path(trend_folder, "LFQ_mean_per_gene_per_group_NO_HEALTHY.xlsx"), overwrite = TRUE)
write.table(
lfq_mean,
file.path(trend_folder, "LFQ_mean_per_gene_per_group_NO_HEALTHY.txt"),
sep = "\t", row.names = FALSE
)

cat("🎯 Exported MEAN LFQ per gene per group for trend:", trend_label, "(NO HEALTHY)\n")
}

# ------------- RUN ------------------
process_trend_clusters("Oncoprotein")
process_trend_clusters("Tumor suppressor protein")

cat("\n✔️ ALL Welch t-test cluster plots (mean version, NO HEALTHY) completed.\n")


################################################################################
################################################################################
#############################  STEP 4 — ANNOTATION  ############################
################################################################################
# Annotate each SOTA cluster with:
#   - GO Biological terms
#   - KEGG pathways
#   - Reactome pathways
#
# Input:
#   SOTA_Cluster/<Trend>/Cluster_X/Cluster_X_genes.txt
#
# Output:
#   Cluster_X_genes_annotated.txt
#   Cluster_X_genes_annotated.xlsx
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(ReactomePA)
library(clusterProfiler)
library(KEGGREST)
library(biomaRt)
library(openxlsx)
})

# Root folder produced by Step 1–3
root_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/SOTA_Cluster"

# Global normalized expression file
norm_file <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics/DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
norm_df <- read.table(norm_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

################################################################################
### Helper functions (all fixed to avoid select() ambiguity)
################################################################################

split_genes <- function(gene_list) {
data.frame(Original = gene_list) %>%
mutate(Gene_Single = strsplit(Original, ";")) %>%
tidyr::unnest(Gene_Single) %>%
mutate(Gene_Single = trimws(Gene_Single))
}

map_to_entrez <- function(symbols) {

mapping <- tryCatch(
AnnotationDbi::select(org.Hs.eg.db,
keys = symbols,
columns = c("ENTREZID","UNIPROT"),
keytype = "SYMBOL"),
error = function(e) data.frame()
)

mapping <- mapping[!duplicated(mapping$SYMBOL), ]

# Fallback: UniProt
if (nrow(mapping)==0 || all(is.na(mapping$ENTREZID))) {
uni_map <- tryCatch(
clusterProfiler::bitr(symbols, fromType="UNIPROT",
toType=c("SYMBOL","ENTREZID"),
OrgDb="org.Hs.eg.db"),
error=function(e) NULL)

if (!is.null(uni_map) && nrow(uni_map)>0)
mapping <- uni_map
}

# Fallback: biomaRt
if (any(is.na(mapping$ENTREZID))) {
missing <- mapping$SYMBOL[is.na(mapping$ENTREZID)]

mart <- tryCatch(useEnsembl("genes","hsapiens_gene_ensembl"),
error=function(e) NULL)

if (!is.null(mart) && length(missing)>0) {
bm <- tryCatch(
getBM(attributes=c("hgnc_symbol","entrezgene_id"),
filters="hgnc_symbol", values=missing, mart=mart),
error=function(e) NULL)

if (!is.null(bm) && nrow(bm)>0) {
for (i in seq_len(nrow(bm))) {
mapping$ENTREZID[mapping$SYMBOL == bm$hgnc_symbol[i]] <-
as.character(bm$entrezgene_id[i])
}
}
}
}

mapping <- mapping[!is.na(mapping$ENTREZID) & mapping$ENTREZID != "", ]
unique(mapping)
}


safe_kegg <- function(entrez_ids, gene_entrez) {
entrez_ids <- as.character(entrez_ids)

# Try KEGG enrichment first
k <- tryCatch(
enrichKEGG(gene=entrez_ids, organism="hsa", keyType="ncbi-geneid",
pvalueCutoff=1, qvalueCutoff=1),
error=function(e) NULL)

if (!is.null(k) && nrow(as.data.frame(k))>0) {
df <- as.data.frame(k)

return(
df %>%
dplyr::select(KEGG_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep="/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by="ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(KEGG_Pathway = paste(unique(KEGG_Pathway), collapse="; "), .groups="drop")
)
}

# KEGGREST fallback
path_list <- KEGGREST::keggLink("pathway","hsa")
path_list <- path_list[gsub("hsa:","", names(path_list)) %in% entrez_ids]

if (length(path_list)==0) return(data.frame())

desc <- KEGGREST::keggList("pathway","hsa")

df <- data.frame(
ENTREZID = sub("hsa:","", names(path_list)),
PATHWAY  = sub("path:","", path_list)
)

df <- merge(df,
data.frame(PATHWAY=sub("path:","",names(desc)),
KEGG_Pathway=desc),
by="PATHWAY")

df %>%
dplyr::left_join(gene_entrez, by="ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(KEGG_Pathway = paste(unique(KEGG_Pathway), collapse="; "),
.groups="drop")
}

################################################################################
### MAIN LOOP — annotate all clusters
################################################################################

trends <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)

for (trend_folder in trends) {

clusters <- list.dirs(trend_folder, recursive = FALSE, full.names = TRUE)
clusters <- clusters[grepl("Cluster_", basename(clusters))]

for (cl_dir in clusters) {

cl_name <- basename(cl_dir)
gene_file <- file.path(cl_dir, paste0(cl_name, "_genes.txt"))

if (!file.exists(gene_file)) {
message("Skipping missing file: ", gene_file)
next
}

message("🔵 Annotating ", trend_folder, " / ", cl_name)

df <- read.table(gene_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

if (!"PG.Genes" %in% colnames(df)) {
message("❌ File does not contain PG.Genes: ", gene_file)
next
}

genes <- unique(na.omit(df$PG.Genes))
map <- split_genes(genes)

symbols <- unique(map$Gene_Single)
gene_entrez <- map_to_entrez(symbols)

entrez <- unique(gene_entrez$ENTREZID)

if (length(entrez) < 2) {
message("⚠️ Too few Entrez IDs for cluster ", cl_name)
next
}

# ----------------------------- GO ----------------------------------------
go <- tryCatch({
go_anno <- AnnotationDbi::select(org.Hs.eg.db,
keys = entrez,
columns = c("GO","ONTOLOGY"),
keytype="ENTREZID")

go_terms <- AnnotationDbi::select(GO.db,
keys = unique(go_anno$GO),
columns="TERM",
keytype="GOID")

gm <- merge(go_anno, go_terms, by.x="GO", by.y="GOID", all.x=TRUE)
gm <- dplyr::left_join(gm, gene_entrez, by="ENTREZID")

gm$GO_Pathway <- paste(gm$GO, gm$ONTOLOGY, gm$TERM, sep=" | ")

gm %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(GO_Pathway = paste(unique(GO_Pathway), collapse="; "), .groups="drop")

}, error=function(e) data.frame())

# ----------------------------- KEGG --------------------------------------
kegg <- safe_kegg(entrez, gene_entrez)

# --------------------------- Reactome ------------------------------------
rea <- tryCatch({
enr <- ReactomePA::enrichPathway(
gene = as.numeric(entrez),
organism = "human",
pvalueCutoff = 1, qvalueCutoff = 1
)

df_r <- as.data.frame(enr)

df_r %>%
dplyr::select(Reactome_Pathway = Description, geneID) %>%
tidyr::separate_rows(geneID, sep="/") %>%
dplyr::rename(ENTREZID = geneID) %>%
dplyr::left_join(gene_entrez, by="ENTREZID") %>%
dplyr::group_by(SYMBOL) %>%
dplyr::summarise(Reactome_Pathway = paste(unique(Reactome_Pathway), collapse="; "),
.groups="drop")

}, error=function(e) data.frame())

# -------------------- MERGE annotations back ------------------------------
ann <- map %>%
dplyr::left_join(go,   by=c("Gene_Single"="SYMBOL")) %>%
dplyr::left_join(kegg, by=c("Gene_Single"="SYMBOL")) %>%
dplyr::left_join(rea,  by=c("Gene_Single"="SYMBOL")) %>%
dplyr::group_by(Original) %>%
dplyr::summarise(
GO_Pathway       = paste(na.omit(unique(GO_Pathway)), collapse="; "),
KEGG_Pathway     = paste(na.omit(unique(KEGG_Pathway)), collapse="; "),
Reactome_Pathway = paste(na.omit(unique(Reactome_Pathway)), collapse="; "),
.groups="drop"
)

annotated <- dplyr::left_join(df, ann, by=c("PG.Genes"="Original"))

# ------------------------------- SAVE -------------------------------------
out_txt  <- file.path(cl_dir, paste0(cl_name, "_genes_annotated.txt"))
out_xlsx <- file.path(cl_dir, paste0(cl_name, "_genes_annotated.xlsx"))

write.table(annotated, out_txt, sep="\t", quote=FALSE, row.names=FALSE)
write.xlsx(annotated, out_xlsx, overwrite = TRUE)

message("💾 Saved: ", out_txt)
}
}

cat("\n🎉 STEP 4 COMPLETED — all clusters annotated.\n")


################################################################################
############################## 🔵 STEP 5 — FISHER TEST #########################
################################################################################
# Fisher Exact Test Enrichment for Step 4 cluster-level annotations
#
# This script:
# - Uses a GLOBAL background file:
#       DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt
# - Processes annotated cluster files:
#       SOTA_Cluster/<Trend>/Cluster_X/Cluster_X_genes_annotated.txt
# - Performs Fisher Exact Test on:
#       GO_Pathway, KEGG_Pathway, Reactome_Pathway
# - Saves (inside each Cluster_X folder):
#       Fisher_GO_Cluster_X_genes_annotated.xlsx
#       Fisher_KEGG_Cluster_X_genes_annotated.xlsx
#       Fisher_Reactome_Cluster_X_genes_annotated.xlsx
#       Merged_Fisher_Cluster_X_genes_annotated.xlsx
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
})

# ==============================================================================
# 1. Define main directories
# ==============================================================================

main_dir  <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
root_dir  <- file.path(main_dir, "SOTA_Cluster")

cat("🏠 Main directory:\n", main_dir, "\n")
cat("📂 SOTA root:\n", root_dir, "\n")

# ==============================================================================
# 2. Load GLOBAL background file (same logic as Step 16)
# ==============================================================================

global_background_file <- file.path(
main_dir,
"DKFZ_HNSCC_patient_EV_filtered_50percent_corrFiltered_quantile_normalized.txt"
)

if (!file.exists(global_background_file)) {
stop("❌ Global normalized file not found: ", global_background_file)
}

global_norm_df <- read.table(global_background_file, header = TRUE, sep = "\t",
quote = "", check.names = FALSE)

if ("PG.Genes" %in% colnames(global_norm_df)) {
background_genes <- unique(trimws(as.character(global_norm_df$PG.Genes)))
} else {
background_genes <- unique(trimws(as.character(global_norm_df$PG.ProteinGroups)))
}

background_genes <- background_genes[background_genes != "" & !is.na(background_genes)]
cat("✅ Using GLOBAL background:", length(background_genes), "proteins/genes\n")

# ==============================================================================
# 3. Fisher Test Function (identical logic to Step 16, reused)
# ==============================================================================

fisher_on_annotation <- function(ann_df, annot_col, background_genes, out_xlsx) {

if (!(annot_col %in% colnames(ann_df))) return(invisible(NULL))
if (!("PG.Genes" %in% colnames(ann_df))) {
cat("❌ Missing PG.Genes column; skipping.\n")
return(invisible(NULL))
}

# Expand annotation column into individual pathway records
df <- ann_df %>%
filter(!is.na(.data[[annot_col]]) & .data[[annot_col]] != "") %>%
transmute(PG.Genes, Pathway = .data[[annot_col]]) %>%
mutate(Pathway = strsplit(as.character(Pathway), ";")) %>%
unnest(Pathway) %>%
mutate(Pathway = str_trim(Pathway)) %>%
filter(!is.na(Pathway) & Pathway != "" & Pathway != "NA")

# Remove pure "NA | NA | NA" style junk
df <- df[!grepl("^NA\\s*\\|\\s*NA\\s*\\|\\s*NA$", df$Pathway, ignore.case = TRUE), ]

# Pathway name simplification for KEGG-type names
df$Pathway_Name <- ifelse(
grepl(" - Homo sapiens", df$Pathway),
sub(" - Homo sapiens.*", "", df$Pathway),
df$Pathway
)

# Foreground (cluster) genes
sig_genes <- unique(trimws(as.character(ann_df$PG.Genes)))
sig_genes <- sig_genes[sig_genes != "" & !is.na(sig_genes)]

terms <- unique(df$Pathway_Name)

if (length(terms) == 0 || length(sig_genes) == 0 || length(background_genes) == 0)
return(invisible(NULL))

out_list <- list()

for (term in terms) {

# Genes in this pathway that are in the background
term_genes <- unique(df$PG.Genes[
df$Pathway_Name == term & df$PG.Genes %in% background_genes
])

a <- length(intersect(term_genes, sig_genes))                             # in term & in cluster
b <- length(setdiff(sig_genes, term_genes))                               # not in term & in cluster
c <- length(setdiff(term_genes, intersect(term_genes, sig_genes)))        # in term & not in cluster
d <- length(setdiff(background_genes, union(sig_genes, term_genes)))      # not in term & not in cluster

if (any(c(a,b,c,d) < 0)) next

ft <- tryCatch(fisher.test(matrix(c(a,b,c,d), nrow = 2)), error = function(e) NULL)
if (is.null(ft)) next

gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
bg_ratio   <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)

out_list[[term]] <- data.frame(
Pathway              = term,
Sig_protein_volcano  = a,
Sig_NotIn_Pathway    = b,
Nonsig_In_Pathway    = c,
Nonsig_NotIn_Pathway = d,
Genes                = paste(intersect(term_genes, sig_genes), collapse = ";"),
p_value              = ft$p.value,
enrichment_factor    = enrichment_factor,
stringsAsFactors     = FALSE
)
}

res <- dplyr::bind_rows(out_list)
if (is.null(res) || nrow(res) == 0) return(invisible(NULL))

res$p_adj <- p.adjust(res$p_value, method = "BH")
res <- res %>% arrange(p_adj)

openxlsx::write.xlsx(res, out_xlsx, overwrite = TRUE)
cat("💾 Saved:", basename(out_xlsx), "(", nrow(res), "terms)\n")

invisible(res)
}

# ==============================================================================
# 4. Main Loop Over SOTA Trend Folders and Cluster Folders
# ==============================================================================

trend_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)

for (trend_folder in trend_dirs) {

cat("\n==================================================================\n")
cat("📂 Trend folder:", trend_folder, "\n")
cat("==================================================================\n")

cluster_dirs <- list.dirs(trend_folder, recursive = FALSE, full.names = TRUE)
cluster_dirs <- cluster_dirs[grepl("Cluster_", basename(cluster_dirs))]

if (!length(cluster_dirs)) {
cat("⚠️ No Cluster_* folders found in", trend_folder, "→ skipping.\n")
next
}

for (cl_dir in cluster_dirs) {

cl_name <- basename(cl_dir)
annot_file <- file.path(cl_dir, paste0(cl_name, "_genes_annotated.txt"))

if (!file.exists(annot_file)) {
cat("⚠️ No annotated file found for", cl_name, "→ expected:", basename(annot_file), "\n")
next
}

cat("🔬 Processing cluster:", cl_name, "\n")
ann_df <- read.table(
annot_file, header = TRUE, sep = "\t",
quote = "", stringsAsFactors = FALSE, check.names = FALSE
)

base <- tools::file_path_sans_ext(basename(annot_file))

# ------------------------ Run Fisher per annotation type ------------------
res_go <- res_kegg <- res_rea <- NULL

if ("GO_Pathway" %in% colnames(ann_df)) {
res_go <- fisher_on_annotation(
ann_df, "GO_Pathway", background_genes,
file.path(cl_dir, paste0("Fisher_GO_", base, ".xlsx"))
)
}

if ("KEGG_Pathway" %in% colnames(ann_df)) {
res_kegg <- fisher_on_annotation(
ann_df, "KEGG_Pathway", background_genes,
file.path(cl_dir, paste0("Fisher_KEGG_", base, ".xlsx"))
)
}

if ("Reactome_Pathway" %in% colnames(ann_df)) {
res_rea <- fisher_on_annotation(
ann_df, "Reactome_Pathway", background_genes,
file.path(cl_dir, paste0("Fisher_Reactome_", base, ".xlsx"))
)
}

# ------------------------ Merge GO + KEGG + Reactome ---------------------
merged_df <- dplyr::bind_rows(
if (!is.null(res_go))   dplyr::mutate(res_go,   Pathway_Type = "GO"),
if (!is.null(res_kegg)) dplyr::mutate(res_kegg, Pathway_Type = "KEGG"),
if (!is.null(res_rea))  dplyr::mutate(res_rea,  Pathway_Type = "Reactome")
)

if (!is.null(merged_df) && nrow(merged_df) > 0) {

merged_df <- merged_df %>%
arrange(p_adj, Pathway_Type, Pathway)

out_merged <- file.path(
cl_dir,
paste0("Merged_Fisher_", base, ".xlsx")
)

openxlsx::write.xlsx(merged_df, out_merged, overwrite = TRUE)
cat("📘 Merged Fisher (GO+KEGG+Reactome) for", cl_name, "→",
basename(out_merged), " (", nrow(merged_df), "rows)\n")
} else {
cat("⚠️ No enriched terms for", cl_name, " (all annotation types empty).\n")
}
}

cat("✅ Fisher enrichment complete for trend folder:", basename(trend_folder), "\n")
}

cat("\n🎉 STEP 5 complete — Fisher results saved inside each Cluster_X folder.\n")
################################################################################
################################################################################
################################################################################


################################################################################
# STEP 6 — CLUSTER FUNCTIONAL RANKING: All clusters (filtered by p-value logic)
# Uses separate files: AllClusters_Mean_SEM.xlsx and AllClusters_Pvalues.xlsx
# - Includes all clusters, marks filters, ranks clusters passing all criteria
# - No "sheet 1/sheet 2" logic, direct files only
################################################################################

suppressPackageStartupMessages({
library(dplyr)
library(openxlsx)
library(tidyr)
library(stringr)
})

main_dir  <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
trend_types <- c("Oncoprotein", "Tumor_suppressor_protein")
pvalue_cutoff <- 0.05
p_adj_cutoff <- 0.01
min_proteins <- 2

clusters_keywords <- list(
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


score_cluster <- function(fisher_df, category_keywords) {
term_names <- tolower(fisher_df$Pathway)
sapply(category_keywords, function(keys) {
sum(sapply(keys, function(k) any(grepl(k, term_names, ignore.case = TRUE))))
})
}

for (trend_label in trend_types) {

trend_dir <- file.path(main_dir, "SOTA_Cluster", trend_label)
mean_sem_file <- file.path(trend_dir, "AllClusters_Mean_SEM.xlsx")
pval_file     <- file.path(trend_dir, "AllClusters_Pvalues.xlsx")

if (!file.exists(mean_sem_file) || !file.exists(pval_file)) {
cat("Missing summary file(s) for:", trend_label, "\n")
next
}

# Read data
mean_sem <- openxlsx::read.xlsx(mean_sem_file)
pvals    <- openxlsx::read.xlsx(pval_file)

cluster_names <- unique(mean_sem$Cluster)
all_stats <- list()
cluster_hits <- list()

for (cl_name in cluster_names) {
# --- Step 1: FILTERING by significance ---
step1_status <- "Passed"
step2_status <- "Passed"
filter_reason <- NA

# Extract p-values for this cluster
cl_pvals <- pvals %>% filter(Cluster == cl_name)
p_HB <- cl_pvals %>% filter(Comparison == "Healthy-Pre_RT") %>% pull(p_adj)
p_BE <- cl_pvals %>% filter(Comparison == "Pre_RT-Ende_RT") %>% pull(p_adj)

# Remove cluster if either comparison is not significant
if (length(p_HB)==0 || is.na(p_HB) || p_HB >= pvalue_cutoff) {
step2_status <- "Removed"
filter_reason <- "Healthy-Pre_RT not significant"
} else if (length(p_BE)==0 || is.na(p_BE) || p_BE >= pvalue_cutoff) {
step2_status <- "Removed"
filter_reason <- "Pre_RT-Ende_RT not significant"
}

# --- RANKING & SCORING ---
cat_score <- NA
ProtWeight <- NA
weighted_score <- NA
TotalScore <- NA
Proteins_in_Pathways <- NA
Pathways_Selected <- NA
Ratio <- NA
Best_Pathway <- NA
Best_p_adj <- NA
prot_vec <- NA
n_prots <- NA
best_idx <- NA

# Fisher file check and scoring if passed
cl_dir <- file.path(trend_dir, cl_name)
fisher_file <- file.path(cl_dir, paste0("Merged_Fisher_", cl_name, "_genes_annotated.xlsx"))
if (!file.exists(fisher_file)) {
step1_status <- "Removed"
filter_reason <- "No Fisher file"
}

# Only score if both steps passed
if (step1_status == "Passed" && step2_status == "Passed") {
fisher_df <- openxlsx::read.xlsx(fisher_file)
filtered <- fisher_df %>%
filter(p_adj < p_adj_cutoff, Sig_protein_volcano >= min_proteins)
if (nrow(filtered) > 0) {
cat_score <- score_cluster(filtered, clusters_keywords)
ProtWeight <- mean(filtered$Sig_protein_volcano)
weighted_score <- cat_score * ProtWeight
TotalScore <- sum(weighted_score)
prot_vec <- unlist(strsplit(filtered$Genes, ";"))
prot_vec <- trimws(prot_vec)
prot_vec <- prot_vec[prot_vec != ""]
cluster_hits[[cl_name]] <- table(prot_vec)
Pathways_Selected <- nrow(filtered)
Proteins_in_Pathways <- length(unique(names(cluster_hits[[cl_name]])))
Ratio <- Proteins_in_Pathways / Pathways_Selected
best_idx <- which.min(filtered$p_adj)
Best_Pathway <- filtered$Pathway[best_idx]
Best_p_adj <- filtered$p_adj[best_idx]
} else {
step1_status <- "Removed"
filter_reason <- "No significant pathways"
}
}

# Number of proteins
gene_file <- file.path(cl_dir, paste0(cl_name, "_genes.txt"))
n_prots <- NA
if (file.exists(gene_file)) {
cluster_genes <- tryCatch(read.table(gene_file, header=TRUE, sep="\t")$PG.Genes, error=function(e) character(0))
n_prots <- length(unique(cluster_genes))
}

all_stats[[cl_name]] <- data.frame(
Cluster = cl_name,
Step1_Filtered = step1_status,
Step2_Filtered = step2_status,
Reason = filter_reason,
Total_Proteins = n_prots,
Proteins_in_Pathways = Proteins_in_Pathways,
Pathways_Selected = Pathways_Selected,
Ratio = Ratio,
Best_Pathway = Best_Pathway,
Best_p_adj = Best_p_adj,
Total_Score = TotalScore,
stringsAsFactors = FALSE
)
}

ranking <- bind_rows(all_stats)
ranking <- ranking %>%
arrange(desc(Total_Score)) %>%
mutate(Rank = ifelse(!is.na(Total_Score), rank(-Total_Score, ties.method = "first"), NA))

best_cluster <- NA
if (any(!is.na(ranking$Rank))) {
best_cluster <- ranking$Cluster[which.min(ranking$Rank)]
cat("\n🏆 BEST CLUSTER for", trend_label, "=", best_cluster, "\n\n")
}

# Export all clusters (filter columns and explanations)
ranking_file <- file.path(trend_dir, paste0(trend_label, "_Cluster_Functional_Ranking_ALL.xlsx"))
openxlsx::write.xlsx(ranking, ranking_file, overwrite = TRUE)

# TXT full report
report_txt <- file.path(trend_dir, paste0("Full_Report_ALL_", trend_label, ".txt"))
lines <- c(
"===============================================================",
paste0(" FULL REPORT (all clusters) — ", trend_label),
"===============================================================",
""
)
for (i in seq_len(nrow(ranking))) {
r <- ranking[i, ]
lines <- c(lines,
paste0(i, ". Cluster: ", r$Cluster),
paste0("   Step1_Filtered: ", r$Step1_Filtered),
paste0("   Step2_Filtered: ", r$Step2_Filtered),
paste0("   Reason: ", ifelse(is.na(r$Reason),"",r$Reason)),
paste0("   Total proteins: ", ifelse(is.na(r$Total_Proteins),"",r$Total_Proteins)),
paste0("   Proteins in pathways: ", ifelse(is.na(r$Proteins_in_Pathways),"",r$Proteins_in_Pathways)),
paste0("   # enriched pathways: ", ifelse(is.na(r$Pathways_Selected),"",r$Pathways_Selected)),
paste0("   Ratio: ", ifelse(is.na(r$Ratio),"",round(as.numeric(r$Ratio), 3))),
paste0("   Best pathway: ", ifelse(is.na(r$Best_Pathway),"",r$Best_Pathway)),
paste0("   Best p_adj: ", ifelse(is.na(r$Best_p_adj),"",signif(as.numeric(r$Best_p_adj), 4))),
paste0("   Total Score: ", ifelse(is.na(r$Total_Score),"",round(as.numeric(r$Total_Score), 2))),
paste0("   Rank: ", ifelse(is.na(r$Rank),"",r$Rank)),
""
)
}
writeLines(lines, con = report_txt)

# Excel report only for clusters that passed both filters
report_xlsx <- file.path(trend_dir, paste0("Full_Report_", trend_label, ".xlsx"))
openxlsx::write.xlsx(ranking %>% filter(Step1_Filtered=="Passed", Step2_Filtered=="Passed"), 
report_xlsx, overwrite = TRUE)

# Top 20 proteins for best cluster (if any)
if (!is.na(best_cluster) && !is.null(cluster_hits[[best_cluster]])) {
best_hits <- sort(cluster_hits[[best_cluster]], decreasing = TRUE)
top20_df <- data.frame(
Protein = names(best_hits)[1:min(20, length(best_hits))],
Count = as.numeric(best_hits)[1:min(20, length(best_hits))]
)
openxlsx::write.xlsx(
top20_df,
file.path(trend_dir, paste0("Top20_Proteins_", best_cluster, ".xlsx")),
overwrite = TRUE
)
cat("⭐ Top 20 proteins saved for:", best_cluster, "\n")
}

# Top 0 proteins for best cluster (if any)
if (!is.na(best_cluster) && !is.null(cluster_hits[[best_cluster]])) {
best_hits <- sort(cluster_hits[[best_cluster]], decreasing = TRUE)
top10_df <- data.frame(
Protein = names(best_hits)[1:min(10, length(best_hits))],
Count = as.numeric(best_hits)[1:min(10, length(best_hits))]
)
openxlsx::write.xlsx(
top10_df,
file.path(trend_dir, paste0("Top10_Proteins_", best_cluster, ".xlsx")),
overwrite = TRUE
)
cat("Top 10 proteins saved for:", best_cluster, "\n")
}
}



cat("\n🎉 STEP 6 COMPLETED — All clusters report with double p-value filter generated.\n")









################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
######################## EV PROTEINS VENN DIAGRAM (FINAL) ######################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
# EV Proteins Venn Diagram — CLEAN FULL VERSION WITH 0Gy + 6Gy FOR 3 CELL LINES
# + Patient SOTA includes Cluster 1 (Oncoprotein) AND Cluster 5 (Tumor suppressor)
################################################################################

base_dir      <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
cellular_root <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics"
output_dir    <- file.path(base_dir, "Venn_Diagrams_individual_cell_line_merge")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
library(openxlsx)
library(VennDiagram)
library(grid)
})

#--------------------------- GENE FILE READER ----------------------------------
read_gene_file <- function(file_path) {
df <- read.table(file_path, sep = "\t", header = TRUE, quote = "",
stringsAsFactors = FALSE, check.names = FALSE)

candidates <- c("Gene", "PG.Genes", "Gene.names", "PG.ProteinGroups")
col_found <- candidates[candidates %in% colnames(df)][1]

if (is.na(col_found)) {
if (ncol(df) == 1) {
col_found <- colnames(df)[1]
message("⚠ No gene column detected, using first column: ", col_found)
} else {
stop("No usable gene column in: ", file_path)
}
}

genes <- as.character(df[[col_found]])
genes <- sub(";.*", "", genes)
unique(trimws(genes[genes != "" & !is.na(genes)]))
}

#--------------------------- LOAD 0 Gy + 6 Gy FOR EACH CELL LINE ----------------
load_condition <- function(cell_line, dose_pattern) {
folder <- file.path(cellular_root, cell_line)
file <- list.files(folder, pattern = paste0(dose_pattern, ".*annotated\\.txt$"),
full.names = TRUE)

if (length(file) == 0)
stop("❌ No ", dose_pattern, " annotated file found in: ", folder)

if (length(file) > 1)
message("⚠ Multiple files detected in ", folder, " — using: ", basename(file[1]))

read_gene_file(file[1])
}

cell_lines <- c("SCC47", "SCC1", "SCC90")

# Load 0Gy + 6Gy per cell line
SCC47_0Gy <- load_condition("SCC47", "0Gy")
SCC47_6Gy <- load_condition("SCC47", "6Gy")

SCC1_0Gy <- load_condition("SCC1", "0Gy")
SCC1_6Gy <- load_condition("SCC1", "6Gy")

SCC90_0Gy <- load_condition("SCC90", "0Gy")
SCC90_6Gy <- load_condition("SCC90", "6Gy")

# Combine both doses per cell line
SCC47_all <- unique(c(SCC47_0Gy, SCC47_6Gy))
SCC1_all  <- unique(c(SCC1_0Gy,  SCC1_6Gy))
SCC90_all <- unique(c(SCC90_0Gy, SCC90_6Gy))

# Global combined dataset
Cellular_significant_EV_proteins <- unique(c(SCC47_all, SCC1_all, SCC90_all))

cat("\n=== CELL LINE DATA LOADED ===\n")
cat("SCC47: ", length(SCC47_all), "genes (0Gy+6Gy)\n")
cat("SCC1 : ", length(SCC1_all),  "genes (0Gy+6Gy)\n")
cat("SCC90: ", length(SCC90_all), "genes (0Gy+6Gy)\n")
cat("Total cellular:", length(Cellular_significant_EV_proteins), "\n\n")

#--------------------------- LOAD PATIENT PROTEIN SETS --------------------------
# Cluster 1 = Oncoprotein
onco <- read_gene_file(file.path(
base_dir, "SOTA_Cluster/Oncoprotein/Cluster_1/Cluster_1_genes.txt"
))

# Cluster 5 = Tumor suppressor (ADDED BACK)
ts <- read_gene_file(file.path(
base_dir, "SOTA_Cluster/Tumor_suppressor_protein/Cluster_4/Cluster_4_genes.txt"
))

# Patient SOTA = Cluster 1 + Cluster 4
Patient_EVSOTA_cluster_EV_proteins <- unique(c(onco, ts))

diff_lists <- c(
"Pre_RT_vs_Healthy/Higher_in_Pre_RT.txt",
"Pre_RT_vs_Post_RT/Higher_in_Pre_RT.txt",
"Pre_RT_vs_Post_RT/Higher_in_Post_RT.txt",
"Healthy_vs_Post_RT/Higher_in_Post_RT.txt"
)

Patient_significant_EV_proteins <- unique(unlist(
lapply(diff_lists, function(p) read_gene_file(file.path(base_dir, p)))
))

cat("Patient oncoproteins (Cluster 1):      ", length(onco), "\n")
cat("Patient tumor suppressors (Cluster 4): ", length(ts), "\n")
cat("Patient SOTA (C1 + C4):                ", length(Patient_EVSOTA_cluster_EV_proteins), "\n")
cat("Patient significant EVs:               ", length(Patient_significant_EV_proteins), "\n")

#--------------------------- VENN FUNCTION -------------------------------------
plot_venn_big <- function(A, B, C, prefix = "Venn") {

venn <- venn.diagram(
x = list(
"Patient SOTA\ncluster EV proteins" = A,
"Patient significant\nEV proteins"  = B,
"Cellular significant\nEV proteins" = C
),
filename = NULL,
fill = c("#B6D7A8", "#F9CB9C", "#A4C2F4"),
alpha = 0.55,
lwd = 0,
cex = 3.2,
fontface = "bold",
cat.cex = 2.2,
cat.fontface = "bold",

main = "Venn diagram",
main.cex = 4.0,
main.fontface = "bold",
main.pos = c(0.5, 0.65),

margin = 1,
disable.logging = TRUE
)

venn_grobs <- Filter(is.grob, venn)
venn_grob  <- grobTree(children = do.call(gList, venn_grobs))

png(file.path(output_dir, paste0(prefix, ".png")),
width = 3000, height = 3000, res = 300)
grid.newpage(); grid.draw(venn_grob); dev.off()

pdf(file.path(output_dir, paste0(prefix, ".pdf")), width = 10, height = 10)
grid.newpage(); grid.draw(venn_grob); dev.off()

#--- Export gene lists
AB  <- intersect(A, B)
AC  <- intersect(A, C)
BC  <- intersect(B, C)
ABC <- Reduce(intersect, list(A, B, C))

wb <- createWorkbook()
addWorksheet(wb, "Shared_ABC"); writeData(wb, 1, ABC)
addWorksheet(wb, "Shared_AB");  writeData(wb, 2, AB)
addWorksheet(wb, "Shared_AC");  writeData(wb, 3, AC)
addWorksheet(wb, "Shared_BC");  writeData(wb, 4, BC)
addWorksheet(wb, "Unique_A");   writeData(wb, 5, setdiff(A, union(B, C)))
addWorksheet(wb, "Unique_B");   writeData(wb, 6, setdiff(B, union(A, C)))
addWorksheet(wb, "Unique_C");   writeData(wb, 7, setdiff(C, union(A, B)))

saveWorkbook(wb, file.path(output_dir, paste0(prefix, "_GeneLists.xlsx")),
overwrite = TRUE)

cat("✔ Saved:", prefix, "\n")
}

#--------------------------- RUN MAIN VENN -------------------------------------
plot_venn_big(
Patient_EVSOTA_cluster_EV_proteins,
Patient_significant_EV_proteins,
Cellular_significant_EV_proteins,
prefix = "EV_Proteins_Venn_0Gy_6Gy_all_cell_lines"
)

cat("\nAll files saved to:\n", output_dir, "\n")

################################################################################
################################################################################
################################################################################

################################################################################
# EV Proteins Venn Diagram — FINAL VERSION (ONLY GLOBAL GROUP COMPARISON DATA)
# + Patient SOTA includes Cluster 1 (Oncoprotein) AND Cluster 4 (Tumor suppressor)
# + Directional margins (L/R/Top/Bottom) via grid viewport padding
################################################################################

base_dir     <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics"
cellular_dir <- "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_Cellular_EV_proteomics/Group_comparison"
output_dir   <- file.path(base_dir, "Venn_Diagrams")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
library(openxlsx)
library(VennDiagram)
library(grid)
})

#--------------------------- GENE FILE READER ----------------------------------
read_gene_file <- function(file_path) {
if (!file.exists(file_path)) stop("File not found: ", file_path)

df <- read.table(
file_path,
sep = "\t",
header = TRUE,
quote = "",
stringsAsFactors = FALSE,
check.names = FALSE
)

candidates <- c("Gene", "PG.Genes", "Gene.names", "PG.ProteinGroups")
col_found <- NULL

for (col in candidates) {
if (col %in% colnames(df)) {
col_found <- col
break
}
}

if (is.null(col_found)) {
if (ncol(df) == 1) {
col_found <- colnames(df)[1]
message(sprintf("⚠ Using first column '%s' in %s", col_found, basename(file_path)))
} else {
stop("No usable gene column in: ", file_path)
}
}

genes <- as.character(df[[col_found]])
genes <- sub(";.*", "", genes)                 # keep first gene if multiple separated by ';'
genes <- trimws(genes)
unique(genes[genes != "" & !is.na(genes)])
}

#----------------------- LOAD ONLY GLOBAL CELLULAR DATA ------------------------
cellular_0Gy <- read_gene_file(file.path(cellular_dir, "Group_comparison_0Gy_high_annotated.txt"))
cellular_6Gy <- read_gene_file(file.path(cellular_dir, "Group_comparison_6Gy_high_annotated.txt"))

Cellular_significant_EV_proteins <- unique(c(cellular_0Gy, cellular_6Gy))

#----------------------- LOAD PATIENT DATASETS ---------------------------------
# Cluster 1 = Oncoprotein
onco <- read_gene_file(file.path(
base_dir, "SOTA_Cluster/Oncoprotein/Cluster_1/Cluster_1_genes.txt"
))

# Cluster 4 = Tumor suppressor
ts <- read_gene_file(file.path(
base_dir, "SOTA_Cluster/Tumor_suppressor_protein/Cluster_4/Cluster_4_genes.txt"
))

# Patient SOTA = Cluster 1 + Cluster 4
Patient_EVSOTA_cluster_EV_proteins <- unique(c(onco, ts))

# Differentially expressed EV proteins across comparisons
Patient_significant_EV_proteins <- unique(c(
read_gene_file(file.path(base_dir, "Pre_RT_vs_Healthy/Higher_in_Pre_RT.txt")),
read_gene_file(file.path(base_dir, "Pre_RT_vs_Post_RT/Higher_in_Pre_RT.txt")),
read_gene_file(file.path(base_dir, "Pre_RT_vs_Post_RT/Higher_in_Post_RT.txt")),
read_gene_file(file.path(base_dir, "Healthy_vs_Post_RT/Higher_in_Post_RT.txt"))
))

#----------------------- PRINT COUNTS ------------------------------------------
cat("n(Cellular)                  =", length(Cellular_significant_EV_proteins), "\n")
cat("n(Patient SOTA - Cluster 1)  =", length(onco), "\n")
cat("n(Patient SOTA - Cluster 4)  =", length(ts), "\n")
cat("n(Patient SOTA - C1 + C4)    =", length(Patient_EVSOTA_cluster_EV_proteins), "\n")
cat("n(Patient signif)            =", length(Patient_significant_EV_proteins), "\n")

#--------------------------- VENN FUNCTION -------------------------------------
# NOTE: VennDiagram's 'margin' is global only.
# To get per-side margins, we draw the grob inside a padded viewport.
plot_venn_big <- function(
A, B, C,
prefix = "Venn"
) {

# ---- SINGLE PLACE TO CONTROL MARGINS (units: "lines") ----
margin_left   <- 10
margin_right  <- 2.5
margin_top    <- 2.5
margin_bottom <- 2.0

venn <- venn.diagram(
x = list(
"Patient SOTA\ncluster EV proteins" = A,
"Patient significant\nEV proteins"  = B,
"Cellular significant\nEV proteins" = C
),
filename = NULL,
fill = c("#B6D7A8", "#F9CB9C", "#A4C2F4"),
alpha = 0.55,
lwd = 0,
cex = 3.2,
fontface = "bold",
cat.cex = 2.2,
cat.fontface = "bold",
main = "Venn diagram",
main.cex = 4.0,
main.fontface = "bold",
main.pos = c(0.5, 1),
disable.logging = TRUE
)

venn_grobs <- Filter(is.grob, venn)
venn_grob  <- grobTree(children = do.call(gList, venn_grobs))

# ---- Directional margins (units: "lines") ----
padded_vp <- viewport(
x = unit(0.5, "npc"),
y = unit(0.5, "npc"),
width  = unit(1, "npc") - unit(margin_left + margin_right, "lines"),
height = unit(1, "npc") - unit(margin_top + margin_bottom, "lines")
)

draw_plot <- function() {
grid.newpage()
pushViewport(padded_vp)
grid.draw(venn_grob)
popViewport()
}

# ---- Save PNG ----
png(file.path(output_dir, paste0(prefix, ".png")),
width = 3000, height = 3000, res = 300)
draw_plot()
dev.off()

# ---- Save PDF ----
pdf(file.path(output_dir, paste0(prefix, ".pdf")), width = 10, height = 10)
draw_plot()
dev.off()

# ---- Export gene lists ----
AB  <- intersect(A, B)
AC  <- intersect(A, C)
BC  <- intersect(B, C)
ABC <- Reduce(intersect, list(A, B, C))

wb <- createWorkbook()

addWorksheet(wb, "Shared_ABC"); writeData(wb, "Shared_ABC", ABC)
addWorksheet(wb, "Shared_AB");  writeData(wb, "Shared_AB",  AB)
addWorksheet(wb, "Shared_AC");  writeData(wb, "Shared_AC",  AC)
addWorksheet(wb, "Shared_BC");  writeData(wb, "Shared_BC",  BC)

addWorksheet(wb, "Unique_A");   writeData(wb, "Unique_A", setdiff(A, union(B, C)))
addWorksheet(wb, "Unique_B");   writeData(wb, "Unique_B", setdiff(B, union(A, C)))
addWorksheet(wb, "Unique_C");   writeData(wb, "Unique_C", setdiff(C, union(A, B)))

saveWorkbook(wb,
file.path(output_dir, paste0(prefix, "_GeneLists.xlsx")),
overwrite = TRUE)

cat("✔ Saved:", prefix, "\n")
}

#--------------------------- RUN SINGLE VENN -----------------------------------
plot_venn_big(
Patient_EVSOTA_cluster_EV_proteins,      # Cluster 1 + Cluster 4
Patient_significant_EV_proteins,
Cellular_significant_EV_proteins,
prefix = "EV_Proteins_Venn_0Gy_6Gy"
)

cat("\nAll files saved to:\n", output_dir, "\n")


################################################################################

################################################################################
# Step 8 (PATIENT, MEAN GROUPS): TARGET HEATMAPS of MEAN LFQ — FULL explicit :: 
# - Targets: CD9, CD81, Alix, Syntenin-1  (CD63 & TSG101 REMOVED)
# - Makes TWO mean-group heatmaps:
#     (1) Healthy vs Pre_RT vs Post_RT
#     (2) Pre_RT vs Post_RT
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
  main_dir = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics",
  
  comparisons = list(
    list(name = "Healthy_PreRT_PostRT", groups = c("Healthy", "Pre_RT", "Post_RT")),
    list(name = "PreRT_PostRT",         groups = c("Pre_RT", "Post_RT"))
  ),
  
  target_genes = c("CD9", "CD81", "PDCD6IP", "SDCBP"),
  gene_labels = c(
    CD9     = "CD9",
    CD81    = "CD81",
    PDCD6IP = "Alix",
    SDCBP   = "Syntenin-1"
  ),
  
  out_subdir = "Step8_TargetProtein_Heatmaps_MeanGroup",
  pattern_step7 = "_corrFiltered\\.txt$",
  pattern_step5 = "_filtered_50percent\\.txt$",
  pattern_step2 = "_filtered_log10LFQ\\.txt$",
  
  heatmap_name = "mean_log10_LFQ",
  na_col = "black",
  colors = c("white", "#4F81BD", "#C00000"),
  
  font_row = 20,
  font_col = 22,
  col_angle = 0,
  title_fontsize = 20,
  
  legend_title_fontsize  = 18,
  legend_labels_fontsize = 18,
  legend_title_position  = "topcenter",
  
  legend_gap_right_mm = 16,
  legend_side = "right",
  
  group_order = c("Healthy", "Pre_RT", "Post_RT"),
  
  width_in = 8,
  height_in = 8,
  dpi = 300,
  
  margins = c(t = 12, r = 12, b = 12, l = 30),
  margin_unit = "mm"
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

pick_input_file_main <- function(main_dir, pattern7, pattern5, pattern2) {
  f7 <- list.files(main_dir, pattern = pattern7, full.names = TRUE)
  if (length(f7) > 0) return(f7[1])
  f5 <- list.files(main_dir, pattern = pattern5, full.names = TRUE)
  if (length(f5) > 0) return(f5[1])
  f2 <- list.files(main_dir, pattern = pattern2, full.names = TRUE)
  if (length(f2) > 0) return(f2[1])
  NA_character_
}

# Parse LFQ_<Group>_<Rep>  -> Group + Rep
parse_lfq_name <- function(x) {
  m <- stringr::str_match(x, "^LFQ_(.+)_(\\d+)$")
  data.frame(
    Sample = x,
    Group  = m[, 2],
    Rep    = suppressWarnings(as.integer(m[, 3])),
    stringsAsFactors = FALSE
  )
}

build_target_mean_matrix <- function(df, lfq_cols, target_genes, gene_labels, group_order) {
  
  info <- parse_lfq_name(lfq_cols)
  
  long <- df %>%
    dplyr::select(PG.Genes, dplyr::all_of(lfq_cols)) %>%
    dplyr::mutate(GeneHit = sapply(PG.Genes, function(x) {
      hits <- base::intersect(split_genes(x), target_genes)  # ✅ conflicted-safe
      if (length(hits) == 0) NA_character_ else hits[1]
    })) %>%
    dplyr::filter(!is.na(GeneHit)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(lfq_cols),
      names_to = "Sample",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Value     = suppressWarnings(as.numeric(Value)),
      GeneLabel = unname(gene_labels[GeneHit])
    ) %>%
    dplyr::left_join(info, by = "Sample") %>%
    dplyr::filter(!is.na(Group))
  
  collapsed_sample <- long %>%
    dplyr::group_by(GeneLabel, Group, Sample) %>%
    dplyr::summarise(
      Value = if (all(is.na(Value))) NA_real_ else base::max(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  collapsed_group <- collapsed_sample %>%
    dplyr::group_by(GeneLabel, Group) %>%
    dplyr::summarise(
      Value = if (all(is.na(Value))) NA_real_ else base::mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  wide <- collapsed_group %>%
    dplyr::mutate(Group = factor(Group, levels = group_order)) %>%
    tidyr::pivot_wider(names_from = Group, values_from = Value)
  
  ordered_labels <- unname(gene_labels[target_genes])
  wide <- tibble::tibble(GeneLabel = ordered_labels) %>%
    dplyr::left_join(wide, by = "GeneLabel")
  
  mat <- as.data.frame(wide)
  rownames(mat) <- mat$GeneLabel
  mat$GeneLabel <- NULL
  
  for (g in group_order) if (!g %in% colnames(mat)) mat[[g]] <- NA_real_
  mat <- mat[, group_order, drop = FALSE]
  
  as.matrix(mat)
}

save_target_heatmap <- function(mat, base_fn, cfg) {
  
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
    heatmap_legend_param = legend_param
  )
  
  pad <- grid::unit(
    c(cfg$margins["t"], cfg$margins["r"], cfg$margins["b"], cfg$margins["l"]),
    cfg$margin_unit
  )
  
  # ✅ IMPORTANT: ht_opt is NOT namespaced
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

# =========================
# Run
# =========================
setwd(cfg$main_dir)

in_file <- pick_input_file_main(cfg$main_dir, cfg$pattern_step7, cfg$pattern_step5, cfg$pattern_step2)
if (is.na(in_file)) stop("❌ No input file found.")

df <- read.table(in_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
if (!"PG.Genes" %in% colnames(df)) stop("❌ PG.Genes column not found in input file.")

lfq_cols_all <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols_all) < 2) stop("❌ Not enough LFQ_ columns found.")

df_t <- df[sapply(df$PG.Genes, row_has_target, targets = cfg$target_genes), , drop = FALSE]

out_dir <- file.path(cfg$main_dir, cfg$out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cmp in cfg$comparisons) {
  
  groups <- cmp$groups
  group_order_cmp <- cfg$group_order[cfg$group_order %in% groups]
  
  info_all <- parse_lfq_name(lfq_cols_all)
  lfq_cols_cmp <- lfq_cols_all[info_all$Group %in% groups]
  
  mat <- build_target_mean_matrix(
    df_t, lfq_cols_cmp,
    cfg$target_genes, cfg$gene_labels,
    group_order_cmp
  )
  
  mat_df <- as.data.frame(mat) %>% tibble::rownames_to_column("Protein")
  
  write.table(
    mat_df,
    file.path(out_dir, paste0("DKFZ_", cmp$name, "_TargetProteins_MeanGroup_log10LFQ_matrix.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  openxlsx::write.xlsx(
    mat_df,
    file.path(out_dir, paste0("DKFZ_", cmp$name, "_TargetProteins_MeanGroup_log10LFQ_matrix.xlsx")),
    overwrite = TRUE
  )
  
  base_fn <- file.path(out_dir, paste0("DKFZ_", cmp$name, "_TargetProteins_MeanGroup_log10LFQ_heatmap"))
  save_target_heatmap(mat, base_fn, cfg)
}

cat("\n✅ Step complete — patient MEAN-group target heatmaps saved (CD63 & TSG101 removed).\n")

################################################################################


################################################################################
# Step 8: TARGET HEATMAPS (simple config-driven version) — PATIENT DATASET
# FULL explicit :: (conflicted-safe + ComplexHeatmap::Heatmap/draw)
# - Column order: LEFT = Healthy -> Pre_RT -> Post_RT; within group: Rep (1,2,3...)
# - Makes TWO heatmaps:
#     (1) Healthy vs Pre_RT vs Post_RT
#     (2) Pre_RT vs Post_RT
# - UPDATED LEGEND:
#     * legend title spacing fixed (newline trick)
#     * legend moved away from heatmap WITHOUT moving heatmap (HEATMAP_LEGEND_PADDING)
# IMPORTANT:
#   - Use base::intersect / base::max / base::mean to avoid conflicted errors
#   - ht_opt is NOT namespaced (must stay ht_opt)
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
  main_dir = "C:/Users/ga53hil/Desktop/DKFZ/02.02.26_DKFZ_patient_EV_proteomics",
  
  # comparisons to generate
  comparisons = list(
    list(name = "Healthy_PreRT_PostRT", groups = c("Healthy", "Pre_RT", "Post_RT")),
    list(name = "PreRT_PostRT",         groups = c("Pre_RT", "Post_RT"))
  ),
  
  target_genes = c("CD9", "CD63", "CD81", "PDCD6IP", "SDCBP", "TSG101"),
  gene_labels = c(
    CD9     = "CD9",
    CD63    = "CD63",
    CD81    = "CD81",
    PDCD6IP = "Alix",
    SDCBP   = "Syntenin-1",
    TSG101  = "TSG101"
  ),
  
  # output + file picking (prefer step7, else step5, else log10)
  out_subdir = "Step8_TargetProtein_Heatmaps",
  pattern_step7 = "_corrFiltered\\.txt$",
  pattern_step5 = "_filtered_50percent\\.txt$",
  pattern_step2 = "_filtered_log10LFQ\\.txt$",
  
  # heatmap look
  heatmap_name = "log10 LFQ",
  na_col = "black",
  colors = c("white", "#4F81BD", "#C00000"),
  
  # text + layout
  font_row = 20,
  font_col = 20,
  col_angle = 90,
  title_fontsize = 20,
  
  # legend text
  legend_title_fontsize  = 18,
  legend_labels_fontsize = 18,
  legend_title_position  = "topcenter",
  
  # distance heatmap -> legend (does NOT move heatmap)
  legend_gap_right_mm = 16,
  legend_side = "right",
  
  # column ordering rules
  group_order = c("Healthy", "Pre_RT", "Post_RT"),
  
  # image export
  width_in = 9,
  height_in = 9,
  dpi = 300,
  
  # margins (t/r/b/l) with units (outer whitespace only)
  margins = c(t = 12, r = 8, b = 14, l = 15),
  margin_unit = "mm"
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

pick_input_file_main <- function(main_dir, pattern7, pattern5, pattern2) {
  f7 <- list.files(main_dir, pattern = pattern7, full.names = TRUE)
  if (length(f7) > 0) return(f7[1])
  f5 <- list.files(main_dir, pattern = pattern5, full.names = TRUE)
  if (length(f5) > 0) return(f5[1])
  f2 <- list.files(main_dir, pattern = pattern2, full.names = TRUE)
  if (length(f2) > 0) return(f2[1])
  NA_character_
}

# Parse: LFQ_<Group>_<Rep> where <Group> CAN include underscores (Pre_RT, Post_RT)
parse_lfq_name <- function(x) {
  m <- stringr::str_match(x, "^LFQ_(.+)_(\\d+)$")
  data.frame(
    Sample = x,
    Group  = m[, 2],
    Rep    = suppressWarnings(as.integer(m[, 3])),
    stringsAsFactors = FALSE
  )
}

# ORDER:
# 1) Group (Healthy -> Pre_RT -> Post_RT)
# 2) Rep (1..n)
order_lfq_cols <- function(cols, group_order = c("Healthy", "Pre_RT", "Post_RT")) {
  info <- parse_lfq_name(cols)
  ok <- !is.na(info$Group) & !is.na(info$Rep)
  
  info$GroupRank <- match(info$Group, group_order)
  info$GroupRank[is.na(info$GroupRank)] <- 999L
  
  idx_ok <- which(ok)
  idx_sorted <- idx_ok[order(info$GroupRank[idx_ok], info$Rep[idx_ok], info$Sample[idx_ok])]
  
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
  
  wide <- collapsed %>% tidyr::pivot_wider(names_from = Sample, values_from = Value)
  
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
  
  # robust legend title spacing (newline trick)
  legend_param <- list(
    title          = paste0(cfg$heatmap_name, "\n"),  # change to "\n\n" if you want more space
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
    column_title = NULL,
    column_title_gp = grid::gpar(fontsize = cfg$title_fontsize, fontface = "bold"),
    heatmap_legend_param = legend_param
  )
  
  pad <- grid::unit(
    c(cfg$margins["t"], cfg$margins["r"], cfg$margins["b"], cfg$margins["l"]),
    cfg$margin_unit
  )
  
  # move legend away from heatmap WITHOUT shifting heatmap
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

# =========================
# Run
# =========================
setwd(cfg$main_dir)

cat("\n=== Step 8: Patient target heatmaps ===\n")

in_file <- pick_input_file_main(cfg$main_dir, cfg$pattern_step7, cfg$pattern_step5, cfg$pattern_step2)
if (is.na(in_file)) stop("❌ No input file found in main_dir (step7/step5/log10 patterns).")

cat("Reading:", in_file, "\n")
df <- read.table(in_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

if (!"PG.Genes" %in% colnames(df)) stop("❌ PG.Genes not found — cannot run Step 8.")

lfq_cols_all <- grep("^LFQ_", colnames(df), value = TRUE)
if (length(lfq_cols_all) < 2) stop("❌ Not enough LFQ columns detected.")

info_all <- parse_lfq_name(lfq_cols_all)
cat("Detected LFQ groups:\n")
print(sort(unique(info_all$Group)))

df_t <- df[sapply(df$PG.Genes, row_has_target, targets = cfg$target_genes), , drop = FALSE]
cat("Matched rows:", nrow(df_t), "\n")
if (nrow(df_t) == 0) stop("❌ No target genes found — stopping.")

out_dir <- file.path(cfg$main_dir, cfg$out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cmp in cfg$comparisons) {
  cmp_name   <- cmp$name
  cmp_groups <- cmp$groups
  
  cat("\n--- Comparison:", cmp_name, "(", paste(cmp_groups, collapse = ", "), ") ---\n")
  
  lfq_cols_cmp <- lfq_cols_all[info_all$Group %in% cmp_groups]
  cat("LFQ columns used:", length(lfq_cols_cmp), "\n")
  
  if (length(lfq_cols_cmp) < 2) {
    cat("⚠️ Not enough LFQ columns for", cmp_name, "— skipping.\n")
    next
  }
  
  mat <- build_target_matrix(df_t, lfq_cols_cmp, cfg$target_genes, cfg$gene_labels)
  
  # order columns: subset of the global order
  col_order <- cfg$group_order[cfg$group_order %in% cmp_groups]
  mat <- mat[, order_lfq_cols(colnames(mat), group_order = col_order), drop = FALSE]
  
  mat_df <- as.data.frame(mat) %>% tibble::rownames_to_column("Protein")
  
  write.table(
    mat_df,
    file.path(out_dir, paste0("DKFZ_", cmp_name, "_TargetProteins_log10LFQ_matrix.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  openxlsx::write.xlsx(
    mat_df,
    file.path(out_dir, paste0("DKFZ_", cmp_name, "_TargetProteins_log10LFQ_matrix.xlsx")),
    overwrite = TRUE
  )
  
  heat_title <- paste0("EV markers (log10 LFQ) – ", paste(cmp_groups, collapse = " vs "))
  base_fn <- file.path(out_dir, paste0("DKFZ_", cmp_name, "_TargetProteins_log10LFQ_heatmap"))
  
  save_target_heatmap(mat, heat_title, base_fn, cfg)
  
  cat("✅ Saved outputs to:", out_dir, "\n")
}

cat("\n✅ Step complete — patient target heatmaps saved.\n")

