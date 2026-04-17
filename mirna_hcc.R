# 0. Load Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c("GEOquery", "dplyr", "stringr", "tidyr", "DESeq2", "limma", 
          "ComplexHeatmap", "circlize", "ggplot2", "tibble", "multiMiR", 
          "RColorBrewer", "ggrepel", "BiocParallel", "grid")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)
}

register(SerialParam())
options(mc.cores = 1)


# 1. Setup and Extract RAW Data
setwd("C:/Users/akhil/comet_downlaod/mirna_2")

dir.create("data", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/DEG", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/QC_PCA", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/Heatmap", showWarnings = FALSE, recursive = TRUE)
dir.create("results/rds", showWarnings = FALSE, recursive = TRUE)

if (!file.exists("data/GSE227378_RAW.tar")) {
  getGEOSuppFiles("GSE227378", baseDir = "data", makeDirectory = FALSE)
}
dir.create("data/GSE227378", showWarnings = FALSE)
untar("data/GSE227378_RAW.tar", exdir = "data/GSE227378")

if (!file.exists("data/GSE76903_RAW.tar")) {
  getGEOSuppFiles("GSE76903", baseDir = "data", makeDirectory = FALSE)
}
dir.create("data/GSE76903", showWarnings = FALSE)
untar("data/GSE76903_RAW.tar", exdir = "data/GSE76903")


# 2. Fetch and Build Metadata
gse227 <- getGEO("GSE227378", GSEMatrix = TRUE)[[1]]
meta227 <- pData(gse227) %>%
  mutate(sample_id = geo_accession,
         condition = ifelse(grepl("normal", title, ignore.case = TRUE), "Normal", "Tumor"),
         batch = "GSE227378") %>%
  dplyr::select(sample_id, condition, batch)

gse76 <- getGEO("GSE76903", GSEMatrix = TRUE)[[1]]
meta76 <- pData(gse76) %>%
  mutate(sample_id = geo_accession,
         condition = case_when(grepl("^N", title) ~ "Normal",
                               grepl("^T", title) ~ "Tumor",
                               TRUE ~ "Exclude"),
         batch = "GSE76903") %>%
  filter(condition != "Exclude") %>%
  dplyr::select(sample_id, condition, batch)

meta <- bind_rows(meta227, meta76)
meta$condition <- factor(meta$condition, levels = c("Normal", "Tumor"))
meta$batch <- factor(meta$batch)
write.csv(meta, "results/tables/combined_metadata.csv", row.names = FALSE)


# 3. Build Combined Count Matrix
read_mirna_counts <- function(gsm_id, batch) {
  folder <- ifelse(batch == "GSE227378", "data/GSE227378", "data/GSE76903")
  files <- list.files(folder, pattern = gsm_id, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  df <- read.table(gzfile(files[1]), header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
  cnames <- tolower(colnames(df))
  id_col <- grep("id|mirna|name", cnames)[1]
  count_col <- grep("count|read", cnames)[1]
  
  if (is.na(id_col) || is.na(count_col)) {
    df <- read.table(gzfile(files[1]), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    id_col <- 1; count_col <- 2
  }
  
  df_clean <- df[, c(id_col, count_col)]
  colnames(df_clean) <- c("miRNA", gsm_id)
  
  df_clean %>%
    filter(grepl("^hsa-miR-|^hsa-let-", miRNA, ignore.case = TRUE)) %>%
    mutate(!!gsm_id := as.integer(.[[2]])) %>%
    filter(!is.na(.[[2]]))
}

count_list <- list()
for (i in seq_len(nrow(meta))) {
  res <- read_mirna_counts(meta$sample_id[i], meta$batch[i])
  if (!is.null(res)) count_list[[meta$sample_id[i]]] <- res
}

long_list <- lapply(names(count_list), function(gsm) {
  df <- count_list[[gsm]]
  colnames(df) <- c("miRNA", "Count")
  df$Sample <- gsm
  return(df)
})

count_matrix <- bind_rows(long_list) %>%
  pivot_wider(names_from = Sample, values_from = Count, values_fill = list(Count = 0), values_fn = sum) %>%
  as.data.frame()

mirna_names <- count_matrix$miRNA
count_matrix$miRNA <- NULL
count_matrix <- as.data.frame(lapply(count_matrix, as.integer))
rownames(count_matrix) <- mirna_names
count_matrix <- count_matrix[, meta$sample_id]

keep <- rowSums(count_matrix >= 10) >= floor(ncol(count_matrix) / 2)
count_filtered <- count_matrix[keep, ]
write.csv(count_filtered, "results/tables/combined_count_matrix.csv", row.names = TRUE)


# 4. Run DESeq2 with Batch Correction
dds <- DESeqDataSetFromMatrix(countData = count_filtered, colData = meta, design = ~ batch + condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)

res_df <- as.data.frame(res) %>%
  rownames_to_column("miRNA") %>%
  arrange(padj) %>%
  mutate(regulation = case_when(
    !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Up",
    !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

write.csv(res_df, "results/tables/DESeq2_results.csv", row.names = FALSE)
vst_data <- assay(varianceStabilizingTransformation(dds, blind = FALSE))


# 5. Volcano Plot
top10_up <- res_df %>% filter(regulation == "Up") %>% arrange(padj) %>% slice_head(n = 10)
top10_dn <- res_df %>% filter(regulation == "Down") %>% arrange(padj) %>% slice_head(n = 10)
top20_labeled <- bind_rows(top10_up, top10_dn)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = regulation)) +
  geom_point(alpha = 0.5, size = 1.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_text_repel(data = top20_labeled, aes(label = miRNA), size = 3.5, fontface = "bold.italic", 
                  min.segment.length = 0, segment.size = 0.4, max.overlaps = Inf, box.padding = 0.8, 
                  point.padding = 0.3, show.legend = FALSE) +
  scale_colour_manual(name = "Regulation", values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
                      labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "NS" = "NS")) +
  labs(x = expression(bold(Log[2]~Fold~Change~(Tumor~vs~Normal))), y = expression(bold(-log[10](italic(p)[adj])))) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
        legend.position = "right", legend.justification = "top", legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.minor = element_blank(), plot.title = element_blank())

ggsave("results/plots/DEG/Volcano_plot.png", volcano_plot, width = 8, height = 7, dpi = 600)
ggsave("results/plots/DEG/Volcano_plot.tiff", volcano_plot, width = 8, height = 7, dpi = 600, compression = "none")


# 6. PCA Plots
if (!inherits(vst_data, "SummarizedExperiment")) {
  vst_obj <- DESeqTransform(SummarizedExperiment(assays = list(vst = vst_data), colData = meta))
} else {
  vst_obj <- vst_data
}

pca_df <- plotPCA(vst_obj, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

pca_cond <- ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  scale_colour_manual(name = "Condition", values = c("Normal" = "#2166AC", "Tumor" = "#D6604D")) +
  scale_fill_manual(name = "Condition", values = c("Normal" = "#2166AC", "Tumor" = "#D6604D")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 13) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
        legend.position = "right", legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.5), panel.grid.minor = element_blank())

ggsave("results/plots/QC_PCA/PCA_by_Condition.png", pca_cond, width = 10, height = 7, dpi = 600)

pca_batch <- ggplot(pca_df, aes(PC1, PC2, colour = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance"), colour = "GEO Study") +
  theme_bw(base_size = 13) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
        legend.position = "right", legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
        legend.box.background = element_rect(colour = "black", linewidth = 0.5), panel.grid.minor = element_blank())

ggsave("results/plots/QC_PCA/PCA_by_Batch.png", pca_batch, width = 10, height = 7, dpi = 600)


# 7. Heatmap
top25_up <- res_df %>% filter(regulation == "Up") %>% arrange(padj) %>% slice_head(n = 25)
top25_dn <- res_df %>% filter(regulation == "Down") %>% arrange(padj) %>% slice_head(n = 25)
top_degs <- bind_rows(top25_dn, top25_up)

vst_matrix <- if(inherits(vst_data, "SummarizedExperiment")) assay(vst_data) else vst_data
design_matrix <- model.matrix(~ condition, data = meta)
vst_matrix_clean <- removeBatchEffect(vst_matrix, batch = meta$batch, design = design_matrix)

col_order <- c(rownames(meta)[meta$condition == "Normal"], rownames(meta)[meta$condition == "Tumor"])
hmat <- vst_matrix_clean[top_degs$miRNA, col_order, drop = FALSE]
rownames(hmat) <- top_degs$miRNA

hmat_z <- t(scale(t(hmat)))
valid <- apply(hmat_z, 1, function(x) !all(is.na(x)) & !any(is.infinite(x)))
hmat_z <- hmat_z[valid, , drop = FALSE]

reg_vec <- top_degs$regulation[match(rownames(hmat_z), top_degs$miRNA)]
row_ord <- c(which(reg_vec == "Down"), which(reg_vec == "Up"))
hmat_z <- hmat_z[row_ord, ]
reg_vec <- reg_vec[row_ord]

col_split <- factor(meta[col_order, "condition"], levels = c("Normal", "Tumor"))
condition_colors <- c("Normal" = "#4575B4", "Tumor" = "#FFD700")
n_batch <- length(levels(meta$batch))
batch_colors <- setNames(brewer.pal(max(3, n_batch), "Set2")[1:n_batch], levels(meta$batch))
col_fun <- colorRamp2(c(-2, 0, 2), c("green3", "black", "red"))

ann_df <- data.frame(Condition = meta[col_order, "condition"], Dataset = meta[col_order, "batch"], row.names = col_order)

ha_top <- HeatmapAnnotation(Condition = ann_df$Condition, Dataset = ann_df$Dataset,
                            col = list(Condition = condition_colors, Dataset = batch_colors),
                            annotation_name_gp = gpar(fontface = "bold", fontsize = 10),
                            show_annotation_name = TRUE, show_legend = FALSE)

ha_left <- rowAnnotation(Regulation = reg_vec,
                         col = list(Regulation = c("Up" = "#8BC34A", "Down" = "#CE93D8")),
                         width = unit(0.4, "cm"), annotation_name_gp = gpar(fontface = "bold", fontsize = 10),
                         annotation_name_side = "top", show_annotation_name = TRUE, show_legend = FALSE)

ht <- Heatmap(hmat_z, name = "Z-score", col = col_fun, top_annotation = ha_top, left_annotation = ha_left,
              cluster_rows = FALSE, cluster_columns = TRUE, cluster_column_slices = FALSE,
              column_split = col_split, column_gap = unit(0, "mm"), rect_gp = gpar(col = NA),
              show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right",
              row_names_gp = gpar(fontface = "bold.italic", fontsize = 9), show_heatmap_legend = FALSE,
              row_split = factor(reg_vec, levels = c("Down", "Up")), row_gap = unit(0, "mm"),
              row_title = NULL, column_title = NULL, border = TRUE)

lgd_condition <- Legend(title = "Condition", at = c("Normal", "Tumor"), legend_gp = gpar(fill = c("#4575B4", "#FFD700")), title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontface = "bold", fontsize = 9))
lgd_regulation <- Legend(title = "Regulation", at = c("Down", "Up"), labels = c("Down-regulated", "Up-regulated"), legend_gp = gpar(fill = c("#CE93D8", "#8BC34A")), title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontface = "bold", fontsize = 9))
lgd_dataset <- Legend(title = "Dataset", at = levels(meta$batch), legend_gp = gpar(fill = batch_colors), title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontface = "bold", fontsize = 9))
lgd_zscore <- Legend(title = "Z-score", col_fun = col_fun, at = c(-2, 0, 2), labels = c("-2", "0", "2"), legend_height = unit(3.5, "cm"), title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontface = "bold", fontsize = 9), direction = "vertical")

all_legends <- packLegend(lgd_condition, lgd_regulation, lgd_dataset, lgd_zscore, direction = "vertical", gap = unit(4, "mm"))

draw_ht <- function() {
  draw(ht, annotation_legend_list = all_legends, heatmap_legend_side = "right", annotation_legend_side = "right", align_annotation_legend = "heatmap_top", merge_legend = FALSE, padding = unit(c(2, 2, 2, 2), "mm"))
}

png("results/plots/Heatmap/Heatmap_top25up_top25dn.png", width = 14, height = 11, units = "in", res = 600, type = "cairo")
draw_ht()
dev.off()

tiff("results/plots/Heatmap/Heatmap_top25up_top25dn.tiff", width = 14, height = 11, units = "in", res = 600, compression = "none")
draw_ht()
dev.off()


# 8. Query multiMiR and Filter Targets
up_mirnas <- res_df %>% filter(regulation == "Up") %>% pull(miRNA)
down_mirnas <- res_df %>% filter(regulation == "Down") %>% pull(miRNA)

res_up_val <- get_multimir(mirna = up_mirnas, org = "hsa", table = "mirtarbase", summary = TRUE)
res_up_pred <- get_multimir(mirna = up_mirnas, org = "hsa", table = "mirdb", summary = TRUE)
res_dn_val <- get_multimir(mirna = down_mirnas, org = "hsa", table = "mirtarbase", summary = TRUE)
res_dn_pred <- get_multimir(mirna = down_mirnas, org = "hsa", table = "mirdb", summary = TRUE)

saveRDS(res_up_val, "results/rds/mirtarbase_up.rds")
saveRDS(res_dn_val, "results/rds/mirtarbase_down.rds")

df_up_val <- res_up_val@data %>% filter(!is.na(target_symbol), target_symbol != "")
df_dn_val <- res_dn_val@data %>% filter(!is.na(target_symbol), target_symbol != "")
df_up_pred <- res_up_pred@data %>% mutate(score = as.numeric(score)) %>% filter(score >= 80, target_symbol != "")
df_dn_pred <- res_dn_pred@data %>% mutate(score = as.numeric(score)) %>% filter(score >= 80, target_symbol != "")

inter_up <- intersect(unique(df_up_val$target_symbol), unique(df_up_pred$target_symbol))
inter_dn <- intersect(unique(df_dn_val$target_symbol), unique(df_dn_pred$target_symbol))

df_up_highconf <- df_up_val %>% filter(target_symbol %in% inter_up)
df_dn_highconf <- df_dn_val %>% filter(target_symbol %in% inter_dn)

filter_targets <- function(df, n_cutoff) {
  result <- aggregate(mature_mirna_id ~ target_symbol, data = df, FUN = function(x) length(unique(x)))
  colnames(result)[2] <- "n_mirnas"
  result %>% filter(n_mirnas >= n_cutoff) %>% arrange(desc(n_mirnas))
}

results_table <- list()
for (n in 1:8) {
  up_n <- filter_targets(df_up_highconf, n)
  dn_n <- filter_targets(df_dn_highconf, n)
  results_table[[as.character(n)]] <- list(up = up_n, down = dn_n)
}

final_up_df <- results_table[["8"]]$up
final_dn_df <- results_table[["4"]]$down

writeLines(final_up_df$target_symbol, "results/tables/STRING_input_UP.txt")
writeLines(final_dn_df$target_symbol, "results/tables/STRING_input_DOWN.txt")
write.csv(final_up_df, "results/tables/targets_UP_degree_counts.csv", row.names = FALSE)
write.csv(final_dn_df, "results/tables/targets_DOWN_degree_counts.csv", row.names = FALSE)


# 9. Format IDs for UALCAN Analysis
top5_up <- res_df %>%
  filter(regulation == "Up") %>%
  arrange(padj) %>%
  head(5) %>%
  dplyr::select(miRNA, log2FoldChange, padj, regulation)

top5_down <- res_df %>%
  filter(regulation == "Down") %>%
  arrange(padj) %>%
  head(5) %>%
  dplyr::select(miRNA, log2FoldChange, padj, regulation)

top10_ualcan <- bind_rows(top5_up, top5_down)
top10_ualcan$ualcan_id <- gsub("-[35]p$", "", tolower(top10_ualcan$miRNA))

print("UP-regulated miRNAs for UALCAN:")
print(top10_ualcan$ualcan_id[top10_ualcan$regulation == "Up"])

print("DOWN-regulated miRNAs for UALCAN:")
print(top10_ualcan$ualcan_id[top10_ualcan$regulation == "Down"])

write.csv(top10_ualcan, "results/tables/Top10_miRNAs_for_UALCAN.csv", row.names = FALSE)





