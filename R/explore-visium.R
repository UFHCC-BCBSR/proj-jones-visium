library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ============================================================================
# PART 1: LOAD ALL SAMPLES WITH METADATA
# ============================================================================

# Define sample information
sample_info <- data.frame(
  sample_id = c("A1", "A2", "A3", "A4", "D1", "D2", "D3", "D4"),
  path = paste0("A", 1:4, "_AM_Spatial_Seq"),
  group = c(rep("A", 4), rep("D", 4)),
  qc_grade = c("Marginal", "Marginal", "Marginal", "Good",
               "Good", "Marginal", "Good", "Good"),
  median_umi = c(554, 134, 556, 1384, 602, 289, 567, 1198),
  median_genes = c(479, 105, 454, 981, 486, 241, 474, 921),
  confident_mapping = c(22.6, 70.8, 83.8, 82.1, 92.3, 91.2, 85.3, 89.9),
  stringsAsFactors = FALSE
)

# Load all samples
cat("Loading samples...\n")
spatial_list <- list()

for(i in 1:nrow(sample_info)) {
  sample_name <- sample_info$sample_id[i]
  sample_path <- sample_info$path[i]

  cat(paste0("Loading ", sample_name, "...\n"))

  tryCatch({
    obj <- Load10X_Spatial(
      data.dir = paste0("/orange/cancercenter-dept/JONES/AM_10x_Visium_v2/SpaceRanger_Output_and_Scripts/",
                        sample_path, "/outs"),
      filename = "filtered_feature_bc_matrix.h5"
    )

    # Add sample metadata
    obj$sample_id <- sample_name
    obj$group <- sample_info$group[i]
    obj$qc_grade <- sample_info$qc_grade[i]
    obj$batch <- ifelse(sample_info$group[i] == "A", "Batch_A", "Batch_D")

    spatial_list[[sample_name]] <- obj
    cat(paste0("  - Loaded ", ncol(obj), " spots\n"))

  }, error = function(e) {
    cat(paste0("  - ERROR loading ", sample_name, ": ", e$message, "\n"))
  })
}

# ============================================================================
# PART 2: CALCULATE QC METRICS FOR ALL SAMPLES
# ============================================================================

cat("\nCalculating QC metrics...\n")

for(sample_name in names(spatial_list)) {
  obj <- spatial_list[[sample_name]]

  # Core QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["log10_umi"]] <- log10(obj$nCount_Spatial + 1)
  obj[["log10_genes"]] <- log10(obj$nFeature_Spatial + 1)

  # Calculate complexity (genes per UMI - indicator of library complexity)
  obj[["complexity"]] <- obj$nFeature_Spatial / obj$nCount_Spatial

  spatial_list[[sample_name]] <- obj
}

# ============================================================================
# PART 3: VISUALIZE PRE-QC METRICS
# ============================================================================

cat("\nGenerating pre-QC visualizations...\n")

# Create pre-QC summary plots for each sample
pdf("00_PreQC_Summary.pdf", width = 16, height = 12)

for(sample_name in names(spatial_list)) {
  obj <- spatial_list[[sample_name]]

  cat(paste0("  Plotting ", sample_name, "\n"))

  # Spatial feature plots
  p1 <- SpatialFeaturePlot(obj, features = "nCount_Spatial") +
    ggtitle(paste0(sample_name, " - UMI Counts"))
  p2 <- SpatialFeaturePlot(obj, features = "nFeature_Spatial") +
    ggtitle(paste0(sample_name, " - Gene Counts"))
  p3 <- SpatialFeaturePlot(obj, features = "percent.mt") +
    ggtitle(paste0(sample_name, " - MT%"))
  p4 <- SpatialFeaturePlot(obj, features = "complexity") +
    ggtitle(paste0(sample_name, " - Complexity"))

  # Violin plots
  p5 <- VlnPlot(obj, features = c("nCount_Spatial", "nFeature_Spatial",
                                  "percent.mt", "complexity"),
                pt.size = 0, ncol = 4) +
    plot_annotation(title = paste0(sample_name, " - Distribution"))

  # Scatter plots to identify outliers
  p6 <- ggplot(obj@meta.data, aes(x = log10_umi, y = log10_genes)) +
    geom_point(aes(color = percent.mt), size = 0.5) +
    scale_color_viridis_c() +
    geom_smooth(method = "lm", color = "red") +
    labs(title = paste0(sample_name, " - UMI vs Genes"),
         x = "Log10(UMI)", y = "Log10(Genes)") +
    theme_minimal()

  print((p1 | p2) / (p3 | p4))
  print(p5)
  print(p6)
}

dev.off()

# ============================================================================
# PART 4: SAMPLE-SPECIFIC ADAPTIVE QC FILTERING
# ============================================================================

cat("\n" %+% "="*80 %+% "\n")
cat("ADAPTIVE QC FILTERING - RATIONALE\n")
cat("="*80 %+% "\n\n")

cat("These FFPE samples show severe quality issues requiring surgical filtering.\n")
cat("We will use ADAPTIVE, SAMPLE-SPECIFIC thresholds because:\n\n")
cat("1. Fixed thresholds would eliminate entire samples\n")
cat("2. FFPE degradation varies by sample\n")
cat("3. We need to preserve biological signal while removing technical noise\n\n")
cat("Strategy: Remove bottom 5-10% of spots by multiple criteria\n\n")

spatial_filtered <- list()
qc_summary <- data.frame()

for(sample_name in names(spatial_list)) {
  obj <- spatial_list[[sample_name]]
  n_spots_initial <- ncol(obj)

  cat(paste0("\n", "─"*80, "\n"))
  cat(paste0("FILTERING: ", sample_name, "\n"))
  cat(paste0("─"*80, "\n"))
  cat(paste0("Initial spots: ", n_spots_initial, "\n"))
  cat(paste0("Known issues: Median UMI = ",
             sample_info$median_umi[sample_info$sample_id == sample_name],
             ", Mapping = ",
             sample_info$confident_mapping[sample_info$sample_id == sample_name], "%\n\n"))

  # -------------------------------------------------------------------------
  # FILTER 1: Remove extreme low-count spots (bottom 5%)
  # -------------------------------------------------------------------------
  # Rationale: These spots likely represent technical failures (no tissue,
  # poor permeabilization, or capture probe failure), not biological signal

  umi_threshold <- quantile(obj$nCount_Spatial, 0.05)
  gene_threshold <- quantile(obj$nFeature_Spatial, 0.05)

  cat(paste0("FILTER 1: Remove extreme low-count spots\n"))
  cat(paste0("  - UMI threshold (5th percentile): ", round(umi_threshold, 0), "\n"))
  cat(paste0("  - Gene threshold (5th percentile): ", round(gene_threshold, 0), "\n"))

  filter1_keep <- obj$nCount_Spatial > umi_threshold &
    obj$nFeature_Spatial > gene_threshold
  cat(paste0("  - Spots removed: ", sum(!filter1_keep),
             " (", round(100*sum(!filter1_keep)/n_spots_initial, 1), "%)\n\n"))

  # -------------------------------------------------------------------------
  # FILTER 2: Remove mitochondrial outliers (if MT% > 95th percentile AND MT% > 20%)
  # -------------------------------------------------------------------------
  # Rationale: Very high MT% indicates either:
  # (a) Dying/dead cells with cytoplasmic RNA depletion
  # (b) Poor RNA quality with preferential nuclear RNA loss
  # (c) Mitochondria-rich cell types (legitimate biology)
  # We only remove spots that are BOTH outliers AND above absolute threshold

  mt_upper <- quantile(obj$percent.mt, 0.95)

  cat(paste0("FILTER 2: Remove mitochondrial outliers\n"))
  cat(paste0("  - MT% 95th percentile: ", round(mt_upper, 1), "%\n"))
  cat(paste0("  - Absolute threshold: 20%\n"))

  if(mt_upper > 20) {
    filter2_keep <- obj$percent.mt < mt_upper
    cat(paste0("  - Filtering applied (95th percentile > 20%)\n"))
  } else {
    filter2_keep <- obj$percent.mt < 20
    cat(paste0("  - Using absolute threshold (95th percentile <= 20%)\n"))
  }

  cat(paste0("  - Spots removed: ", sum(!filter2_keep),
             " (", round(100*sum(!filter2_keep)/n_spots_initial, 1), "%)\n\n"))

  # -------------------------------------------------------------------------
  # FILTER 3: Remove low-complexity outliers (bottom 5%)
  # -------------------------------------------------------------------------
  # Rationale: Complexity (genes/UMI) identifies spots where few genes
  # dominate the transcriptome. This can indicate:
  # (a) Ambient RNA contamination (reads dominated by abundant transcripts)
  # (b) Technical artifacts
  # (c) Spots at tissue edge with poor capture
  # Low complexity with low counts = technical noise, not biology

  complexity_threshold <- quantile(obj$complexity, 0.05)

  cat(paste0("FILTER 3: Remove low-complexity outliers\n"))
  cat(paste0("  - Complexity threshold (5th percentile): ",
             round(complexity_threshold, 3), "\n"))

  filter3_keep <- obj$complexity > complexity_threshold
  cat(paste0("  - Spots removed: ", sum(!filter3_keep),
             " (", round(100*sum(!filter3_keep)/n_spots_initial, 1), "%)\n\n"))

  # -------------------------------------------------------------------------
  # FILTER 4: For critically low-mapping samples, be MORE stringent
  # -------------------------------------------------------------------------
  # Rationale: Sample A1 has 22.6% confident mapping, indicating severe
  # contamination or technical failure. For such samples, we need to be
  # more aggressive and remove spots in bottom 10% instead of 5%

  mapping_rate <- sample_info$confident_mapping[sample_info$sample_id == sample_name]

  if(mapping_rate < 50) {
    cat(paste0("FILTER 4: CRITICAL SAMPLE - Enhanced filtering\n"))
    cat(paste0("  - This sample has critically low mapping (", mapping_rate, "%)\n"))
    cat(paste0("  - Applying enhanced filtering (bottom 10% instead of 5%)\n"))

    umi_threshold_strict <- quantile(obj$nCount_Spatial, 0.10)
    gene_threshold_strict <- quantile(obj$nFeature_Spatial, 0.10)

    filter4_keep <- obj$nCount_Spatial > umi_threshold_strict &
      obj$nFeature_Spatial > gene_threshold_strict

    cat(paste0("  - Additional spots removed: ", sum(!filter4_keep),
               " (", round(100*sum(!filter4_keep)/n_spots_initial, 1), "%)\n\n"))
  } else {
    filter4_keep <- rep(TRUE, ncol(obj))
    cat(paste0("FILTER 4: Mapping acceptable - no enhanced filtering\n\n"))
  }

  # -------------------------------------------------------------------------
  # FILTER 5: Remove spatial outliers (isolated low-quality spots)
  # -------------------------------------------------------------------------
  # Rationale: True tissue regions should have spatially coherent quality.
  # Isolated low-quality spots surrounded by high-quality spots likely
  # represent technical artifacts (dust, tissue folding, bubbles)

  cat(paste0("FILTER 5: Spatial outlier detection\n"))
  cat(paste0("  - Identifying isolated low-quality spots\n"))

  # Get spatial coordinates
  coords <- GetTissueCoordinates(obj)

  # For each spot, check if neighboring spots are also low quality
  # If a spot is low quality but all neighbors are high quality, remove it
  filter5_keep <- rep(TRUE, ncol(obj))

  # Simple approach: flag spots in bottom 20% that are isolated
  low_quality_spots <- obj$nCount_Spatial < quantile(obj$nCount_Spatial, 0.20)

  # Note: Full spatial neighbor analysis would require more complex code
  # For now, we'll use this as a placeholder and rely on other filters

  cat(paste0("  - Spatial filtering: Placeholder (manual review recommended)\n\n"))

  # -------------------------------------------------------------------------
  # COMBINE ALL FILTERS
  # -------------------------------------------------------------------------

  final_keep <- filter1_keep & filter2_keep & filter3_keep & filter4_keep & filter5_keep

  n_spots_removed <- sum(!final_keep)
  n_spots_retained <- sum(final_keep)
  percent_removed <- round(100 * n_spots_removed / n_spots_initial, 1)
  percent_retained <- round(100 * n_spots_retained / n_spots_initial, 1)

  cat(paste0("FINAL FILTERING SUMMARY:\n"))
  cat(paste0("  - Initial spots: ", n_spots_initial, "\n"))
  cat(paste0("  - Spots removed: ", n_spots_removed, " (", percent_removed, "%)\n"))
  cat(paste0("  - Spots retained: ", n_spots_retained, " (", percent_retained, "%)\n\n"))

  # Apply filtering
  obj_filtered <- subset(obj, cells = colnames(obj)[final_keep])

  # Store QC summary
  qc_summary <- rbind(qc_summary, data.frame(
    sample = sample_name,
    group = sample_info$group[sample_info$sample_id == sample_name],
    initial_spots = n_spots_initial,
    retained_spots = n_spots_retained,
    removed_spots = n_spots_removed,
    percent_retained = percent_retained,
    median_umi_before = median(obj$nCount_Spatial),
    median_umi_after = median(obj_filtered$nCount_Spatial),
    median_genes_before = median(obj$nFeature_Spatial),
    median_genes_after = median(obj_filtered$nFeature_Spatial),
    mean_mt_before = mean(obj$percent.mt),
    mean_mt_after = mean(obj_filtered$percent.mt)
  ))

  spatial_filtered[[sample_name]] <- obj_filtered
}

# Print QC summary table
cat("\n" %+% "="*80 %+% "\n")
cat("FILTERING SUMMARY - ALL SAMPLES\n")
cat("="*80 %+% "\n\n")
print(qc_summary)

# Save QC summary
write.csv(qc_summary, "QC_Filtering_Summary.csv", row.names = FALSE)

# ============================================================================
# PART 5: POST-QC VISUALIZATION
# ============================================================================

cat("\nGenerating post-QC visualizations...\n")

pdf("01_PostQC_Summary.pdf", width = 16, height = 12)

for(sample_name in names(spatial_filtered)) {
  obj <- spatial_filtered[[sample_name]]

  # Spatial plots
  p1 <- SpatialFeaturePlot(obj, features = "nCount_Spatial") +
    ggtitle(paste0(sample_name, " POST-QC - UMI Counts"))
  p2 <- SpatialFeaturePlot(obj, features = "nFeature_Spatial") +
    ggtitle(paste0(sample_name, " POST-QC - Genes"))
  p3 <- SpatialFeaturePlot(obj, features = "percent.mt") +
    ggtitle(paste0(sample_name, " POST-QC - MT%"))
  p4 <- SpatialFeaturePlot(obj, features = "complexity") +
    ggtitle(paste0(sample_name, " POST-QC - Complexity"))

  print((p1 | p2) / (p3 | p4))
}

dev.off()

# ============================================================================
# PART 6: NORMALIZE AND PROCESS EACH SAMPLE
# ============================================================================

cat("\nNormalizing and processing filtered samples...\n")

for(sample_name in names(spatial_filtered)) {
  cat(paste0("  Processing ", sample_name, "...\n"))

  obj <- spatial_filtered[[sample_name]]

  # SCTransform normalization (recommended for spatial data)
  # This accounts for sequencing depth and stabilizes variance
  obj <- SCTransform(obj, assay = "Spatial",
                     vars.to.regress = "percent.mt",
                     verbose = FALSE)

  # PCA
  obj <- RunPCA(obj, assay = "SCT", npcs = 50, verbose = FALSE)

  # UMAP
  obj <- RunUMAP(obj, reduction = "PCA", dims = 1:30, verbose = FALSE)

  # Clustering
  obj <- FindNeighbors(obj, reduction = "PCA", dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  spatial_filtered[[sample_name]] <- obj
}

# ============================================================================
# PART 7: MERGE SAMPLES BY GROUP
# ============================================================================

cat("\n" %+% "="*80 %+% "\n")
cat("MERGING SAMPLES BY GROUP\n")
cat("="*80 %+% "\n\n")

# Group A samples
group_A_samples <- spatial_filtered[grep("^A", names(spatial_filtered))]
cat(paste0("Group A samples: ", paste(names(group_A_samples), collapse = ", "), "\n"))

# Group D samples
group_D_samples <- spatial_filtered[grep("^D", names(spatial_filtered))]
cat(paste0("Group D samples: ", paste(names(group_D_samples), collapse = ", "), "\n\n"))

# -------------------------------------------------------------------------
# MERGE GROUP A
# -------------------------------------------------------------------------

cat("Merging Group A...\n")

if(length(group_A_samples) > 0) {
  # Merge using standard Seurat merge
  group_A_merged <- merge(x = group_A_samples[[1]],
                          y = group_A_samples[-1],
                          add.cell.ids = names(group_A_samples),
                          project = "Group_A")

  cat(paste0("  - Total spots in Group A: ", ncol(group_A_merged), "\n"))
  cat(paste0("  - Total genes: ", nrow(group_A_merged), "\n\n"))

  # Re-process merged object
  cat("  Re-processing merged Group A...\n")
  group_A_merged <- SCTransform(group_A_merged, assay = "Spatial",
                                vars.to.regress = c("percent.mt", "sample_id"),
                                verbose = FALSE)
  group_A_merged <- RunPCA(group_A_merged, assay = "SCT", npcs = 50, verbose = FALSE)
  group_A_merged <- RunUMAP(group_A_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
  group_A_merged <- FindNeighbors(group_A_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
  group_A_merged <- FindClusters(group_A_merged, resolution = 0.5, verbose = FALSE)

  cat("  Done.\n\n")
}

# -------------------------------------------------------------------------
# MERGE GROUP D
# -------------------------------------------------------------------------

cat("Merging Group D...\n")

if(length(group_D_samples) > 0) {
  group_D_merged <- merge(x = group_D_samples[[1]],
                          y = group_D_samples[-1],
                          add.cell.ids = names(group_D_samples),
                          project = "Group_D")

  cat(paste0("  - Total spots in Group D: ", ncol(group_D_merged), "\n"))
  cat(paste0("  - Total genes: ", nrow(group_D_merged), "\n\n"))

  cat("  Re-processing merged Group D...\n")
  group_D_merged <- SCTransform(group_D_merged, assay = "Spatial",
                                vars.to.regress = c("percent.mt", "sample_id"),
                                verbose = FALSE)
  group_D_merged <- RunPCA(group_D_merged, assay = "SCT", npcs = 50, verbose = FALSE)
  group_D_merged <- RunUMAP(group_D_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
  group_D_merged <- FindNeighbors(group_D_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
  group_D_merged <- FindClusters(group_D_merged, resolution = 0.5, verbose = FALSE)

  cat("  Done.\n\n")
}

# ============================================================================
# PART 8: SAVE PROCESSED OBJECTS
# ============================================================================

cat("Saving processed objects...\n")

# Save individual filtered samples
for(sample_name in names(spatial_filtered)) {
  saveRDS(spatial_filtered[[sample_name]],
          file = paste0("Filtered_", sample_name, "_Spatial.rds"))
  cat(paste0("  - Saved ", sample_name, "\n"))
}

# Save merged objects
saveRDS(group_A_merged, file = "Group_A_Merged_Spatial.rds")
saveRDS(group_D_merged, file = "Group_D_Merged_Spatial.rds")

cat("All objects saved.\n\n")

# ============================================================================
# PART 9: MERGED OBJECT VISUALIZATIONS
# ============================================================================

cat("="*80, "\n")
cat("MERGED OBJECT ANALYSIS\n")
cat("="*80, "\n\n")

pdf("02_Merged_Groups_Analysis.pdf", width = 16, height = 10)

# -------------------------------------------------------------------------
# GROUP A VISUALIZATIONS
# -------------------------------------------------------------------------

cat("Generating Group A visualizations...\n")

# UMAP colored by sample
p1 <- DimPlot(group_A_merged, reduction = "umap", group.by = "sample_id") +
  ggtitle("Group A - UMAP by Sample") +
  theme_minimal()

# UMAP colored by cluster
p2 <- DimPlot(group_A_merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Group A - UMAP by Cluster") +
  theme_minimal()

# QC metrics on UMAP
p3 <- FeaturePlot(group_A_merged, features = c("nCount_Spatial", "nFeature_Spatial"),
                  reduction = "umap", ncol = 2)

print(p1 | p2)
print(p3)

# Batch effect assessment
p4 <- VlnPlot(group_A_merged, features = c("nCount_Spatial", "nFeature_Spatial"),
              group.by = "sample_id", pt.size = 0) +
  ggtitle("Group A - QC by Sample (check for batch effects)")

print(p4)

# Cluster composition by sample
cluster_composition_A <- table(group_A_merged$sample_id, group_A_merged$seurat_clusters)
cluster_composition_A_pct <- prop.table(cluster_composition_A, margin = 1) * 100

# Convert to data frame for plotting
cluster_comp_df_A <- as.data.frame(cluster_composition_A_pct)
colnames(cluster_comp_df_A) <- c("Sample", "Cluster", "Percentage")

p5 <- ggplot(cluster_comp_df_A, aes(x = Sample, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Group A - Cluster Composition by Sample",
       y = "Percentage of Spots") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# -------------------------------------------------------------------------
# GROUP D VISUALIZATIONS
# -------------------------------------------------------------------------

cat("Generating Group D visualizations...\n")

# UMAP colored by sample
p6 <- DimPlot(group_D_merged, reduction = "umap", group.by = "sample_id") +
  ggtitle("Group D - UMAP by Sample") +
  theme_minimal()

# UMAP colored by cluster
p7 <- DimPlot(group_D_merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Group D - UMAP by Cluster") +
  theme_minimal()

# QC metrics on UMAP
p8 <- FeaturePlot(group_D_merged, features = c("nCount_Spatial", "nFeature_Spatial"),
                  reduction = "umap", ncol = 2)

print(p6 | p7)
print(p8)

# Batch effect assessment
p9 <- VlnPlot(group_D_merged, features = c("nCount_Spatial", "nFeature_Spatial"),
              group.by = "sample_id", pt.size = 0) +
  ggtitle("Group D - QC by Sample (check for batch effects)")

print(p9)

# Cluster composition by sample
cluster_composition_D <- table(group_D_merged$sample_id, group_D_merged$seurat_clusters)
cluster_composition_D_pct <- prop.table(cluster_composition_D, margin = 1) * 100

cluster_comp_df_D <- as.data.frame(cluster_composition_D_pct)
colnames(cluster_comp_df_D) <- c("Sample", "Cluster", "Percentage")

p10 <- ggplot(cluster_comp_df_D, aes(x = Sample, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Group D - Cluster Composition by Sample",
       y = "Percentage of Spots") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p10)

dev.off()

# ============================================================================
# PART 10: BATCH EFFECT ASSESSMENT
# ============================================================================

cat("\n", "="*80, "\n")
cat("BATCH EFFECT ASSESSMENT\n")
cat("="*80, "\n\n")

cat("Assessing potential batch effects and sample-specific biases...\n\n")

# Function to assess batch effects
assess_batch_effects <- function(seurat_obj, group_name) {

  cat(paste0("Analyzing ", group_name, ":\n"))

  # 1. Check if samples separate in PC space
  pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")[, 1:10]
  sample_ids <- seurat_obj$sample_id

  # Run simple ANOVA on first 10 PCs
  pc_pvalues <- sapply(1:10, function(pc) {
    summary(aov(pca_embeddings[, pc] ~ sample_ids))[[1]][["Pr(>F)"]][1]
  })

  cat(paste0("  PC1-10 association with sample (p-values):\n"))
  for(i in 1:10) {
    sig <- ifelse(pc_pvalues[i] < 0.001, " ***",
                  ifelse(pc_pvalues[i] < 0.01, " **",
                         ifelse(pc_pvalues[i] < 0.05, " *", "")))
    cat(paste0("    PC", i, ": ", formatC(pc_pvalues[i], format = "e", digits = 2), sig, "\n"))
  }

  if(sum(pc_pvalues < 0.05) >= 5) {
    cat("\n  WARNING: Strong batch effects detected!\n")
    cat("  Consider:\n")
    cat("    1. Harmony integration\n")
    cat("    2. Regressing out sample_id in SCTransform\n")
    cat("    3. ComBat batch correction\n\n")
  } else if(sum(pc_pvalues < 0.05) >= 2) {
    cat("\n  CAUTION: Moderate batch effects detected.\n")
    cat("  Monitor downstream analyses for sample-specific clusters.\n\n")
  } else {
    cat("\n  GOOD: Minimal batch effects detected.\n\n")
  }

  # 2. Check QC metric differences between samples
  qc_by_sample <- seurat_obj@meta.data %>%
    group_by(sample_id) %>%
    summarize(
      mean_umi = mean(nCount_Spatial),
      mean_genes = mean(nFeature_Spatial),
      mean_mt = mean(percent.mt),
      n_spots = n()
    )

  cat("  QC metrics by sample:\n")
  print(qc_by_sample)
  cat("\n")

  # Test for significant differences
  umi_test <- kruskal.test(nCount_Spatial ~ sample_id, data = seurat_obj@meta.data)
  genes_test <- kruskal.test(nFeature_Spatial ~ sample_id, data = seurat_obj@meta.data)

  cat(paste0("  Kruskal-Wallis test for UMI differences: p = ",
             formatC(umi_test$p.value, format = "e", digits = 2), "\n"))
  cat(paste0("  Kruskal-Wallis test for Gene differences: p = ",
             formatC(genes_test$p.value, format = "e", digits = 2), "\n\n"))

  if(umi_test$p.value < 0.001 | genes_test$p.value < 0.001) {
    cat("  WARNING: Significant QC differences between samples!\n")
    cat("  This may confound biological comparisons.\n\n")
  }

  return(list(pc_pvalues = pc_pvalues,
              qc_summary = qc_by_sample))
}

# Assess both groups
batch_assessment_A <- assess_batch_effects(group_A_merged, "Group A")
batch_assessment_D <- assess_batch_effects(group_D_merged, "Group D")

# ============================================================================
# PART 11: OPTIONAL - HARMONY INTEGRATION (if batch effects detected)
# ============================================================================

# Uncomment and run if severe batch effects are detected

# library(harmony)
#
# cat("Running Harmony integration to correct batch effects...\n")
#
# # Group A
# group_A_harmony <- RunHarmony(group_A_merged,
#                               group.by.vars = "sample_id",
#                               reduction = "pca",
#                               dims.use = 1:30,
#                               verbose = FALSE)
#
# group_A_harmony <- RunUMAP(group_A_harmony, reduction = "harmony", dims = 1:30)
# group_A_harmony <- FindNeighbors(group_A_harmony, reduction = "harmony", dims = 1:30)
# group_A_harmony <- FindClusters(group_A_harmony, resolution = 0.5)
#
# # Group D
# group_D_harmony <- RunHarmony(group_D_merged,
#                               group.by.vars = "sample_id",
#                               reduction = "pca",
#                               dims.use = 1:30,
#                               verbose = FALSE)
#
# group_D_harmony <- RunUMAP(group_D_harmony, reduction = "harmony", dims = 1:30)
# group_D_harmony <- FindNeighbors(group_D_harmony, reduction = "harmony", dims = 1:30)
# group_D_harmony <- FindClusters(group_D_harmony, resolution = 0.5)
#
# # Save harmony-corrected objects
# saveRDS(group_A_harmony, "Group_A_Harmony_Spatial.rds")
# saveRDS(group_D_harmony, "Group_D_Harmony_Spatial.rds")

# ============================================================================
# PART 12: FIND CLUSTER MARKERS
# ============================================================================

cat("="*80, "\n")
cat("FINDING CLUSTER MARKERS\n")
cat("="*80, "\n\n")

# Find markers for Group A
cat("Finding markers for Group A clusters...\n")
Idents(group_A_merged) <- "seurat_clusters"

markers_A <- FindAllMarkers(group_A_merged,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.5,
                            test.use = "wilcox",
                            verbose = FALSE)

# Get top 10 markers per cluster
top10_A <- markers_A %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(markers_A, "Group_A_Cluster_Markers.csv", row.names = FALSE)
write.csv(top10_A, "Group_A_Top10_Markers.csv", row.names = FALSE)

cat(paste0("  Found markers for ", length(unique(markers_A$cluster)), " clusters\n"))
cat("  Saved to Group_A_Cluster_Markers.csv\n\n")

# Find markers for Group D
cat("Finding markers for Group D clusters...\n")
Idents(group_D_merged) <- "seurat_clusters"

markers_D <- FindAllMarkers(group_D_merged,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.5,
                            test.use = "wilcox",
                            verbose = FALSE)

top10_D <- markers_D %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(markers_D, "Group_D_Cluster_Markers.csv", row.names = FALSE)
write.csv(top10_D, "Group_D_Top10_Markers.csv", row.names = FALSE)

cat(paste0("  Found markers for ", length(unique(markers_D$cluster)), " clusters\n"))
cat("  Saved to Group_D_Cluster_Markers.csv\n\n")

# ============================================================================
# PART 13: VISUALIZE TOP MARKERS
# ============================================================================

cat("Generating marker visualization plots...\n")

pdf("03_Cluster_Markers.pdf", width = 20, height = 14)

# Group A - top 3 markers per cluster
cat("  Plotting Group A markers...\n")
top3_A <- markers_A %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

for(clust in unique(top3_A$cluster)) {
  genes <- top3_A %>% filter(cluster == clust) %>% pull(gene)

  # Spatial plots
  p_spatial <- SpatialFeaturePlot(group_A_merged,
                                  features = genes,
                                  ncol = 3,
                                  pt.size.factor = 1.6)

  # Violin plots
  p_vln <- VlnPlot(group_A_merged, features = genes, ncol = 3, pt.size = 0)

  print(p_spatial)
  print(p_vln)
}

# Group D - top 3 markers per cluster
cat("  Plotting Group D markers...\n")
top3_D <- markers_D %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

for(clust in unique(top3_D$cluster)) {
  genes <- top3_D %>% filter(cluster == clust) %>% pull(gene)

  # Spatial plots
  p_spatial <- SpatialFeaturePlot(group_D_merged,
                                  features = genes,
                                  ncol = 3,
                                  pt.size.factor = 1.6)

  # Violin plots
  p_vln <- VlnPlot(group_D_merged, features = genes, ncol = 3, pt.size = 0)

  print(p_spatial)
  print(p_vln)
}

dev.off()

# ============================================================================
# PART 14: DIFFERENTIAL EXPRESSION BETWEEN GROUPS
# ============================================================================

cat("\n", "="*80, "\n")
cat("DIFFERENTIAL EXPRESSION: GROUP A vs GROUP D\n")
cat("="*80, "\n\n")

cat("NOTE: This comparison assumes matched tissue regions/cell types.\n")
cat("If tissue composition differs, results may be confounded.\n\n")

# For meaningful group comparison, we need to merge ALL samples
cat("Merging all samples for group comparison...\n")

all_samples_merged <- merge(x = spatial_filtered[[1]],
                            y = spatial_filtered[-1],
                            add.cell.ids = names(spatial_filtered),
                            project = "All_Samples")

# Add group metadata
all_samples_merged$group <- ifelse(grepl("^A", all_samples_merged$sample_id), "Group_A", "Group_D")

cat("Re-processing merged dataset...\n")
all_samples_merged <- SCTransform(all_samples_merged,
                                  assay = "Spatial",
                                  vars.to.regress = c("percent.mt", "sample_id"),
                                  verbose = FALSE)

all_samples_merged <- RunPCA(all_samples_merged, assay = "SCT", npcs = 50, verbose = FALSE)
all_samples_merged <- RunUMAP(all_samples_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
all_samples_merged <- FindNeighbors(all_samples_merged, reduction = "PCA", dims = 1:30, verbose = FALSE)
all_samples_merged <- FindClusters(all_samples_merged, resolution = 0.5, verbose = FALSE)

# Find DEGs between groups
cat("Finding differentially expressed genes between groups...\n")
Idents(all_samples_merged) <- "group"

group_degs <- FindMarkers(all_samples_merged,
                          ident.1 = "Group_A",
                          ident.2 = "Group_D",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          verbose = FALSE)

# Add gene names and filter significant
group_degs$gene <- rownames(group_degs)
group_degs_sig <- group_degs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

write.csv(group_degs_sig, "Group_A_vs_D_DEGs_significant.csv", row.names = FALSE)

cat(paste0("  Found ", nrow(group_degs_sig), " significant DEGs (FDR < 0.05)\n"))
cat(paste0("  Upregulated in Group A: ", sum(group_degs_sig$avg_log2FC > 0), "\n"))
cat(paste0("  Upregulated in Group D: ", sum(group_degs_sig$avg_log2FC < 0), "\n"))
cat("  Saved to Group_A_vs_D_DEGs_significant.csv\n\n")

# Visualize top DEGs
pdf("04_Group_Comparison.pdf", width = 16, height = 10)

# Top 20 DEGs
top_degs <- group_degs_sig %>%
  top_n(n = 20, wt = abs(avg_log2FC))

# Heatmap
p_heatmap <- DoHeatmap(all_samples_merged,
                       features = top_degs$gene,
                       group.by = "group") +
  ggtitle("Top 20 DEGs: Group A vs D")

print(p_heatmap)

# Volcano plot
p_volcano <- ggplot(group_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(title = "Differential Expression: Group A vs D",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)") +
  theme_minimal()

print(p_volcano)

# Feature plots for top DEGs
top5_A <- group_degs_sig %>% filter(avg_log2FC > 0) %>% top_n(n = 5, wt = avg_log2FC)
top5_D <- group_degs_sig %>% filter(avg_log2FC < 0) %>% top_n(n = 5, wt = -avg_log2FC)

p_features <- FeaturePlot(all_samples_merged,
                          features = c(head(top5_A$gene, 4), head(top5_D$gene, 4)),
                          reduction = "umap",
                          ncol = 4)
print(p_features)

# Violin plots
p_vln_compare <- VlnPlot(all_samples_merged,
                         features = c(head(top5_A$gene, 3), head(top5_D$gene, 3)),
                         group.by = "group",
                         ncol = 3,
                         pt.size = 0)
print(p_vln_compare)

dev.off()

# Save final merged object
saveRDS(all_samples_merged, "All_Samples_Merged_Spatial.rds")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n\n", "="*80, "\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("="*80, "\n\n")

cat("SAVED OBJECTS:\n")
cat("  - Individual filtered samples: Filtered_*_Spatial.rds\n")
cat("  - Group A merged: Group_A_Merged_Spatial.rds\n")
cat("  - Group D merged: Group_D_Merged_Spatial.rds\n")
cat("  - All samples merged: All_Samples_Merged_Spatial.rds\n\n")

cat("SAVED RESULTS:\n")
cat("  - QC filtering summary: QC_Filtering_Summary.csv\n")
cat("  - Group A markers: Group_A_Cluster_Markers.csv\n")
cat("  - Group D markers: Group_D_Cluster_Markers.csv\n")
cat("  - Group comparison DEGs: Group_A_vs_D_DEGs_significant.csv\n\n")

cat("GENERATED PLOTS:\n")
cat("  - 00_PreQC_Summary.pdf\n")
cat("  - 01_PostQC_Summary.pdf\n")
cat("  - 02_Merged_Groups_Analysis.pdf\n")
cat("  - 03_Cluster_Markers.pdf\n")
cat("  - 04_Group_Comparison.pdf\n\n")

cat("NEXT STEPS:\n")
cat("  1. Review QC filtering results and adjust thresholds if needed\n")
cat("  2. Check for batch effects - consider Harmony if needed\n")
cat("  3. Annotate clusters using marker genes\n")
cat("  4. Validate key findings with spatial feature plots\n")
cat("  5. Consider pathway enrichment analysis on DEGs\n\n")

cat("="*80, "\n\n")
