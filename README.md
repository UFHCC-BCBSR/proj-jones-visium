# 10x Visium FFPE Spatial Transcriptomics Analysis Pipeline

## Project Overview

This repository contains R Markdown-based analysis pipelines for processing and analyzing 10x Genomics Visium FFPE spatial transcriptomics data from placental tissue samples. The pipeline implements standard Seurat workflows with FFPE-specific quality control considerations and pseudobulk differential expression analysis.

## Goals

1. **Quality Assessment:** Generate comprehensive QC reports with FFPE-specific metrics and thresholds
2. **Sample Filtering:** Apply QC filtering to identify high-quality samples and spots
3. **Batch Correction:** Integrate multiple samples within groups using Harmony
4. **Cluster Analysis:** Identify spatially-resolved cell populations and their marker genes
5. **Differential Expression:** Perform statistically appropriate pseudobulk DE analysis between experimental groups

## Input Data Requirements

### Raw Data

These files are hidden from the repo. Contact `hkates@ufl.edu` for access.

- 10x Genomics Visium Space Ranger output directories for each sample
- Expected structure per sample:
  ```
  {Sample}_AM_Spatial_Seq/
  └── outs/
      ├── filtered_feature_bc_matrix.h5
      ├── spatial/
      └── metrics_summary.csv
  ```

### Metadata

These files are hidden from the repo. Contact `hkates@ufl.edu` for access.

- `data/combined_metrics_summary.csv` - Space Ranger QC metrics for all samples
- `data/web-summaries/` - Space Ranger HTML reports (optional, for reference)

## Analysis Workflow

### 1. QC Report Generation (`R/qc-summary.Rmd`)

Generates detailed quality control report with:
- FFPE-specific QC metric definitions and interpretations
- Technical/molecular/biological explanations for each metric
- Sample-level quality assessment and recommendations
- Evidence-based thresholds from 10x Genomics documentation

**Key Metrics Evaluated:**
- Median UMI counts per spot
- Median genes per spot  
- Reads confidently mapped to probe set
- Split-mapped reads
- Sequencing saturation
- Mean reads per spot
- Spots under tissue
- Total genes detected

### 2. Spatial Analysis Pipeline (`R/explore-visium.Rmd`)

Comprehensive analysis including:

#### Part 1-4: Quality Control
- Load samples and calculate QC metrics
- Apply standard Seurat filtering (>200 genes, <25% MT, outlier removal)
- Visualize pre/post-QC metrics spatially and in feature space

#### Part 5-6: Normalization and Integration
- SCTransform normalization with MT% regression
- Sample merging by experimental group
- Harmony batch correction for technical variation

#### Part 7-11: Clustering and Annotation
- Identify spatial clusters using Louvain algorithm
- Find cluster-specific marker genes
- Annotate with known placental cell type markers
- Visualize markers in spatial and UMAP space

#### Part 12: Differential Expression
- Pseudobulk aggregation by sample (proper biological replicates)
- edgeR-based differential expression for biologically meaningful clusters
- Focus on clusters with interpretable marker expression

## Expected Outputs

These files are hidden from the repo. Contact `hkates@ufl.edu` for access.

### QC Reports
- `results/qc-summary.html` - Interactive HTML report with detailed QC assessments
- `results/QC_Filtering_Summary.csv` - Filtering statistics per sample

### Processed Seurat Objects
- `results/Filtered_{Sample}_Spatial.rds` - QC-filtered individual samples
- `results/Group_A_Merged_Spatial.rds` - Integrated Group A samples
- `results/Group_D_Merged_Spatial.rds` - Integrated Group D samples  
- `results/All_Samples_Merged_Spatial.rds` - All samples combined

### Marker Gene Lists
- `results/Group_A_Cluster_Markers.csv` - All cluster markers for Group A
- `results/Group_A_Top10_Markers.csv` - Top 10 markers per cluster
- `results/Group_D_Cluster_Markers.csv` - All cluster markers for Group D
- `results/Group_D_Top10_Markers.csv` - Top 10 markers per cluster
- `results/Cluster_Markers_Visualization.pdf` - Spatial/violin plots of marker genes

### Differential Expression Results
- `results/Pseudobulk_Cluster3_SignificantDEGs.csv` - Significant DEGs for cluster 3 (FDR < 0.05)
- `results/Pseudobulk_Cluster3_AllGenes.csv` - Full results for cluster 3
- `results/Pseudobulk_Cluster5_SignificantDEGs.csv` - Significant DEGs for cluster 5 (FDR < 0.05)
- `results/Pseudobulk_Cluster5_AllGenes.csv` - Full results for cluster 5

## Usage

### Prerequisites

**R version:** >= 4.3.0

**Required R packages:**
```r
# Core analysis
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "tidyverse"))

# Report generation
install.packages(c("knitr", "kableExtra", "rmarkdown"))

# Batch correction
install.packages("harmony")

# Differential expression
install.packages("edgeR")
```

### Running the Analysis

1. **Generate QC Report:**
   ```bash
   cd R/
   Rscript -e "rmarkdown::render('qc-summary.Rmd')"
   ```

2. **Run Spatial Analysis Pipeline:**
   ```bash
   Rscript -e "rmarkdown::render('explore-visium.Rmd')"
   ```

**Note:** The pipeline uses conditional loading - if processed `.rds` objects exist, they will be loaded rather than reprocessed. Delete specific `.rds` files in `results/` to force reprocessing of those steps.

## Pipeline Design Decisions

### Why Standard QC Filtering?
Initial attempts at complex adaptive filtering showed minimal improvement over standard thresholds. We use simple, transparent, community-standard filters that are reproducible and well-documented.

### Why Exclude Certain Samples?
Samples A2 and D2 show extreme quality issues (median UMI <300, critically low mapping rates) that cannot be salvaged by filtering. Sample A1 has critically low confident mapping (22.6%) suggesting severe technical failure. These are excluded from group-level analyses to prevent technical artifacts from dominating biological signal.

### Why Pseudobulk DE?
Pseudobulk aggregation provides:
- Proper statistical inference at the biological replicate (sample) level
- Control for within-sample spot correlation
- Compatibility with established bulk RNA-seq DE frameworks (edgeR)
- Follows current best practices (Squair et al. 2021, Nature Communications)

### Why Focus on Clusters 3 and 5?
Comprehensive marker analysis revealed only these clusters show biologically interpretable, spatially coherent gene expression patterns. Other clusters appear driven by technical variation (sample quality differences, degradation gradients) rather than distinct cell populations.

## Repository Structure

```
.
├── R/
│   ├── qc-summary.Rmd              # QC report generation
│   ├── explore-visium.Rmd          # Main analysis pipeline
│   └── combine-qc-metrics.bash     # Helper script to aggregate metrics
├── data/                            # Not tracked (contact for access)
├── results/                         # Not tracked (generated by pipeline)
└── README.md                        # This file
```

## Key References

- **Seurat v5:** Hao et al. 2024, Nature Biotechnology
- **SCTransform:** Hafemeister & Satija 2019, Genome Biology  
- **Harmony:** Korsunsky et al. 2019, Nature Methods
- **Pseudobulk DE:** Squair et al. 2021, Nature Communications
- **10x Visium FFPE:** 10x Genomics Technical Notes CG000408, CG000495

## Contact

For questions or to request access to raw data paths, contact `hkates@ufl.edu`.

## License

This analysis pipeline is provided as-is for research purposes.
