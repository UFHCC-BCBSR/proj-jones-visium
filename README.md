# 10x Visium FFPE Spatial Transcriptomics Exploration

## Project Overview

This repository contains R Markdown-based analysis pipelines for
processing and analyzing 10x Genomics Visium FFPE spatial
transcriptomics data from placental tissue samples. The pipeline
implements standard Seurat workflows with FFPE-specific quality control
considerations and pseudobulk differential expression analysis.

## Goals

1.  **Quality Assessment:** Generate comprehensive QC reports with
    FFPE-specific metrics and thresholds
2.  **Sample Filtering:** Apply QC filtering to identify high-quality
    samples and spots
3.  **Batch Correction:** Integrate multiple samples within groups using
    Harmony
4.  **Cluster Analysis:** Identify spatially-resolved cell populations
    and their marker genes

## Input Data Requirements

### Raw Data

These files are hidden from the repo. Contact `hkates@ufl.edu` for
access.

-   10x Genomics Visium Space Ranger output directories for each sample

-   Expected structure per sample:

    ```         
    {Sample}_AM_Spatial_Seq/
    └── outs/
        ├── filtered_feature_bc_matrix.h5
        ├── spatial/
        └── metrics_summary.csv
    ```

### Metadata

These files are hidden from the repo. Contact `hkates@ufl.edu` for
access.

-   `data/combined_metrics_summary.csv` - Space Ranger QC metrics for
    all samples
-   `data/web-summaries/` - Space Ranger HTML reports (optional, for
    reference)

## Analysis Workflow

### 1. QC Report Generation (`R/qc-summary.Rmd`)

Generates detailed quality control report with: - FFPE-specific QC
metric definitions and interpretations - Technical/molecular/biological
explanations for each metric - Sample-level quality assessment and
recommendations - Evidence-based thresholds from 10x Genomics
documentation

**Key Metrics Evaluated:** - Median UMI counts per spot - Median genes
per spot\
- Reads confidently mapped to probe set - Split-mapped reads -
Sequencing saturation - Mean reads per spot - Spots under tissue - Total
genes detected

### 2. Spatial Analysis Pipeline (`R/visium-minimal.Rmd`)

Exploratory analysis including:

#### Quality Control

#### Normalization and Integration

#### Clustering and Visualization

## Expected Outputs

These files are hidden from the repo. Contact `hkates@ufl.edu` for
access.

### Reports

-   `results/qc-summary.html` - HTML report with detailed QC assessments
-   `results/visium-minimal.html` - HTML report with basic normalization, integration, and visualizations

### Processed Seurat Objects

-   `results/merged-seurat.rds` - Integrated samples R object to load into R for further analysis
-   `results/integrated_spatial_analysis.cloupe` - cloupe file of
    integrated sample data for exploration in 10X Loupe Browser

## Usage

### Prerequisites

**R version:** \>= 4.3.0

**Required R packages:**

``` r
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

1. Ensure you have access to files in data/ that are hidden from this repo.

2.  **Generate QC Report:**

    ``` bash
    cd R/
    Rscript -e "rmarkdown::render('qc-summary.Rmd')"
    ```

3.  **Run Spatial Analysis Pipeline:**

    ``` bash
    Rscript -e "rmarkdown::render('visium-minimal.Rmd')"
    ```

## Repository Structure

```         
.
├── data
│ ├── combined_metrics_summary.csv
│ └── web-summaries
│     ├── A1_web_summary.html
│     ├── A2_web_summary.html
│     ├── A3_web_summary.html
│     ├── A4_web_summary.html
│     ├── D1_web_summary.html
│     ├── D2_web_summary.html
│     ├── D3_web_summary.html
│     └── D4_web_summary.html
├── R
│ ├── combine-qc-metrics.bash
│ ├── QC_Filtering_Summary.csv
│ ├── qc-summary.Rmd
│ └── visium-minimal.Rmd
├── README.md
└── results
    ├── integrated_spatial_analysis.cloupe
    ├── merged-seurat.Rds
    ├── qc-summary.html
    └── visium-minimal.html
```

## Key References

-   **Seurat v5:** Hao et al. 2024, Nature Biotechnology
-   **SCTransform:** Hafemeister & Satija 2019, Genome Biology\
-   **Harmony:** Korsunsky et al. 2019, Nature Methods
-   **Pseudobulk DE:** Squair et al. 2021, Nature Communications
-   **10x Visium FFPE:** 10x Genomics Technical Notes CG000408, CG000495

## Contact

For questions or to request access to raw data paths, contact
`hkates@ufl.edu`.

## License

This analysis pipeline is provided as-is for research purposes.
