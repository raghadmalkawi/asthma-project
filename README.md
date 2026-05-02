
# Biomarker Discovery: Transcriptomic Analysis of Asthmatics vs. Healthy Controls

A bioinformatics project identifying differentially expressed genes (DEGs) in nasal epithelial samples from asthmatic patients versus healthy controls, using RNA-seq data and the DESeq2 pipeline in R.


> **Authors:** Raghad Malkawi, Rebekah Lanning
> 
> **Course:** Indiana University Luddy School of Informatics, Computation and Engineering - Taught by Dr. Sarath Janga, PhD

---

## Background

Asthma is a chronic inflammatory respiratory disease affecting approximately 300 million people worldwide. More than 100 genes have been linked to asthma risk, yet the underlying molecular mechanisms remain incompletely understood. This project performs transcriptomic analysis to identify differentially expressed genes in nasal epithelial samples, with the goal of uncovering possible biomarkers relevant to asthma diagnosis and treatment.

---

## Dataset

- **Source:** [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) — Accession: **GSE97668**
- **Samples used:** 10 nasal epithelial samples (5 asthmatics, 5 healthy controls), all collected **7 days before** Rhinovirus HRV-A16 inoculation
- **Files needed:**
  - `GSE97668_raw_counts_GRCh38.p13_NCBI.tsv.gz` — Raw gene counts
  - `Human.GRCh38.p13.annot.tsv.gz` — Gene annotation (GeneID → Symbol mapping)
  - `GSE97668_series_matrix.txt.gz` — Sample metadata

All three files can be downloaded directly from the GEO accession page linked above.

---

## Requirements

**R packages:**

| Package | Source | Purpose |
|---|---|---|
| `DESeq2` | Bioconductor | Differential expression analysis |
| `GEOquery` | Bioconductor | GEO data access |
| `clusterProfiler` | Bioconductor | KEGG pathway enrichment |
| `org.Hs.eg.db` | Bioconductor | Human gene ID mapping |
| `KEGGREST` | Bioconductor | KEGG database queries |
| `pheatmap` | CRAN | Heatmap visualization |
| `ggplot2` | CRAN | General plotting |
| `dplyr` | CRAN | Data wrangling |
| `tibble` | CRAN | Data frame utilities |
| `readr` | CRAN | File reading |

Install Bioconductor packages with:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GEOquery", "clusterProfiler", "org.Hs.eg.db", "KEGGREST"))
```

Install CRAN packages with:
```r
install.packages(c("ggplot2", "pheatmap", "dplyr", "tibble", "readr"))
```

---

## Project Structure

```
asthma-project/
├── advanced-seminar project.R   # Main analysis script
├── README.md                    # This file
└── outputs/                     # Generated output files (after running)
    ├── transformed_metadata_file.csv
    ├── KEGG_Enrichment_Results.csv
    └── plots/
        ├── MA_plot.png
        ├── volcano_plot.png
        ├── heatmap.png
        ├── bar_plot_DEGs.png
        ├── top10_genes_bar.png
        ├── KEGG_barplot.png
        ├── KEGG_dotplot.png
        └── PCA_plot.png
```

---

## Pipeline Overview

### 1. Data Loading & Preprocessing
- Load raw count data and gene annotation file
- Map GeneID to gene symbols via inner join
- Sum counts for duplicate gene symbols
- Filter out low-expression genes (minimum count of 10 in at least half of samples)

### 2. Metadata Processing
- Parse GEO series matrix file to extract sample metadata
- Select the 10 relevant pre-inoculation samples
- Assign group labels: `Asthma` or `control`

### 3. Differential Expression Analysis (DESeq2)
- Build a `DESeqDataSet` object using raw counts and sample metadata
- Set `control` as the reference level
- Run DESeq2 normalization and statistical testing
- Filter results: adjusted p-value < 0.05 and |log2FoldChange| > 1

### 4. Visualizations
| Plot | Purpose |
|---|---|
| **MA Plot** | Log2 fold change vs. mean expression across all genes |
| **Volcano Plot** | Significance vs. magnitude of expression change |
| **Heatmap** | Expression patterns of all 37 DEGs across 10 samples |
| **Bar Plot (DEGs)** | Proportion of upregulated vs. downregulated genes |
| **Top 10 Genes Bar** | Log2 fold change of the 10 most significant DEGs |
| **PCA Plot** | Sample clustering based on overall expression profiles |

### 5. KEGG Pathway Enrichment Analysis
- Map significant DEG symbols to Entrez IDs via `bitr()`
- Run `enrichKEGG()` with adjusted p-value cutoff of 0.05
- Visualize enriched pathways as a bar plot and dot plot
- Export results to `KEGG_Enrichment_Results.csv`

---

## Results Summary

- **37 DEGs** identified between asthmatic and control nasal epithelial samples
  - **35 upregulated** (94.6%)
  - **2 downregulated** (5.4%)

- **Top biomarker candidates:** `TPSB2` and `TPSAB1` — both previously identified as hub genes for asthma in independent studies

- **Enriched KEGG Pathways:**
  - **Osteoclast Differentiation Pathway** (4 DEGs: *LILBR2, FCGR2C, FOSB, SIRPA*)
  - **Complement and Coagulation Cascades** (3 DEGs: *ITGAX, ITGAM, THBD*)

- PC1 and PC2 together account for **62.7%** of the variance in the dataset, with partial separation between asthma and control groups visible in the PCA

---

## How to Run

1. Clone this repository:
   ```bash
   git clone https://github.com/raghadmalkawi/asthma-project.git
   cd asthma-project
   ```

2. Download the required data files from [GEO accession GSE97668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97668) and place them in your working directory.

3. Update the file paths at the top of `advanced-seminar project.R` to match your local file locations:
   ```r
   raw_data <- read.table(gzfile("path/to/GSE97668_raw_counts_GRCh38.p13_NCBI.tsv.gz"), ...)
   file_path <- "path/to/Human.GRCh38.p13.annot.tsv.gz"
   lines <- readLines(gzfile("path/to/GSE97668_series_matrix.txt.gz"))
   ```

4. Open `advanced-seminar project.R` in RStudio and run the script section by section, or run the entire script at once.

5. Output CSV files and plots will be saved to your working directory (use `getwd()` to check the location).

---

## Contact

| Author | Email |
|---|---|
| Raghad Malkawi | ragmalka@iu.edu |
| Rebekah Lanning | rebelann@iu.edu |
