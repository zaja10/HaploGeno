# HaploGeno

**HaploGeno** implements a high-performance framework for haplotype-based genomic prediction and QTL discovery. Designed to address the computational challenges of large-scale genomic datasets, it leverages file-backed matrices (via `bigstatsr`) to enable efficient analysis on standard hardware.

## Motivation

Traditional genomic prediction models often rely on single-marker effects, potentially missing the complex epistatic interactions captured by haplotypes. However, constructing and analyzing haplotypes for dense marker panels can be computationally prohibitive. `HaploGeno` bridges this gap by providing a scalable, memory-efficient pipeline for:

1.  **Constructing haplotypes** based on Linkage Disequilibrium (LD) or fixed windows.
2.  **Modeling complex traits** using Kernel Ridge Regression (KRR) on haplotype alleles.
3.  **Identifying significant regions** through local Genomic Estimated Breeding Values (GEBV).

## Features

### 1. Data Handling

- **File-Backed Matrices**: Efficiently handle large genotype matrices that exceed RAM using `bigstatsr`.
- **Flexible Import**: Load genotypes, marker maps, and phenotypes from various formats.

### 2. Haplotype Construction

- **Fixed Window Blocking**: Define haplotype blocks based on a fixed number of markers.
- **LD-Based Blocking**: Define blocks using Linkage Disequilibrium (LD) scan with customizable $r^2$ threshold and tolerance (Shaffer et al. 2025).
- **Haplotype Encoding**: Map unique haplotype alleles within blocks to integer IDs.

### 3. Genomic Prediction (Global)

- **Kernel Ridge Regression (KRR)**:
  - Computes a **Haplotype Relationship Matrix (HRM)** based on shared haplotype alleles.
  - Fits a KRR model to predict phenotypes with statistical rigor.
  - Solves the dual problem: $\alpha = (K + \lambda I)^{-1} y$.

### 4. QTL Discovery (Local)

- **localGEBV Methodology** (Shaffer et al. 2025):
  - **Marker Effects**: Estimates SNP effects using Ridge Regression.
  - **Local GEBV**: Calculates Genomic Estimated Breeding Values for specific haploblocks.
  - **Significance Testing**: Implements a Scaled Inverse Chi-Squared test to assess block significance with `test_significance()`.

### 5. Statistical Engine (New)

- **Factor Analysis (Haplo-FA)**:
  - Performs internal Factor Analysis on local GEBVs to extract latent genomic gradients.
  - **Percent Variance Explained (PVE)**: Calculates phenotypic variance explained by each block.
  - **Structural Drivers**: Analyzes the relationship between Communality and Position.

### 6. Analysis & Visualization (Base R)

- **Visualization**:
  - **Biplot**: Genomic Architecture Biplot (`plot_haplo_biplot`).
  - **Factor Genome**: Stacked "Manhattan-like" plots for latent factors (`plot_fa_genome`).
  - **Communality**: Visualize structural constraints (`plot_communality`).
  - **Scree Plot**: Variance explained by extracted factors (`plot_scree`).
  - **Manhattan Plot**: Standard QTL visualization (`plot_manhattan`).

### 7. Selection & Profiling (New)

- **Haplotype Extremes**: Identify superior/inferior genotypes for crossing (`get_haplo_extremes`).
- **Haplotype Profile**: Mosaic heatmap of genetic formulas (`plot_haplo_profile`).

### 6. Demo Dataset

- **Quick Start**: Use `load_demo_data()` to load a pre-processed dataset for testing visualization and analysis functions immediately.

## Quick Start

```r
library(HaploGeno)

# Option A: Load Demo Data (Ready for Visualization)
haplo <- load_demo_data()
haplo$plot_manhattan()

# Option B: Start from Scratch
# 1. Initialize
haplo <- HaploObject$new(tempfile())
haplo$import_genotypes(geno_matrix)
haplo$load_map(map_data)
haplo$load_pheno(phenotypes)

# 2. Pre-processing
haplo$filter_monomorphic() # Remove zero-variance markers
haplo$impute_genotypes(method = "mean") # Impute missing values

# 3. Define Blocks (LD-based)
haplo$define_blocks_ld(r2_threshold = 0.5, tolerance = 2)

# --- Workflow A: Global Prediction (KRR) ---
haplo$encode_haplotypes()
haplo$compute_hrm()
alpha <- haplo$fit_krr(lambda = 0.1)

# --- Workflow B: QTL Discovery (localGEBV) ---
haplo$estimate_marker_effects(lambda = 0.1)
haplo$calculate_local_gebv()
results <- haplo$test_significance()

# 5. Visualize & Analyze
haplo$plot_manhattan()
haplo$plot_haplo_biplot(top_n = 10)
haplo$plot_fa_genome()

# 6. Factor Analysis & Structure
haplo$calculate_pve()
haplo$analyze_block_structure(factors = 3)
haplo$plot_scree()
haplo$plot_communality()

# 7. Selection
extremes <- haplo$get_haplo_extremes(top_n = 10)
haplo$plot_haplo_profile(top_n_blocks = 20)

# 7. Save Project
haplo$save_project("my_analysis.rds")
```

## Future Work

- **Bayesian Models**: Implementation of BayesA/B/C for variable selection.
- **Multi-Trait Models**: Extension to multi-trait genomic prediction.
- **Interactive Visualization**: Integration with `plotly` for interactive exploration of results.

## References

Shaffer, et al. (2025). _Local genomic estimates provide a powerful framework for haplotype discovery_. biorxiv.
