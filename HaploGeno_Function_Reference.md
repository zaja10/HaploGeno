# **HaploGeno Function Reference**

This document provides a comprehensive list of all functions available in the **HaploGeno** package, categorized by their role in the workflow.

## **1\. Standalone Utility Functions**

These functions are called directly from the R console to initialize projects or load data.

| Function | Purpose |
| :---- | :---- |
| load\_demo\_data() | Loads a pre-packaged demo dataset (50 samples, 500 markers) to test the package features immediately. |
| load\_haplo\_project(path) | Loads a previously saved HaploObject from an .rds file and automatically reconnects the file-backed genotype matrix. |
| init\_from\_met(met\_object, ...) | A bridge function that initializes a HaploObject directly from Multi-Environment Trial (MET) or Factor Analytic (FA) model outputs, handling ID alignment automatically. |

## **2\. HaploObject Class Methods**

These functions are methods of the HaploObject R6 class. They must be called on an instantiated object (e.g., obj$method\_name()).

### **Data Management & Setup**

| Method | Purpose |
| :---- | :---- |
| initialize(backing\_file\_path) | Constructor method. Creates a new object and sets up the file-backed storage path. |
| import\_genotypes(matrix\_or\_path) | Imports genotype data from PLINK (.bed), VCF, CSV/TXT, or R matrices into the internal file-backed matrix (FBM). |
| load\_map(map\_data) | Loads marker map information (Chromosome, Position, Alleles). |
| load\_pheno(pheno\_data) | Loads phenotypic data corresponding to the genotypes. |
| impute\_genotypes(method) | Fills missing genotype data using "expectation" (preserves variance) or "mean" (rounds to integer) methods. |
| filter\_monomorphic() | Removes markers with zero variance to ensure statistical stability. |
| save\_project(path) | Saves the current state of the object to an .rds file, preserving backing file links. |

### **Core Genomic Analysis**

| Method | Purpose |
| :---- | :---- |
| define\_haploblocks(method, ...) | Segments the genome into haplotype blocks using an LD scan or fixed windows. |
| estimate\_marker\_effects(lambda) | Estimates individual marker effects using Ridge Regression (optional auto-lambda via CV). |
| calculate\_local\_gebv(n\_cores) | Computes Local Genomic Estimated Breeding Values (GEBVs) for every defined block. |
| test\_significance() | Runs an empirical permutation test to determine the p-value of block variances. |
| calculate\_pve(adjust) | Calculates the Phenotypic Variance Explained (PVE) by each block, with optional de-regression adjustment. |
| encode\_haplotypes(n\_cores) | Converts raw genotype blocks into integer-encoded haplotype IDs. |

### **Structural & Factor Analysis**

| Method | Purpose |
| :---- | :---- |
| analyze\_block\_structure(top\_n, ...) | Performs Factor Analysis (Latent Regression) on top blocks to find underlying genomic gradients. |
| get\_model\_correlation() | Reconstructs the genetic correlation matrix implied by the Factor Analysis results. |
| analyze\_structural\_drivers() | Regresses factor communality against variance and genomic position to identify structural constraints (e.g., centromeres). |

### **Haplotype Discovery & Selection**

| Method | Purpose |
| :---- | :---- |
| analyze\_hoi(block\_id) | "Haplotype of Interest" analysis. Analyzes effect distributions within a block to identify distinct haplotype groups. |
| identify\_superior\_haplotypes(top\_n) | Identifies the specific haplotype allele with the highest positive effect in top blocks. |
| get\_haplo\_extremes(top\_n) | Returns lists of genotypes carrying the most extreme (positive/negative) haplotypes for top blocks. |
| score\_stacking(superior\_haplos) | Calculates a score for individuals based on the count of superior haplotypes they possess. |

### **Global Prediction (Genomic Prediction)**

| Method | Purpose |
| :---- | :---- |
| compute\_hrm(n\_cores) | Computes the Haplotype Relationship Matrix (HRM) based on shared haplotype blocks. |
| fit\_krr(lambda, ...) | Fits a Kernel Ridge Regression (KRR) model using the HRM for phenotype prediction. |
| cross\_validate(k, lambdas) | Performs k-fold cross-validation to optimize hyperparameters for the KRR model. |

### **Visualization Methods**

| Method | Purpose |
| :---- | :---- |
| plot\_manhattan(threshold) | Generates a standard Base R Manhattan plot of block significance or PVE. |
| plot\_manhattan\_gg(threshold) | Generates a ggplot2 version of the Manhattan plot. |
| plot\_volcano() | Plots Block Variance vs. Significance (-log10 P) to identify high-impact regions. |
| plot\_factor\_heatmap() | Visualizes the genetic correlation matrix of top blocks derived from Factor Analysis. |
| plot\_fa\_genome() | Displays a "stacked" Manhattan plot showing latent factor loadings across the genome. |
| plot\_scree() | Scree plot showing the variance explained by extracted latent factors. |
| plot\_communality() | plots communality vs. genomic position. |
| plot\_haplo\_biplot(top\_n, ...) | A Genomic Architecture Biplot overlaying high-variance block vectors on HRM PCA. |
| plot\_haplo\_profile(genotypes) | A mosaic heatmap showing the haplotype composition of selected lines across top blocks. |
| plot\_hoi\_distribution(block\_id) | Density plot of scaled effects for a specific block, highlighting separation. |
| plot\_gebv\_image() | Heatmap of local GEBV values (Individuals x Blocks). |
| plot\_hrm() | Heatmap of the Haplotype Relationship Matrix. |
| plot\_pca() | PCA plot derived from the HRM. |
| plot\_block\_sizes() | Histogram showing the distribution of block sizes (in number of markers). |
| plot\_stacking\_trend() | Regression plot showing the relationship between stacking scores and phenotypes. |

## **3\. S3 Dispatch Methods**

Standard R functions that have been extended to work with HaploObject.

| Function | Purpose |
| :---- | :---- |
| print(x) | Prints a concise status summary of the object. |
| summary(x) | Prints a detailed "Breeding Report" including top driver blocks and analysis status. |
| plot(x, type=...) | Generic plotter. Supports types: 'manhattan', 'scree', 'profile', 'stacking'. |

## **4\. Internal / Rcpp Functions**

Optimized C++ helper functions not typically called by the user.

| Function | Purpose |
| :---- | :---- |
| encode\_block\_fast(mat) | Rapidly encodes a matrix of genotypes into integer haplotype IDs. |
| ridge\_solver\_cpp(X, y, lambda) | C++ implementation of a Ridge Regression solver. |

