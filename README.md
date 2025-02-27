
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RAmP (R package for the Analysis of \[modified\] Proteomes)

<!-- badges: start -->
<!-- badges: end -->

The goal of RAmP is to perform a one-stop proteomic and PTM-proteomic
analysis. It also includes differential expression analysis using limma,
imputation and filtering options (from DEP), batch corrections, PCA and
other quality control options. It also includes post-analysis
visualization tools such as volcano plots, heatmaps, GO term analysis
(both ORA and GSEA), KEGG analysis and corresponding visaulizations and
motif search for PTM-based proteomics.

## Installation

You can install the development version of RAmP from
[GitHub](https://github.com/anuroopv/RAmP) with:

``` r
# install.packages("devtools")
devtools::install_github("anuroopv/RAmP", dependencies = TRUE)
library(RAmP)
```

Please note that RAmP requires an R version higher than 4.3.1. Please
intsall the latest R version to avoid any compatibility issues with
dependent packages.

## Accessimg the data

Example proteome and acetylome data from Venkatasubramani AV, et al.,
2023 (PMID: 37724628) can be accessed using:

``` r
data("sampleTable")
data("fasta")
data("protData")
data("enrichedData")
```

## Further information

For more information about RAmP, check the website:
<https://github.com/anuroopv/RAmP> or refer to the help options after
installation.
