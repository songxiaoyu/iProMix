
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `iProMix`

## Cell-Type Specific Association Analysis for Bulk Profiling Data

Human tissues are comprised of multiple cell types with varying
compositions. The bulk-tissue profiles of multiple -omic data types
(e.g. DNA methylation, mRNA, proteomics) are impacted by the cell-type
composition heterogeneity, as their levels in different cell types may
be different. iProMix decomposes data from single/multiple -omic data
types (e.g. DNA methylation, mRNA, proteomics) and evaluates their
cell-type specific dependences. A major difference of iProMix from
previous studies is that it allows association analysis on two data
types that are both affected by cell types (e.g. mRNA vs. protein). It
builds in features to improve cell type composition estimation if
existing estimates are not satisfactory. It also takes into
consideration the effects of decomposition and biased input on
hypothesis tests, and generates valid inference in non-asymptotic
settings.

### Installation

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/iProMix")
```

the most recent officially-released version from CRAN with

``` r
install.packages("iProMix")
```
