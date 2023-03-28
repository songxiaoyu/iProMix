
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

### Load package

``` r
library(iProMix)
```

### Example 1: Correlation between one gene X and one gene Y in a cell type.

Notably, the estimate and inference might be inaccurate if the sample
size is small or the input cell-type composition estimation is
imprecise. If multiple gene pairs are under consideration, iProMix can
leverage the parallel nature of these genes to calculate the robust
eFDR.

``` r
set.seed(123)
n=500
pi=runif(n) # proportion of the cell type of interest
x=rnorm(n) # X level in the tissue. 
y=rnorm(n) # Y level in the tissue. 

# Estimate X-Y correlation in the cell type  of interest
est <- iProMix(y = y, x = x, pi = pi, cov=NULL)
est$cor.score1[2,1] 
# 0.03697112 

# Likelihood Ratio (LR) test on whether X-Y correlation in the cell type of interest is significantly different from zero.
t1 <- iProMix.LRT(y = y, x = x, pi = pi, cov=NULL, reduce1=c(1,2))
t1$LRT.pvalue
# [1] 0.6458332
```

### Example 2: Identify if any of the P genes in Y is associated with one gene X in a cell type.

In this example, iProMix uses the parallel nature of P genes to infer
eFDR, so P needs to be large. Both `iProMix.eFDR.PermReplace` and
`iProMix.eFDR.PermAdd` provide valid results under correctly specified
models, but when the input cell-type composition is inaccurately
estimated, `iProMix.eFDR.PermReplace` is more robust than the
`iProMix.eFDR.PermAdd` strategy.

``` r
library(iProMix)
set.seed(123)
n=500
P=10
pi=runif(n) # proportion of cell type 1
y=matrix(rnorm(n*P), ncol=P) #  Y is a n by P matrix in the tissue. 
x=rnorm(n) # X level in the tissue. 


ft1 <- iProMix.eFDR.PermReplace(yMatrix = y, x = x, pi = pi, cov=NULL, 
                                cl=T, reduce1 = c(2, 1), B=5, FDR=0.1)
# Obtain eFDR 
ft1$gene_eFDR

# Obtain estimated correlation
r1=sapply(1:P, function(f) ft1$ft_data[[f]]$cor.score1[2,1])
```

### Example 3: Identify if a pathway captured by Y is associated with one gene X in a cell type.

Identifying pathways is simpler than identifying genes in iProMix, where
permutation is not required.

``` r
library(iProMix)
set.seed(123)
n=500
P=10
pi=runif(n) # proportion of cell type 1
y=matrix(rnorm(n*P), ncol=P) #  Y is a n by P matrix in the tissue. 
x=rnorm(n) # X level in the tissue. 
path1=seq(1:5)

ft2 <- iProMix.LRT.matrix(yMatrix = y, x = x, pi = pi, cov=NULL, reduce1=c(1,2), cl=T)

# Obtain estimated correlation
r1=sapply(1:P, function(f) ft2[[f]]$cor.score1[2,1])

# Obtain pathway enrichment analysis results
sign_lr=sapply(1:P, function(f) sign(ft2[[f]]$cor.score1[2,1])*ft2[[f]]$LRT) # signed likelihood ratio statistics
Q_in= rank(sign_lr)[path1] # rank of the sign_lr within a pathway 
Q_out= rank(sign_lr)[-path1] # rank of the sign_lr outside of a pathway 
wilcox.test(Q_in, Q_out) # test p-value 
mean(Q_in) - mean(Q_out) # direction (if >0, positively enriched in a pathway)
```
