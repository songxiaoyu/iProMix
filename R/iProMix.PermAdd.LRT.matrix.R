#' For Y matrix, calculate the difference in likelihood ratio test statistics for x and x_tilde
#'
#' @param yMatrix The quantitative measure (e.g. protein/expression) of a gene in the matrix form (Row: Genes; Columns: Samples)
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param x_tilde The permuatated x
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param CellType The cell type to be estimated. Default is 1.
#' @param B The number of permutation
#' @param seed Seed for the permutation
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates.
#' @param cl True of False. If true, parallel computing is used; need the library(doRNG). Default is FALSE
#'
#' @return list with number of elements equal to the number of genes for the yMatrix. Each gene contains a sublist with 2 elements. It contains
#' \item{ll:}{The estimated log likelihood function}
#' \item{M:}{The M value}
#' @export iProMix.PermAdd.LRT.matrix
#' @import doRNG
#' @import foreach
#' @import importFrom(doRNG,"%dorng%")
#' @examples
#' \donttest{library(iProMix)
#' set.seed(111)
#' y <- matrix(rnorm(1000,10,1), ncol=10)
#' x <- rnorm(100,10,1)
#' pi <- runif(100)
#' result <- iProMix.PermAdd.LRT.matrix(yMatrix=y, x=x, pi=pi, B=1, CellType=1)
#' }
iProMix.PermAdd.LRT.matrix=function(yMatrix, x,  x_tilde=NULL, cov=NULL, pi, CellType=1, B=1, seed=NULL, tuningPar=1e-6, cl=F) {

  yMatrix=as.matrix(yMatrix);# -- coef1 and coef2 are not on original scale
  x=as.matrix(x);
  #cov=as.matrix(cov)
  pi=as.matrix(pi)
  g=ncol(yMatrix)
  n=nrow(yMatrix)

  ft.all=list()
  if (cl==F) {
    for (p in 1:g) {
      print(p)
      ft= iProMix.PermAdd.LRT(y=yMatrix[,p], x=x, x_tilde=x_tilde, cov=cov, pi=pi, CellType=CellType, seed=seed, B=B,tuningPar=tuningPar )
      ft.all=append(ft.all, list(ft))
    }
    return(ft.all=ft.all)
  } else {
    ft.all=foreach::foreach (p = 1:g) %dorng% {
      ft= iProMix.PermAdd.LRT(y=yMatrix[,p], x=x, x_tilde=x_tilde, cov=cov, pi=pi, CellType=CellType, seed=seed, B=B,tuningPar=tuningPar )
      ft
    }
  }

}
