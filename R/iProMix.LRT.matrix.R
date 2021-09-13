#' For Y matrix with multiple genes/proteins/variables, calculate the likelihood ratio test statistics (as well as chisq p-value) for iProMix
#'
#' @param yMatrix The quantitative measure (e.g. protein/expression) of a gene in the matrix form (Row: Genes; Columns: Samples)
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates.
#' @param cl True of False. If true, parallel computing is used; need the library(doRNG). Default is FALSE
#'
#' @return list with number of elements equal to the number of genes for the yMatrix. Each gene contains a sublist with 10 elements. It contains
#' \item{var1:}{The estimated (y,x) variance of cell type 1}
#' \item{var2:}{The estimated (y,x) variance of cell type 2}
#' \item{mu1:}{The estimated (y,x) mean function of cell type 1}
#' \item{mu2:}{The estimated (y,x) mean function of cell type 2}
#' \item{cor.score1:}{The estimated X-Y correlation in cell type 1}
#' \item{cor.score2:}{The estimated X-Y correlation in cell type 2}
#' \item{full.ll:}{The estimated log likelihood function from the full model}
#' \item{reduced.ll:}{The estimated log likelihood function from the reduced model}
#' \item{LRT:}{The likelihood ratio test statistics for iProMix}
#' \item{LRT.pvalue:}{The chisq p-value for iProMix}
#' @import foreach
#' @import importFrom(doRNG,"%dorng%")
#' @export iProMix.LRT.matrix
#'
#' @examples
#' \donttest{library(iProMix)
#' set.seed(111)
#' y <- matrix(rnorm(100,10,1), ncol=10)
#' x <- rnorm(10,10,1)
#' pi <- runif(10)
#' iProMix_result <- iProMix.LRT.matrix(yMatrix = y, x = x, pi = pi, reduce1=c(2,1), reduce2=NULL)
#' }
iProMix.LRT.matrix=function(yMatrix, x, cov=NULL, pi, reduce1=c(2,1), reduce2=NULL, tuningPar=1e-6, cl=F) {

  yMatrix=as.matrix(yMatrix);
  x=as.matrix(x);
  #cov=as.matrix(cov)
  pi=as.matrix(pi)
  g=ncol(yMatrix)
  n=nrow(yMatrix)

  ft.all=list()
  if (cl==F) { # parallel computing
    for (p in 1:g) {
      print(p)
      ft= iProMix.LRT(y=yMatrix[,p], x=x, cov=cov, pi=pi, reduce1=reduce1, reduce2=reduce2, tuningPar=tuningPar)
      ft.all=append(ft.all, list(ft))
    }
    return(ft.all=ft.all)
  } else {
    ft.all=foreach::foreach (p = 1:g) %dopar% {
      ft= iProMix.LRT(y=yMatrix[,p], x=x, cov=cov, pi=pi, reduce1=reduce1, reduce2=reduce2, tuningPar=tuningPar)
      ft
    }
  }

}
