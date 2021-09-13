#' Calculate the likelihood ratio test statistics (as well as chisq p-value) for iProMix
#'
#' @param y The quantitative measure (e.g. protein/expression) of a gene
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates.
#' @param diffNum Default 0.0001. The convergency criterion for EM algorithm.
#' @param numitersNum Default 200. The number of iterations in EM algorithm.
#'
#' @return list with 10 elements. It contains
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
#' @export iProMix.LRT
#'
#' @examples
#' library(iProMix)
#' set.seed(111)
#' y <- rnorm(100,10,1)
#' x <- rnorm(100,10,1)
#' pi <- runif(100)
#' iProMix.LRT_result <- iProMix.LRT(y = y, x = x, pi = pi, reduce1=c(2,1), reduce2=NULL)
iProMix.LRT=function(y, x, cov=NULL, pi, reduce1=c(2,1), reduce2=NULL, tuningPar=1e-6, diffNum=0.0001,numitersNum=200 ) {


  # Under the null model: no dependency in either cell type
  ft.reduced1=iProMix(y=y, x=x,cov=cov, pi=pi, reduce1=c(2,1), reduce2=c(2,1), tuningPar=tuningPar)
  # under the specified model, using parameters from the null model
  ft.reduced=iProMix(y=y, x=x,cov=cov, pi=pi, reduce1=reduce1, reduce2=reduce2, tuningPar=tuningPar,
                     inital.mu1=ft.reduced1$mu1, inital.mu2=ft.reduced1$mu2, inital.var1=ft.reduced1$var1, inital.var2=ft.reduced1$var2)
  # Under the full model: under estimated dependence, using parameters from the null model
  ft.full=iProMix(y=y, x=x, cov=cov, pi=pi, reduce1=NULL, reduce2=NULL, tuningPar=tuningPar,
                  inital.mu1=ft.reduced1$mu1, inital.mu2=ft.reduced1$mu2, inital.var1=ft.reduced1$var1, inital.var2=ft.reduced1$var2)

  LRT=-2*(ft.reduced$ll-ft.full$ll)
  pvalue=stats::pchisq(LRT, df=1, lower.tail = F)
  pvalue
  return(list(var1=ft.full$var1, var2=ft.full$var2, mu1=ft.full$mu1, mu2=ft.full$mu2,
              cor.score1=ft.full$cor.score1, cor.score2=ft.full$cor.score2,
              full.ll=ft.full$ll, reduced.ll=ft.reduced$ll,LRT=LRT, LRT.pvalue=pvalue))

}
