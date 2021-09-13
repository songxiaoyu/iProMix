#' M step of the EM algorithm: Update the parameter estimates given the expectation of Q functions
#'
#' @param Eu1 The expected value in cell type 1 given the observed (y, x) at the tissue level
#' @param Eu2 The expected value in cell type 2 given the observed (y, x) at the tissue level
#' @param Evar1 The expected variance-covariance matrix in cell type 1 given the observed (y, x) at the tissue level
#' @param Evar2 The expected variance-covariance matrix in cell type 2 given the observed (y, x) at the tissue level
#' @param cov The covariates under consideration for its effects on the mean function of (y, x)
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates.
#'
#' @return A list of the 5 elements
#' \item{var1:}{The estimated (y,x) variance of cell type 1}
#' \item{var2:}{The estimated (y,x) variance of cell type 1}
#' \item{mu1:}{The estimated (y,x) mean function of cell type 1}
#' \item{mu2:}{The estimated (y,x) mean function of cell type 2}
#' \item{coef1:}{The estimated covariate effects on (y,x) in cell type 1}
#' \item{coef2:}{The estimated covariate effects on (y,x) in cell type 2}
MStep=function(Eu1, Eu2, Evar1, Evar2, cov, reduce1, reduce2, tuningPar=1e-8) {
  n=nrow(Eu1)
  if (is.null(cov)) {design=matrix(1, nrow=n, ncol=1)} else {design=cbind(1, cov)}
  # print(design[1:3,])
  coef1=stats::coef(stats::lm(Eu1~design-1))
  mu1= design%*% coef1
  coef2=stats::coef(stats::lm(Eu2~design-1))
  mu2= design%*% coef2

  s1=t(Eu1-mu1)%*%(Eu1-mu1)/n+Evar1 # var(u)= E(var(u|v)) + var(E(u|v))
  s2=t(Eu2-mu2)%*%(Eu2-mu2)/n+Evar2
  w1=glasso::glasso(s1, rho=tuningPar, zero=reduce1, penalize.diagonal = F)$w
  w2=glasso::glasso(s2, rho=tuningPar, zero=reduce2, penalize.diagonal = F)$w
  return(list(mu1=mu1, mu2=mu2, var1=w1, var2=w2, coef1=coef1, coef2=coef2))
}
