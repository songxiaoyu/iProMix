#' A decomposition model
#'
#' @param y The quantitative measure (e.g. protein/expression) of a gene
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates
#' @param diffNum Default 0.0001. The convergency criterion for EM algorithm.
#' @param numitersNum Default 200. The number of iterations in EM algorithm.
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param inital.mu1 The initial mean value of (Y,X) for cell type 1, starting in Y and then in X. At default inital.mu1 =NULL, tissue-level mean is used.
#' @param inital.mu2 The initial mean value of (Y,X) for cell type 2, starting in Y and then in X. At default inital.mu2 =NULL, tissue-level mean is used.
#' @param inital.var1 The initial variance-covariance matrix of (Y,X) for cell type 1. At default inital.var1 =NULL, tissue-level variance-covariance matrix is used.
#' @param inital.var2 The initial variance-covariance matrix of (Y,X) for cell type 2. At default inital.var2 =NULL, tissue-level variance-covariance matrix is used.
#'
#' @return list with 9 elements. It contains
#' \item{var1:}{The estimated (y,x) variance of cell type 1}
#' \item{var2:}{The estimated (y,x) variance of cell type 2}
#' \item{mu1:}{The estimated (y,x) mean function of cell type 1}
#' \item{mu2:}{The estimated (y,x) mean function of cell type 2}
#' \item{cor.score1:}{The estimated X-Y correlation in cell type 1}
#' \item{cor.score2:}{The estimated X-Y correlation in cell type 2}
#' \item{ll:}{The estimated log likelihood function}
#' \item{coef1:}{The estimated covariate effects on (y,x) in cell type 1}
#' \item{coef2:}{The estimated covariate effects on (y,x) in cell type 2}
#' @export iProMix
#'
#' @examples
#' library(iProMix)
#'set.seed(111)
#' y <- rnorm(100,10,1)
#' x <- rnorm(100,10,1)
#' pi <- runif(100)
#' iProMix_result <- iProMix(y = y, x = x, pi = pi, reduce1=c(2,1), reduce2=NULL)
iProMix= function(y, x, cov=NULL, pi,
                  tuningPar=1e-8, diffNum=0.001,numitersNum=200,
                  reduce1=NULL, reduce2=NULL,
                  inital.mu1=NULL, inital.mu2=NULL,
                  inital.var1=NULL, inital.var2=NULL) {
  # reduce=c(2,1) # index of the row and column to reduce to zero
  # diffNum=0.001;numitersNum=200; inital.mu1=NULL; inital.mu2=NULL; inital.var1=NULL; inital.var2=NULL
  y=as.matrix(y);  x=as.matrix(x);  pi=as.matrix(pi)
  u=cbind(y,x)
  u=u[stats::complete.cases(u),]
  n=nrow(y);  p=ncol(u)

  # initial
  # Eu1=Eu2=u
  if (is.null(inital.mu1))  {mu1=matrix(apply(u, 2, mean), nrow=n, ncol=p, byrow=T)} else{mu1=inital.mu1}
  if (is.null(inital.mu2))  {mu2=matrix(apply(u, 2, mean), nrow=n, ncol=p, byrow=T)} else{mu2=inital.mu2}
  if (is.null(inital.var1))  {var1=diag(diag(stats::var(u)))} else{var1=inital.var1}
  if (is.null(inital.var2))  {var2=diag(diag(stats::var(u)))} else{var2=inital.var2}


  diff<-100
  numiters<-1
  ll=  LL(u=u, pi=pi, mu1=mu1, mu2=mu2, var1=var1, var2=var2)

  while (diff > diffNum & numiters < numitersNum) {
    # print(numiters)
    numiters=numiters+1
    # Estep - give parameters; calcualte x1, x2, y1, y2
    e_output=EStep(u=u, pi=pi, mu1=mu1, var1=var1, mu2=mu2, var2=var2)
    Eu1=e_output$Eu1;Eu2=e_output$Eu2

    # Mstep: Given X1, Y1, X2, Y2, get mu_1, mu_2, sigma_1, sigma_2, rho
    m_output=MStep(Eu1=Eu1, Eu2=Eu2,
                   Evar1=e_output$Evar1, Evar2=e_output$Evar2,
                   reduce1=reduce1, reduce2=reduce2,
                   cov=cov, tuningPar=tuningPar)
    mu1=m_output$mu1;mu2=m_output$mu2;  var1=m_output$var1; var2=m_output$var2
    #print(var)
    # calculate LL
    ll= c(ll, LL(u=u, pi=pi, mu1=mu1, mu2=mu2, var1=var1, var2=var2))

    diff<-ll[numiters]-ll[numiters-1]
    #ll
  }
  # estimate - likelihood
  r1=stats::cov2cor(var1)
  r2=stats::cov2cor(var2)
  r1;r2

  if (is.null(cov)==F) {
    ft1=summary(lm(Eu1~cov));
    ft2=summary(lm(Eu1~cov))
  } else{
    ft1=ft2=NULL
  }


  return(list(var1=var1, mu1=mu1,var2=var2, mu2=mu2,
              cor.score1=r1, cor.score2=r2, ll=ll[numiters],
              coef1=m_output$coef1, coef2=m_output$coef2, ft1=ft1, ft2=ft2))
}
