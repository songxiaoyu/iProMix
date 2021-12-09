#' iProMix.PermAdd.LRT
#'
#' @param y The quantitative measure (e.g. protein/expression) of a gene
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param x_tilde Permutated x
#' @param cov he covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param B The number of permutation
#' @param CellType The cell type to be estimated. Default is 1.
#' @param seed Seed for the permutation
#' @param tuningPar Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates
#' @param diffNum diffNum
#' @param numitersNum Number of Iteration
#'
#' @return A list
#'
iProMix.PermAdd.LRT=function(y, x, x_tilde=NULL, cov=NULL, pi,  B=1, CellType=1,
                             seed=NULL, tuningPar=1e-8, diffNum=0.0001,numitersNum=200) {

  set.seed(seed)
  y=as.matrix(y);  x=as.matrix(x);  pi=as.matrix(pi);
  n=nrow(x)
  if (is.null(x_tilde)) { x_tilde=lapply(1:B, function(f) sample(x))};

  M=lls=NULL
  for (b in 1:B) {

    xx=cbind(x, x_tilde[[b]])

    ft.zero=iProMix(y=y, x=xx, cov=cov, pi=pi,  reduce1=rbind(c(2,1), c(3,1)), reduce2=rbind(c(2,1), c(3,1)))

    # full model
    ft.full=iProMix(y=y, x=xx, cov=cov, pi=pi,  reduce1=NULL, reduce2=NULL,
                    inital.mu1=ft.zero$mu1, inital.mu2=ft.zero$mu2, 
                    inital.var1=ft.zero$var1, inital.var2=ft.zero$var2,
                    tuningPar=tuningPar, diffNum=diffNum,numitersNum=numitersNum)
    if (CellType==1) {
      # reduce x,y
      ft.dat=iProMix(y=y, x=xx,cov=cov,  pi=pi, reduce1=c(2,1), reduce2=NULL,
                     inital.mu1=ft.zero$mu1, inital.mu2=ft.zero$mu2, 
                     inital.var1=ft.zero$var1, inital.var2=ft.zero$var2,
                     tuningPar=tuningPar, diffNum=diffNum,numitersNum=numitersNum)
      # reduce x_tilde, y
      ft.noise=iProMix(y=y, x=xx,cov=cov,  pi=pi, reduce1=c(3,1), reduce2=NULL,
                       inital.mu1=ft.zero$mu1, inital.mu2=ft.zero$mu2, inital.var1=ft.zero$var1, inital.var2=ft.zero$var2)
    }
    if (CellType==2) {
      # reduce x,y
      ft.dat=iProMix(y=y, x=xx,cov=cov,  pi=pi, reduce1=NULL, reduce2=c(2,1),
                     inital.mu1=ft.zero$mu1, inital.mu2=ft.zero$mu2, 
                     inital.var1=ft.zero$var1, inital.var2=ft.zero$var2,
                     tuningPar=tuningPar, diffNum=diffNum,numitersNum=numitersNum)
      # reduce x_tilde, y
      ft.noise=iProMix(y=y, x=xx,cov=cov,  pi=pi, reduce1=NULL, reduce2=c(3,1),
                       inital.mu1=ft.zero$mu1, inital.mu2=ft.zero$mu2, 
                       inital.var1=ft.zero$var1, inital.var2=ft.zero$var2,
                       tuningPar=tuningPar, diffNum=diffNum,numitersNum=numitersNum)

    }

    LRT.dat=-2*(ft.dat$ll-ft.full$ll)
    LRT.noise=-2*(ft.noise$ll-ft.full$ll)

    lls=rbind(lls, c(full=ft.full$ll, data=ft.dat$ll, noise=ft.noise$ll))
    # diff in LRT
    M=c(M, LRT.dat-LRT.noise)
  }

  return(list(lls=lls, M=M))

}
