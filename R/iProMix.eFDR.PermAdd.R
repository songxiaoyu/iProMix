#' FDR controlled iProMix identification with permutation-add procedure
#'
#' @param yMatrix The quantitative measure (e.g. protein/expression) of a gene in the matrix form (Row: Genes; Columns: Samples)
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param CellType The cell type to be estimated. Default is 1.
#' @param cl True of False. If true, parallel computing is used; need the library(doRNG). Default is FALSE
#' @param B The number of permutation
#' @param seed Seed for the permutation
#' @param FDR FDR cutoff; Default 0.1
#'
#' @return list with 10 elements. It contains
#' \item{NoSigGene:}{Number of significant genes at the pre-determined FDR cutoff, using consensus voting procedure if B>1}
#' \item{IdxSigGene:}{The index of significant genes at the pre-determined FDR cutoff, using consensus voting procedure if B>1}
#' \item{IdxSigGeneEach:}{The index of significant genes at each permutation}
#' \item{ft_add:}{The iProMix output (e.g. estimated parameters, log likelihoods) using combined original and permutation data (y, x, x_tilde).}
#' @export iProMix.eFDR.PermAdd
#'
#' @examples
#' \donttest{library(iProMix)
#' set.seed(111)
#' y <- matrix(rnorm(1000,10,1), ncol=10)
#' x <- rnorm(100,10,1)
#' pi <- runif(100)
#' result <- iProMix.eFDR.PermAdd(yMatrix=y, x=x, pi=pi,reduce1=c(2,1), B=1, seed=1132, FDR=0.1)
#' }
iProMix.eFDR.PermAdd=function(yMatrix, x, cov=NULL, pi,
                              reduce1=c(2,1), reduce2=NULL, CellType=1, cl=F, B=1, seed=NULL, FDR=0.1) {
  if (is.null(seed)==F) {set.seed(seed)}
  x_tilde=lapply(1:B, function(f) sample(x))

  # add model likelihood
  ft3= iProMix.PermAdd.LRT.matrix(yMatrix=yMatrix, x=x, x_tilde=x_tilde, cov=cov,CellType=CellType,
                                  pi=pi,  B=B, seed=NULL, tuningPar=1e-6, cl=cl)
  print("Permutation finished")
  M=sapply(1:length(ft3), function(f) ft3[[f]]$M)
  if (B==1) {
    cut2=knockoff::knockoff.threshold(W=M, fdr = FDR, offset = 1);
    M_vote=M_identified <- which(M>cut2)
  }
  if (B>1) {
    cut2=sapply(1:B, function(f) knockoff::knockoff.threshold(W=M[,f], fdr = FDR, offset = 1));
    M_identified <- lapply(1:B, function(f) which(M[,f]>cut2[f]))
    vote=table(unlist(M_identified))
    M_vote=names(vote)[(vote>B/2)]
  }

  summary=list(NoSigGene=length(M_vote), IdxSigGene=M_vote, IdxSigGeneEach=M_identified,  ft_add = ft3)
  return(summary)
}
