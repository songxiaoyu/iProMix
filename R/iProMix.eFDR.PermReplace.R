#' FDR controlled iProMix identification with permutation-replace procedure
#'
#' @param yMatrix The quantitative measure (e.g. protein/expression) of a gene in the matrix form (Row: Genes; Columns: Samples)
#' @param x The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y
#' @param cov The covariates for adjustment. Their impact on the mean value of X and Y are adjusted
#' @param pi The proportion of cell type 1
#' @param reduce1 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 1.  Default is NULL.
#' @param reduce2 A index of the row and column of the variance-covariance matrix that should be forced to be zero in cell type 2.  Default is NULL.
#' @param cl True of False. If true, parallel computing is used; need the library(doRNG). Default is FALSE
#' @param B The number of permutation
#' @param seed Seed for the permutation
#' @param FDR FDR cutoff; Default 0.1
#' @param verbose logic variable. If TRUE then will print additional information in permutation of iProMix
#'
#' @return list with 4 elements. It contains
#' \item{NoSigGene:}{Number of significant genes at the pre-determined FDR cutoff}
#' \item{IdxSigGene:}{The index of significant genes}
#' \item{ft_data:}{The iProMix output (e.g. estimated parameters, log likelihoods) using original data}
#' \item{ft_permutated:}{The iProMix output (e.g. estimated parameters, log likelihoods) using permuted data}
#' @export iProMix.eFDR.PermReplace
#'
#' @examples
#' \donttest{library(iProMix)
#' set.seed(111)
#' y <- matrix(rnorm(100,10,1), ncol=10)
#' x <- rnorm(10,10,1)
#' pi <- runif(10)
#' result <- iProMix.eFDR.PermReplace(yMatrix=y, x=x, pi=pi,reduce1=c(2,1), B=1, seed=1132, FDR=0.1)
#' }
iProMix.eFDR.PermReplace=function(yMatrix, x, cov=NULL, pi, x_tilde=NULL,
                                  reduce1=c(2,1), reduce2=NULL, cl=FALSE, 
                                  B=1, seed=NULL, FDR=0.1, verbose = FALSE) {
  if (is.null(seed)==FALSE) {set.seed(seed)}
  if (is.null(x_tilde)) {  x_tilde=lapply(1:B, function(f) sample(x))}

  #print(x_tilde[[1]][1:3])
  # calculate likelihood of the original data
  ft1=iProMix.LRT.matrix(yMatrix=yMatrix, x=x, cov=cov, pi=pi,
                         reduce1=reduce1, reduce2=reduce2, cl=cl)
  if (verbose == TRUE){
  print("Finished likelihood")
  }
  # calculate likelihood of the permuted data
  result_ft2 <- vector("list", B)
  for (b in 1:B){
    print(b)
    ft2=iProMix.LRT.matrix(yMatrix=yMatrix, x=x_tilde[[b]], cov=cov, pi=pi,
                           reduce1=reduce1, reduce2=reduce2, cl=cl)
    result_ft2[[b]] <- ft2
    if (verbose == TRUE){
    print(paste0("Finished likelihood permutation ", b))
    }
  }

  #  eFDR
  LRT_d <- sapply(1:ncol(yMatrix), function(f) ft1[[f]]$LRT)
  LRT_p <- sapply(1:B, function(f2) sapply(1:ncol(yMatrix), function(f1) result_ft2[[f2]][[f1]]$LRT))
  cut=seq(1, max(LRT_d), by =0.1)
  idx=which(sapply(cut, function(f) (mean(LRT_p>f)+1)/mean(LRT_d>f))<FDR)[1]
  NoSigGene=ifelse(is.na(idx), 0, sum(LRT_d>cut[idx]))
  
  # gene-specific eFDR
  gene_eFDR=sapply(1:ncol(yMatrix), function(f) mean(LRT_p>LRT_d[f])/(mean(LRT_d>LRT_d[f])+1e-12) )
  summary=list(NoSigGene=NoSigGene, IdxSigGene=  which(LRT_d>cut[idx]), gene_eFDR=gene_eFDR,
               ft_data = ft1, ft_permuted = result_ft2)
  return(summary)
}
