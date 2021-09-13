#' Calculate the log likelihood function given parameters
#'
#' @param u The combined data of (y, x)
#' @param pi The proportion of cell type 1 in all cells
#' @param mu1 The estimated (y,x) mean function of cell type 1
#' @param mu2 The estimated (y,x) mean function of cell type 2
#' @param var1 The estimated (y,x) variance of cell type 1
#' @param var2 The estimated (y,x) variance of cell type 2
#'
#' @return log likelihood function
#'
LL<-function(u, pi, mu1, mu2, var1, var2)
{
  n=nrow(u); p=ncol(u)

  sum(sapply(1:n,function(i)
  {
    mu_u= pi[i,] %*% mu1[i,] + (1-pi[i,]) %*% mu2[i,]
    var_u=pi[i,]^2 * var1  + (1-pi[i,])^2 *var2
    emdbook::dmvnorm(u[i,], mu_u, var_u, log = TRUE) # log =T returns log likelihood function
  }))
}
