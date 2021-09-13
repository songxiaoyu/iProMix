#' E step of the EM algorithm: Calculate the expectation of Q function given parameters
#'
#' @param u The combined data of (y, x)
#' @param pi The proportion of cell type 1 in all cells
#' @param mu1 The estimated (y,x) mean function of cell type 1
#' @param mu2 The estimated (y,x) mean function of cell type 2
#' @param var1 The estimated (y,x) variance of cell type 1
#' @param var2 The estimated (y,x) variance of cell type 2
#'
#' @return A list of the 4 elements
#' \item{Eu1:}{The expected value in cell type 1 given the observed (y, x) at the tissue level given input parameters}
#' \item{Eu2:}{The expected value in cell type 2 given the observed (y, x) at the tissue level given input parameters}
#' \item{Evar1:}{The expected variance-covariance matrix in cell type 1 given the observed (y, x) at the tissue level given input parameters}
#' \item{Evar2:}{The expected variance-covariance matrix in cell type 2 given the observed (y, x) at the tissue level given input parameters}
EStep=function(u, pi, mu1, mu2, var1, var2){
  n=nrow(u); p=ncol(u);

  Eu1=matrix(NA,n, p);   Eu2=matrix(NA,n, p)
  Evar1=matrix(0,p,p);   Evar2=matrix(0,p,p);
  for (i in 1:n) {
    # print(i)
    mu_u= pi[i,] %*% mu1[i,] + (1-pi[i,]) %*% mu2[i,]
    var_u=pi[i,]^2 * var1  + (1-pi[i,])^2 *var2
    inv_var_u=MASS::ginv(var_u)
    # mean
    cov_u1u=pi[i,]*var1
    mu_u1_u=mu1[i,] + t( t(cov_u1u) %*% inv_var_u %*%  matrix(u[i,] - mu_u))
    cov_u2u=(1-pi[i,])*var2
    # mu_u2_u=mu2[i,] + t( t(cov_u2u) %*% inv_var_u%*%  matrix(u[i,] - mu_u))
    mu_u2_u=(u[i,]-pi[i,]*mu_u1_u)/(1-pi[i,])

    # var
    var_u1_u=var1- t(cov_u1u) %*% inv_var_u %*%  cov_u1u
    var_u2_u=var2- t(cov_u2u) %*% inv_var_u %*%  cov_u2u


    Eu1[i,]= mu_u1_u;
    Eu2[i,]= mu_u2_u;
    Evar1=Evar1+ var_u1_u;
    Evar2=Evar2+ var_u2_u;
    # print(var_u2_u)
  }
  Evar1=Evar1/n; Evar1<-(Evar1+t(Evar1))/2
  Evar2=Evar2/n; Evar2<-(Evar2+t(Evar2))/2

  # get mean
  return(list(Eu1=Eu1,Eu2=Eu2, Evar1=Evar1,  Evar2=Evar2))
}
