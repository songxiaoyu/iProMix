% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iProMix.PermAdd.LRT.R
\name{iProMix.PermAdd.LRT}
\alias{iProMix.PermAdd.LRT}
\title{iProMix.PermAdd.LRT}
\usage{
iProMix.PermAdd.LRT(
  y,
  x,
  x_tilde = NULL,
  cov = NULL,
  pi,
  B = 1,
  CellType = 1,
  seed = NULL,
  tuningPar = 1e-08,
  diffNum = 1e-04,
  numitersNum = 200
)
}
\arguments{
\item{y}{The quantitative measure (e.g. protein/expression) of a gene}

\item{x}{The quantitative measure of anther gene  (e.g. ACE2 protein levels) that we would like to know their cell-type specific dependency with Y}

\item{x_tilde}{Permutated x}

\item{cov}{he covariates for adjustment. Their impact on the mean value of X and Y are adjusted}

\item{pi}{The proportion of cell type 1}

\item{B}{The number of permutation}

\item{CellType}{The cell type to be estimated. Default is 1.}

\item{seed}{Seed for the permutation}

\item{tuningPar}{Default is 1e-8. It is used in the embedded graphic lasso procedure for estimating correlation. A larger tuningPar can be selected if one is interested in penalized estimates}

\item{diffNum}{diffNum}

\item{numitersNum}{Number of Iteration}
}
\value{
A list
}
\description{
iProMix.PermAdd.LRT
}
