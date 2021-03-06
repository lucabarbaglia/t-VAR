% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Large.tVAR.R
\name{Large.tVAR}
\alias{Large.tVAR}
\title{Large.tVAR}
\usage{
Large.tVAR(
  Data_ECM,
  P,
  type_lasso = "Lasso",
  nu_init = NULL,
  maxit.ECM = 25,
  tol.ECM = 0.01,
  tol.nu = 0.001,
  maxit.both = 50,
  tol.both = 0.01,
  maxit.nu = 100,
  lambda1_min = NULL,
  lambda1_max = NULL,
  lambda1_steps = NULL,
  gamma1_min = NULL,
  gamma1_max = NULL,
  gamma1_steps = NULL,
  lambda1_OPT = NULL,
  gamma1_OPT = NULL
)
}
\arguments{
\item{Data_ECM}{a NxJ matrix of log-volatilities. J: number of time series. N: time series length.}

\item{P}{VAR order}

\item{type_lasso}{type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".}

\item{nu_init}{degrees-of-freedom initial value. Default is 1000.}

\item{maxit.ECM}{maximum iterations for ECM algorithm. Default is 25.}

\item{tol.ECM}{tolerance ECM algorithm. Default is 0.01.}

\item{maxit.both}{maximum iterations for gaussian lasso algorithm. Default is 50.}

\item{tol.both}{tolerance gaussian lasso algorithm. Default is 0.01.}

\item{maxit.nu}{maximum iteration for estimation of the degrees-of-freedom. Default is 1000.}

\item{lambda1_min}{minimum value of the regularization parameter on Beta. Default is NULL.}

\item{lambda1_max}{maximum value of the regularization parameter on Beta.  Default is NULL.}

\item{lambda1_steps}{number of steps in the lambda grid.  Default is NULL.}

\item{gamma1_min}{minimum value of the regularization parameter on Omega.  Default is NULL.}

\item{gamma1_max}{maximum value of the regularization parameter on Omega.  Default is NULL.}

\item{gamma1_steps}{number of steps in the gamma grid.  Default is NULL.}

\item{lambda1_OPT}{optimal value of the regularization parameter on Beta.  Default is NULL.}

\item{gamma1_OPT}{optimal value of the regularization parameter on Omega.  Default is NULL.}
}
\value{
A list containing:
\item{"Beta_new"}{a vector containing the estimated Beta.}
\item{"Beta_arr"}{a JxJxP array containing the estimated Beta.}
\item{"innov"}{a (N-P)xJ containing the estimated VAR residuals.}
\item{"Omega_new"}{a JxJ matrix containing the estimated Omega.}
\item{"tau_new"}{a vector of length N-P containing the estimated gamma variable tau.}
\item{"nu_new"}{estimated degrees-of-freedom of the multivaraite t-distribution of the VAR innovations.} 
\item{"Obj_ECM"}{objective function.}
\item{"iter_ECM"}{number of iterations of the ECM algorithm.}
\item{"iter_vec"}{number of iteration of the Gaussian Lasso algorithm for each ECM iteration.}
\item{"lambda1_opt"}{selected value of the regularization parameter on Beta.}
\item{"gamma1_opt"}{selected value of the regularization parameter on Omega.}
}
\description{
ECM algorithm for the Large VAR with errors following a
multivariate t-distribution with estimation of the degrees of freedom
}
\references{
Barbaglia, L., Croux, C., & Wilms, I. (2020). Volatility spillovers in commodity markets: A large t-vector autoregressive approach. Energy Economics, 85, 104555.
}
\author{
Luca Barbaglia \url{https://lucabarbaglia.github.io/}
}
