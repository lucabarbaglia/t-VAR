#' VARbyGLS 
#' 
#' Function for the classic unpenalized VAR esitmation, Generalized Least Squares
#' 
#' @param Y mulitvariate time series.
#' @param p lag length.
#' 
#' @return A list containing:
#' \item{"Beta"}{a list of betas.}
#' \item{"Sigma"}{var-cov matrix.}
#' \item{"Beta.vec"}{vector of betas [b^1_11,b^2_11,..,b^P_11,b^1_21,b^2_21,..,b^P_21,.., b^P_JJ].}
#' 
#' @export

VARbyGLS <- function(Y,p){
  
  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  Response <- c(Y[-(1:p),])
  Xmatrix <- VARdesign(Y,p)
  Xmatrix.kron <- cbind(1,kronecker(diag(1,q),Xmatrix))
  
  # Obtain OLS fit residuals
  OLSresid <- matrix(resid(lm(Response~Xmatrix.kron)),ncol=q)
  Sigma <- cov(OLSresid)
  # Do GLS
  Betahat <- solve(t(Xmatrix.kron)%*%kronecker(solve(Sigma),diag(1,n-p))%*%Xmatrix.kron)%*%t(Xmatrix.kron)%*%kronecker(solve(Sigma),diag(1,n-p))%*%Response
  Betahat <- Betahat[-1]
  Betahatt <- list()
  for(i in 1:p){
    Betahatt[[i]] <- matrix(Betahat[seq(i,(p*q*q),by=p)],ncol=q,byrow=T)
  }
  out <- list(Beta=Betahatt,Sigma=Sigma, Beta.vec=Betahat)
}
