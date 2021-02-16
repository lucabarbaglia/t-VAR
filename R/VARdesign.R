#' VARdesign
#' 
#' VARdesign is a function that construct the design matrix 
#' of the VAR model for a multivariate time series Y
#' 
#' @param Y mulitvariate time series.
#' @param p lag length.
#' @return a design matrix.
#' 
#' @keywords internal
#' @noRd


VARdesign <- function(Y,p){
  

  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  designmatrix <- matrix(rep(NA,n-p),ncol=1)
  for(iq in 1:q){
    for(ip in 1:p){
      designmatrix <- cbind(designmatrix,Y[((p:1)[ip]):(n-ip),iq])
    }
  }
  designmatrix <- designmatrix[,-1]
  out <- designmatrix
}
