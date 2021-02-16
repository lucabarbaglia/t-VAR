#' BICfunction_glasso
#' 
#' Function for the choice of the regularization parameters 
#' for the inverse error covariance matrix based on BIC
#' 
#' @param sbic var-cov matrix
#' @param gamma1 regularization parameter for the inverse error covariance matrix
#' @param ndata  number of data points
#' @return BIC_values_glasso BIC values
#'  
#' @keywords internal
#' @noRd

BICfunction_glasso <- function(gamma1, sbic, ndata){
  
  #########
  # INPUT #
  #########
  # sbic : var-cov matrix 
  # gamma1 : regularization parameter for the inverse error covariance matrix 
  # ndata : number od data points
  
  ##########
  # OUTPUT #
  ##########
  # BIC_values_glasso: BIC values
  
  # Estimate the autoregressive coefficients using SPG algorithm
  FIT_glasso<-glasso::glasso(s = sbic,  rho=gamma1)
  FIT_glasso$wi
  
  # Compute the BIC
  NLog_lik_glasso<-(- FIT_glasso$loglik) # negative log-lik
  df_glasso<-length(which(FIT_glasso$wi!=0))
  BIC_value_glasso<-(2*NLog_lik_glasso+log(ndata)*df_glasso)
  
  return(BIC_value_glasso)
  
}
