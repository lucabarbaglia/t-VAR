#' BICfunction_grplasso
#' 
#' Function for the choice of the regularization parameters 
#' for the autoregressive coeffcients based on BIC
#' 
#' @param Ybic vector of dimension NJ x 1 of dependent variables 
#' @param Xbic matrix of dimension NJ x J x P of indipendent variables 
#' @param lambda1 regularization parameter for the autoregressive coeffcients
#' @param indexbic group structure
#' @return BIC_values_glasso BIC values
#'  
#' @keywords internal
#' @noRd

BICfunction_grplasso <- function(lambda1, Ybic, Xbic, indexbic){
  
  #########
  # INPUT #
  #########
  # Ybic : vector of dimension NJ x 1 of dependent variables 
  # Xbic : matrix of dimension NJ x J?P of indipendent variables 
  # lambda1 : regularization parameter for the autoregressive coeffcients
  # indexbic : group structure
  
  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values
  
  # Estimate the autoregressive coefficients using SPG algorithm
  FIT_grplasso<-grplasso::grplasso(x = Xbic, y = Ybic, index=indexbic, lambda = lambda1, 
                         model = grplasso::LinReg(), center = F,  control = grplasso::grpl.control(trace=0))
  
  # Compute the BIC
  NLog_lik_grplasso<-FIT_grplasso$nloglik # negative log-lik
  df_grplasso<-length(which(FIT_grplasso$coefficients!=0))
  BIC_value_grplasso<-(2*NLog_lik_grplasso+log(length(Ybic))*df_grplasso)
  
  return(BIC_value_grplasso)
  
}
