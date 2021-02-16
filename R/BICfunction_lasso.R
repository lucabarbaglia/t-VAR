#' BICfunction_lasso
#' 
#' Function for the choice of the regularization parameters 
#' for the autoregressive coeffcients based on BIC
#' 
#' @param Ybic vector of dimension NJ x 1 of dependent variables 
#' @param Xbic matrix of dimension NJ x J x P of indipendent variables 
#' @param lambda1 regularization parameter for the autoregressive coeffcients
#' @return BIC_values_glasso BIC values
#'  
#' @keywords internal
#' @noRd

BICfunction_lasso <- function(lambda1, Ybic, Xbic){
  
  #########
  # INPUT #
  #########
  # Ybic : vector of dimension NJ x 1 of dependent variables 
  # Xbic : matrix of dimension NJ x J?P of indipendent variables 
  # lambda1 : regularization parameter for the autoregressive coeffcients
  
  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values
  
  # Estimate the autoregressive coefficients using SPG algorithm
  FIT_lasso <- glmnet::glmnet(x = Xbic,  y = Ybic, family = "gaussian", lambda = lambda1)
  # Compute the BIC
  nj<-length(Ybic)
  bbb<-as.matrix(FIT_lasso$beta)
  Log_lik_lasso<- (1/nj) *t((Ybic - Xbic %*% bbb))%*%(Ybic - Xbic %*% bbb)
  df_lasso<-length(which(bbb != 0))
  BIC_value_lasso<-( - 2*Log_lik_lasso + log(nj)*df_lasso)
  
  return(BIC_value_lasso)
  
}
