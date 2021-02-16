#' Spillovers
#' 
#' Function to compute volatility spillovers and 
#' its index based on forecast error variance decomposition
#' 
#' @param fit t-VAR fit obtained from the function Large.tVAR
#' @param dec.type decomposition for forecast error varaince decomposition: "Cholesky", "Spectral", "Generalized", "GeneralizedIRF" or "Generalized_Lanne". Default "Generalized".
#' @param h forecast horizon for variance decomposition. Default 1.
#' 
#' @return A list containing:
#' \item{"spill"}{matrix of volatility spillovers.}
#' \item{"spill_index"}{volatility spillover index.}
#' 
#' @export


Spillovers<-function(fit, dec.type="Generalized", h=1){
  
  fevd_mat<-FEVDec(B_arr=fit$Beta_arr, Sig=solve(fit$Omega_new), lag=h, dec.type = dec.type)
  
  spill<-100*fevd_mat
  
  # Remove the own-effects (main diagonal)
  spill_cross<-spill
  diag(spill_cross)<-0
  
  # Global Spillovers (exclude the own-effects)
  spill_gl<-sum(spill_cross)
  
  Spillovers<-list("spill"=spill, "spill_index"=spill_gl)
}
