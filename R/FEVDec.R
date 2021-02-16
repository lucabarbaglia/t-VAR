#' FEVDec
#' 
#' Obtain Forecast Error variance Decomposition 
#' with Spectral, Cholesky or Generealized decomposition
#' 
#' @param B_arr arrays of autoregressive coeff 
#' @param Sig var-cov matrix
#' @param lag forecast horizon
#' @param dec.type  "Cholesky", "Spectral", "Generalized", "GeneralizedIRF" or "Generalized_Lanne". Default "Generalized".
#' 
#' @return fevd: forecast error variance decompositoin matrix
#' 
#' @references Barbaglia, L., Croux, C., & Wilms, I. (2020). Volatility spillovers in commodity markets: A large t-vector autoregressive approach. Energy Economics, 85, 104555.
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export



FEVDec<-function (B_arr, Sig, lag, dec.type="Generalized") {

  # Function
  spectral_dec<-function(Omega){
    decomp<-eigen(Omega, symmetric = TRUE)  
    CC <- decomp$vectors
    Lambda <- diag(sqrt(decomp$values))
    P <- CC %*% Lambda %*% t(CC)
  } 
  
  # Dimensions
  J<-dim(B_arr)[1]
  P<-dim(B_arr)[3]
  nt<-lag+P # number of MA coefficients
  
  # MA coefficients
  T_arr<-array(0, c(J,J,nt)) # VMA coeff. Take 0 for ntheta>P
  T_arr[,,1]<-diag(J) # Initial MA is the Identity
  
  for (ima in (1 : (P))){
    T_arr[,,ima+1]<-T_arr[,,ima]%*%B_arr[,,ima]
  }
  
  if (dec.type=="Spectral"){
    # Orthogonal MA
    Sig_dec<-spectral_dec(Sig) # Spectral decomposition of Sigma
    T_arr_ort<-array(NA, c(J,J,nt))
    for (ima in (1:nt)){
      T_arr_ort[,,ima]<-T_arr[,,ima]%*%Sig_dec
    }
    # FEVD
    fevd<-matrix(NA, nrow=J, ncol=J)
    mse<-matrix(NA, nrow=J, ncol=1)
    for (irow in (1:J)){
      mse[irow,]<-sum((T_arr_ort[irow,,])^2)
      for (icol in 1:J){
        fevd[irow,icol]<-sum((T_arr_ort[irow,icol,])^2)/mse[irow,]
      }
    }
  } # end of Spectral
  
  if (dec.type=="Cholesky"){
    # Orthogonal MA
    Sig_dec<-chol(Sig) # Spectral decomposition of Sigma
    T_arr_ort<-array(NA, c(J,J,nt))
    for (ima in (1:nt)){
      T_arr_ort[,,ima]<-T_arr[,,ima]%*%Sig_dec
    }
    # FEVD
    fevd<-matrix(NA, nrow=J, ncol=J)
    mse<-matrix(NA, nrow=J, ncol=1)
    for (irow in (1:J)){
      mse[irow,]<-sum((T_arr_ort[irow,,])^2)
      for (icol in 1:J){
        fevd[irow,icol]<-sum((T_arr_ort[irow,icol,])^2)/mse[irow,]
      }
    }
  } # end of Cholesky
  
  if (dec.type=="Generalized"){
    # FEVD
    fevd1<-matrix(NA, nrow=J, ncol=J)
    for (irow in (1:J)){
      e_irow<-matrix(0, nrow=J, ncol=1)
      e_irow[irow,]<-1
      mse_int<-array(NA,c(1,1,nt))
      for (int in 1:nt){
        mse_int[,,int]<-t(e_irow)%*%T_arr[,,int]%*%Sig%*%t(T_arr[,,int])%*%e_irow 
      }
      mse<-sum(mse_int)
      for (icol in 1:J){
        e_icol<-matrix(0, nrow=J, ncol=1)
        e_icol[icol,]<-1
        num_int<-array(NA,c(1,1,nt))
        for (int in 1:nt){
          num_int[,,int]<-(t(e_irow)%*%T_arr[,,int]%*%Sig%*%e_icol)^2
        }
        num<-sum(num_int)/Sig[icol,icol]
        fevd1[irow,icol]<-num/mse
      }
    }
    # Normalize to have the sum of each row =1 (Diebold Yilmaz, 2012)
    fevd<-matrix(NA, nrow=J, ncol=J)
    for (irow in 1:J){
      sum_row<-sum(fevd1[irow,])
      for (icol in 1:J){
        fevd[irow,icol]<-fevd1[irow, icol]/sum_row 
      }
    }
  } # end of Generalized
  
  
  if (dec.type=="GeneralizedIRF"){
    # Generealized IRF
    fevd<-matrix(NA, nrow=J, ncol=J)
    for (irow in (1:J)){
      e_irow<-matrix(0, nrow=J, ncol=1)
      e_irow[irow,]<-1
      for (icol in 1:J){
        e_icol<-matrix(0, nrow=J, ncol=1)
        e_icol[icol,]<-1
        num_int<-array(NA,c(1,1,nt))
        for (int in 1:nt){
          num_int[,,int]<-(t(e_irow)%*%T_arr[,,int]%*%Sig%*%e_icol)^2
        }
        fevd[irow,icol]<-sum(num_int)/Sig[icol,icol]
      }
    }
  } # end of GeneralizedIRF
  
  if (dec.type=="Generalized_Lanne"){
    # FEVD
    fevd1<-matrix(NA, nrow=J, ncol=J)
    fevd<-matrix(NA, nrow=J, ncol=J)
    for (irow in (1:J)){
      e_irow<-matrix(0, nrow=J, ncol=1)
      e_irow[irow,]<-1
      for (icol in 1:J){
        e_icol<-matrix(0, nrow=J, ncol=1)
        e_icol[icol,]<-1
        num_int<-array(NA,c(1,1,nt))
        for (int in 1:nt){
          num_int[,,int]<-(t(e_irow)%*%T_arr[,,int]%*%Sig%*%e_icol)^2
        }
        num<-sum(num_int)/Sig[icol,icol]
        fevd1[irow,icol]<-num
      }
      for (icol in 1 :J){
        fevd[irow,icol]<-fevd1[irow,icol]/sum(fevd1[irow,])  
      }
    }
  } # end of Generalized_Lanne
  
  
  return(fevd)
}
