#' stack_Xbig
#' 
#' Build the matrix of regressors
#'  
#' @keywords internal
#' @noRd

stack_Xbig <- function(A,burnin=NULL,P,N){ 
  
  #########
  # INPUT #
  #########
  # A: data.
  # burnin: burnin to remove at the beginnnig. Deafault NULL.
  # P : lag order.
  # N : time series length.
  
  ##########
  # OUTPUT #
  ##########
  # a matrix X.
  
  if(is.null(burnin)){
    burnin = 1
  }
  J<-ncol(A)
  X0_arr1<-array(NA, c((N-P-burnin+1), J, P)) # X0 for class k in array form
  for (ip in P:1){ 
    A1<-matrix(A[(burnin+ip-1):(N-P+ip-1),], nrow=(N-P-burnin+1) , ncol=J )
    X0_arr1[1:(N-P-burnin+1),(1:J),ip]<-A1 # [,,1]: elements at lag p=1, J columns and N rows going from P-p to N-p. [,,P] elements at lag p=P
  }
  X0_list<-plyr::alply(X0_arr1, 3,.dims=TRUE) # Convert the array into a list
  X0<-do.call(cbind,X0_list) # Matrix; the first J columns are at lag p=1, the last J columns are at lag p=P. The rows go from P-p to N-p.
  Xbig_mat<-kronecker(diag(J), X0)
  return(Xbig_mat)
}
