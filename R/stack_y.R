#' stack_y
#' 
#' Build the vector of responses
#'  
#' @keywords internal
#' @noRd

stack_y <- function(A,burnin=NULL,P,N){ 
  
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
  # a vector of responses.
  
  if(is.null(burnin)){
    burnin = 1
  }
  a<-data.frame(A)
  a<-a[(burnin+P):N,]
  a_stacked<-stack(a)
  return(a_stacked[,1])
}
