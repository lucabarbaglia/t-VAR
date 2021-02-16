#' EM_VAR
#' 
#' EM algorithm for the Robust Sparse VAR with errors following a 
#' multivariate t-distribution (without estimation of df)
#' 
#' @param Data_EM a NxJ matrix of log-volatilities. J: number of time series. N: time series length
#' @param P VAR order
#' @param type_lasso type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
#' @param nu degrees-of-freedom of the multivaraite t-distribution of the VAR innovations.
#' @param maxit.EM maximum iterations for EM algorithm. Default is 25.
#' @param tol.EM tolerance EM algorithm. Default is 0.01.
#' @param maxit.both maximum iterations for gaussian lasso algorithm. Default is 50.
#' @param tol.both tolerance gaussian lasso algorithm. Default is 0.01.
#' @param lambda1_min minimum value of the regularization parameter on Beta. Default is NULL. 
#' @param lambda1_max maximum value of the regularization parameter on Beta.  Default is NULL.
#' @param lambda1_steps number of steps in the lambda grid.  Default is NULL.
#' @param gamma1_min minimum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_max maximum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_steps number of steps in the gamma grid.  Default is NULL.
#' @param gamma1_steps number of steps in the gamma grid.  Default is NULL.
#' @param lambda1_OPT optimal value of the regularization parameter on Beta.  Default is NULL.
#' @param gamma1_OPT optimal value of the regularization parameter on Omega.  Default is NULL.
#' 
#' @return A list containing:
#' \item{"Beta_new"}{a vector containing the estimated Beta.}
#' \item{"Beta_arr"}{a JxJxP array containing the estimated Beta.}
#' \item{"innov"}{a (N-P)xJ containing the estimated VAR residuals.}
#' \item{"Omega_new"}{a JxJ matrix containing the estimated Omega.}
#' \item{"tau_new"}{a vector of length N-P containing the estimated gamma variable tau.}
#' \item{"Obj_ECM"}{objective function.}
#' \item{"iter_ECM"}{number of iterations of the ECM algorithm.}
#' \item{"iter_vec"}{number of iteration of the Gaussian Lasso algorithm for each ECM iteration.}
#' \item{"lambda1_opt"}{selected value of the regularization parameter on Beta.}
#' \item{"gamma1_opt"}{selected value of the regularization parameter on Omega.}
#' 
#' @references Barbaglia, L., Croux, C., & Wilms, I. (2020). Volatility spillovers in commodity markets: A large t-vector autoregressive approach. Energy Economics, 85, 104555.
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export



EM_VAR<-function(Data_EM, P, type_lasso="Lasso", 
                 nu, maxit.EM=25, tol.EM=0.01,maxit.both=50, tol.both=0.01, 
                 lambda1_min=NULL, lambda1_max=NULL, lambda1_steps=NULL, 
                 gamma1_min=NULL, gamma1_max=NULL, gamma1_steps=NULL,
                 lambda1_OPT=NULL, gamma1_OPT=NULL){
  
  

  ########################
  # Additional functions #
  ########################
  
  # Opposite of null
  is.not.null <- function(x) ! is.null(x)
  
  # Opposite of is.matrix
  is.not.matrix <- function(x) ! is.matrix(x)
  
  
  #####################
  # Preliminary steps #
  #####################
  
  if(is.not.matrix(Data_EM)){stop("Data must be in matrix form.")}
  
  if (is.not.null(lambda1_OPT) & (is.not.null(lambda1_min) | is.not.null(lambda1_max) | is.not.null(lambda1_steps))){
    stop("Either set the optimal value of lambda, either the values for the grid search")
  }
  if (is.not.null(gamma1_OPT) & (is.not.null(gamma1_min) | is.not.null(gamma1_max) | is.not.null(gamma1_steps))){
    stop("Either set the optimal value of gamma, either the values for the grid search")
  }
  
  if (is.null(lambda1_OPT) & (is.null(lambda1_min) | is.null(lambda1_max) | is.null(lambda1_steps))){
    stop("Set all the values of the lambda grid or an optimal value for lambda")
  }
  if (is.null(gamma1_OPT) & (is.null(gamma1_min) | is.null(gamma1_max) | is.null(gamma1_steps))){
    stop("Set all the values of the gamma grid or an optimal value for gamma")
  }
  
  
  
  #########
  # START #
  #########
  
  # VAR design
  N<-nrow(Data_EM)
  J<-ncol(Data_EM)
  n<-N-P
  Y_EM<-as.matrix(stack_y(A = Data_EM, burnin = NULL, P=P, N=N)) # Responses
  X_EM<-stack_Xbig(A = Data_EM, burnin = 1, P=P, N=N) # Inputs
  
  # Data preliminaries
  Obj_EM=rep(1,maxit.EM)
  Obj_Conv_EM=rep(1,maxit.EM)
  
  beta_init_arr<-array(NA, c(J, J, P))
  for (ip in 1:P){
    beta_init_arr[,,ip]<-diag(J)      # Take an identity as initial value for each lag 
  }
  beta_EM<-stack(data.frame(beta_init_arr))[,1]
  
  omega_EM<-diag(J)                   # Take an identity as initial value
  tau_EM<-matrix(NA, ncol=1, nrow=n)
  delta_EM<-matrix(NA, ncol=1, nrow=n)
  
  
  iter_vec<-matrix(NA, nrow=maxit.EM, ncol=1)
  lambda1_opt_vec<-matrix(NA, nrow=maxit.EM, ncol=1)
  gamma1_opt_vec<-matrix(NA, nrow=maxit.EM, ncol=1)
  
  
  
  iter_EM<-2
  
  
  while (Obj_Conv_EM[iter_EM]>tol.EM & iter_EM<maxit.EM) {
    
    
    # Data prelimeinaries
    
    error_EM<-(Y_EM-X_EM%*%beta_EM)
    e_EM<-matrix(error_EM, nrow = n, ncol = J) 
    
    ##########
    # E-step #
    ##########
    
    delta_EM<-matrix(diag(e_EM%*%omega_EM%*%t(e_EM)), nrow=n, ncol=1)
    for (i.n in 1:n){
      tau_EM[i.n,]<-(nu + J)/(nu + delta_EM[i.n,])
    }
    
    
    ##########
    # M-step #
    ##########
    
    # Weighted responses
    Y_mat<-matrix(Y_EM, nrow = n, ncol=J) # ok 
    Y_weighted_mat<-matrix(NA, nrow = n, ncol=J)
    for (i in 1:n){
      Y_weighted_mat[i,]<-sqrt(tau_EM[i,1]) * Y_mat[i,] # ok
    }
    Y_weighted<-stack(data.frame(Y_weighted_mat))[,1]
    
    # Weighted inputs
    tau_EM_X<-matrix(rep(tau_EM, J), nrow = n*J, ncol = 1)
    X_weighted<-matrix(NA, nrow = nrow(X_EM),ncol = ncol(X_EM))
    for (i in 1:(n*J)){
      X_weighted[i,]<- sqrt(tau_EM_X[i,1]) * X_EM[i,] # ok
    }
    
    # Joint Group Lasso
    jgrlasso_fit<-JGrLasso(Y = Y_weighted, X = X_weighted, P = P, n = n, J = J, type_lasso = type_lasso,
                           maxiter.both = maxit.both, tol.both = tol.both,  
                           lambda1.min = lambda1_min, lambda1.max = lambda1_max, lambda1.steps = lambda1_steps,
                           gamma1.min = gamma1_min, gamma1.max = gamma1_max, gamma1.steps = gamma1_steps,
                           lambda1.OPT=lambda1_OPT, gamma1.OPT=gamma1_OPT)
    
    
    # Output of the Joint Group Lasso
    
    iter_vec[iter_EM,]<-jgrlasso_fit$iter
    lambda1_opt_vec[iter_EM,]<-jgrlasso_fit$lambda
    gamma1_opt_vec[iter_EM,]<-jgrlasso_fit$gamma
    
    # Coefficients
    beta_EM_NEW<-jgrlasso_fit$beta.new
    omega_EM_NEW<-jgrlasso_fit$omega.new
    
    
    ###############
    #    Check    #
    # Convergence #
    ###############
    
    # Built the objective function
    ee_EM <- (Y_weighted - X_weighted %*% beta_EM_NEW)
    ee_EM_mat<-matrix(ee_EM, nrow = n, ncol = J)
    ll_EM <- sum(diag(omega_EM_NEW %*% t(ee_EM_mat) %*% ee_EM_mat)) - (n*J*log(det(omega_EM_NEW)))
    
    Obj_EM[iter_EM] <- ll_EM
    Obj_Conv_EM[iter_EM+1]<-(abs(Obj_EM[iter_EM] - Obj_EM[iter_EM - 1]) / Obj_EM[iter_EM - 1])
    
    
    # Create the new INPUTS
    iter_EM=iter_EM + 1
    beta_EM <- beta_EM_NEW
    omega_EM <- omega_EM_NEW
    
  } # end while loop
  
  # Additional output
  beta_arr<-array(beta_EM, c(J,J,P))
  error_EM<-(Y_EM - X_EM %*% beta_EM)
  e_EM<-matrix(error_EM, nrow = n, ncol = J)
  lambda1_opt<-lambda1_opt_vec[max(which(!is.na(lambda1_opt_vec))),] 
  gamma1_opt<-gamma1_opt_vec[max(which(!is.na(gamma1_opt_vec))),] 
  
  EM_VAR <- list(Beta_new=beta_EM, Beta_arr=beta_arr, innov=e_EM,
                 Omega_new=omega_EM, tau_new=tau_EM, 
                 Obj_EM=Obj_EM, iter_EM=iter_EM, iter_vec=iter_vec, 
                 lambda1_opt=lambda1_opt, gamma1_opt=gamma1_opt)
  
  
} # end function EM_VAR
