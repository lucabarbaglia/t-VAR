#' Large.tVAR
#' 
#' ECM algorithm for the Large VAR with errors following a
#' multivariate t-distribution with estimation of the degrees of freedom
#' 
#' @param Data_ECM a NxJ matrix of log-volatilities. J: number of time series. N: time series length.
#' @param P VAR order
#' @param type_lasso type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
#' @param nu_init degrees-of-freedom initial value. Default is 1000.
#' @param maxit.ECM maximum iterations for ECM algorithm. Default is 25.
#' @param tol.ECM tolerance ECM algorithm. Default is 0.01.
#' @param maxit.both maximum iterations for gaussian lasso algorithm. Default is 50.
#' @param tol.both tolerance gaussian lasso algorithm. Default is 0.01.
#' @param maxit.nu maximum iteration for estimation of the degrees-of-freedom. Default is 1000.
#' @param lambda1_min minimum value of the regularization parameter on Beta. Default is NULL. 
#' @param lambda1_max maximum value of the regularization parameter on Beta.  Default is NULL.
#' @param lambda1_steps number of steps in the lambda grid.  Default is NULL.
#' @param gamma1_min minimum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_max maximum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_steps  number of steps in the gamma grid.  Default is NULL.
#' @param lambda1_OPT optimal value of the regularization parameter on Beta.  Default is NULL.
#' @param gamma1_OPT optimal value of the regularization parameter on Omega.  Default is NULL.
#' 
#' @return A list containing:
#' \item{"Beta_new"}{a vector containing the estimated Beta.}
#' \item{"Beta_arr"}{a JxJxP array containing the estimated Beta.}
#' \item{"innov"}{a (N-P)xJ containing the estimated VAR residuals.}
#' \item{"Omega_new"}{a JxJ matrix containing the estimated Omega.}
#' \item{"tau_new"}{a vector of length N-P containing the estimated gamma variable tau.}
#' \item{"nu_new"}{estimated degrees-of-freedom of the multivaraite t-distribution of the VAR innovations.} 
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


Large.tVAR<-function(Data_ECM, P, type_lasso="Lasso", nu_init=NULL,
                     maxit.ECM=25, tol.ECM=0.01, tol.nu=0.001, maxit.both=50, tol.both=0.01, maxit.nu=100,
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
  
  if(is.null(nu_init)){nu_init=1000}
  
  if(is.not.matrix(Data_ECM)){stop("Data must be in matrix form.")}
  
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
  N<-nrow(Data_ECM)
  J<-ncol(Data_ECM)
  n<-N-P
  Y_ECM<-as.matrix(stack_y(A = Data_ECM, burnin = NULL, P=P, N=N)) # Responses
  X_ECM<-stack_Xbig(A = Data_ECM, burnin = NULL, P=P, N=N) # Inputs
  
  # Data preliminaries
  Obj_ECM=rep(1,maxit.ECM)
  Obj_Conv_ECM=rep(1,maxit.ECM)
  
  beta_init_arr<-array(NA, c(J, J, P))
  for (ip in 1:P){
    beta_init_arr[,,ip]<-diag(J)      # Take an identity as initial value for each lag 
  }
  beta_ECM<-stack(data.frame(beta_init_arr))[,1]
  
  omega_ECM<-diag(J)                   # Take an identity as initial value
  tau_ECM<-matrix(NA, ncol=1, nrow=n)
  delta_ECM<-matrix(NA, ncol=1, nrow=n)
  
  nu_ECM<-matrix(NA, ncol=1, nrow = maxit.ECM)
  nu_ECM[1,1]<-nu_init
  
  iter_vec<-matrix(NA, nrow=maxit.ECM, ncol=1)
  lambda1_opt_vec<-matrix(NA, nrow=maxit.ECM, ncol=1)
  gamma1_opt_vec<-matrix(NA, nrow=maxit.ECM, ncol=1)
  
  
  
  iter_ECM<-2
  
  
  while (Obj_Conv_ECM[iter_ECM]>tol.ECM & iter_ECM<maxit.ECM) {
    
    
    # Data prelimeinaries
    
    error_ECM<-(Y_ECM-X_ECM%*%beta_ECM)
    e_ECM<-matrix(error_ECM, nrow = n, ncol = J) 
    
    ############
    # E-step 1 #
    ############
    
    delta_ECM<-matrix(diag(e_ECM%*%omega_ECM%*%t(e_ECM)), nrow=n, ncol=1)
    for (i.n in 1:n){
      tau_ECM[i.n,]<-(nu_ECM[(iter_ECM-1),1] + J)/(nu_ECM[(iter_ECM-1),1] + delta_ECM[i.n,])
    }
    
    
    #############
    # CM-step 1 #
    #############
    
    # Weighted responses
    Y_mat<-matrix(Y_ECM, nrow = n, ncol=J) # ok 
    Y_weighted_mat<-matrix(NA, nrow = n, ncol=J)
    for (i in 1:n){
      Y_weighted_mat[i,]<-sqrt(tau_ECM[i,1]) * Y_mat[i,] # ok
    }
    Y_weighted<-stack(data.frame(Y_weighted_mat))[,1]
    
    # Weighted inputs
    tau_ECM_X<-matrix(rep(tau_ECM, J), nrow = n*J, ncol = 1)
    X_weighted<-matrix(NA, nrow = nrow(X_ECM),ncol = ncol(X_ECM))
    for (i in 1:(n*J)){
      X_weighted[i,]<- sqrt(tau_ECM_X[i,1]) * X_ECM[i,] # ok
    }
    
    # Joint Group Lasso
    jgrlasso_fit<-JGrLasso(Y = Y_weighted, X = X_weighted, P = P, n = n, J = J, type_lasso = type_lasso,
                           maxiter.both = maxit.both, tol.both = tol.both,  
                           lambda1.min = lambda1_min, lambda1.max = lambda1_max, lambda1.steps = lambda1_steps,
                           gamma1.min = gamma1_min, gamma1.max = gamma1_max, gamma1.steps = gamma1_steps,
                           lambda1.OPT=lambda1_OPT, gamma1.OPT=gamma1_OPT)
    
    # Output of the Joint Group Lasso
    iter_vec[iter_ECM,]<-jgrlasso_fit$iter
    lambda1_opt_vec[iter_ECM,]<-jgrlasso_fit$lambda
    gamma1_opt_vec[iter_ECM,]<-jgrlasso_fit$gamma
    
    # Coefficients
    beta_ECM_NEW<-jgrlasso_fit$beta.new
    omega_ECM_NEW<-jgrlasso_fit$omega.new
    
    
    ############
    # E-step 2 #
    ############
    
    error_ECM_update<-(Y_ECM - X_ECM %*% beta_ECM_NEW)
    e_ECM_update<-matrix(error_ECM_update, nrow = n, ncol = J)
    
    delta_ECM_update <- matrix(diag(e_ECM_update %*% omega_ECM_NEW %*% t(e_ECM_update)), nrow=n, ncol=1)
    for (i.n in 1:n){
      tau_ECM[i.n,]<-(nu_ECM[(iter_ECM-1),1] + J)/(nu_ECM[(iter_ECM-1),1] + delta_ECM_update[i.n,])
    }
    
    #############
    # CM-step 2 #
    #############
    
    # One dimensional search
    df_grad<-function(nu, delta, J){
      -digamma(nu/2) + log(nu/2) + (sum(log((nu+J)/(delta+nu)) - ((nu+J)/(delta+nu)))/(length(delta))) + 
        1 +digamma((nu+J)/2) - log((nu+J)/2)
    }
    fit_nu<-stats::uniroot(df_grad, c(0.1, 100), extendInt = "yes", tol = tol.nu, maxiter = maxit.nu, delta=delta_ECM_update, J=J)
    
    if (iter_ECM==2){
      if(fit_nu$root > nu_init) {
        nu_ECM[iter_ECM,1]<-nu_ECM[iter_ECM-1,1]
      } else {
        nu_ECM[iter_ECM,1]<-fit_nu$root
      }
    }
    if (iter_ECM>2){
      if(fit_nu$root > nu_init) {
        nu_ECM[iter_ECM,1]<-nu_ECM[iter_ECM-2,1]
      } else {
        nu_ECM[iter_ECM,1]<-fit_nu$root
      }
    }
    
    
    ###############
    #    Check    #
    # Convergence #
    ###############
    
    # Built the objective function
    ee_ECM <- (Y_weighted - X_weighted %*% beta_ECM_NEW)
    ee_ECM_mat<-matrix(ee_ECM, nrow = n, ncol = J)
    ll_ECM <- -1/2*sum(diag(omega_ECM_NEW %*% t(ee_ECM_mat) %*% ee_ECM_mat)) - (n*log(det(omega_ECM_NEW))/2) 
    
    
    Obj_ECM[iter_ECM] <- ll_ECM
    Obj_Conv_ECM[iter_ECM+1]<-(abs(Obj_ECM[iter_ECM] - Obj_ECM[iter_ECM - 1]) / Obj_ECM[iter_ECM - 1])
    
    
    # Create the new INPUTS
    iter_ECM=iter_ECM + 1
    beta_ECM <- beta_ECM_NEW
    omega_ECM <- omega_ECM_NEW
    
  } # end while loop
  
  # Additional output
  beta_arr<-array(beta_ECM, c(J,J,P))
  error_ECM<-(Y_ECM - X_ECM %*% beta_ECM)
  e_ECM<-matrix(error_ECM, nrow = n, ncol = J)
  nu_new<-nu_ECM[max(which(!is.na(nu_ECM))),] 
  lambda1_opt<-lambda1_opt_vec[max(which(!is.na(lambda1_opt_vec))),] 
  gamma1_opt<-gamma1_opt_vec[max(which(!is.na(gamma1_opt_vec))),] 
  
  
  Large.tVAR <- list(Beta_new=beta_ECM, Beta_arr=beta_arr, innov=e_ECM,
                     Omega_new=omega_ECM, tau_new=tau_ECM, nu_new=nu_new, 
                     Obj_ECM=Obj_ECM, iter_ECM=iter_ECM, iter_vec=iter_vec, 
                     lambda1_opt=lambda1_opt, gamma1_opt=gamma1_opt)
  
  
}