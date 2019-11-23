###################################################
### Large_tVAR.R: functions to estimate the VAR ###
###################################################


##################
### Large.tVAR ###
##################

# ECM algorithm for the Large VAR with errors following a multivariate t-distribution with estimation of the degrees of freedom

Large.tVAR<-function(Data_ECM, P, type_lasso="Lasso", nu_init=NULL,
                    maxit.ECM=25, tol.ECM=0.01, tol.nu=0.001, maxit.both=50, tol.both=0.01, maxit.nu=100,
                    lambda1_min=NULL, lambda1_max=NULL, lambda1_steps=NULL, 
                    gamma1_min=NULL, gamma1_max=NULL, gamma1_steps=NULL,
                    lambda1_OPT=NULL, gamma1_OPT=NULL){
  
  ##########
  # INPUTS #
  ##########
  # Data_ECM: a NxJ matrix of log-volatilities:
  # - J: number of time series
  # - N: time series length
  # P: VAR order
  # type_lasso: type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
  # nu_init: degrees-of-freedom initial value. Default is 1000.
  # maxit.ECM: maximum iterations for ECM algorithm. Default is 25.
  # tol.ECM: tolerance ECM algorithm. Default is 0.01.
  # maxit.both: maximum iterations for gaussian lasso algorithm. Default is 50.
  # tol.both: tolerance gaussian lasso algorithm. Default is 0.01.
  # mxit.nu: maximum iteration for estimation of the degrees-of-freedom. Default is 1000.
  # lambda1_min : minimum value of the regularization parameter on Beta. Default is NULL. 
  # lambda1_max : maximum value of the regularization parameter on Beta.  Default is NULL.
  # lambda1_steps : number of steps in the lambda grid.  Default is NULL.
  # gamma1_min : minimum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_max : maximum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_steps : number of steps in the gamma grid.  Default is NULL.
  # lambda1_OPT: optimal value of the regularization parameter on Beta.  Default is NULL.
  # gamma1_OPT: optimal value of the regularization parameter on Omega.  Default is NULL.
  
  ##########
  # OUTPUT #
  ##########
  # Beta_new: a vector cointaining the estimated Beta.
  # Beta_arr: a JxJxP array cointaining the estimated Beta.
  # innov: a (N-P)xJ containing the estimated VAR residuals.
  # Omega_new: a JxJ matrix cointaining the estimated Omega.
  # tau_new: a vector of length N-P cointaining the estimated gamma variable tau.
  # nu_new: estimated degrees-of-freedom of the multivaraite t-distribution of the VAR innovations.
  # Obj_ECM: objective function.
  # iter_ECM: number of iterations of the ECM algorithm.
  # iter_vec: number of iteration of the Gaussian Lasso algoithm for each ECM iteration.
  # lambda1_opt: selected value of the regularization parameter on Beta. 
  # gamma1_opt: selected value of the regularization parameter on Omega. 
  
  
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
    fit_nu<-uniroot(df_grad, c(0.1, 100),extendInt = "yes", tol = tol.nu, maxiter = maxit.nu, delta=delta_ECM_update, J=J)
    
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
  
  
} # end function Large.tVAR



########################
### EM algorithm VAR ###
########################

# EM algorithm for the Robus Sparse VAR with errors following a multivariate t-distribution (without estimation of df)


EM_VAR<-function(Data_EM, P, type_lasso="Lasso", 
                 nu, maxit.EM=25, tol.EM=0.01,maxit.both=50, tol.both=0.01, 
                 lambda1_min=NULL, lambda1_max=NULL, lambda1_steps=NULL, 
                 gamma1_min=NULL, gamma1_max=NULL, gamma1_steps=NULL,
                 lambda1_OPT=NULL, gamma1_OPT=NULL){
  
  
  ##########
  # INPUTS #
  ##########
  # Data_EM: a NxJ matrix of log-volatilities:
  # - J: number of time series
  # - N: time series length
  # P: VAR order
  # type_lasso: type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
  # nu: degrees-of-freedom of the multivaraite t-distribution of the VAR innovations.
  # maxit.EM: maximum iterations for EM algorithm. Default is 25.
  # tol.EM: tolerance EM algorithm. Default is 0.01.
  # maxit.both: maximum iterations for gaussian lasso algorithm. Default is 50.
  # tol.both: tolerance gaussian lasso algorithm. Default is 0.01.
  # lambda1_min : minimum value of the regularization parameter on Beta. Default is NULL. 
  # lambda1_max : maximum value of the regularization parameter on Beta.  Default is NULL.
  # lambda1_steps : number of steps in the lambda grid.  Default is NULL.
  # gamma1_min : minimum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_max : maximum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_steps : number of steps in the gamma grid.  Default is NULL.
  # lambda1_OPT: optimal value of the regularization parameter on Beta.  Default is NULL.
  # gamma1_OPT: optimal value of the regularization parameter on Omega.  Default is NULL.
  
  ##########
  # OUTPUT #
  ##########
  # Beta_new: a vector cointaining the estimated Beta.
  # Beta_arr: a JxJxP array cointaining the estimated Beta.
  # innov: a (N-P)xJ containing the estimated VAR residuals.
  # Omega_new: a JxJ matrix cointaining the estimated Omega.
  # tau_new: a vector of length N-P cointaining the estimated gamma variable tau.
  # Obj_ECM: objective function.
  # iter_ECM: number of iterations of the ECM algorithm.
  # iter_vec: number of iteration of the Gaussian Lasso algoithm for each ECM iteration.
  # lambda1_opt: selected value of the regularization parameter on Beta. 
  # gamma1_opt: selected value of the regularization parameter on Omega. 
  
  
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



#################
### JGrLasso ####
#################

# Function for the joint estimation of of beta and omega 

JGrLasso <- function(Y, X, P, n, J, type_lasso="Lasso",
                     maxiter.both=50,  tol.both=0.01, 
                     lambda1.min=NULL, lambda1.max=NULL, lambda1.steps=NULL, 
                     gamma1.min=NULL, gamma1.max=NULL, gamma1.steps=NULL,
                     lambda1.OPT=NULL, gamma1.OPT=NULL){
  
  ##########
  # INPUTS #
  ##########
  # Y: a vector of length J(N-P) cointaining the responses. Obtained with funxtion "stack_y".
  # - J: number of time series
  # - N: time series length
  # X: a J(N-P)xPJ^2 matrix cointainig the regressor. Obtained with function "stack_Xbig".
  # P: VAR order.
  # type_lasso: type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
  # maxit.both: maximum iterations for gaussian lasso algorithm. Default is 50.
  # tol.both: tolerance gaussian lasso algorithm. Default is 0.01.
  # lambda1_min : minimum value of the regularization parameter on Beta. Default is NULL. 
  # lambda1_max : maximum value of the regularization parameter on Beta.  Default is NULL.
  # lambda1_steps : number of steps in the lambda grid.  Default is NULL.
  # gamma1_min : minimum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_max : maximum value of the regularization parameter on Omega.  Default is NULL.
  # gamma1_steps : number of steps in the gamma grid.  Default is NULL.
  # lambda1_OPT: optimal value of the regularization parameter on Beta.  Default is NULL.
  # gamma1_OPT: optimal value of the regularization parameter on Omega.  Default is NULL.
  
  ##########
  # OUTPUT #
  ##########
  # beta.new: a vector cointaining the estimated Beta.
  # beta.arr: a JxJxP array cointaining the estimated Beta.
  # omega.new: a JxJ matrix cointaining the estimated Omega.
  # Obj_JGrL: objective function.
  # iter: number of iterations.
  # lambda: selected value of the regularization parameter on Beta. 
  # gamma: selected value of the regularization parameter on Omega. 
  

  ########################
  # Additional Functions #
  ########################

  # Spectral decomposition
  spectral_dec<-function(Omega){
    decomp<-eigen(Omega, symmetric = TRUE)  
    CC <- decomp$vectors
    Lambda <- diag(sqrt(decomp$values))
    P <- CC %*% Lambda %*% t(CC)
  } 
  
  # Opposite of null
  is.not.null <- function(x) ! is.null(x)

  
  #####################
  # Preliminary steps #
  #####################
  
  if (is.not.null(lambda1.OPT) & (is.not.null(lambda1.min) | is.not.null(lambda1.max) | is.not.null(lambda1.steps))){
    stop("Either set the optimal value of lambda, either the values for the grid search")
  }
  if (is.not.null(gamma1.OPT) & (is.not.null(gamma1.min) | is.not.null(gamma1.max) | is.not.null(gamma1.steps))){
    stop("Either set the optimal value of gamma, either the values for the grid search")
  }
  
  if (is.null(lambda1.OPT) & (is.null(lambda1.min) | is.null(lambda1.max) | is.null(lambda1.steps))){
    stop("Set all the values of the lambda grid or an optimal value for lambda")
  }
  if (is.null(gamma1.OPT) & (is.null(gamma1.min) | is.null(gamma1.max) | is.null(gamma1.steps))){
    stop("Set all the values of the gamma grid or an optimal value for gamma")
  }
  
  
  #########
  # Start #
  #########
  Obj_JGrL = rep(1,maxiter.both)
  Obj_Conv_JGrL=rep(1,maxiter.both)
  iter<-2
  
  Y.new<-Y
  X.new<-X
  
  while (Obj_Conv_JGrL[iter]>tol.both & iter<maxiter.both) {
    
    ############
    ### BETA ###
    ############
    
    if (type_lasso == "Group"){    # Use the Group Lasso to estimate the autoregressive coefficients
      
      index_arr<-array(NA, c(1,J^2,P))
      for(ip in 1:P){
        index_arr[,,ip]<-seq(1:(J^2))
      }
      index_list<-alply(index_arr, 3, .dims = TRUE)
      index<-as.vector(unlist(index_list))
      
      # Regularization parameters
      if (is.not.null(lambda1.min) & is.not.null(lambda1.max) & is.not.null(lambda1.steps)){
        lambda1_grid<-seq(from=lambda1.min, to=lambda1.max, length=lambda1.steps)
        BIC_values_beta<-lapply(lambda1_grid, BICfunction_grplasso, 
                                Xbic = X.new, Ybic = Y.new, indexbic = index)
        lambda1_opt<-lambda1_grid[which.min(BIC_values_beta)]                
      }
      if (is.not.null(lambda1.OPT)){
        lambda1_opt<-lambda1.OPT
      }
      
      
      # Fit
      grplasso_fit<-grplasso(x = X.new, y = Y.new, index=index, lambda = lambda1_opt, 
                             model = LinReg(), center = F,  control = grpl.control(trace=0))
      beta.new<-grplasso_fit$coefficients
      
    } 
    
    if (type_lasso == "Lasso") {                # Use the simple the non-grouped lasso
      
      # Grid search for parameters
      index_arr<-array(seq(1:(P*J^2)), c(1,J^2,P))
      index_list<-alply(index_arr, 3, .dims = TRUE)
      index<-as.vector(unlist(index_list))
      
      # Regularization parameters
      if (is.not.null(lambda1.min) & is.not.null(lambda1.max) & is.not.null(lambda1.steps)){
        lambda1_grid<-seq(from=lambda1.min, to=lambda1.max, length=lambda1.steps)
        BIC_values_beta<-lapply(lambda1_grid, BICfunction_grplasso, 
                                Xbic = X.new, Ybic = Y.new, indexbic = index)
        lambda1_opt<-lambda1_grid[which.min(BIC_values_beta)]                
      }
      if (is.not.null(lambda1.OPT)){
        lambda1_opt<-lambda1.OPT
      }
      
      # Fit
      grplasso_fit<-grplasso(x = X.new, y = Y.new, index=index, lambda = lambda1_opt, 
                             model = LinReg(), center = F,  control = grpl.control(trace=0))
      beta.new<-grplasso_fit$coefficients
      
    }  
    

    
    #########
    # OMEGA #
    #########
    
    # Consider the residuals 
    e_grplasso<-(Y.new - X.new %*% beta.new) # Residuals
    e_grplasso_mat<-matrix(e_grplasso, nrow = n, ncol = J)
    s_grplasso<-cov(e_grplasso_mat)                    # Var-Cov matrix 
    
    # Regularuzation parameters
    if (is.not.null(gamma1.min) & is.not.null(gamma1.max) & is.not.null(gamma1.steps)){
      gamma1_grid<-seq(from=gamma1.min, to=gamma1.max, length=gamma1.steps)
      BIC_values_omega<-lapply(gamma1_grid, BICfunction_glasso, sbic = s_grplasso, ndata = n*J)
      gamma1_opt<-gamma1_grid[which.min(BIC_values_omega)]  
    }
    if (is.not.null(gamma1.OPT)){
      gamma1_opt<-gamma1.OPT
    }
    
    # Fit
    glasso_fit<-glasso(s=s_grplasso,rho=gamma1_opt)
    omega.new<-glasso_fit$wi
    
    
    ###############
    #    Check    #
    # Convergence #
    ###############
    
    ee <- (Y.new - X.new %*% beta.new)
    omega.new_tilde<-kronecker(diag(n), omega.new)
    ll<- ((t(ee) %*% omega.new_tilde %*% ee) - (n*J*log(det(omega.new))))
    
    Obj_JGrL[iter] <- ll
    Obj_Conv_JGrL[iter+1]<-(abs(Obj_JGrL[iter]-Obj_JGrL[iter-1])/Obj_JGrL[iter-1])
    

    ##############
    # NEW INPUTS #
    ##############
    
    # Spectral decomposition
    P_dec<-spectral_dec(omega.new) # Spectral decomposition of omega
    Y.new_mat<-matrix(Y.new, nrow = n, ncol = J)
    Y_new_mat<-Y.new_mat %*% P_dec
    Y_new<-stack(data.frame(Y_new_mat))[,1]
    
    P_decX<-kronecker(diag(P), kronecker(diag(J), P_dec))
    X_new<-X.new %*% P_decX
    
    # New inputs
    iter<-iter+1
    Y.new<-Y_new
    X.new<-X_new
    
    
  } # end of while loop
  
  beta.arr<-array(beta.new, c(J,J,P))
  
  JGrLasso <- list(beta.new = beta.new, beta.arr=beta.arr, omega.new = omega.new, 
                   iter = iter, Obj_JGrL = Obj_JGrL, 
                   lambda = lambda1_opt, gamma = gamma1_opt)
  
}


#################
### VARdesign ###
#################

# VARdesign is a function that construct the design matrix of the VAR model for a multivariate time series Y

VARdesign <- function(Y,p){


  #########
  # INPUT #
  #########
  # Y: mulitvariate time series.
  # p: lag length.
 
  ##########
  # OUTPUT #
  ##########
  # a design matrix.
  
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


################
### VARbyGLS ###
################

# Function for the classic unpenalized VAR esitmation, Generalized Least Squares

VARbyGLS <- function(Y,p){

  #########
  # INPUT #
  #########
  # Y: mulitvariate time series.
  # p: lag length.

  ##########
  # OUTPUT #
  ##########
  # Beta: a list of betas.
  # Sigma: var-cov matrix.
  # Beta.vec=vector of betas [b^1_11,b^2_11,..,b^P_11,b^1_21,b^2_21,..,b^P_21,.., b^P_JJ].
  
  
  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  Response <- c(Y[-(1:p),])
  Xmatrix <- VARdesign(Y,p)
  Xmatrix.kron <- cbind(1,kronecker(diag(1,q),Xmatrix))
  
  # Obtain OLS fit residuals
  OLSresid <- matrix(resid(lm(Response~Xmatrix.kron)),ncol=q)
  Sigma <- cov(OLSresid)
  # Do GLS
  Betahat <- solve(t(Xmatrix.kron)%*%kronecker(solve(Sigma),diag(1,n-p))%*%Xmatrix.kron)%*%t(Xmatrix.kron)%*%kronecker(solve(Sigma),diag(1,n-p))%*%Response
  Betahat <- Betahat[-1]
  Betahatt <- list()
  for(i in 1:p){
    Betahatt[[i]] <- matrix(Betahat[seq(i,(p*q*q),by=p)],ncol=q,byrow=T)
  }
  out <- list(Beta=Betahatt,Sigma=Sigma, Beta.vec=Betahat)
}



###########
# stack_y #
###########

# Build the vector of responses

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

##############
# stack_Xbig #
##############

# Build the matrix of regressors

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
  X0_list<-alply(X0_arr1, 3,.dims=TRUE) # Convert the array into a list
  X0<-do.call(cbind,X0_list) # Matrix; the first J columns are at lag p=1, the last J columns are at lag p=P. The rows go from P-p to N-p.
  Xbig_mat<-kronecker(diag(J), X0)
  return(Xbig_mat)
}



########################
# BICfunction_grplasso #
########################

# Function for the choice of the regularization parameters for the autoregressive coeffcients based on BIC

"BICfunction_grplasso" = function(lambda1, Ybic, Xbic, indexbic){
  
  #########
  # INPUT #
  #########
  # Ybic : vector of dimension NJ x 1 of dependent variables 
  # Xbic : matrix of dimension NJ x J²P of indipendent variables 
  # lambda1 : regularization parameter for the autoregressive coeffcients
  # indexbic : group structure
  
  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values
  
  # Estimate the autoregressive coefficients using SPG algorithm
  FIT_grplasso<-grplasso(x = Xbic, y = Ybic, index=indexbic, lambda = lambda1, 
                         model = LinReg(), center = F,  control = grpl.control(trace=0))
  
  # Compute the BIC
  NLog_lik_grplasso<-FIT_grplasso$nloglik # negative log-lik
  df_grplasso<-length(which(FIT_grplasso$coefficients!=0))
  BIC_value_grplasso<-(2*NLog_lik_grplasso+log(length(Ybic))*df_grplasso)
  
  return(BIC_value_grplasso)
  
}

#####################
# BICfunction_lasso #
#####################

# Function for the choice of the regularization parameters for the autoregressive coeffcients based on BIC

"BICfunction_lasso" = function(lambda1, Ybic, Xbic){
  
  #########
  # INPUT #
  #########
  # Ybic : vector of dimension NJ x 1 of dependent variables 
  # Xbic : matrix of dimension NJ x J²P of indipendent variables 
  # lambda1 : regularization parameter for the autoregressive coeffcients

  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values
  
  # Estimate the autoregressive coefficients using SPG algorithm
  FIT_lasso <- glmnet(x = Xbic,  y = Ybic, family = "gaussian", lambda = lambda1)
  # Compute the BIC
  nj<-length(Ybic)
  bbb<-as.matrix(FIT_lasso$beta)
  Log_lik_lasso<- (1/nj) *t((Ybic - Xbic %*% bbb))%*%(Ybic - Xbic %*% bbb)
  df_lasso<-length(which(bbb != 0))
  BIC_value_lasso<-( - 2*Log_lik_lasso + log(nj)*df_lasso)
  
  return(BIC_value_lasso)
  
}


######################
# BICfunction_glasso #
######################

# Function for the choice of the regularization parameters for the inverse error covariance matrix based on BIC

"BICfunction_glasso" = function(gamma1, sbic, ndata){
  
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
  FIT_glasso<-glasso(s = sbic,  rho=gamma1)
  FIT_glasso$wi
  
  # Compute the BIC
  NLog_lik_glasso<-(- FIT_glasso$loglik) # negative log-lik
  df_glasso<-length(which(FIT_glasso$wi!=0))
  BIC_value_glasso<-(2*NLog_lik_glasso+log(ndata)*df_glasso)
  
  return(BIC_value_glasso)
  
}


##############
# Spillovers #
##############

# Function to compute volatility spillovers and its index based on forecast error variance decomposition

Spillovers<-function(fit, dec.type="Generalized", h=1){
  
  #########
  # INPUT #
  #########
  # fit: t-VAR fit obtained from the function Large.tVAR
  # dec.type : decomposition for forecast error varaince decomposition: "Cholesky", "Spectral", "Generalized", "GeneralizedIRF" or "Generalized_Lanne". Default "Generalized".
  # h: forecast horizon for variance decomposition. Default 1.
  
  ##########
  # OUTPUT #
  ##########
  # spill: matrix of volatility spillovers
  # spill_index: volatility spillover index
  
  
  fevd_mat<-FEVDec(B_arr=fit$Beta_arr, Sig=solve(fit$Omega_new), lag=h, dec.type = dec.type)
  
  spill<-100*fevd_mat
  
  # Remove the own-effects (main diagonal)
  spill_cross<-spill
  diag(spill_cross)<-0
  
  # Global Spillovers (exclude the own-effects)
  spill_gl<-sum(spill_cross)
  
  Spillovers<-list("spill"=spill, "spill_index"=spill_gl)
}






########
# FEVD #
########

# Obtain Forecast Error variance Decomposition with Spectral, Cholesky or Generealized decomposition

FEVDec<-function (B_arr, Sig, lag, dec.type="Generalized") {
  
  
  #########
  # INPUT #
  #########
  # B_arr: arrays of autoregressive coeff 
  # Sig: var-cov matrix
  # lag: forecast horizon
  # dec.type : "Cholesky", "Spectral", "Generalized", "GeneralizedIRF" or "Generalized_Lanne". Default "Generalized".
  
  ##########
  # OUTPUT #
  ##########
  # fevd: forecast error variance decompositoin matrix
  
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


##############
# network.vs #
##############

# Function to obtain network of volatility spillovers

network.vs<-function(spill, main.label=NULL, names=NULL){
  
  
  #########
  # INPUT #
  #########
  # spill: a JxJ matrix of volatility spillovers obtained with the Spillovers functions.
  # main.label: network names. Must be of length K.  Default NULL.
  # names: node names. Default NULL.
  
  ##########
  # OUTPUT #
  ##########
  # network of volatility spillovers
  
  # Dimension
  J<-ncol(spill)
  
  # Nodes labels
  if (!is.null(names)){
    colnames(spill)<-rownames(spill)<-names
  }
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  lab.locs <- radian.rescale(x=1:J, direction=-1, start=0) 
  
  GRAPH<-graph.adjacency(adjmatrix=t(spill), mode="directed", diag=FALSE, weighted=T)
  
  E(GRAPH)$width <- abs(E(GRAPH)$weight)
  E(GRAPH)$color<-"darkblue"
  
  loc<-layout.circle(GRAPH) # location of the nodes
  
  plot(margin=rep(10^-10,4), GRAPH, layout=layout.circle(GRAPH), 
       vertex.label.dist=1, 
       vertex.color="white", 
       vertex.frame.color="white",
       vertex.shape="circle", vertex.label.color="black"
       , vertex.label.degree=lab.locs
       , edge.curved=seq(-0.3, 0.3, length = ecount(GRAPH))
       , edge.arrow.width=1.3
  )
  
  if(!is.null(main.label)){
    title(main=main.label, cex.main=2)
  }
  
  
} # end function



