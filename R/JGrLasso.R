#' JGrLasso
#' 
#' Function for the joint estimation of of beta and omega 
#' 
#' @param Y a vector of length J(N-P) cointaining the responses. Obtained with funxtion "stack_y". J: number of time series. N: time series length.
#' @param X a J(N-P)xPJ^2 matrix cointainig the regressor. Obtained with function "stack_Xbig".
#' @param P VAR order.
#' @param type_lasso type of lasso penalty. "Lasso" for standard lasso, "Group" for group lasso. Default is "Lasso".
#' @param maxit.both maximum iterations for gaussian lasso algorithm. Default is 50.
#' @param tol.both tolerance gaussian lasso algorithm. Default is 0.01.
#' @param lambda1_min minimum value of the regularization parameter on Beta. Default is NULL. 
#' @param lambda1_max maximum value of the regularization parameter on Beta.  Default is NULL.
#' @param lambda1_steps number of steps in the lambda grid.  Default is NULL.
#' @param gamma1_min minimum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_max maximum value of the regularization parameter on Omega.  Default is NULL.
#' @param gamma1_steps number of steps in the gamma grid.  Default is NULL.
#' @param lambda1_OPT optimal value of the regularization parameter on Beta.  Default is NULL.
#' @param gamma1_OPT optimal value of the regularization parameter on Omega.  Default is NULL.
#' 
#' @return A list containing two objects:
#' \item{"beta.new"}{a vector containing the estimated Beta.}
#' \item{"beta.arr"}{a JxJxP array containing the estimated Beta.}
#' \item{"omega.new"}{a JxJ matrix containing the estimated Omega.}
#' \item{"Obj_JGrL"}{objective function.}
#' \item{"iter"}{number of iterations.}
#' \item{"lambda"}{selected value of the regularization parameter on Beta.}
#' \item{"gamma"}{selected value of the regularization parameter on Omega.}
#' 
#' @references Barbaglia, L., Croux, C., & Wilms, I. (2020). Volatility spillovers in commodity markets: A large t-vector autoregressive approach. Energy Economics, 85, 104555.
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export



JGrLasso <- function(Y, X, P, n, J, type_lasso="Lasso",
                     maxiter.both=50,  tol.both=0.01, 
                     lambda1.min=NULL, lambda1.max=NULL, lambda1.steps=NULL, 
                     gamma1.min=NULL, gamma1.max=NULL, gamma1.steps=NULL,
                     lambda1.OPT=NULL, gamma1.OPT=NULL){
  
  
  
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
      index_list<-plyr::alply(index_arr, 3, .dims = TRUE)
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
      grplasso_fit<-grplasso::grplasso(x = X.new, y = Y.new, index=index, lambda = lambda1_opt, 
                             model = grplasso::LinReg(), center = F,  control = grplasso::grpl.control(trace=0))
      beta.new<-grplasso_fit$coefficients
      
    } 
    
    if (type_lasso == "Lasso") {                # Use the simple the non-grouped lasso
      
      # Grid search for parameters
      index_arr<-array(seq(1:(P*J^2)), c(1,J^2,P))
      index_list<-plyr::alply(index_arr, 3, .dims = TRUE)
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
      grplasso_fit<-grplasso::grplasso(x = X.new, y = Y.new, index=index, lambda = lambda1_opt, 
                             model = grplasso::LinReg(), center = F,  control = grplasso::grpl.control(trace=0))
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
    glasso_fit<-glasso::glasso(s=s_grplasso,rho=gamma1_opt)
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
