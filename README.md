
# Read-me Large t-VAR

Barbaglia, Croux, Wilms (2019) "Volatility Spillovers and Heavy Tails: A Large t-Vector AutoRegressive Approach", Energy Economics.
https://www.sciencedirect.com/science/article/pii/S0140988319303500?via%3Dihub

This file contains a brief illustration of the code in order to implement the t-Lasso estimator, carry on the volatility spillover analysis and plot the volatility spillover network.


Author's website: https://lucabarbaglia.github.io/


``` r 
###############
# PRELIMINARY #
###############

# Preliminary step: check if all necessary packages are installed 
checkpackage<-function(U){
  if((U %in% rownames(installed.packages()))==F){
    install.packages(U)
  }
}
packagelist<-list("plyr", "grplasso", "glasso","zoo","Matrix","MASS","xtable","CADFtest","plm","foreach","doSNOW","QRM","vars","igraph","ggplot2","scales","rARPACK","ghyp","magic","MVN","glmnet")
lapply(packagelist,checkpackage)
# Load packages
suppressMessages(suppressWarnings(packages <- lapply(packagelist, FUN = function(x) {
  library(x, character.only = TRUE)
})))

###########
# t-Lasso #
###########

# 1.1 t-Lasso with estimation of the degrees-of-freedom:
fit<-Large.tVAR(Data=DATA, P=P, lambda1_OPT = 5, gamma1_OPT = 0.2) 

# 1.2 t-Lasso with estimation of the degrees-of-freedom:
fit<-EM_VAR(Data=DATA, P=P, lambda1_OPT = 5, gamma1_OPT = 0.2, nu=3) 

# where:
# - DATA is a NxJ matrix of log-volatilities, where:
#	    - N is the time series length
#            - J is the number of time series
# - P is the VAR order.
# - lambda1_OPT: value of the regularization parameter lambda.
# - gamma1_OPT: value of the regularization parameter gamma.
# - nu: degrees-of-freedom of the multivaraite t-distribution for the VAR innnovations.

# Alternatively, lambda and gamma can be selected over a grid search.


#########################
# Volatility Spillovers #
#########################

# 2. Volatility Spillovers:
vs<-Spillovers(fit = fit)
vs$spill 		# volatility spillover matrix
vs$spill_index 		# volatility spillover index

###########
# Network #
###########

# 3. Volatility spillover network:
network.vs(spill=vs$spill)

```





