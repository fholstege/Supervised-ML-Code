
# Put here the packages required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(stargazer,
               lmtest,
               matlib, # install dependencies for dsmle
               SVMMaj,
               ISLR,
               plotrix,
               e1071)

# get local package, dsmle
install.packages("dsmle_1.0-4.tar.gz", repos = NULL, type="source")
library(dsmle)



# get dataset 
load("Airline.Rdata")

# define independent and dependent variable
mX = as.matrix(Airline[,-4])
mX_scaled = scale(mX)
mY = as.matrix(Airline[,4])


# creates rbf kernel from independent variables
create_RBF_kernel = function(mX, gamma, n){
  
  # calculate distance matrix
  XtX <- tcrossprod(mX)
  XX <- matrix(1, n) %*% diag(XtX)
  D <- XX - 2*XtX + t(XX) # distance matrix
  
  # calculate rbf-kernel matrix
  K <- exp(gamma* D)
  
  return(K)
  
}

# loss function for kernel ridge regression
loss_krr = function(q_hat, mXXt, mY, mBeta_zero, lambda, J){
  
  # define the SSE for the intercept and the other beta's
  SSE_intercept = sum(abs(mY - mBeta_zero[[1]])^2)
  SSE_independent = sum(abs((J %*% mY) - q_hat ))
  
  # penalty term
  penalty = lambda * (t(q_hat) %*% inv(mXXt) %*% q_hat)
  
  # combine, return loss
  loss = SSE_intercept + SSE_independent + penalty 
  return(loss)
  
}

# executes kernel ridge regression
kernel_ridge <- function(mX, mY, lambda, type_kernel = "Linear"){
  
  # define n, p
  n = nrow(mX)
  p = ncol(mX)
  
  # initial XX^T, X^TX
  mXXt = mX %*% t(mX)
  mXtX = t(mX) %*% mX
  
  # define J, using vector of 1s
  v_1 =  matrix(1, 1, 90)
  J = diag(n)- ((1/n)*  t(v_1) %*% v_1)[[1]] 
  
  # find optimal intercept, and q
  optimal_mBeta_zero = (1/n) * v_1 %*% mY
  # use moore-penrose inverse here 
  optimal_q_hat = inv(diag(n) + (lambda * Ginv(mXXt))) %*% J %*% mY

  # from q's, get other beta's
  beta_from_q_hat = inv(mXtX) %*% t(mX) %*% optimal_q_hat 
  
  # define beta's together 
  mBeta_Intercept = rbind(optimal_mBeta_zero, beta_from_q_hat)
  
  # find predicted values, sums of squared errors
  yhat = cbind(1,mX) %*% mBeta_Intercept
  SSE = sum((mY - yhat)^2)
  
  # return results as list
  result = list(Beta = mBeta_Intercept,
                yhat = yhat,
                SSE = SSE,
                K = mXXt
    
  )
  
  return(result)
  
  

}

# first, get dsmle result
r_dsmle = dsmle_result = krr(mY, mX, kernel.type = "linear", lambda =1000)
r_dsmle$const
sum((r_dsmle$y - r_dsmle$yhat)^2)

# observe here that the two kernels are the same!
mXXt_dsmle = r_dsmle$K
mXXt_ours = mX_scaled %*% t(mX_scaled)

# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours = kernel_ridge(mX_scaled, mY, 1000)
r_ours$Beta
r_ours$SSE





