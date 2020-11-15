

# Put here the packages required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(stargazer,
               lmtest,
               matlib, # install dependencies for dsmle
               SVMMaj,
               ISLR,
               plotrix,
               e1071, 
               ggplot2,
               reshape2)

# get local package, dsmle
install.packages("dsmle_1.0-4.tar.gz", repos = NULL, type="source")
library(dsmle)

# get dataset 
setwd("KRR")
load("Airline.Rdata")

# define independent and dependent variable
mX = as.matrix(Airline[,-4])
mX_scaled = scale(mX)
mY = as.matrix(Airline[,4])


# creates rbf kernel from independent variables
create_RBF_kernel = function(mXtX, gamma, n){
  
  # calculate distance matrix
  XX <- matrix(1, n) %*% diag(mXtX)
  D <- XX - 2*XtX + t(XX) # distance matrix
  
  # calculate rbf-kernel matrix
  K <- exp(gamma* D)
  
  return(K)
  
}

create_nonhompolynom_kernel = function(mXXt, degree){(1 +  mXXt)^degree}



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
kernel_ridge <- function(mX, mY, lambda, type_kernel = "Linear", degree =1, gamma=1){
  
  # define n, p
  n = nrow(mX)
  p = ncol(mX)
  
  # initial XX^T, X^TX
  mXXt = mX %*% t(mX)
  mXtX = t(mX) %*% mX
  
  if(type_kernel == "Linear"){
    
    # experiment with K
    K = mXXt
    KtK = t(K) %*% K
    
    
  }else if(type_kernel == "nonhompolynom"){
    
    K = create_nonhompolynom_kernel(mXXt,degree)
    
  }else if(type_kernel == "RBF"){
    
    K = create_RBF_kernel(mXtX, gamma, n)
  }
  
  
  # define J, using vector of 1s
  v_1 =  matrix(1, 1, 90)
  J = diag(n)- ((1/n)*  t(v_1) %*% v_1)[[1]] 
  D = diag(eigen(K)$values^-2)
  U = eigen(K)$vectors
  
  # find optimal intercept, and q
  optimal_mBeta_zero = (1/n) * v_1 %*% mY
  
  # use eigenvalues to quickly compute optimal Q
  optimal_q_hat = U %*% inv(diag(n) + lambda* D) %*% t(U) %*% (J %*% mY)

  # from q's, get other beta's
  beta_from_q_hat = inv(mXtX) %*% t(mX) %*% optimal_q_hat 
  
  # define beta's together 
  mBeta_Intercept = rbind(optimal_mBeta_zero, beta_from_q_hat)
  
  # find predicted values, sums of squared errors
  #yhat = cbind(1,mX) %*% mBeta_Intercept
  yhat = optimal_q_hat + optimal_mBeta_zero[[1]]
  SSE = sum((mY - yhat)^2)
  
  # return results as list
  result = list(Beta = mBeta_Intercept,
                yhat = yhat,
                SSE = SSE,
                K = K
    
  )
  
  return(result)
  
  

}


### linear

# first, get dsmle result
r_dsmle_linear = krr(y=mY, X= mX_scaled, kernel.type = "linear", lambda =1/1000, scale=FALSE, center=FALSE)


# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours_linear = kernel_ridge(mX_scaled, mY, lambda=1000)
r_ours$Beta
r_ours$SSE

# observe here that the two kernels are the same!
mXXt_dsmle = r_dsmle$K
mXXt_dsmle
mXXt_ours =r_ours_linear$K
mXXt_ours



# create df to compare predictions
df_compare = data.frame(our_pred = r_ours$yhat,
                        dsmle_pred = r_dsmle$yhat,
                        actual = Airline$output,
                        index = 1:length(r_ours$yhat))

df_compare = melt(df_compare, id.vars = "index")

# show in plot
ggplot(data = df_compare, aes(x = index, y = value, col = variable)) + 
  geom_line()



### RBF

# first, get dsmle result
r_dsmle_rbf = krr(y=mY, X= mX, kernel.type = "RBF", kernel.RBF.sigma = ncol(mX)/2, lambda = 1000)

# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours_rbf = kernel_ridge(mX_scaled, mY, lambda=1000, type_kernel = "RBF", gamma = 1/ncol(mX))


# observe here that if two kernels are the same!
mXXt_dsmle_rbf = r_dsmle_rbf$K
mXXt_dsmle_rbf
mXXt_ours_rbf = 
mXXt_ours


# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours_rbf = kernel_ridge(mX_scaled, mY, lambda=1000)
r_ours$Beta
r_ours$SSE

# create df to compare predictions
df_compare = data.frame(our_pred = r_ours$yhat,
                        dsmle_pred = r_dsmle$yhat,
                        actual = Airline$output,
                        index = 1:length(r_ours$yhat))

df_compare = melt(df_compare, id.vars = "index")

# show in plot
ggplot(data = df_compare, aes(x = index, y = value, col = variable)) + 
  geom_line()





