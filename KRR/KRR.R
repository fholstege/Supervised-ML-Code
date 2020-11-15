

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
               reshape2,
               tidyverse)

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
create_RBF_kernel = function(mX, gamma){exp(-as.matrix(dist(mX)^2) * gamma)}

# non homogeneous polynomial kernel
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

calcRMSE = function(mX, mY, est_beta, n){
 
  error <- mY - mX %*% est_beta
  rsme <- sqrt((1/n) * (t(error)%*%error))
}

cv_kernel_ridge <- function(mX, mY, lambda, type_kernel,folds,param){
  
  # initial value for total rmse, min and max
  total_RMSE <- 0
  
  #Perform k fold cross validation
  for(i in 1:length(folds)){
    
    #Split the data according to the folds
    test = folds[[i]]
    train = -folds[[i]]
    
    # define train and test set for y and x
    y_train <- mY[train,]
    X_train <- mX[train,]
    y_test <- mY[test,]
    X_test <- mX[test,]
    
    # get Beta's from training set
    result = kernel_ridge(X_train,y_train,type_kernel,param)
    mBeta = result$Beta
    
    # use the beta's from training set to test on test set
    error <- y_test - cbind(1,X_test) %*% mBeta
    fold_RMSE <- sqrt((1/nrow(X_test)) * (t(error)%*%error))
    
    total_RMSE = total_RMSE + fold_RMSE
    
  }
  
  avg_RMSE = total_RMSE/length(folds)
  
  other_hyperparam = NA
  
  if(type_kernel == "RBF"){
    other_hyperparam = gamma
  }else if (type_kernel == "nonhomopolynom"){
    other_hyperparam = degree
  }
  
  # returns result
  result = list(avg_RMSE = avg_RMSE,
                param = param
                
  )
  
  return(result)
  
  
}

# executes kernel ridge regression
kernel_ridge <- function(mX, mY, type_kernel = "Linear", param){
  
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
    
    
  }else if(type_kernel == "nonhomopolynom"){
    
    K = create_nonhompolynom_kernel(mXXt,param$degree)
    
  }else if(type_kernel == "RBF"){
    
    K = create_RBF_kernel(mX, param$gamma)
    
  }
  
  
  # define J, using vector of 1s
  v_1 =  matrix(1, 1, n)
  J = diag(n)- ((1/n)*  t(v_1) %*% v_1)[[1]] 
  #D = diag(eigen(K)$values^-2)
  D = diag(1/eigen(K)$values)
 
  U = eigen(K)$vectors
  
  # find optimal intercept, and q
  optimal_mBeta_zero = (1/n) * v_1 %*% mY
  
  #lEig <- eigen(mK, symmetric = TRUE)
  #mD <- diag(1 / lEig$values)
  #mU <- lEig$vectors
  mShrunk = eigen(K)$values / (eigen(K)$values + outer(rep(1,n), param$lambda))
  optimal_q_hat = U %*% (mShrunk * as.vector(t(U) %*% mY))
  
  #mShrunk <- lEig$values / (lEig$values + outer(rep(1,iN), dLambda))  
  #vQTilde <- lEig$vectors %*% (mShrunk * as.vector(t(mU) %*% vY))
  
  # use eigenvalues to quickly compute optimal Q
  #optimal_q_hat = U %*% Ginv(diag(n) + lambda* D) %*% t(U) %*% (J %*% mY)

  # from q's, get other beta's
  beta_from_q_hat = inv(mXtX) %*% t(mX) %*% optimal_q_hat 
  
  # define beta's together 
  mBeta_Intercept = rbind(optimal_mBeta_zero, beta_from_q_hat)
  
  # find predicted values, sums of squared errors
  yhat = optimal_q_hat + optimal_mBeta_zero[[1]]
  error = mY - yhat
  SSE = sum(t(error)%*%error)
  RMSE = sqrt((1/n) * (t(error)%*%error))
  
  # return results as list
  result = list(Beta = mBeta_Intercept,
                yhat = yhat,
                SSE = SSE,
                RMSE = RMSE,
                K = K

  )
  
  return(result)
  
  

}


hyperparam_search = function(mX, mY, k,type_kernel = "Linear", vlambda, vdegree = NA, vgamma = NA){
  
  if(type_kernel == "Linear"){
    
    paramGrid <- expand.grid(vlambda)
    colnames(paramGrid) <- c("lambda")
    
  }else if(type_kernel == "RBF"){
    
    paramGrid <- expand.grid(vlambda, vgamma)
    colnames(paramGrid) <- c("lambda", "gamma")
    
  }else if(type_kernel == "nonhomopolynom"){
    
    paramGrid <- expand.grid(vlambda, vdegree)
    colnames(paramGrid) <- c("lambda", "degree")
  }
  

  # create k equally size folds
  folds = createFolds(mY, k = k, list = TRUE, returnTrain = FALSE)
  
  paramGrid$avg_RMSE = NA

  # iterate over the grid
  for(i in 1:nrow(paramGrid)){
    
    param <- paramGrid[i,]

    cv_result <- cv_kernel_ridge(mX, mY, type_kernel = type_kernel, folds = folds, param = param)

    paramGrid$avg_RSME[i] <- cv_result$avg_RMSE

    
  }
  
  return(paramGrid)
}


lambda = 10^seq(-5, 5, length.out = 100)
degree = seq(1,5,1)
gamma = ncol(mX_scaled)^seq(-2,3,length.out =5)


# grid search
grid_result_linear <- hyperparam_search(mX_scaled, mY, 5, type_kernel = "Linear", vlambda = lambda)
grid_result_nonhomopolynom <- hyperparam_search(mX_scaled, mY, 5, type_kernel = "nonhomopolynom", vlambda = lambda, vdegree = degree)
grid_result_RBF <- hyperparam_search(mX_scaled, mY, 5, type_kernel = "RBF", vlambda = lambda, vgamma = gamma)

# get optimal parameters
optimal_param_linear <- grid_result_linear[which.min(grid_result_linear$avg_resume),]
optimal_param_nonhomopolynom <- grid_result_nonhomopolynom[which.min(grid_result_nonhomopolynom$avg_resume),]
optimal_param_RBF <- grid_result_RBF[which.min(grid_result_RBF$avg_resume),]












# old code, need to rewrite
### linear

# first, get dsmle result
r_dsmle_linear = krr(y=mY, X= mX, kernel.type = "linear", lambda =1000)


# our result
r_ours_linear = kernel_ridge(mX_scaled, mY, lambda=1000)

# observe here that the two kernels are the same!
mXXt_dsmle = r_dsmle_linear$K
mXXt_dsmle
mXXt_ours =r_ours_linear$K
mXXt_ours

# create df to compare predictions
df_compare = data.frame(our_pred = r_ours_linear$yhat,
                        dsmle_pred = r_dsmle_linear$yhat,
                        actual = Airline$output,
                        index = 1:length(r_ours_linear$yhat))

df_compare = melt(df_compare, id.vars = "index")

# show in plot
ggplot(data = df_compare, aes(x = index, y = value, col = variable)) + 
  geom_line()



### RBF

help(krr)

# first, get dsmle result
r_dsmle_rbf = krr(y=mY, X= mX, kernel.type = "RBF", kernel.RBF.sigma = ncol(mX)/2, lambda = 1000)

# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours_rbf = kernel_ridge(mX_scaled, mY, lambda=1000, type_kernel = "RBF", gamma = 1/ncol(mX))

# The two kernels are different! maybe due to difference in sigma/gamma
mXXt_dsmle_rbf = r_dsmle_rbf$K
mXXt_ours_rbf = r_ours_rbf$K


# create df to compare predictions
df_compare_rbf = data.frame(our_pred = r_ours_rbf$yhat,
                        dsmle_pred = r_dsmle_rbf$yhat,
                        actual = Airline$output,
                        index = 1:length(r_ours_rbf$yhat))

df_compare_rbf = melt(df_compare_rbf, id.vars = "index")

# show in plot
ggplot(data = df_compare_rbf, aes(x = index, y = value, col = variable)) + 
  geom_line()


### non-homogeneous polynomial


# first, get dsmle result
r_dsmle_nonhomopolynom = krr(y=mY, X= mX, kernel.type = "nonhompolynom", kernel.degree =2, lambda = 1/1000)

# but our sse is much lower...different predicted values, most likely due to different beta's
r_ours_nonhomopolynom = kernel_ridge(mX_scaled, mY, lambda=1000, type_kernel = "nonhompolynom", degree = 2)

# The two kernels are different! maybe due to difference in sigma/gamma
mXXt_dsmle_nonhomopolynom= r_dsmle_nonhomopolynom$K
mXXt_ours_nonhomopolynom = r_ours_nonhomopolynom$K

# create df to compare predictions
df_compare_nonhomopolynom = data.frame(our_pred = r_ours_nonhomopolynom$yhat,
                            dsmle_pred = r_dsmle_nonhomopolynom$yhat,
                            actual = Airline$output,
                            index = 1:length(r_ours_rbf$yhat))

df_compare_nonhomopolynom = melt(df_compare_nonhomopolynom, id.vars = "index")

# show in plot
ggplot(data = df_compare_nonhomopolynom, aes(x = index, y = value, col = variable)) + 
  geom_line()



