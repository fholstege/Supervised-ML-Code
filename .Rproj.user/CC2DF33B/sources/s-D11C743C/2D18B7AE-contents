

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

# set seed for reproduction
set.seed(123)

# define independent and dependent variable
mX = as.matrix(Airline[,-4])
mX_scaled = scale(mX)
mY = as.matrix(Airline[,4])

# creates rbf kernel from independent variables
create_RBF_kernel = function(mX, gamma){exp(-as.matrix(dist(mX)^2) * gamma)}

# creates non homogeneous polynomial kernel
create_nonhompolynom_kernel = function(mXXt, degree){((1 + mXXt)^degree)-1}

# implements k-fold cross validation for kernel ridge regression
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
    error <- y_test - (cbind(1,X_test) %*% mBeta)
    fold_RMSE <- sqrt(1/nrow(X_test)) * (t(error)%*%error)

    # add to total RMSE to later average out
    total_RMSE = total_RMSE + fold_RMSE
    
  }
  
  # compute average RMSE
  avg_RMSE = total_RMSE/length(folds)
  
  # define results, including parameters
  result = list(avg_RMSE = avg_RMSE,
                param = param)
  
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
  
  ## create kernel based on type chosen
  if(type_kernel == "Linear"){
  
    K = mXXt
    KtK = t(K) %*% K
    
  }else if(type_kernel == "nonhompolynom"){
    
    K = create_nonhompolynom_kernel(mXXt,param$degree)
    
  }else if(type_kernel == "RBF"){
    
    K = create_RBF_kernel(mX, param$gamma)
    
  }
  
  
  # define J, using vector of 1s
  v_1 =  matrix(1, 1, n)
  J = diag(n)- ((1/n)*  t(v_1) %*% v_1)[[1]] 
  
  # define U and d
  D = diag(1/eigen(K)$values)
  U = eigen(K)$vectors
  
  # find optimal intercept
  optimal_mBeta_zero = (1/n) * v_1 %*% mY
  
  # define shrunk matrix first, then find optimal intercept
  mShrunk = eigen(K)$values / (eigen(K)$values + outer(rep(1,n), param$lambda))
  optimal_q_hat = U %*% (mShrunk * as.vector(t(U) %*% mY))

  # from q's, get other beta's
  beta_from_q_hat = inv(mXtX) %*% t(mX) %*% optimal_q_hat 
  
  # add Beta's together with constant 
  mBeta_Intercept = rbind(optimal_mBeta_zero, beta_from_q_hat)
  
  # predict Y-hat
  yhat = optimal_q_hat + optimal_mBeta_zero[[1]]
  
  
  ## to check - is this appropriate?
  if (type_kernel == "RBF"){
    yhat = optimal_q_hat
  }
  
  # get error and RMSE   
  error = mY - yhat
  RMSE = sqrt((1/n) * (t(error)%*%error))
  
  # return results as list
  result = list(Beta = mBeta_Intercept,
                yhat = yhat,
                X = mX,
                y = mY,
                RMSE = RMSE,
                K = K,
                lambda = param$lambda,
                param = param
  )
  
  return(result)
  
  

}

# search for best lambda and other parameters 
hyperparam_search = function(mX, mY, k,type_kernel, vlambda, vdegree = NA, vgamma = NA){
  
  # define parameters based on type of kernel
  if(type_kernel == "Linear"){
    
    paramGrid <- expand.grid(vlambda)
    colnames(paramGrid) <- c("lambda")
    
  }else if(type_kernel == "RBF"){
    
    paramGrid <- expand.grid(vlambda, vgamma)
    colnames(paramGrid) <- c("lambda", "gamma")
    
  }else if(type_kernel == "nonhompolynom"){

    paramGrid <- expand.grid(vlambda, vdegree)
    colnames(paramGrid) <- c("lambda", "degree")
  }
  
  # create k equally size folds
  folds = createFolds(mY, k = k, list = TRUE, returnTrain = FALSE)
  
  # empty column, later to be filled
  paramGrid$avg_RSME <- NA

  # iterate over the grid
  for(i in 1:nrow(paramGrid)){
    
    # select parameters from the grid
    param <- paramGrid[i,]

    # test these with k-fold
    cv_result <- cv_kernel_ridge(mX, mY, type_kernel = type_kernel, folds = folds, param = param)

    #save the result
    paramGrid$avg_RSME[i] <- cv_result$avg_RMSE

    
  }
  
  return(paramGrid)
}


# try out these different parameters
vlambda = 10^seq(-5, 5, length.out = 50)
vdegree = seq(1,5,1)
vgamma = ncol(mX_scaled)^seq(-2,3,length.out =5)
k = 5

# grid search
grid_result_linear <- hyperparam_search(mX_scaled, mY, k, type_kernel = "Linear", vlambda = vlambda)
grid_result_nonhompolynom <- hyperparam_search(mX_scaled, mY, k, type_kernel = "nonhompolynom", vlambda = vlambda, vdegree = vdegree)
grid_result_RBF <- hyperparam_search(mX_scaled, mY, k, type_kernel = "RBF", vlambda = vlambda, vgamma = vgamma)

# get optimal parameters
optimal_param_linear <- grid_result_linear[which.min(grid_result_linear$avg_RSME),]
optimal_param_nonhompolynom <- grid_result_nonhompolynom[which.min(grid_result_nonhompolynom$avg_RSME),]
optimal_param_RBF <- grid_result_RBF[which.min(grid_result_RBF$avg_RSME),]

# check if dsmle finds same optimal parameters
cv_dsmle_linear = cv.krr(y=as.vector(mY), X= mX, kernel.type = "linear", lambda = list_lambda, k.folds = k)
cv_dsmle_linear$lambda.min
cv_dsmle_linear$rmse
grid_result_linear$avg_RSME

cv.krr


# compare for optimal to dsmle
result_linear = kernel_ridge(mX_scaled, mY, type_kernel = "Linear", param = test_param)
result_dsmle_linear = krr(y=mY, X= mX, kernel.type = "linear", lambda = test_param$lambda)
predict(result_dsmle_linear, mX)




error = result_dsmle_linear$yhat - result_dsmle_linear$y
RMSE = sqrt((1/nrow(mX)) * (t(error)%*%error))



# compare for optimal to dsmle
result_linear = kernel_ridge(mX_scaled, mY, type_kernel = "Linear", param = optimal_param_linear)
result_nonhompolynom = kernel_ridge(mX_scaled, mY, type_kernel = "nonhompolynom", param = optimal_param_nonhompolynom)
result_RBF = kernel_ridge(mX_scaled, mY, type_kernel = "RBF", param = optimal_param_RBF)

# for the linear kernel, exactly the same predicted values as the dsmle package (only a difference of 15 decimals)
result_dsmle_linear = krr(y=mY, X= mX, kernel.type = "linear", lambda = optimal_param_linear$lambda)
result_dsmle_linear$yhat - result_linear$yhat

# for the non homogeneous polynomial kernel, exactly the same predicted values as the dsmle package (only a difference of 15 decimals)
result_dsmle_nonhompolynom = krr(y=mY, X= mX, kernel.type = "nonhompolynom", 
                                 lambda = optimal_param_nonhompolynom$lambda, 
                                 kernel.degree = optimal_param_nonhompolynom$degree)
result_dsmle_nonhompolynom$yhat - (result_nonhompolynom$yhat)

# result for the RBF - different, due to different kernel. Small difference
result_dsmle_rbf =  krr(y=mY, X= mX, kernel.type = "RBF", 
                        lambda = optimal_param_RBF$lambda, 
                        kernel.RBF.sigma  =  2* optimal_param_RBF$gamma)

result_dsmle_rbf$yhat - (result_RBF$yhat)

# compare the kernels
result_linear$RMSE
result_nonhompolynom$RMSE
result_RBF$RMSE

# create dataframe for plot
df_compare <- data.frame(index = 1:length(result_linear$yhat),
                         rbf_yhat = result_RBF$yhat,
                         nonhompolynom_yhat = result_nonhompolynom$yhat,
                         linear_yhat = result_linear$yhat, 
                         dsmle_rbf_yhat = result_dsmle_rbf$yhat,
                         dsmle_nonhompolynom_yhat = result_dsmle_nonhompolynom$yhat,
                         dsmle_linear_yhat = result_dsmle_linear$yhat,
                         actual = mY)

df_compare_melted = melt(df_compare, id.vars = "index") %>%
  mutate(type = ifelse(str_detect(variable, "rbf"), "RBF",
                       ifelse(str_detect(variable, "linear"),"linear",
                              ifelse(str_detect(variable, "nonhompolynom"), "non-homogeneous polynomial", "Actual"))))


# plots to show differences
ggplot(data = df_compare_melted %>%
         filter(type %in% c("non-homogeneous polynomial", "Actual")), 
       aes(x = index, y = value, col = variable)) + 
  geom_line()
