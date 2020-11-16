# Goal: Implementation of kernel ridge regression, comparision with dsmle package, including k-fold cv
# Author: Floris Holstege, Ruoying Dai, Ruben Eschauzier, Joyce Zhi
#

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

# get local package, dsmle, and add
install.packages("dsmle_1.0-4.tar.gz", repos = NULL, type="source")
library(dsmle)


create_kernel = function(type_kernel, mXXt, mX, param){
  
  ##create kernel based on type chosen
  if(type_kernel == "linear"){
    
    K = mXXt
    
  }else if(type_kernel == "nonhompolynom"){
    
    K = ((1 + mXXt)^param$degree)-1
    
  }else if(type_kernel == "RBF"){
    
    K = exp(-as.matrix(dist(mX)^2) * param$gamma)
    
  }
  return(K)
}
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
    mY_train <- mY[train,]
    mX_train <- mX[train,]
    mY_test <- mY[test,]
    mX_test <- mX[test,]
    
    # get result from training set
    result = kernel_ridge(mX_train,mY_train,type_kernel,param)

    # get beta's and predict on the test set
    # Likely something going wrong here - get different RMSE than dsmle
    mBeta = rbind(result$intercept, result$Beta)
    yhat_test = result$intercept  + cbind(1,mX_test) %*% mBeta
    
    # calculate RMSE for this fold
    error <- mY_test - yhat_test
    fold_RMSE <- sqrt(1/nrow(mX_test)) * (t(error)%*%error)

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
kernel_ridge <- function(mX, mY, type_kernel = "linear", param){
  
  # define n, p
  n = nrow(mX)
  p = ncol(mX)
  
  # initial XX^T, X^TX
  mXXt = mX %*% t(mX)
  mXtX = t(mX) %*% mX
  
  # create kernel
  K = create_kernel(type_kernel, mXXt, mX, param)
  
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
  beta_from_q_hat =  inv(mXtX) %*% t(mX) %*% optimal_q_hat 
  
  # predict Y-hat
  yhat = optimal_q_hat + optimal_mBeta_zero[[1]]
  
  ## adaptation to ensure correct Y-hat with RBF kernel 
  if (type_kernel != "linear"){
    yhat = optimal_q_hat
  }
  
  # get error and RMSE   
  error = mY - yhat
  RMSE = sqrt((1/n) * (t(error)%*%error))
  
  # return results as list
  result = list(Beta = beta_from_q_hat,
                intercept = optimal_mBeta_zero[[1]],
                qhat = optimal_q_hat,
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

# search for best lambda and other parameters, using grid search
hyperparam_search = function(mX, mY, k,type_kernel, vlambda, vdegree = NA, vgamma = NA){
  
  # define parameter grid based on type of kernel
  if(type_kernel == "linear"){
    
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

### section 1: data preparation

# get dataset 
setwd("KRR")
load("Airline.Rdata")

# set seed for reproduction
set.seed(123)

# define independent and dependent variables
mX = as.matrix(Airline[,-4])
mX_scaled = scale(mX)
mY = as.matrix(Airline[,4])

# define the parameters to try out in the gridsearch
vlambda = 10^seq(-5, 5, length.out = 50)
vdegree = seq(1,5,1)
vgamma = ncol(mX_scaled)^seq(-2,3,length.out =5)
k = 5

### section 2: finding the optimal parameters

# apply grid search for all the types of kernels
grid_result_linear <- hyperparam_search(mX_scaled, mY, k, type_kernel = "linear", vlambda = vlambda)
grid_result_nonhompolynom <- hyperparam_search(mX_scaled, mY, k, type_kernel = "nonhompolynom", vlambda = vlambda, vdegree = vdegree)
grid_result_RBF <- hyperparam_search(mX_scaled, mY, k, type_kernel = "RBF", vlambda = vlambda, vgamma = vgamma)

# get optimal parameters for each kernel
optimal_param_linear <- grid_result_linear[which.min(grid_result_linear$avg_RSME),]
optimal_param_nonhompolynom <- grid_result_nonhompolynom[which.min(grid_result_nonhompolynom$avg_RSME),]
optimal_param_RBF <- grid_result_RBF[which.min(grid_result_RBF$avg_RSME),]

# train the data with optimal parameters
result_linear = kernel_ridge(mX_scaled, mY, type_kernel = "linear", param = optimal_param_linear)
result_nonhompolynom = kernel_ridge(mX_scaled, mY, type_kernel = "nonhompolynom", param = optimal_param_nonhompolynom)
result_RBF = kernel_ridge(mX_scaled, mY, type_kernel = "RBF", param = optimal_param_RBF)


### section 3: applying optimal parameters, comparing our kernel_ridge function to krr()

## linear kernel
# for the linear kernel, we get a very similar result to the dsmle package (only a difference smaller than 9/10 decimals)
result_dsmle_linear = krr(y=mY, X= mX, kernel.type = "linear", lambda = optimal_param_linear$lambda)
result_dsmle_linear$yhat - result_linear$yhat

# for the non homogeneous polynomial kernel, exactly the same predicted values as the dsmle package (only a difference smaller 5 decimals)
result_dsmle_nonhompolynom = krr(y=mY, X= mX, kernel.type = "nonhompolynom", 
                                 lambda = optimal_param_nonhompolynom$lambda, 
                                 kernel.degree = optimal_param_nonhompolynom$degree)
result_dsmle_nonhompolynom$yhat - (result_nonhompolynom$yhat)

# result for the RBF - slightly different predicted values. This is likely due to a different kernel
result_dsmle_rbf =  krr(y=mY, X= mX, kernel.type = "RBF", 
                        lambda = optimal_param_RBF$lambda, 
                        kernel.RBF.sigma  =  2* optimal_param_RBF$gamma)

result_dsmle_rbf$yhat - (result_RBF$yhat)

# compare the kernels for RBF
## it is odd that the dsmle_rbf has the diagonal not at 1...
result_RBF$K
result_dsmle_rbf$K


## compare the kernels visually
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


# plots to show differences - sometimes the lines overlap exactly (linear, nonhompolynom)
df_linear = ggplot(data = df_compare_melted %>%
         filter(type %in% c("linear", "Actual")), # in this line, one can fill in "RBF" or "linear" instead of "non-homogeneous polynomial"
       aes(x = index, y = value, col = variable)) + 
  geom_line()+ 
  labs(col = "Model", y = "Output", title = "Linear Kernel") + 
  theme_classic()

# lines of DSMLE and ours overlap
df_linear

# plots to show differences - sometimes the lines overlap exactly (linear, nonhompolynom)
df_RBF = ggplot(data = df_compare_melted %>%
                     filter(type %in% c("RBF", "Actual")), # in this line, one can fill in "RBF" or "linear" instead of "non-homogeneous polynomial"
                   aes(x = index, y = value, col = variable)) + 
  geom_line()+ 
  labs(col = "Model", y = "Output", title = "RBF") + 
  theme_classic()

# comparison of RBF - it seems like the lambda is too low for the RBF, its severely overfitting for DSMLE_rbf
df_RBF

# plots to show differences - sometimes the lines overlap exactly (linear, nonhompolynom)
df_nonhompolynom = ggplot(data = df_compare_melted %>%
                  filter(type %in% c("non-homogeneous polynomial", "Actual")), # in this line, one can fill in "RBF" or "linear" instead of "non-homogeneous polynomial"
                aes(x = index, y = value, col = variable)) + 
  geom_line()+ 
  labs(col = "Model", y = "Output", title = "non-homogeneous polynomial") + 
  theme_classic()

# again, seems like there is something wrong with the optimal parameters, since overfitting
df_nonhompolynom


## Section 4: compare k-fold 

# check if dsmle finds same optimal parameters for linear kernel
cv_dsmle_linear = cv.krr(y=as.vector(mY), X= mX, kernel.type = "linear", lambda = vlambda, k.folds = k)

# dsmle finds different optimal lambda's and avg rmse

#compare lambda
optimal_param_linear$lambda
cv_dsmle_linear$lambda.min

# compare rmse
cv_dsmle_linear$rmse
grid_result_linear$avg_RSME 

## the mistake in cross validation likely leads to the parameters that in earlier plots show overfitting
# this is either due to 
# 1) a wrong calculation of rmse in our implementation
# 2) a wrong method for out of sample prediction (see comment in function)






