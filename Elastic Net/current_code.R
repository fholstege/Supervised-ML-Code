# load libraries
library(MASS)
library(matlib)
library(caret)
library(glmnet)


# calc_mD: Calculates diagonal matrix D
# Parameters
#
#

calc_mD = function(mBeta, p, epsilon){
  
  # Create the diagonal vector that will be filled with diagonal elements of D
  to_diagonalise = vector(length=p)
  
  # get diagonal elements of D (max of mBeta_i, epsilon)
  for (i in 1:p) {
    
    to_diagonalise[i] = 1 / max(abs(mBeta[i]), epsilon)
  }
  
  # Multiply the vector containing diagonals with Identity matrix to get D
  return(to_diagonalise * diag(p))
}


elasticLoss = function(mBeta, mA, mXtY, mYtY, alpha, lambda, n){
  
  transposed_mBeta = t(mBeta)
  # Compute constant
  constant = (1/2*n) %*% mYtY + (1/2) * alpha * lambda *sum(abs(mBeta))
  
  return(1/2 * (transposed_mBeta%*%mA%*%mBeta)-(1/n)*(transposed_mBeta%*%mXtY) + constant)
}


# calcTypeNet: calculates a often, used variable in subsequent computations λ(1 − α)I + λα.
# This indicates the division between the lasso (a = 1) and ridge method (a = 0)
#  
# 

calc_typeNet = function(lambda, alpha, mD,p){
  lambda *(1 - alpha)*diag(p) + (lambda * alpha  * mD)
  }

calcRMSE = function(X, y, est_beta, n){
  error = y - X %*% est_beta
  rsme = sqrt((1/n) * (t(error)%*%error))
}


ElasticNetEst = function(mX, mY, beta_init, lambda, alpha, tolerance, epsilon, max_iter = 100000){
  
  # Set iterations and improvement
  k = 1
  improvement = 0
  
  # Define number of regressors and datapoints
  n = nrow(mX)
  p = ncol(mX)
  
  # Pre-compute constants
  mXtX = crossprod(mX,mX)
  mXtY = crossprod(mX,mY)
  mYtY = crossprod(mY,mY)
  scaled_I = lambda * (1-alpha) * diag(p)
  
  
  # get initial values for mD, mA
  mD = calc_mD(beta_init, p, epsilon)
  typeNetInit = calcTypeNet(lambda, alpha, mD, p)
  mA = 1/n * mXtX + typeNetInit
  
  
  Beta_prev = beta_init
  
  
  while (k == 1 | k < max_iter && (improvement > tolerance)) {
    
    # Increase number steps k
    k = k + 1
    
    # calculate mD, MA
    mD = calc_mD(Beta_prev, p, epsilon)
    typeNet = calc_typeNet(lambda, alpha, mD, p)
    mA = ((1/n) * mXtX) + typeNet
    
    # get new set of Beta's
    Beta_current = solve(mA, 1/n * mXtY)
    
    
    
    loss_current = elasticLoss(Beta_current,mA, mXtY, mYtY, alpha, lambda, n)
    
    loss_prev = elasticLoss(Beta_prev, mA, mXtY, mYtY, alpha, lambda, n)
    
    improvement = (loss_prev - loss_current)/loss_prev
    
    Beta_prev = Beta_current
    
  }
  return(Beta_current)
}

crossValidation = function (df, k, beta_init, lambda, alpha, tolerance) {
  
  total_rmse = 0
  n = nrow(df)
  
  y = scale(as.matrix(df[,1]))
  X = scale(as.matrix(df[,-1]))
  
  df = scale(df)
  min_rsme = Inf 
  max_rsme = 0
  
  #Create k equally size folds
  folds = createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  
  #Perform k fold cross validation
  for(i in 1:length(folds)){
    
    
    #Split the data according to the folds
    test = df[folds[[i]],]
    train = df[-folds[[i]],]
    
    y_train = as.matrix(train[,1])
    X_train = as.matrix(train[,-1])
    
    y_test = as.matrix(test[,1])
    X_test = as.matrix(test[,-1])
    
    Beta_est = ElasticNetEst(X_train, y_train, beta_init, lambda, alpha, tolerance, epsilon)
    
    rmse = calcRMSE(X_test, y_test, as.matrix(Beta_est), nrow(X_test))
    total_rmse = total_rmse + rmse
    
    if(rmse > max_rsme){
      max_rsme = rmse
    }else if (rmse < min_rsme){
      min_rsme = rmse
    }
    
  }
  
  
  avg_rmse = total_rmse / length(folds)
  
  # returns results
  result = list(alpha = alpha,
                lambda = lambda,
                avg_rmse = avg_rmse,
                min_rsme = min_rsme,
                max_rsme = max_rsme
  )
  
  return(result)
  
}

HyperSearch = function(df, k, grid, beta_init, tolerance){
  
  results <- data.frame(Lambda= numeric(), 
                        Alpha= numeric(),
                        avg_rmse = numeric(), 
                        min_rsme = numeric(),
                        max_rsme = numeric())
  
  
  for(i in 1:nrow(grid)){
    
    lambda = as.numeric(grid[i,][1])
    alpha = as.numeric(grid[i,][2])
    
    
    result_cv = crossValidation(df, k, beta_init, lambda, alpha, tolerance)
    result_row <- c(lambda, 
                    alpha,
                    result_cv$avg_rmse,
                    result_cv$min_rsme,
                    result_cv$max_rsme)
    results[i,] <- result_row
    
  }
  
  return(results)
  
  
}


# load the data
load("supermarket1996.RData")
df = subset(supermarket1996, select = -c(STORE, CITY, ZIP, GROCCOUP_sum, SHPINDX) )  


y = as.matrix(df[,1])
X = as.matrix(df[,-1])

Beta_init = as.matrix(runif(ncol(df)-1, min=-5, max=5))


tolerance = 1e-12
epsilon = 1e-12






listLambda <- 10^seq(-2, 10, length.out = 50)
listAlpha <- seq(0,1,0.1)
paramGrid <- expand.grid(listLambda, listAlpha)

gridSearch <- HyperSearch(df, 20, paramGrid, Beta_init, tolerance)
gridSearch_best <- 


# plot results per lambda
plot(gridSearch$Lambda, gridSearch$avg_rmse, log = "x", col = "red", type = "p", pch = 20,
     xlab = expression(lambda), ylab = "RMSE", las = 1)


ggplot() + 
  geom_point(data = gridSearch, aes( x = log(Lambda), y = avg_rmse, col = 'Red')) + 
  geom_point(aes(x = log(result.cv$lambda), y = result.cv$cvm, col = 'Blue')) +
  theme_minimal()



ggplot(data = gridSearch, aes( x = log(Lambda), y = avg_rmse)) + 
  geom_point() + 
  geom_ribbon(aes(ymin=min_rsme, ymax=max_rsme), alpha=0.4) + 
  theme_minimal()


ggplot() + 
  geom_point(aes(x = log(result.cv$lambda), y = result.cv$cvm, col = 'Red')) +
  geom_ribbon(aes(x = log(result.cv$lambda), ymin=result.cv$cvup, ymax=result.cv$cvlo), alpha=0.4) 


result.cv <- cv.glmnet(scale(X), scale(y), alpha = 0.5,
                       lambda = 10^seq(-2, 10, length.out = 50), nfolds = 10)

plot(result.cv$lambda, result.cv$cvm, log = "x", col = "red", type = "p", pch = 20,
     xlab = expression(lambda), ylab = "RMSE", las = 1)
