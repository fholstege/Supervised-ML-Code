# load libraries
library(MASS)
library(matlib)
library(caret)
library(glmnet)
library(tidyverse)
library(reshape2)
library(stargazer)


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
  typeNetInit = calc_typeNet(lambda, alpha, mD, p)
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

crossValidation = function (df, k, beta_init, lambda, alpha, tolerance, folds) {
  
  total_rmse = 0
  n = nrow(df)
  
  y = scale(as.matrix(df[,1]))
  X = scale(as.matrix(df[,-1]))
  
  df = scale(df)
  min_rsme = Inf 
  max_rsme = 0
  
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
  
  #Create k equally size folds
  folds = createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  
  
  for(i in 1:nrow(grid)){
    
    lambda = as.numeric(grid[i,][1])
    alpha = as.numeric(grid[i,][2])
    
    
    result_cv = crossValidation(df, k, beta_init, lambda, alpha, tolerance, folds)
    result_row <- c(lambda, 
                    alpha,
                    result_cv$avg_rmse,
                    result_cv$min_rsme,
                    result_cv$max_rsme)
    results[i,] <- result_row
    
  }
  
  return(results)
  
  
}


# set seed to ensure stability of results
set.seed(0)

# load the data
load("supermarket1996.RData")

# create dataframe with dependent and independent variables
df = subset(supermarket1996, select = -c(STORE, CITY, ZIP, GROCCOUP_sum, SHPINDX) )
y = scale(as.matrix(df[,1]))
X = scale(as.matrix(df[,-1]))

# define the initial set of beta's
Beta_init = as.matrix(runif(ncol(df)-1, min=-1, max=1))

# define the parameters 
tolerance = 0.0000000000001
epsilon = 0.0000000000001

# create grid of lambda and alpha combinations
listLambda <- 10^seq(-2, 10, length.out = 50)
listAlpha <- seq(0.0,1,0.1)
paramGrid <- expand.grid(listLambda, listAlpha)

# find the results of gridsearch
search_result <- HyperSearch(df, 10, paramGrid, Beta_init, tolerance)

# show best set of parameters
best_param <- search_result[search_result$avg_rmse==min(search_result$avg_rmse),]
best_param
best_lambda <- round(best_param$Lambda,2)
best_alpha <- best_param$Alpha

# Beta's for our estimate
BetaEst <- ElasticNetEst(X, y, Beta_init, lambda = best_lambda, alpha = best_alpha, tolerance, epsilon)
top_Beta <- data.frame(BetaEst) %>%
  filter(abs(BetaEst) > 0.01)%>%
  arrange(BetaEst)
stargazer(top_Beta,
          digits=2,
          summary = FALSE)



# Beta's for glm.net estimate, same param as our ideal
result.cv.ideal <- glmnet(scale(X), scale(y), alpha = best_alpha,
                              lambda =best_lambda, nfolds = 10)
glm.net_Beta <- as.matrix(result.cv.ideal$beta)
colnames(glm.net_Beta) <- c("BetaEst")

glm.net_top_Beta <- data.frame(glm.net_Beta) %>%
  filter(abs(BetaEst) > 0.01)%>%
  arrange(BetaEst)

stargazer(glm.net_top_Beta,
          digits=2,
          summary = FALSE)


# Beta's for glm.net estimate
result.cv.lambda <- cv.glmnet(scale(X), scale(y), alpha = 0,
                       lambda =listLambda, nfolds = 10)



# compare convergence
paramGrid_compare <- expand.grid(listLambda, best_alpha)

search_compare <- search_result[search_result$Alpha == best_alpha,]


df_compare = data.frame(lambda = log(listLambda), 
                        MM = search_compare$avg_rmse, 
                        GLM.NET = rev(result.cv.lambda$cvm))
df_compare = melt(df_compare, id.vars = 'lambda', variable.name = 'series') 



# plot compare to glm.net
ggplot() + 
  geom_point(data = df_compare, aes( x = lambda, y = value, col=series)) + 
  theme_minimal() +
  labs(y = "Average RMSE", 
       x = "Log Lambda",
       title = "At alpha = 0.2") + 
  theme_classic() + 
  theme(legend.title = element_blank())






####

# 
# dfBeta <- data.frame()
# 
# for(lambda in listLambda){
#   
#   y = scale(as.matrix(df[,1]))
#   x = scale(as.matrix(df[,-1]))
#   
#   BetaEst <- ElasticNetEst(x, y, Beta_init, lambda, 0.3, tolerance, epsilon)
#   
#   
#   dfBeta <- rbind(dfBeta, t(BetaEst))
#   ==
# }
# 
# dfBeta <- cbind(dfBeta, listLambda)
# 
# dfBetaViz <- melt(dfBeta ,  id.vars = 'listLambda', variable.name = 'series') %>%
#   filter(abs(value) >0.01)
# 
# ggplot(dfBetaViz, aes(log(listLambda),value)) + geom_line(aes(colour = series))
# 





#####


# 
# ggplot(data = gridSearch, aes( x = log(Lambda), y = avg_rmse)) + 
#   geom_point() + 
#   geom_ribbon(aes(ymin=min_rsme, ymax=max_rsme), alpha=0.4) + 
#   theme_minimal()
# 
# 
# ggplot() + 
#   geom_point(aes(x = log(result.cv$lambda), y = result.cv$cvm, col = 'Red')) +
#   geom_ribbon(aes(x = log(result.cv$lambda), ymin=result.cv$cvup, ymax=result.cv$cvlo), alpha=0.4) 
# 
# 
# plot(result.cv$lambda, result.cv$cvm, log = "x", col = "red", type = "p", pch = 20,
#      xlab = expression(lambda), ylab = "RMSE", las = 1)
# 
# cor(scale(X))
