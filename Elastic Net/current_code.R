# first, we define several functions that we later use for the analysis

# Calculates diagonal matrix D
calc_mD = function(mBeta, p, epsilon){
  
  # Create the diagonal vector that will be filled with diagonal elements of D
  to_diagonalise <- vector(length=p)
  
  # Get diagonal elements of D (max of mBeta_i, epsilon)
  for (i in 1:p) {
    
    to_diagonalise[i] <- 1 / max(abs(mBeta[i]), epsilon)
  }
  # Multiply the vector containing diagonals with Identity matrix to get D
  return(to_diagonalise * diag(p))
}


# Calculates the loss function for the elastic
elasticLoss = function(mBeta, mA, mXtY, mYtY, alpha, lambda, n){
  
  # Calculate transposed matrix of Beta's
  transposed_mBeta = t(mBeta)
  
  # Compute constant
  constant <- (1/2*n) %*% mYtY + (1/2) * alpha * lambda *sum(abs(mBeta))
  
  # Return loss function
  return(1/2 * (transposed_mBeta%*%mA%*%mBeta)-(1/n)*(transposed_mBeta%*%mXtY) + constant)
}

# Calculates a often, used variable in subsequent computations λ(1 − α)I + λα.
# This indicates the division between the lasso (a = 1) and ridge method (a = 0)
calc_typeNet = function(lambda, alpha, mD,p){
  lambda *(1 - alpha)*diag(p) + (lambda * alpha  * mD)
  }

# Calculates the root mean squared error (RMSE)
calcRMSE = function(mX, mY, est_beta, n){
  error <- mY - mX %*% est_beta
  rsme <- sqrt((1/n) * (t(error)%*%error))
}

# calculates estimate for the elastic net, given lambda and alpha, using MM algorithm
ElasticNetEst = function(mX, mY, beta_init, lambda, alpha, tolerance, epsilon, max_iter = 100000){
  
  # Set iterations and improvement
  k <- 1
  improvement <- 0
  
  # Define number of predictor variables and datapoints
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Pre-compute constants
  mXtX <- crossprod(mX,mX)
  mXtY <- crossprod(mX,mY)
  mYtY <- crossprod(mY,mY)
  scaled_I <- lambda * (1-alpha) * diag(p)

  # get initial values for mD, mA, and Beta's
  mD <- calc_mD(beta_init, p, epsilon)
  typeNetInit <- calc_typeNet(lambda, alpha, mD, p)
  mA <- 1/n * mXtX + typeNetInit
  Beta_prev <- beta_init
  
  # start stepwise improvement of Beta's
  while (k == 1 | k < max_iter && (improvement > tolerance)) {
    
    # Increase number steps k
    k <- k + 1
    
    # calculate mD, MA
    mD <- calc_mD(Beta_prev, p, epsilon)
    typeNet <- calc_typeNet(lambda, alpha, mD, p)
    mA <- ((1/n) * mXtX) + typeNet
    
    # get new set of Beta's
    Beta_current <- solve(mA, 1/n * mXtY)
    
    # calculate loss function for previous, current Beta's - and the improvement
    loss_current <- elasticLoss(Beta_current,mA, mXtY, mYtY, alpha, lambda, n)
    loss_prev <- elasticLoss(Beta_prev, mA, mXtY, mYtY, alpha, lambda, n)
    improvement <- (loss_prev - loss_current)/loss_prev
    
    # set the previous beta's to current beta's for next step
    Beta_prev <- Beta_current
    
  }
  
  # return est. Beta's
  return(Beta_current)
}

# k-fold crossvalidation of the elastic net
crossValidation = function (df, k, beta_init, lambda, alpha, tolerance, folds) {
  
  # initial value for total rmse, min and max
  total_rmse <- 0
  min_rsme <- Inf 
  max_rsme <- 0
  
  # save the n of observations
  n <- nrow(df)

  #Perform k fold cross validation
  for(i in 1:length(folds)){
    
    #Split the data according to the folds
    test = df[folds[[i]],]
    train = df[-folds[[i]],]
    
    # define train and test set for y and x
    y_train <- as.matrix(train[,1])
    X_train <- as.matrix(train[,-1])
    y_test <- as.matrix(test[,1])
    X_test <- as.matrix(test[,-1])
    
    # get est. Beta's from the elastic net
    Beta_est <- ElasticNetEst(X_train, y_train, beta_init, lambda, alpha, tolerance, epsilon)
    
    # define rmse for this set of lambda, alpha
    rmse <- calcRMSE(X_test, y_test, as.matrix(Beta_est), nrow(X_test))
    
    # add current rms to total, to tlater take average
    total_rmse <- total_rmse + rmse
    
    # save min and max of rmse
    if(rmse > max_rsme){
      max_rsme = rmse
    }else if (rmse < min_rsme){
      min_rsme = rmse
    }
    
  }
  # calculate the avg. rmse across the folds
  avg_rmse <- total_rmse / length(folds)
  
  # returns results
  result = list(alpha = alpha,
                lambda = lambda,
                avg_rmse = avg_rmse,
                min_rsme = min_rsme,
                max_rsme = max_rsme
  )
  
  return(result)
  
}

# search the  hyperparameters lambda, alpha that minimize rmse in k-fold
HyperSearch = function(df, k, grid, beta_init, tolerance){
  
  # scale both the dependent and independent
  df <- scale(df)
  
  # empty dataframe 
  results <- data.frame(Lambda= numeric(), 
                        Alpha= numeric(),
                        avg_rmse = numeric(), 
                        min_rsme = numeric(),
                        max_rsme = numeric())
  
  # create k equally size folds
  folds = createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  
  # iterate over the grid
  for(i in 1:nrow(grid)){
    
    # get current lambda & alpha
    lambda <- as.numeric(grid[i,][1])
    alpha <- as.numeric(grid[i,][2])
    
    # get result of cross validation for lambda, alpha
    result_cv <- crossValidation(df, k, beta_init, lambda, alpha, tolerance, folds)
    
    # define row to add to dataframe with results
    result_row <- c(lambda, 
                    alpha,
                    result_cv$avg_rmse,
                    result_cv$min_rsme,
                    result_cv$max_rsme)
    results[i,] <- result_row
    
  }
  
  return(results)
  
  
}


# Here, we perform the analysis used in the report, using the functions above

# load libraries
library(MASS)
library(matlib)
library(caret)
library(glmnet)
library(tidyverse)
library(reshape2)
library(stargazer)

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

# set best set of parameters
best_param <- search_result[search_result$avg_rmse==min(search_result$avg_rmse),]
best_lambda <- round(best_param$Lambda,2)
best_alpha <- best_param$Alpha

# find Beta's for our estimate
BetaEst <- ElasticNetEst(X, y, Beta_init, lambda = best_lambda, alpha = best_alpha, tolerance, epsilon)

# select the top beta's in terms of absolute value
top_Beta <- data.frame(BetaEst) %>%
  filter(abs(BetaEst) > 0.01)%>%
  arrange(BetaEst)

# create overview table
stargazer(top_Beta,
          digits=2,
          summary = FALSE)

# Beta's for glm.net estimate, same param as from our grid search
result.cv.ideal <- glmnet(scale(X), scale(y), alpha = best_alpha,
                              lambda =best_lambda, nfolds = 10)
glm.net_Beta <- as.matrix(result.cv.ideal$beta)
colnames(glm.net_Beta) <- c("BetaEst")

# find top beta's for glm.net in terms of absolute value
glm.net_top_Beta <- data.frame(glm.net_Beta) %>%
  filter(abs(BetaEst) >= 0.01)%>%
  arrange(BetaEst)

# create overview table
stargazer(glm.net_top_Beta,
          digits=2,
          summary = FALSE)

# Beta's for glm.net estimate - for all values of lambda evaluated
result.cv.lambda <- cv.glmnet(scale(X), scale(y), alpha = 0,
                       lambda =listLambda, nfolds = 10)



# compare convergence of the two methods
# set variables for our method
paramGrid_compare <- expand.grid(listLambda, best_alpha)
search_compare <- search_result[search_result$Alpha == best_alpha,]


# create dataframe for the visualization
df_compare = data.frame(lambda = log(listLambda), 
                        MM = search_compare$avg_rmse, 
                        GLM.NET = rev(result.cv.lambda$cvm))
df_compare = melt(df_compare, id.vars = 'lambda', variable.name = 'series') 



# plot to compare convergence to glm.net
ggplot() + 
  geom_point(data = df_compare, 
             aes( x = lambda, y = value, col=series),
             size=4) + 
  theme_minimal() +
  labs(y = "Average RMSE", 
       x = "Log Lambda") + 
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("At alpha = 0.2") + 
  theme(
    plot.title = element_text(color="black", size=30)) + 
  theme_classic()+
  theme(legend.title = element_blank())




# create correlation matrix
corr <- cor(X)

#prepare to drop duplicates and correlations of 1     
corr[lower.tri(corr,diag=TRUE)] <- NA 
#drop perfect correlations
corr[corr == 1] <- NA 
#turn into a 3-column table
corr <- as.data.frame(as.table(corr))
#remove the NA values from above 
corr <- na.omit(corr) 
# only show high correlations
large_corr <- corr %>%
  filter(abs(Freq)>0.7) %>%
  distinct(Var1, .keep_all = TRUE) %>%
  arrange(-Freq)

# turn into table
stargazer(large_corr, summary = FALSE)
