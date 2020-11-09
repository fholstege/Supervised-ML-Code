# Author: Floris Holstege
# Purpose: implements the elastic net
# Date: 27/10/2020

library(MASS)
library(matlib)
library(caret)
library(glmnet)


calcD = function(prevBeta, p, e){
  
  # Create the diagonal vector that will create diagonal matrix D
  Diagonalise = vector(length=p)

  for (i in 1:p) {
      Diagonalise[i] <- 1 / max(abs(prevBeta[i]), e)
  }
  # Multiply the vector containing diagonals with Identity matrix to get D
  return(Diagonalise * diag(p))
}


# calcTypeNet: calculates a often, used variable in subsequent computations λ(1 − α)I + λα.
# This indicates the division between the lasso (a = 1) and ridge method (a = 0)
#  
# 

calcTypeNet <- function(lambda, alpha, D,p){lambda *(1 - alpha)*diag(p) + (lambda * alpha  * D) }

#calcA: calculates the term A -which, which multiplies XtX with the TypeNet term  inv(n) * t(X)X + λ(1 − α)I + λα
#
#

calcA <- function(mXtX,n, TypeNet){ (1/n) * mXtX + TypeNet}

# calcLoss: calculates the result of the loss function for a set of Beta's
#
#

calcLoss <- function(mBeta, mYtY, n, alpha, lambda, mA){

  c = (((2*n)^-1) %*% mYtY) + ((1/2) * lambda *alpha * sum(abs(mBeta)))
  mXtY <- t(mX) %*% mY
  
  lossResult <- ((1/2) * (t(mBeta) %*% mA %*% mBeta)) - ((1/n) * t(mBeta) %*% mXtY) + c

  return(lossResult)
  
  
}


# calcRSME

calcRSME <- function(mX, mY, finalBeta, n){
  
  mYpredicted = mX %*% finalBeta
  
  error <- mY - mYpredicted
  
  MSE <- (1/n) * (t(error) %*% error)
  
  return(sqrt(RMSE))
  
  
  
}
  
  



# calcElasticNet: calculates the Elastic net results for a given lambda, alpha


calcElasticNet <- function(mX,mY,lambda,alpha,e){
  
  
  # get number of observations
  n <- nrow(mX)
  p <- ncol(mX)
  
  
  # set the previous beta variable to initial, random beta's
  prevBeta <- runif(p, min=-1000, max=1000)
  
  print("Initial beta's")
  print(prevBeta)
  
  # calculate X'X, y'y
  mXtX <- t(mX) %*% mX
  mYtY <- t(mY) %*% mY
  mXtY <- t(mX) %*% mY
  
  # set initial stepscore to 0, k to 1. 
  StepScore <- 0
  k <- 1
  
  # run while, either if k is equal to 1, or the improvement between k-1th and kth set of beta's is smaller than the parameter e
  while (k == 1 | StepScore > e ){
    


    # step to next k
    k <- k + 1
    
    mD <- calcD(prevBeta, p, e)

    
    TypeNet <- calcTypeNet(lambda,alpha, mD,p)
    mA <- calcA(mXtX, n, TypeNet)


    print("WE ARE HERE")
    currentBeta =  solve(mA, (1/n) * mXtY)
    
    prevScore <- calcLoss(prevBeta, mYtY, n, alpha, lambda, mA)
    newScore <- calcLoss(currentBeta, mYtY, n, alpha, lambda, mA)

    StepScore <- (prevScore - newScore)/prevScore
    print("STEPSCORE: ")
    print(StepScore)
    
    prevBeta <- currentBeta
    
  }
  
  RMSE <- calcRSME(mX, mY, prevBeta, n)
  

  dfBeta <- data.frame(prevBeta)
  rownames(dfBeta) <- colnames(mX)
  
  
  result = list(alpha = alpha,
                lambda = lambda,
                RMSE = RMSE,
                Beta = dfBeta
    
  )
  
  return(result)

  
}


# kfoldEval: evaluate the hyperparameters using k-fold cross validation 


kfoldEval <- function(mX, mY,lambda, alpha, e, k){
  
  # create k folds
  folds <- createFolds(mY, k =k, list = TRUE, returnTrain = FALSE)
  totalRSME = 0
  
  
  # test the model on k-1 folds, k iterations
  for(i in seq(1, k)){
    
    mXfold <- mX[-folds[[i]],]
    mYfold <- mY[-folds[[i]],]
    
    mXtest <- mX[folds[[i]],]
    mYtest <- mY[folds[[i]],]
    
    result <- calcElasticNet(mXfold, mYfold, lambda, alpha, e)
    
    BetaFold <- result$Beta
    
    foldRSME <- calcRSME(mXtest, mYtest, BetaFold, nrow(mX))
    
    
    totalRSME  = totalRSME + foldRSME


  }
  
  # average RSME of all the tests 
  AvgRSME = totalRSME / k
    
  
  # returns results
  result = list(alpha = alpha,
                lambda = lambda,
                AvgRSME = AvgRSME
  )
  
  return(result)
}


findHyperParam <- function(mX, mY, e, k, ParamCombinations){
  
  results <- data.frame(Lambda= numeric(), 
                        Alpha= numeric(),
                        AvgRSME = numeric())
  

  for(i in seq(1, length(ParamCombinations[[1]]))){
    
    print(paste0("Iteration: ", i))
    print(paste0("Lambda: ", ParamCombinations$Var1[i]))
    print(paste0("alpha: ", ParamCombinations$Var2[i]))
    
    resultKfold <- kfoldEval(mX, mY, ParamCombinations$Var1[i],ParamCombinations$Var2[i], e, k)
    
    print(paste0("Avg RMSE: ", resultKfold$AvgRSME))
    
    
    resultRow <- c(ParamCombinations$Var1[i], ParamCombinations$Var2[i],resultKfold$AvgRSME)
    results[i,] <- resultRow

  }
  
  return(results)
 
  
}


getwd()
setwd("Elastic Net")

# load supermarket data
load("supermarket1996.RData")

# pick dependent and indpendent variable
mY = as.matrix(supermarket1996$GROCERY_sum)
mX = as.matrix(supermarket1996[,-5:-1])

# scale the independent
mXscaled <- scale(mX)

# list of all possible lambda's, and all possible alpha's, create grid of these
listLambda <- 10^seq(-2, 10, length.out = 10)
listAlpha <- 0
ParamCombinations <- expand.grid(listLambda, listAlpha)

r <- calcElasticNet(mXscaled, scale(mY), 1000, alpha = 0, e=0.000000000001)
r$Beta

result <- glmnet(mXscaled, mY, alpha = 0, lambda = 1000,
                 standardize = FALSE)
result$beta


glm_RMSE <- calcRSME(mXscaled, mY, result$beta, nrow(mXscaled))

r$RMSE[[1]] - glm_RMSE
