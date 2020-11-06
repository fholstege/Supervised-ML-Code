# Author: Floris Holstege
# Purpose: implements the elastic net
# Date: 27/10/2020

require(MASS)
require(matlib)
require(caret)




calcLossterm <- function(lambda, alpha, D,p){lambda *(1 - alpha)*diag(p) + (lambda * alpha  * D) }

calcA <- function(mXtX,n, lossterm){ (1/n) * mXtX + lossterm}

calcLoss <- function(mX,mXtX, mBeta, mY,mYtY,lossterm, n, alpha){

  c = (1/2*n) %*% mYtY + (1/2) * alpha * sum(abs(mBeta))
  
  
  S = (t(mBeta) %*% t(mX) %*% mY)
  
  U = ((1/n) * mXtX %*% lossterm)
  V = (1/n)* S
  
  
  Z = t(mBeta) %*% U %*% mBeta
  

  result = (1/2)* Z - V + c
  
  return(result)
  
  
}



calcElasticNet <- function(mX,mY,lambda,alpha,e){
  
  
  # get number of observations
  n <- nrow(mX)
  p <- ncol(mX)
  
  
  # set the previous beta variable to initial, random beta's
  prevBeta <- runif(p, min=0, max=10)
  
  # calculate X'X, y'y
  mXtX <- t(mX) %*% mX
  mYtY <- t(mY) %*% mY
  

  # set initial stepscore to 0, k to 1. 
  StepScore <- 0
  k <- 1
  
  # run while, either if k is equal to 1, or the improvement between k-1th and kth set of beta's is smaller than the parameter e
  while (k == 1 | StepScore > e ){
    
    # step to next k
    k <- k + 1
    
    
    D <- 1/max(abs(prevBeta), e) * diag(p)
    lossterm <- calcLossterm(lambda,alpha, D,p)
    A <- calcA(mXtX, n, lossterm)
    

    currentBeta = (1/n) * inv(A) %*% t(mX) %*% mY
    

    prevScore <- calcLoss(mX, mXtX, prevBeta, mY, mYtY, lossterm, n, alpha)
    newscore <- calcLoss(mX, mXtX, currentBeta, mY, mYtY, lossterm, n, alpha)
    
    StepScore <- (prevScore - newscore)/prevScore
    
    prevBeta <- currentBeta
    
  }
  
  finalBeta <- prevBeta
  
  mYpredicted = mX %*% finalBeta
  
  error <- mY - mYpredicted
  
  RMSE <- (1/n) * (t(error) %*% error)
  
  
  result = list(alpha = alpha,
                lambda = lambda,
                RMSE = RMSE,
                Beta = finalBeta
    
  )
  
  return(result)

  
}



kfoldEval <- function(mX, mY,lambda, alpha, e, k){
  
  folds <- createFolds(mY, k =k, list = TRUE, returnTrain = FALSE)

  totalRSME = 0
  
  
  for(i in seq(1, k)){
    
    mXfold <- mX[-folds[[i]],]
    mYfold <- mY[-folds[[i]],]
    

    result <- calcElasticNet(mXfold, mYfold, lambda, alpha, e)
    
    totalRSME  = totalRSME + result$RMSE


  }
  AvgRSME = totalRSME / k
    
  
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
    print(results)
    
  }
 
  
}

load("supermarket1996.RData")
summary(supermarket1996)

mY = as.matrix(supermarket1996$GROCERY_sum)
mX = as.matrix(supermarket1996[,-5:-1])

listLambda <- 10^seq(-2, 10, length.out = 50)
listAlpha <- seq(0,1,0.1)

ParamCombinations <- expand.grid(listLambda, listAlpha)
ParamCombinations$Var1[1]

findHyperParam(mX, mY, e=0.000001, k=5, ParamCombinations)


results <- data.frame(Lambda= numeric(), 
                      Alpha= numeric(),
                      AvgRSME = numeric())
