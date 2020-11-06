# Author: Floris Holstege
# Purpose: implements the elastic net
# Date: 27/10/2020


load("supermarket1996.RData")
summary(supermarket1996)

mY = as.matrix(supermarket1996$GROCERY_sum)
mX = as.matrix(supermarket1996[,-5:-1])

library(MASS)
library(matlib)




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
    
    print(k)
    
    
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
  
  finalRSS <- prevScore
  
  mYpredicted = mX %*% finalBeta
  
  error <- mY - mYpredicted
  
  RMSE <- (1/n) * (t(error) %*% error)
  print(RMSE)
  
  
}

listLambda = 10^(seq(-2,100,1))
listLambda


for(i in seq(1,length(listLambda),1)){
  
  lambda = listLambda[i]
  calcElasticNet(mX, mY, lambda, 0.5,0.0001)


}



