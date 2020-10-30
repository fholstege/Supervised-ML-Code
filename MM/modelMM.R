
# Author: Floris Holstege
# Purpose: implements the maximize by minorization algorithm for a linear model
# Date: 27/10/2020



# calcRSS
# Calculates the residual squared errors for a multiple regression of the form Y = XBeta + e
# 
# Parameters: 
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mY: Matrix of n x 1 dependent variables (n = observations)
#   mBeta: Column Matrix of p x 1 with coefficients of said model
#   
# Output:
#   ESquared: float, residual squared errors
# 

calcRSS <- function(mX, mY,mBeta){

  # calculate the errors
  mE <-  mY - mX %*% mBeta
  
  # get errors squared
  ESquared <- t(mE) %*% mE
  
  # return the residual sum of squared errors
  return(ESquared)
  
}

# calcLargestEigen
# Calculates the largest eigenvalue of an matrix of independent variables
# 
# Parameters: 
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   
# Output:
#   LargestEigenval: float, largest eigenvalue of said matrix
#   

calcLargestEigen <- function(X){
  
  # get the eigenvalues of X 
  EigenValX <- eigen(X)$values
  
  # from these eigenvalues, get the largest one
  LargestEigenVal <- max(EigenValX, na.rm = TRUE)
  
  return(LargestEigenVal)
  
}

#CalcBetaK
# Calculates the kth beta in an MM algorithm
# 
# Parameters:
#   prevBeta: double, k-1th beta
#   Lambda: double, largest eigenvalue of X (independent variables) squared
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mY: Matrix of n x 1 dependent variables (n = observations)
# 
# Output: 
#   BetaK; matrix of new set of coefficients

calcBetaK <- function(prevBeta,Lambda, mX, mY){
  
  # get matrix of squared X
  XSquared <- t(mX) %*% mX
  
  # calculate the Kth 
  BetaK = prevBeta - ((1/Lambda) *  XSquared %*% prevBeta) + ((1/Lambda) * t(mX) %*% mY )
  
  return(BetaK)
  
  
  
}

# CalcStepScore
# Calculates the % improvement between the k-1th and kth set of beta's
# 
# Parameters:
#   prevBeta: double, k-1th beta
#   currbeta: double, kth beta
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
# 
# Output: 
#   StepScore; double, % improvement between the RSS of the two sets of beta's

calcStepScore <- function(X,Y, prevBeta, currBeta){
  
  # difference in RSS between previous and current set of beta's
  diffRSS <- (calcRSS(X,Y,prevBeta) - calcRSS(X,Y,currBeta))
  
  # divide difference with previous score to get % change
  StepScore <- diffRSS /calcRSS(X,Y,prevBeta)
  
  return(StepScore)
  
}

# getB0
# Set Beta's, randomly, between 0 and 1
# 
# Parameters:
#   mX: Matrix of n x p (n = observations, p = independent variables)
# 
# Output: 
#   B0: matrix, set of Beta's between 0 and 1

getB0 <- function(mX){
  
  # determine the number of independent variables, generate as many random beta's
  nIndVar = ncol(mX)
  Beta0 <- runif(nIndVar, min=0, max=1)
  
  # turn to matrix format and returns
  return(as.matrix(Beta0))
  
}

# calcYest
# Calculates the predicted Y, based on the X and est. Beta's of a linear model
# 
# Parameters:
#   X: matrix of n x p (n = observations, p = independent variables)
#   BetaEst: Matrix, Estimated Beta's, px1 column, 
#
# Output:
#   Yestdf: matrix, predicted Y
#

calcYest <- function(mX,mBetaEst){
  
  
  # multiply X with Beta (est.) to get predicted Y
  Yest <- mX %*% mBetaEst
  
  # turn into dataframe
  Yestdf <- as.data.frame(Yest)
  colnames(Yestdf) <- c("Yest")
  
  return(as.matrix(Yestdf))
  
}

# calcRsquared
# Calculates the r-squared
#
# Parameters:
#   Y: matrix, the true dependent variable   
#   Yest: matrix, the predicted dependent variable
#   (optional) adjusted: if True, return adjusted r squared
#   (optional) p: if adjusted is calculated, add number of variables
# 
# Output:
#   Rsquared: double, the Rsquared or adjusted Rsquared for a linear model

calcRsquared <- function(mY, mYest, adjusted = FALSE, p=0){
  
  # standardize Y, and Yest (mean of 0)
  mStandY = mY - mean(mY)
  mStandYest = mYest - mean(mYest)
  
  # calculate Rsquared
  numerator <- (t(mStandY) %*% mStandYest)^2
  denominator <- (t(mStandY) %*% mY) %*% (t(mStandYest) %*% mStandYest)
  resultRsquared <- (numerator/denominator)
  
  if(adjusted){
    
    n <- nrow(mY)
    
    adjRsquared = 1 - (((1-resultRsquared)*(n - 1))/(n-p-1))
    
    resultRsquared <- adjRsquared
    
  }
  
  
  
  return(resultRsquared)
  
  
}


# calcModelMM
# Calculates a linear model, using the majorization in minimization (MM) algorithm
#
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
#   e: epsilon, parameter for threshold of improvement after which the algorithm should halt
#   nBeta: number of variables one wants to use
#
# Output:
#   model: dataframe, with the following attributes
#       - Beta: dataframe, the calculated Beta's
#       - RSS: double, Sum of squared residuals
#       - Yest: dataframe, the predicted Y
#       - Rsquared: double, R^2 for the predicted Y
#       - Residuals: dataframe, Y - Yest.
#


calcModelMM <- function(mX,mY,e, nBeta){
  
  print("Number of variables:")
  print(nBeta)
  
  if(nBeta > ncol(mX) + 1){
    
    stop("You want to use more variables than there are in the dataset of independent variables")
    
  }
  
  
  # set the previous beta to initial, random beta's
  prevBeta <- getB0(mX)

  # calculate X'X
  mXtX <- t(mX) %*% mX
  
  # get largest eigenvalue for the square of independent variables
  Lambda <- calcLargestEigen(mXtX)
  
  # set initial stepscore to 0, k to 1. 
  StepScore <- 0
  k <- 1
  
  
  # run while, either if k is equal to 1, or the improvement between k-1th and kth set of beta's is smaller than the parameter e
  while (k == 1 | StepScore > e ){
    
    # step to next k
    k <- k + 1
    
    # calculate beta's for this k
    BetaK <- calcBetaK(prevBeta, Lambda, mX,mY)

    # sort the beta's based on absolute value, remove the smallest ones to keep m 
    absBetaKOrdered <- order(abs(BetaK[,1]), decreasing = T)
    BetaK[!BetaK %in% BetaK[absBetaKOrdered,][1:nBeta]] <- 0

    # new stepscore, % difference in RSS between new Beta's and previous beta's
    StepScore <- calcStepScore(mX,mY,prevBeta,BetaK)

    # assign current beta's to prevBeta variable for next iteration
    prevBeta <- BetaK

    
  }
  
  # calculate several attributes of the linear model, put in dataframes or doubles
  BetaFinal <- as.matrix((BetaK))
  RSSBetaK <- calcRSS(mX,mY, BetaK)
  mYest <- calcYest(mX, BetaFinal)
  Rsquared <- calcRsquared(mY, mYest)
  adjRsquared <- calcRsquared(mY,mYest, adjusted = T, p = nBeta-2)
  Resi <- data.frame(residuals = mY - mYest)
  

  # add these attributes together as a list to make it easily accessible
  results <- list(Beta = BetaFinal, RSS = RSSBetaK, Yest = mYest, Rsquared = Rsquared, adjRsquared = adjRsquared, Residuals = Resi)
  

  return(results)
  
}

# findModelMM
# finds the best linear model, using the MM algorithm, by testing model with 1, 2...up to all variables in X
#
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)

findModelMM <- function(mX, mY, e){
  
  nIndVar = ncol(mX)
  M = 1
  
  results <- list()
  
  
  while(M <= nIndVar){
    
    M <- M + 1
    
    resultM <- calcModelMM(mX, mY, e, M)
    
    strSave <- paste0("Model with ", M-1, " variable(s)")
    results[[strSave]] <- resultM

  }
  
  return(results)
  
  
  
}

# make sure working directory is correct
setwd("C:/Users/flori/OneDrive/Documents/GitHub/Supervised ML Code/MM")


# load the air quality data
load("Data/Airq_numeric.Rdata")

# set to dataframe
dfAirQ <- data.frame(Airq)

# select dependent variable of air quality
Yair = dfAirQ$airq

# select all other variables as independent variables
Xair = dfAirQ[,-1]

# scale the independent variables, and add an intercept to these
XairScaled <- scale(Xair)
XairIntercept <- cbind(intercept = 1, XairScaled)

# set the data to matrix format
mYair <- as.matrix(Yair)
mXairIntercept <- as.matrix(XairIntercept)



# set seed to ensure stability of results
set.seed(1)

# set e small
e <- 0.000001


# select the number of beta's you want to use in the model
nBeta <- ncol(mXairIntercept)

# calculate the model using the MM algorithm, using the max (6) variables
modelMM <- calcModelMM(mXairIntercept, mYair, e, nBeta)


# calculate the model with MM, for 1-6 variables
compareModelMM <- findModelMM(mXairIntercept, mYair, e)









