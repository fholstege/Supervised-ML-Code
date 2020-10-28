
# Author: Floris Holstege
# Purpose: implements the maximize by minorization algorithm for a linear model
# Date: 27/10/2020
#testingtesting




# calcRSS
# Calculates the residual squared errors for a multiple regression of the form Y = XBeta + e
# 
# Parameters: 
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
#   Beta: p x 1 double with coefficients of said model
#   
# Output:
#   ESquared: float, residual squared errors
# 

calcRSS <- function(X, Y, Beta){
  
  # set the dataframes to matrices
  mX = as.matrix(X)
  mY = as.matrix(Y)
  mBeta = as.matrix(Beta)
  

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
  
  # get matrix of squared X
  mX = as.matrix(X)
  XSquared <- t(mX) %*% mX
  
  # get the eigenvalues of X squared
  EigenValXSquared <- eigen(XSquared)$values
  
  # from these eigenvalues, get the largest one
  LargestEigenVal <- max(EigenValXSquared, na.rm = TRUE)
  
  return(LargestEigenVal)
  
}

#CalcBetaK
# Calculates the kth beta in an MM algorithm
# 
# Parameters:
#   prevBeta: double, k-1th beta
#   Lambda: double, largest eigenvalue of X (independent variables) squared
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
# 
# Output: 
#   BetaK; matrix of new set of coefficients

calcBetaK <- function(prevBeta,Lambda, X, Y){
  
  # get matrix of squared X
  mX = as.matrix(X)
  XSquared <- t(mX) %*% mX
  
  # turn Y into matrix
  mY = as.matrix(Y)
  
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
# Set the initial set of Beta's, randomly, between 0 and 1
# 
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
# 
# Output: 
#   B0: double, set of Beta's between 0 and 1

getB0 <- function(X){
  
  # determine the number of independent variables, generate as many random beta's
  nIndVar = ncol(X)
  Beta0 <- runif(nIndVar, min=0, max=1)
  
  return(Beta0)
  
}

# calcYest
# Calculates the predicted Y, based on the X and est. Beta's of a linear model
# 
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   BetaEst: Estimated Beta's, px1 vector, 
#
# Output:
#   Yestdf: Dataframe, predicted Y
#

calcYest <- function(X,BetaEst){
  
  # turn X and Beta's (est.) into matrix
  mBetaEst <- as.matrix(BetaEst)
  mX <- as.matrix(X)
  
  # multiply X with Beta (est.) to get predicted Y
  Yest <- mX %*% mBetaEst
  
  # turn into dataframe
  Yestdf <- as.data.frame(Yest)
  colnames(Yestdf) <- c("Yest")
  
  return(Yestdf)
  
}

# calcRsquared
# Calculates the r-squared
#
# Parameters:
#   Y: dataframe, the true dependent variable   
#   Yest: dataframe, the predicted dependent variable
# 
# Output:
#   Rsquared: double, the Rsquared for a linear model

calcRsquared <- function(Y, Yest){
  
  # standardize Y, and Yest (mean of 0)
  standardY = Y - mean(Y)
  standardYest = Yest - mean(Yest$Yest)

  # turn into matrix to perform multiplication
  mY <- as.matrix(standardY)
  mYest <- as.matrix(standardYest)
  
  
  # calculate Rsquared
  numerator <- (t(mY) %*% mYest)^2
  denominator <- (t(mY) %*% Y) %*% (t(mYest) %*% mYest)
  Rsquared <- (numerator/denominator)
  
  return(Rsquared)
  
  
}


# calcModelMM
# Calculates a linear model, using the majorization in minimization (MM) algorithm
#
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
#
#
# Output:
#   model: dataframe, with the following attributes
#       - Beta: dataframe, the calculated Beta's
#       - RSS: double, Sum of squared residuals
#       - Yest: dataframe, the predicted Y
#       - Rsquared: double, R^2 for the predicted Y
#       - Residuals: dataframe, Y - Yest.
#


calcModelMM <- function(X,Y,e){
  
  # set the previous beta to initial, random beta's
  prevBeta <- getB0(X)
  
  # get largest eigenvalue for the square of independent variables
  Lambda <- calcLargestEigen(X)
  
  # set initial stepscore to 0, k to 1. 
  StepScore <- 0
  k <- 1
  
  
  # run while, either if k is equal to 1, or the improvement between k-1th and kth set of beta's is smaller than the parameter e
  while (k == 1 | StepScore > e ){
    
    # step to next k
    k <- k + 1
    
    # calculate beta's for this k
    BetaK <- calcBetaK(prevBeta, Lambda, X,Y)

    # new stepscore, % difference in RSS between new Beta's and previous beta's
    StepScore <- calcStepScore(X,Y,prevBeta,BetaK)

    # assign current beta's to prevBeta variable for next iteration
    prevBeta <- BetaK
    
  }
  
 
  # calculate several attributes of the linear model, put in dataframes or doubles
  BetaFinal <- data.frame(BetaK)
  RSSBetaK <- calcRSS(X,Y, BetaK)
  Yest <- calcYest(X, BetaFinal)
  Rsquared <- calcRsquared(Y, Yest)
  Resi <- data.frame(residuals = Y - Yest$Yest)

  # add these attributes together as a list to make it easily accessible
  results <- list(Beta = BetaFinal, RSS = RSSBetaK, Yest = Yest, Rsquared = Rsquared, Residuals = Resi)
  

  return(results)
  
}






