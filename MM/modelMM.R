
# Author: Floris Holstege
# Purpose: implements the maximize by minorization algorithm for a linear model, including better subset selection
# Date: 27/10/2020


# calcRSS: Calculates the residual squared errors for a multiple regression of the form Y = XBeta + e
# 
# Parameters: 
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mY: Column matrix of n x 1 dependent variables (n = observations)
#   mBeta: Column Matrix of p x 1 coefficients
#   
# Output:
#   ESquared: double, residual squared errors

calcRSS <- function(mX, mY,mBeta){
  
  # calculate the errors
  mE <-  mY - mX %*% mBeta
  
  # get errors squared
  ESquared <- t(mE) %*% mE
  
  # return the residual sum of squared errors
  return(ESquared[1,1])
  
}

# calcCovar: Calculates the covariance matrix 
#
# Parameters: 
#   RSS: Residual squared errors
#   mXtX: pxp matrix, created from independent variables (X), multiplied with itself
#   n: double, number of observations
#   p: double, number of variables
#
# Output:
#   Covar: matrix, covariance matrix

calcCovar <- function(RSS, mXtX,n, p){
  
  # est. for sigma squared
  SigmaSquared <- (RSS) / (n - p -1)
  
  Covar <- SigmaSquared * as.matrix(inv(mXtX))
  
  return(Covar)
  
}

# calcSignificance: Calculates the statistical significance of a set of beta's
#
# Parameters: 
#   RSS: Residual squared errors
#   mXtX: pxp matrix, created from independent variables (X), multiplied with itself
#   n: double, number of observations
#   p: double, number of variables
#   mBetaEst: matrix of estimated Beta's
#
# Output:
#   dfSignificance: dataframe, containing the results on statistical signficance

calcSignificance <- function(RSS, mXtX, n,p, mBetaEst){
  
  # get covariance matrix
  mCovar <- calcCovar(RSS,mXtX,n,p)
  
  # calculate the standard deviations
  stdev <- sqrt(diag(mCovar))
  
  # define t, which is t-distributed with n-p-1 degrees of freedom 
  t <- mBetaEst/stdev
  pval <- 2*pt(-abs(t),df=n-p-1)
  
  dfSignificance <- data.frame(BetaEst = mBetaEst, 
                               stdev = stdev, 
                               t = t, 
                               pval = pval)
  
  return(dfSignificance)
}


# calcLargestEigen: Calculates the largest eigenvalue of an matrix of independent variables
# 
# Parameters: 
#   mX: Dataframe of n x p (n = observations, p = independent variables)
#   
# Output:
#   LargestEigenval: float, largest eigenvalue of said matrix

calcLargestEigen <- function(mX){
  
  # get the eigenvalues of X 
  EigenValX <- eigen(mX)$values
  
  # from these eigenvalues, get the largest one
  LargestEigenVal <- max(EigenValX, na.rm = TRUE)
  
  return(LargestEigenVal)
  
}


# CalcStepScore: Calculates the % improvement between the k-1th and kth set of beta's
# 
# Parameters:
#   prevBeta: double, k-1th beta
#   currbeta: double, kth beta
#   mX: Dataframe of n x p (n = observations, p = independent variables)
# 
# Output: 
#   StepScore; double, % improvement between the RSS of the two sets of beta's

calcStepScore <- function(mX,mY, prevBeta, currBeta){
  
  # difference in RSS between previous and current set of beta's
  diffRSS <- (calcRSS(mX,mY,prevBeta) - calcRSS(mX,mY,currBeta))
  
  # divide difference with previous score to get % change
  StepScore <- diffRSS /calcRSS(mX,mY,prevBeta)
  
  return(StepScore)
  
}

# calcRsquared:  Calculates the r-squared
#
# Parameters:
#   Y: matrix, the true dependent variable   
#   Yest: matrix, the predicted dependent variable
#   (optional) adjusted: if True, return adjusted r squared
#   (optional) p: if adjusted is calculated, add number of variables
# 
# Output:
#   Rsquared: double, the Rsquared or adjusted Rsquared for a linear model

calcRsquared <- function(mY, mYest, adjusted = FALSE, p=0, n=0){
  
  # standardize Y, and Yest (mean of 0)
  mStandY = mY - mean(mY)
  mStandYest = mYest - mean(mYest)
  
  # calculate Rsquared
  numerator <- (t(mStandY) %*% mStandYest)^2
  denominator <- (t(mStandY) %*% mStandY) %*% (t(mStandYest) %*% mStandYest)
  resultRsquared <- (numerator/denominator)
  
  # if want adjusted R squared, 
  if(adjusted){
    
    adjRsquared = 1 - (((1-resultRsquared)*(n - 1))/(n-p-1))
    
    resultRsquared <- adjRsquared
    
  }
  
  return(resultRsquared)
  
}


# calcModelMM:  Calculates a linear model, using the majorization in minimization (MM) algorithm
#
# Parameters:
#   X: Dataframe of n x p (n = observations, p = independent variables)
#   Y: Dataframe of n x 1 dependent variables (n = observations)
#   e: epsilon, parameter for threshold of improvement after which the algorithm should halt
#   nBeta: number of variables one wants to use
#
# Output:
#   result: dataframe with attributes of the model: 
#       - Beta: dataframe, the calculated Beta's
#       - RSS: double, Sum of squared residuals
#       - Yest: dataframe, the predicted Y
#       - Rsquared: double, R^2 for the predicted Y
#       - AdjRsquared: Adjusted Rsquared
#       - Significance results: dataframe with significance results on the beta's
#       - Residuals: dataframe, Y - Yest.
#

calcModelMM <- function(mX,mY,e, p){
  
  
  # get number of observations
  n <- nrow(mX)
  
  
  # check the user has filled in an appropriate amount of beta's
  if(p > ncol(mX) - 1){
    
    stop("You want to use more variables than there are in the dataset of independent variables")
    
  }
  
  # set the previous beta variable to initial, random beta's
  prevBeta <- runif(ncol(mX), min=0, max=10)
  
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
    BetaK <- prevBeta - ((1/Lambda) *  mXtX %*% prevBeta) + ((1/Lambda) * t(mX) %*% mY )
    
    # sort the beta's based on absolute value, remove the smallest ones to keep m 
    absBetaKOrdered <- order(abs(BetaK[,1]), decreasing = T)
    BetaK[!BetaK %in% BetaK[absBetaKOrdered,][0:p+1]] <- 0
    
    # new stepscore, % difference in RSS between new Beta's and previous beta's
    StepScore <- calcStepScore(mX,mY,prevBeta,BetaK)
    
    # assign current beta's to prevBeta variable for next iteration
    prevBeta <- BetaK
    
    
  }
  
  ## Calculate several attributes of the linear model, put in dataframes or doubles
  BetaFinal <- as.matrix(BetaK)
  
  # calculate the RSS of this final est.
  RSSBetaK <- calcRSS(mX,mY, BetaK)
  
  # get the est. dependent variables
  mYest <- mX %*% BetaFinal
  
  # get the r2 and adjusted r2
  Rsquared <- calcRsquared(mY, mYest)
  adjRsquared <- calcRsquared(mY,mYest, adjusted = T, p, n)
  
  # get the residuals
  Resi <- mY - mYest
  
  # get the results on significance
  dfSignificance <- calcSignificance(RSSBetaK, mXtX, n, p, BetaFinal)
  
  # add these attributes together as a list to make it easily accessible
  result <- list(Beta = BetaFinal,
                 RSS = RSSBetaK, 
                 Yest = mYest,
                 Rsquared = Rsquared, 
                 adjRsquared = adjRsquared, 
                 SignificanceResults = dfSignificance,
                 Residuals = Resi, 
                 n = n,
                 p = p)
  
  return(result)
  
}

# findModelMM: finds the best linear model, using the MM algorithm, by testing model with 1, 2...up to all variables in X
#
# Parameters:
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mY: Matrix of n x 1 dependent variables (n = observations)
#
# Output:
#   results: list with the results for each model version

findModelMM <- function(mX, mY, e){
  
  # get the number of independent variables used
  nIndVar = ncol(mX) - 1
  
  # start at m = 1, create empty list to be filled with results
  M = 1
  results <- list()
  
  # for each m, check the best model and save the results
  while(M <= nIndVar){
    
    resultM <- calcModelMM(mX, mY, e, M)
    
    strSave <- paste0("Model with ", M, " variable(s)")
    results[[strSave]] <- resultM
    
    M <- M + 1
    
  }
  
  return(results)
  
  
  
}



# load the air quality data
load("Data/Airq_numeric.Rdata")

# set to dataframe
dfAirQ <- data.frame(Airq)

# select dependent variable of air quality
Yair = dfAirQ$airq

# select all other variables as independent variables
Xair = dfAirQ[,-1]

# summary stats & correlation matrix
stargazer(dfAirQ)

corMatrix <- cor(Xair)


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

# calculate the model with MM, for 1-5 variables. This contains all the values shown in the paper 
compareModelMM <- findModelMM(mXairIntercept, mYair, e)

# rescale to check 
BetaCoastInModel2 <- compareModelMM$`Model with 2 variable(s)`$Beta[,1][4]
BetaCoestInModel2_rescaled <- (mean(dfAirQ$coasyes) + BetaCoastInModel2 )/sd(dfAirQ$coasyes)

# check heteroskedasticity
preferredModel <- compareModelMM$`Model with 2 variable(s)`
plot(preferredModel$Yest,preferredModel$Residuals, 
     ylab = "Residuals", 
     xlab = "Y est.",
     pch=19)
