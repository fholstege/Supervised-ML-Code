# getB0: Set Beta's, randomly, between 0 and 1
# 
# Parameters:
#   mX: Matrix of n x p (n = observations, p = independent variables)
# 
# Output: 
#   B0: matrix, set of Beta's between 0 and 1

getB0 <- function(mXm, p){
  
  Beta0 <- runif(p + 1, min=0, max=1)
  
  # turn to matrix format and return
  return(as.matrix(Beta0))
  
}



# calcBetaK: Calculates the kth beta in an MM algorithm
# 
# Parameters:
#   prevBeta: double, k-1th beta
#   Lambda: double, largest eigenvalue of X (independent variables) squared
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mY: Matrix of n x 1 dependent variables (n = observations)
# 
# Output: 
#   BetaK; matrix of new set of coefficients

calcBetaK <- function(prevBeta,Lambda, mX, mXtX, mY){
  
  # calculate the Kth 
  BetaK = prevBeta - ((1/Lambda) *  mXtX %*% prevBeta) + ((1/Lambda) * t(mX) %*% mY )
  
  return(BetaK)
  
  
}



# calcYest: Calculates the predicted Y, based on the X and est. Beta's of a linear model
# 
# Parameters:
#   mX: Matrix of n x p (n = observations, p = independent variables)
#   mBetaEst: Matrix, Estimated Beta's, px1 column, 
#
# Output:
#   Yestdf: matrix, predicted Y
#

calcYest <- function(mX,mBetaEst){
  
  # multiply X with Beta (est.) to get predicted Y
  Yest <- mX %*% mBetaEst
  
  return(as.matrix(Yest))
  
}