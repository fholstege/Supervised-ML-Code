
# Functions
#-------------------------------------------------------------------------------

# Residual Sum of Squares (RSS)
# For calculating the residual sum of squares 
#
# Parameters: 
#   y: numeric vector of n x 1
#   X: numeric matrix of n x p
#   beta: numeric vector of p x 1
#   
# Output:
#   RSS: numeric scalar
# 
RSS <- function(y,X,beta) {
  u <- (y - X %*% beta)
  RSS <- t(u) %*% u
  return(RSS[1])
}

# Fast Ortogonalizing Subset Screening ((F)OSS)
# Implementation of the (F)OSS algorithm as presented in
# Xiong, S. (2014). Better subset regression.Biometrika,101(1), 71-84.
# 
# Parameters: 
#   y: numeric vector of n x 1
#   X: numeric matrix of n x p
#   m: positive integer (1<=m<=p)
#   intercept: bool
#   eps: numeric
#   beta: numeric vector of p x 1
#   verbose_at: positive integer
#   convergence_limit: positive integer
#   
# Output:
#   list

OSS <- function(y, X, m, intercept=FALSE, eps=1e-8, fast=FALSE, 
                verbose_at=0, converge_limit=2e6) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  # Error handling
  if (intercept) {
    X <- cbind(c=1,X)
    p = p + 1
    m = m + 1
  }
  
  if (m > p | m < 1) {
    print("ERROR: m should be greater than 0 and less than or equal to p!")
    return(list("beta"=NULL, "RSS"=NULL,
                "k"=NULL, "converged"=NULL, "Ok"=FALSE))
  }
  
  if (nrow(y) != nrow(X)) {
    print("ERROR: Number of rows in y should be equal to that of X!")
    return(list("beta"=NULL, "RSS"=NULL, "k"=NULL, "converged"=NULL))
  }
  
  # calculate X'X and X'y
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  
  # retrieve the largest eigenvalue of X'X
  lambda <- max(eigen(XtX)$values)
  
  # initial estimated solution (ols)
  beta_k <- solve(XtX, Xty)
  
  # initial RSS
  RSS_k <- RSS(y, X, beta_k)
  RSS_k1 <- RSS_k - 1e6
  
  k <- 0
  
  # to store the selected regressors 
  X_A <- NULL
  
  # iteratively update beta_k by the steps in Xiong (2014) page 12
  while (k < converge_limit & abs(RSS_k - RSS_k1) > eps) {
    
    k <- k + 1
    RSS_k <- RSS_k1
    
    varphi <- beta_k - (XtX %*% beta_k - Xty) / lambda
    
    U <- order(abs(varphi), decreasing = TRUE)
    
    # use the m largest values of U to update beta_k (and X_A if FOSS)
    if (fast) {
      X_A = X[,U[1:m]]
      beta_k[U[1:m]] <- solve(t(X_A) %% X_A, t(X_A) %% y)
    } else {
      beta_k[U[1:m]] <- varphi[U[1:m]]
    }
    
    # set smallest p-m values of beta to 0
    if (m < p) {
      beta_k[U[m+1:p]] <- 0
    }
    
    # obtain RSS of new estimate
    RSS_k1 <- RSS(y,X,beta_k)
    
    if (verbose_at & k %% verbose_at == 0){
      cat("k = ",k," - RSS = ", RSS_k1,sep="",fill=TRUE)
    } 
  }
  
  names(beta_k) <- colnames(X)
  
  if (k < converge_limit) {
    converged = TRUE
  } else {
    converged = FALSE
  }
  
  # if OSS algorithm is used select the subset of regressors obtained from
  # the above iterative procedure
  if (!fast) {
    X_A <- X[,U[1:m]]
  }
  
  
  if (m == 1) {
    X_A <- data.matrix(X_A)
    colnames(X_A) <- colnames(X)[U[1:m]]
    
    XAtXA_inv <- 1/(t(X_A) %*% X_A)
  } else {
    
    # obtain inverse of subset
    XAtXA_inv <- solve(t(X_A) %*% X_A)
  }
  
  
  # degrees of freedom
  df <- (n - m)
  
  # error
  e <- y - X %*% beta_k
  
  # variance of the error
  sigma2 <- RSS_k1 / df
  
  # obtain total sum of squares
  RSS_tot <- t(y - mean(y)) %*% (y - mean(y))
  
  # calculate (adjusted) R-sqaured
  R2 <- 1 - RSS_k1 / RSS_tot
  
  if (m != 1) {
    R2_adjust <- 1 - ( (1 - R2) * (n - 1) ) / df
  } else{
    R2_adjust <- 1 - (1 - R2) * n / (n - 1)
  }
  
  #F-stat and p-value
  
  if (m == 1) {
    q <- m
  } else {
    q <- m-1
  }
  F_stat <- ( (RSS_tot - RSS_k1) / RSS_k1 )* ( df / q )
  F_p <- pf(F_stat, df1=q, df2=df, lower.tail=FALSE)
  
  beta_k <- data.frame(t(beta_k))
  
  # obtain standard deviations
  beta_std <- beta_k * 0
  
  if (m == 1) {
    beta_std[colnames(X_A)] <- sqrt(sigma2 * XAtXA_inv)
  }
  beta_std[colnames(X_A)] <- sqrt(sigma2 * diag(XAtXA_inv))
  
  # obtain t-values
  beta_t <- beta_k * 0
  beta_t[colnames(X_A)] <- ( beta_k[colnames(X_A)] / beta_std[colnames(X_A)])
  
  # obtain p-values
  beta_p <- beta_k * 0
  beta_p[colnames(X_A)] <- 2 * pt(abs(data.matrix(beta_t[colnames(X_A)]))
                                  ,df=df, lower.tail=FALSE)
  
  data <- rbind("beta"=beta_k,"beta_std"=beta_std,
                "beta_t"=beta_t,"beta_p"=beta_p)
  
  return(list("results"=t(data), "R2"=R2[1], "R2_adjust"=R2_adjust[1], 
              "Fstat"=F_stat, "F_p_value"=F_p, "RSS"=RSS_k1, "X_A"=X_A,
              "err"=e, "k"=k, "converged"=converged))
}


# Read data
#-------------------------------------------------------------------------------

data(Airq)


X_names_numeric <- c("rain","dens","medi","vala")

colnames(Airq)

# split data in dependent (y) and independent (X) variables
y <- data.matrix(Airq)
y <- y - mean(y)
X <- Airq[,-which(colnames(Airq) == "airq")]


# make a dummy variable
lookup <- c("No" = 0,"Yes" = 1)
X$coas <- lookup[X[,3]]

# convert other to numeric variables
X_numeric <- X[X_names_numeric]
X_numeric <- sapply(X_numeric,as.numeric)

# scale numeric variables
X_numeric_scaled <- scale(X_numeric)

# combine dummy and scaled variables
X_full <- cbind(X_numeric_scaled,coas = X[,3] )

# total numer of regressors
p <- ncol(X_full)
p

########
# run this section if you want to compare our output with that of the lm package
#bs <- OSS(y, X_full, 1, intercept=TRUE, fast=FALSE, verbose_at=1e5)
#bs

#X_df <- data.frame(cbind(y,bs$X_A))
#X_df$c <- NULL
#model <- lm(V1~.,X_df)
#summary(model)
########

# obtain variance inflation factors
VIF <- vector(mode="list", length=p+1)
for (M in 1:p) {
  X_full2 <- X_full
  
  X_target <- data.matrix(X_full2[,M])
  X_ <- X_full2[,-M]
  
  bs <- OSS(X_target, X_, ncol(X_), intercept=TRUE, fast=TRUE)
  VIF[[colnames(X_full2)[M]]] <- 1/(1-bs$R2)
  
}
print(VIF)

# obtain the FOSS results for different values of M
FOSS_I <- vector(mode="list", length=p+1)
for (M in 0:p) {
  bs <- OSS(y, X_full, M, intercept=TRUE, fast=TRUE)
  FOSS_I[[M+1]] <- bs
  
  print(M+1)
  # print(summary(model))
  #print(bs$F_p_value)
}

# select the model with 3 regressors
FOSS_3 <- FOSS_I[[3]]

g <- FOSS_3$err

# plots
h <- hist(g, density = 10,breaks=10,
          col = "lightgray", xlab = "Error", main = " ") 
xfit <- seq(-80, 80, length = 100) 
yfit <- dnorm(xfit, mean = mean(g), 
              sd = sqrt(FOSS_3$RSS/(nrow(FOSS_3$X_A)-ncol(FOSS_3$X_A))) )
yfit <- yfit * diff(h$mids[1:2]) * length(g) 

lines(xfit, yfit, col = "black", lwd = 2)

plot(FOSS_3$err^2,xlab='Observation',ylab='Squared error')

# shapiro
shapiro.test(FOSS_3$err)