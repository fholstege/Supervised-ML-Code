
setwd("C:/Users/flori/OneDrive/Documents/GitHub/Supervised ML Code/SVM")


# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies,
               SVMMaj,
               dplyr,
               caret)


# load bank data
load("bank.RData")



# get dependent variable and transform
mY <- bank[,ncol(bank)]
mY_num <- as.matrix(ifelse(mY == "yes", 1, -1))

# get independent variables, scale them, add intercept
mX = bank[,-ncol(bank)]

# create dummy variables 
df_toDummy = as.matrix(data.frame(mX$job,
                        mX$marital, 
                        mX$education, 
                        mX$default, 
                        mX$housing, 
                        mX$loan, 
                        mX$contact,
                        mX$month, 
                        mX$day_of_week,
                        mX$poutcome))

dummy_vars <-  dummy_columns(df_toDummy)
dummy_vars <- dummy_vars[,(1+ncol(df_toDummy)):ncol(dummy_vars)]

# select numeric variables, scale these
numeric_vars <- as.matrix(data.frame(mX$age, 
                           mX$duration,
                           mX$campaign, 
                           mX$previous, 
                           mX$cons.price.idx, 
                           mX$cons.conf.idx, 
                           mX$nr.employed))



# combine numeric and dummy variables
mX_var <- as.matrix(cbind(1, scale(numeric_vars),dummy_vars))

getHinge

create_hinge_param <- function(mY, z, hinge = "absolute", epsilon = 1e-08, delta =1){
  
  if(hinge == "absolute"){
  
    m <- abs(1 - z)
    m <- m * (m > epsilon) + epsilon * (m <= epsilon)
  
    a <- .25 * m^-1
    b <- mY * (a + 0.25)
    
    return(result = list(a = a,
                          b = b))

  }else if (hinge == "quadratic"){
    
    m <- z* (z > 1) + (1 >= z)
    b <- m * mY
    
    return(list(b = b))
    
  }else if(hinge == "huber"){
    
    m <- (z <- delta) * (z + delta)
  }
  
}


calc_loss <- function(z, mW, lambda, hinge = "absolute"){
  
  mWtW = t(mW) %*% mW
  penalty = lambda * mWtW
  
  if(hinge == "absolute"){
    
    vloss <- (1 > z) * (1 - z)
    
  }else if(hinge == "quadratic"){
    
    
    vloss <- (1 > z) * (1 - z)^2
  }

  return(sum(vloss) + penalty)
}




svm_mm <- function(mY, mX, hinge = "absolute", lambda = 10, epsilon = 1e-08){
  
  # set initial weights and constant. From these, create initial v = c + w
  vW = runif(ncol(mX)-1, -1, 1)
  fConstant = 0.0
  vV_previous = as.matrix(c(fConstant, vW))
  
  # get n of observations in the sample, and p variables
  n = nrow(mX)
  p = ncol(mX)

  # define matrix p
  P = diag(p)
  P[1,1]<-0
  
  # set initial values; k, stepscore 
  k = 1
  stepscore = 0 
  
  if(hinge == "quadratic"){
    Z = inv(t(mX) %*% mX + lambda * P) %*% t(mX)
  }
  
  
  while(k ==1 || stepscore > epsilon){
    
    # get current prediction (q) and z
    mCurrent_q = mX %*% vV_previous
    fCurrent_z = mY * mCurrent_q
    
    # get parameters given the y and z
    lhinge_param = create_hinge_param(mY, fCurrent_z, hinge, epsilon)
    
    b = lhinge_param$b
    
    if(hinge == "quadratic"){
      
      vV_update = Z %*% b
      
    }else{
      
      # define A and b
      A = diag(as.vector(lhinge_param$a))
      
      # get updated v
      vV_update = solve(t(mX) %*% A %*% mX + lambda * P, t(mX) %*% b)
      
    }
    
    mW_previous = tail(vV_previous,-1)
    mW_update = tail(vV_update, -1)


    # get new prediction (q) and z
    mNew_q = mX %*% vV_update
    fNew_z = mY * mNew_q
    
    # calculate new, and previous loss
    fCurrent_loss = calc_loss(fCurrent_z, mW_previous, lambda, hinge)
    fNew_loss = calc_loss(fNew_z,mW_update, lambda, hinge)
    
    # calculate improvement
    stepscore <- (fCurrent_loss - fNew_loss)/fCurrent_loss
    

    # check: if all predicted correctly, turns to NaN since divided by zero
    if (is.na(stepscore)){
      stepscore = 0
    }
    
    # move to next iteration
    k = k + 1
    vV_previous = vV_update
  
    
    
  }
  
  mYhat = sign(mNew_q)
  
  resultOverview <- confusionMatrix(as.factor(mY), as.factor(mYhat))
  mConfusionTable <- resultOverview$table
  fAccuracy <- resultOverview$overall[1]

  
  lResults = list(v = vV_update,
                  mY = mY,
                  mX = mX, 
                  loss = fCurrent_loss,
                  q = mNew_q,
                  yhat = mYhat,
                  Accuracy = fAccuracy,
                  ConfusionTable = mConfusionTable)
  
  return(lResults)

  
  
}




# ensure that we can reproduce the code
set.seed(1)

# pick a random sample of 1000
sample_id <- sample(4000, 1000)
sample_x = mX_var[sample_id,]
sample_y = mY_num[sample_id]


# compare our algorithm with SVM Maj - we get the same predictions (see confusion matrix), and exactly the same loss
# but still different beta's? most likely something to do with the transformation of the data...
result = svm_mm(sample_y, sample_x, lambda = 100, hinge = "quadratic")
result$v
result$loss
result$ConfusionTable
result$Accuracy


result_svmmaj <- svmmaj(sample_x,sample_y, lambda = 100, scale = "none",hinge = "quadratic")
result_svmmaj$beta 
result_svmmaj$loss
result_svmmaj$q





create_hinge_param()

