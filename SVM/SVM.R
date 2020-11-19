
setwd("C:/Users/flori/OneDrive/Documents/GitHub/Supervised ML Code/SVM")


# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies,
               SVMMaj,
               dplyr)


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


create_hinge_param <- function(mY, z, type_hinge = "absolute", epsilon = 1e-08){
  
  m <- abs(1 - z)
  m <- m * (m > epsilon) + epsilon * (m <= epsilon)

  a <- .25 * m^-1
  b <- mY * (a + 0.25)
  c <- a + 0.5 + 0.25 * m
  
  result = list(a = a,
                b = b,
                c = c)
}

calc_loss <- function(z, lambda, mWtW, type_loss = "absolute"){
  
  #penalty = lambda * mWtW
  
  vloss <- (1 > z) * (1 - z)# + penalty[[1]]
  
  return(sum(vloss))
}

getHinge


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
  
  
  while(k ==1 || stepscore > epsilon){
    
    # get current prediction (q) and z
    mCurrent_q = mX %*% vV_previous
    fCurrent_z = mY * mCurrent_q
    
    # get parameters given the y and z
    lhinge_param = create_hinge_param(mY, fCurrent_z, epsilon)
    
    # define A and b
    A = diag(as.vector(lhinge_param$a))
    b = lhinge_param$b
    
    # get updated v
    vV_update = solve(t(mX) %*% A %*% mX + lambda * P, t(mX) %*% b)
    
    mW = tail(vV_update,-1)
    mWtW = t(mW) %*% mW

    # get new prediction (q) and z
    mNew_q = mX %*% vV_update
    fNew_z = mY * mNew_q
    
    # calculate new, and previous loss
    fCurrent_loss = calc_loss(fCurrent_z, lambda, mWtW)
    fNew_loss = calc_loss(fNew_z, lambda, mWtW)
    
    # calculate improvement
    stepscore <- (fCurrent_loss - fNew_loss)/fCurrent_loss
    
    print("improvement")
    print(stepscore)
    
    # check: if all predicted correctly, turns to NaN since divided by zero
    if (is.na(stepscore)){
      stepscore = 0
    }
    
    # move to next iteration
    k = k + 1
    vV_previous = vV_update
    
    
    
    
  }
  
  mYhat = sign(mNew_q)
  
  #vPredict_outcome = ifelse(
       #  mYhat == mY == 1, "true positive",
       #  mYhat == mY == -1, "true negative",
       #  myhat == 1, my == -1, "false positive", 
       #  myhat == -1, my == 1, "false negative")
  
  
  #dfPredict_overview = table(vPredict_outcome)
  
  
  
  
  lResults = list(v = vV_update,
                  mY = mY,
                  mX = mX, 
                  loss = fCurrent_loss,
                  q = mNew_q,
                  yhat = mYhat)
  
  return(lResults)

  
  
}


# ensure that we can reproduce the code
set.seed(1)

sample_id <- sample(4000, 1000)


sample_x = mX_var[sample_id,]
sample_y = mY_num[sample_id]

result = svm_mm(sample_y, sample_x)
result$v
result$loss

mYhat = result$yhat
sample_y <- as.matrix(sample_y)
mYhat
vPredict_outcome = ifelse(
  (mYhat == sample_y & sample_y == 1), "true positive",
  ifelse((mYhat == sample_y & sample_y == -1), "true negative",
  ifelse((mYhat == 1 & sample_y == -1), "false positive", 
  ifelse((mYhat == -1 & sample_y == 1), "false negative", NA))))




dfPredict_overview = table(vPredict_outcome)
dfPredict_overview


result_svmmaj <- svmmaj(sample_x,sample_y, lambda = 10, scale = "none",hinge = "absolute")
result_svmmaj$beta 
result_svmmaj$loss
result_svmmaj$q



# check score

# get current prediction (q) and z
mCurrent_q = sample_x %*% result_svmmaj$beta 
fCurrent_z = sample_y* mCurrent_q

calc_loss(fCurrent_z)
