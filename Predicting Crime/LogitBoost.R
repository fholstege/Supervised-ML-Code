

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(gbm, 
               caret,
               caTools)

# get smoking data
load("../trees/smoking.Rdata")

# preprocessing to smoke data
smoking$intention_to_smoke <- ifelse(smoking$intention_to_smoke == "yes", 1,0)
smoking$lied_to_parents <- ifelse(smoking$lied_to_parents == "yes", 1, 0)
smoking$friends_smoke <- ifelse(smoking$friends_smoke == "one or more", 1,0)


# set seed, create training and test data
set.seed(2)
folds <- createDataPartition(y = smoking$intention_to_smoke, p = 0.7)
smoke_train <- smoking[folds$Resample1, ]
smoke_test <- smoking[-folds$Resample1, ]

##################################
# Own implementation
##################################



BoostLogit <- function(formula, df, iM_classifiers, iZmax = 3, sBaseLearner = "WLS"){
  
  # coerce to data.frame
  df <- as.data.frame(df)
  
  # handle formula
  formula <- terms.formula(formula, data = df)
  
  # get the independent variable matrix
  mX <- model.matrix(formula, df)[,-1]

  # get the dependent variable, and turn into dummy
  vY_raw <- df[, as.character(formula)[2]]
  vY <- dummy_cols(vY_raw)[-1]
  colnames(vY) <- sort(unique(vY_raw), decreasing = FALSE)

  # define the sample size
  iN = nrow(mX)
  iK = ncol(vY)

  # define the initial probabilities per class 
  mP_overall <- matrix(1 / iK, nrow = iN, ncol = iK)
  mf <- matrix(0, nrow = iN, ncol = iK)
  mF <- matrix(0, nrow = iN, ncol = iK)
  
  # create vector where coefficients are stored
  if(sBaseLearner == "WLS"){
    vCoefficients <- matrix(0,nrow = ncol(mX) + 1, ncol = iK)
  }
  
  # update integer with each weak classifier added
  iClassifier = 0
  
  # minimum for weights
  iEps = 1e-24
  
  # learning rate
  iLearnRate <-  (iK - 1) / iK
  
  # likelihood 
  fll_previous =  sum(vY * log(mP_overall)) / iN
  
  # Continue updating until as many weak classifiers as defined in iM_classifiers
  while(iClassifier < iM_classifiers){
    
    iClassifier = iClassifier+ 1
    
    vCategoryRange <- 1:iK
    
    calcProb_category <- function(iCategory, mP_overall,vY, iEps, mf,mF, ...){
      # get current probabilities
      vP_category = mP_overall[,iCategory]
      
      # get dependent for this category
      vY_category = vY[,iCategory]
      
      # determine current weights
      vWeights = pmax(vP_category * (1-vP_category),iEps)
      
      # determine the vZ - ensure it does not surpass min/max defined
      vZ_category <- ifelse(vY_category == 1, 1/vP_category, ifelse(vY_category == 0, -1/(1-vP_category), (vY_category - vP_category)/vWeights))
      vZ_category <- pmin(vZ_category, iZmax)
      vZ_category <- pmax(vZ_category, -iZmax)
      
      # create weak classifier with WLS
      vWeakClassifier <- lm(vZ_category ~ mX, weights = vWeights) 
      
      print(vWeakClassifier$coefficients)
      
      # save coefficients of weak classifier
      vCoefficients[,iCategory] <- vWeakClassifier$coefficients
      
      # get probabilities from weak classifier, and add to matrix to save result
      vWeakClassifier_p <- vWeakClassifier$fitted.values
      mf[,iCategory] <- vWeakClassifier_p
      
      # determine overall probabilities 
      mF[, iCategory] <- mF[, iCategory] + iLearnRate * (mf[, iCategory] - rowSums(mf) / iK)
      
    }
    
    mclapply(vCategoryRange, calcProb_category, mP_overall = mP_overall, vY = vY, iEps = iEps, mf=mf, mF=mF,vCoefficients=vCoefficients, mc.cores = detectCores())

   
    
    # update the probabilies after all categories are calculated
    for (iCategory in 1:iK){
      
      # update current probabilities 
      mP_overall[, iCategory]<- exp(mF[, iCategory]) / (rowSums(exp(mF)))

    }
    
    # calculate the new log likelihood
    fll_new <- sum(vY*log(mP_overall))/iN
    
    # minimize log likelihood, break while loop when done
    if (fll_new < fll_previous) {
      break
    }
    
    # update
    fll_previous <- fll_new
  }
  
  # get dataframe of probabilities per category
  dfP_overall <- data.frame(mP_overall)
  colnames(dfP_overall) <- colnames(vY)

  # create vector with most likely category per observation
  vPredicted_category <- colnames(dfP_overall)[apply(dfP_overall,1,which.max)]
  
  # turn that vector numeric if the categories are numeric
  if (any(!is.na(as.numeric(test_r)))){
    vPredicted_category <- as.numeric(vPredicted_category)
  }
  
  print(vCoefficients)
}



LB <- function(y_train, X_train, M = 1000, v = 1, eps = 1e-4) {
  
  N <- dim(y_train)[1]
  K <- dim(y_train)[2]
  
  f <- matrix(0, nrow = N, ncol = K)
  z <- matrix(0, nrow = N, ncol = K)
  p <- matrix(1 / K, nrow = N, ncol = K)
  F_ <- matrix(0, nrow = N, ncol = K)
  
  ll_prev <- sum(y_train * log(p)) / N
  err <- ll_prev + 1e4
  gamma <- v * (K - 1) / K
  
  coefs = list()
  for (m in 1:M) {
    
    for (k in 1:K) {
      wk <- p[, k] * (1 - p[, k]) + 1e-5
      
      z[y_train[, k] == 1, k] <- pmin(1 / p[y_train[, k] == 1, k], 3)
      z[y_train[, k] == 0, k] <- pmax(-1 / (1 - p[y_train[, k] == 0, k]), -3)
      
      wls <- lm(z[,k] ~ X_train, weights = wk) 
      coefs[[k]] <- wls$coefficients
      
      print(coefs[[k]])
      f[, k] <- wls$fitted.values
      F_[, k] <- F_[, k] + gamma * (f[, k] - rowSums(f) / K)
    }
    
    for (k in 1:K) {
      p[, k]<- exp(F_[, k]) / (rowSums(exp(F_)))
    }
    
    ll_new <- sum(y_train*log(p))/N
    err <- abs(ll_new - ll_prev)
    print(paste0('Iteration: ',m, ' - Loglik: ', ll_new,sep=''))
    
    if (ll_new < ll_prev) {
      break
    }
    
    ll_prev <- ll_new
  }
  
  function(X_test) {
    N_test <- dim(X_test)[1]
    X_test <- cbind(rep(1, N_test), X_test)
    P <- matrix(0, nrow = N_test, ncol = K)
    y_hat <- matrix(0, nrow = N_test, ncol = K)
    
    for (k in 1:K) {
      P[,k] <- X_test %*% coefs[[k]]
      print(coefs[[k]])
    }
    
    for (i in 1:N_test) {
      idx <- which.max(P[i,])
      y_hat[i, idx] <- 1
      y_hat[i, -idx] <- 0
    }
    return(list("preds" = y_hat,"probs" = P))
  }
}

  
  



y_train <- dummy_cols(smoke_train$intention_to_smoke)[,-1]
x_train <- smoke_train[,-5]

t <- LB(y_train, as.matrix(x_train), M = 10)
test_r <- BoostLogit(intention_to_smoke ~.,  df = smoke_train , iM_classifiers = 10)




# compare to mboost
ctrl <- boost_control(mstop = 10)
Mboost_result = gamboost(as.factor(intention_to_smoke) ~ ., data = smoke_train, family = Binomial(), control = ctrl, baselearner = "bols")





# go through each category
for(iCategory in 1:iK){
  
  # get current probabilities
  vP_category = mP_overall[,iCategory]
  
  # get dependent for this category
  vY_category = vY[,iCategory]
  
  # determine current weights
  vWeights = pmax(vP_category * (1-vP_category),iEps)
  
  # determine the vZ - ensure it does not surpass min/max defined
  vZ_category <- ifelse(vY_category == 1, 1/vP_category, ifelse(vY_category == 0, -1/(1-vP_category), (vY_category - vP_category)/vWeights))
  vZ_category <- pmin(vZ_category, iZmax)
  vZ_category <- pmax(vZ_category, -iZmax)
  
  # create weak classifier with WLS
  vWeakClassifier <- lm(vZ_category ~ mX, weights = vWeights) 
  
  print(vWeakClassifier$coefficients)
  
  # save coefficients of weak classifier
  vCoefficients[,iCategory] <- vWeakClassifier$coefficients
  
  # get probabilities from weak classifier, and add to matrix to save result
  vWeakClassifier_p <- vWeakClassifier$fitted.values
  mf[,iCategory] <- vWeakClassifier_p
  
  # determine overall probabilities 
  mF[, iCategory] <- mF[, iCategory] + iLearnRate * (mf[, iCategory] - rowSums(mf) / iK)
  
}
