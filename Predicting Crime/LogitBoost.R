

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(gbm, 
               caret,
               caTools,
               mclust)

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

BoostLogit_train <- function(formula, df, iM_classifiers, iZmax = 3, sBaseLearner = "WLS"){
  
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
  
  # range from 1 to number of categories
  vCategoryRange <- 1:iK
  
  # Continue updating until as many weak classifiers as defined in iM_classifiers
  while(iClassifier < iM_classifiers){
    
    # update until max base learners reached
    iClassifier = iClassifier+ 1

    for(iCategory in vCategoryRange){
      
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

      # save coefficients of weak classifier
      vCoefficients[,iCategory] <- vWeakClassifier$coefficients

      # get probabilities from weak classifier, and add to matrix to save result
      vWeakClassifier_p <- vWeakClassifier$fitted.values
      mf[,iCategory] <- vWeakClassifier_p
      
      # determine overall probabilities 
      mF[, iCategory] <- mF[, iCategory] + iLearnRate * (mf[, iCategory] - rowSums(mf) / iK)

    }

    
    # update the probabilities after all categories are calculated
    for (iCategory in 1:iK){
      
      # update current probabilities 
      mP_overall[, iCategory]<- exp(mF[, iCategory]) / (rowSums(exp(mF)))

    }
    
    # calculate the new log likelihood
    fll_new <- sum(vY*log(mP_overall))/iN
    
    print(paste0('Iteration: ',iClassifier, ' - Loglik: ', fll_new,sep=''))
    
    # minimize log likelihood, break while loop when done
    if (fll_new < fll_previous) {
      break
    }
    
    # update
    fll_previous <- fll_new
  }
  
  if (sBaseLearner == "WLS"){
    modelParam = list(coefficients = vCoefficients, 
                    BaseLearner = sBaseLearner,
                    Categories = colnames(vY),
                    iK = iK)
  }

  return(modelParam)
  
}

BoostLogit_test <- function(modelParam, mX_test){
  
  mP_predicted <- matrix(0, nrow = nrow(mX_test), ncol = modelParam$iK)
  
  mX_test <- cbind(1, mX_test)
  
  
  for (iCategory in 1:modelParam$iK) {
    mP_predicted[,iCategory] <- mX_test %*% modelParam$coefficients[,iCategory]
  }
  

  # get dataframe of probabilities per category
  dfP_predicted <- data.frame(mP_predicted)
  colnames(dfP_predicted) <- modelParam$Categories
  
  # create vector with most likely category per observation
  vPredicted_category <- colnames(dfP_predicted)[apply(dfP_predicted,1,which.max)]

  # turn that vector numeric if the categories are numeric
  if (any(!is.na(as.numeric(test_r)))){
    vPredicted_category <- as.numeric(vPredicted_category)
  }
  
  lResult_test <- list(predicted_category = vPredicted_category, probabilities = dfP_predicted)
  
  return(lResult_test)
}

BoostLogit <- function(formula, dfTrain, dfTest, iM_classifiers, iZmax = 3, sBaseLearner = "WLS"){
  

  # force to data.frame
  dfTrain <- as.data.frame(dfTrain)
  dfTest <- as.data.frame(dfTest)
  
  # get the independent and dependent for test 
  mX_test <-  model.matrix(formula, dfTest)[,-1]
  vY_test_raw <- dfTest[, as.character(formula)[2]]
  
  # train a model
  Trained_Model_param <- BoostLogit_train(formula,dfTrain, iM_classifiers, iZmax, sBaseLearner)
  
  # get predictions on test set
  lResult_test <- BoostLogit_test(Trained_Model_param, mX_test)
  
  mConfusionTable <- table(lResult_test$predicted_category, vY_test_raw)
  iARI <- adjustedRandIndex(vY_test_raw,lResult_test$predicted_category)
  
  lResult_test$ConfusionTable <- table(lResult_test$predicted_category, vY_test_raw)
  lResult_test$ARI <- iARI
  
  return(lResult_test)

}




smoke_train <- smoking[folds$Resample1, ]
smoke_test <- smoking[-folds$Resample1, ]



Result <- BoostLogit(intention_to_smoke ~., dfTrain = smoke_train, dfTest = smoke_test, iM_classifiers = 5, iZmax = 3, sBaseLearner = "WLS")
Result
BoostLogit_train(intention_to_smoke ~., df = smoke_test, iM_classifiers = 5, iZmax = 3, sBaseLearner = "WLS")


t <- LB(dummy_cols(as.matrix(smoke_train$intention_to_smoke))[-1], as.matrix(smoke_train[,-5]), M=10)






# compare to mboost
ctrl <- boost_control(mstop = 10)
Mboost_result = gamboost(as.factor(intention_to_smoke) ~ ., data = smoke_train, family = Binomial(), control = ctrl, baselearner = "bols")



