

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies,
               readxl,
               nortest,
               dsmle,
               tidyverse,
               ggplot2, 
               parallel,
               devtools,
               gbm,
               mboost,
               gtools, 
               splitstackshape,
               rpart,
               rpart.plot)

# add dsmle package
install.packages("../KRR/dsmle_1.0-4.tar.gz", repos = NULL, type="source")

# add repository to make mclapply() run in parallel (only necessary on windows)
install_github('nathanvan/parallelsugar')
library(parallelsugar)

# add helper functions
source("helperfunctions.R")

# load in dataframe, add variable names
dfCrime <- read.table("data/crime_unnormalized.txt", sep = ",")
dfVarnames <- read_xlsx("data/Varnames_unnormalized.xlsx")
colnames(dfCrime) <- dfVarnames$Variable 

# remove variables that we will not consider for the regression, such as the county number
vRemoveVar <- c("communityname", "state", "countyCode","communityCode", "fold")
dfCrime_removeIrrelevant <- dfCrime %>% select(-vRemoveVar)
  
# remove the other potential dependent variables
vIndex_OtherDependent <- c((ncol(dfCrime_removeIrrelevant) -  18):(ncol(dfCrime_removeIrrelevant)-2))
dfCrime_removeDependent <- dfCrime_removeIrrelevant[,-vIndex_OtherDependent]

# turn all to numeric
dfCrime_num <- data.frame(apply(dfCrime_removeDependent, 2, as.numeric))

# remove rows with na in dependent (this only removes a handful of rows)
dfCrime_cleanDependent <- dfCrime_num[!is.na(dfCrime_num$ViolentCrimesPerPop),]
dfCrime_clean <- dfCrime_cleanDependent[,colSums(is.na(dfCrime_cleanDependent)) == 0]

# Log the dependent variable
vY <- as.matrix((dfCrime_clean$ViolentCrimesPerPop))
vY_logged <- log(1 + vY)

# define independent variables, check which ones are highly correlated
mX <- as.matrix(dfCrime_clean %>% select(-ViolentCrimesPerPop))
dfCorr_mX <- calc_mostCorrelated(mX, 0.7)

# scale the independent variables
mX_scaled <- scale(mX)

# check histograms
hist(vY, breaks = 30)
hist(vY_logged, breaks = 30)


###########################
#  Kernel ridge regression
###########################

train_krr_fold <- function(iFold_index, lFold,vY, mX,sKernel_type, lParam){
  
  #Split the data according to the folds
  vTest_indexes = lFold[[iFold_index]]
  vTrain_indexes = -lFold[[iFold_index]]
  
  # define train and test set for y and x
  vY_train <- vY[vTrain_indexes]
  mX_train <- mX[vTrain_indexes,]
  vY_test <- vY[vTest_indexes]
  mX_test <- mX[vTest_indexes,]
  
  # define parameter grid based on type of kernel
  if(sKernel_type == "linear"){
    
    # run the kernel ridge regression for specified kernel and parameters on the fold
    lmResult <- krr(vY_train, mX_train, kernel.type = sKernel_type, lambda = lParam$lambda, scale = FALSE, center = FALSE)
    
  }else if(sKernel_type == "RBF"){
    
    # run the kernel ridge regression for specified kernel and parameters on the fold
    lmResult <- krr(vY_train, mX_train, kernel.type = sKernel_type, lambda = lParam$lambda, kernel.RBF.sigma = lParam$sigma, scale = FALSE, center = FALSE)
    
  }else if(sKernel_type == "nonhompolynom"){
    
    # run the kernel ridge regression for specified kernel and parameters on the fold
    lmResult <- krr(vY_train, mX_train, kernel.type = sKernel_type, lambda = lParam$lambda, kernel.degree = lParam$degree, scale = FALSE, center = FALSE)
  }
  
  # run the kernel ridge regression for specified kernel and parameters on the fold
  lmResult <- krr(vY_train, mX_train, kernel.type = sKernel_type, lambda = lParam$lambda)
  
  # get the predictions for this fold
  vYhat_fold <- predict(lmResult, mX_test)

  return(data.frame(yhat = vYhat_fold, indexes = vTest_indexes))

}


gridSearch_krr <- function(vY, mX, sKernel_type, vLambda, iK, vSigma = NA, vDegree = NA){
  
  # define parameter grid based on type of kernel
  if(sKernel_type == "linear"){
    
    dfParamGrid <- expand.grid(vLambda)
    colnames(dfParamGrid) <- c("lambda")
    
  }else if(sKernel_type == "RBF"){
    
    dfParamGrid <- expand.grid(vLambda, vSigma)
    colnames(dfParamGrid) <- c("lambda", "gamma")
    
  }else if(sKernel_type == "nonhompolynom"){
    
    dfParamGrid <- expand.grid(vLambda, vDegree)
    colnames(dfParamGrid) <- c("lambda", "degree")
  }
  
  folds = createFolds(vY, k = iK, list = TRUE, returnTrain = FALSE)
  
  # create list with all paramters
  ldfParamGrid =  split(dfParamGrid, seq(nrow(dfParamGrid)))
  
  lResult_grid <- lapply(ldfParamGrid, cv_krr, vY = vY, mX = mX, sKernel_type = sKernel_type, lFold = folds)
  lResult_grid_sorted <- do.call("rbind", lResult_grid)
  
  return(lResult_grid_sorted)
  
}


cv_krr <- function(lParam, vY, mX,sKernel_type, lFold){
  
  # get result of parallel 
  lResult_cv = mclapply(1:length(lFold), train_krr_fold,lFold = lFold,sKernel_type = "linear", lParam = lParam, mX = mX_scaled, vY = vY_logged, mc.cores = detectCores())
  
  # turn list of lists to dataframe
  dfResult_cv <- do.call("rbind", lResult_cv) %>% arrange(indexes)
  
  # get average RMSE
  residualsSquared <- (dfResult_cv$yhat - vY)^2
  fAvgMSE = (sum(residualsSquared)/length(vY))
  
  lParam$avg_mse = fAvgMSE
  
  print("Done with checking the following parameters: ")
  print(lParam)

  return(lParam)
  
}


calc_PermuationTest <- function(sVarName,sModelType,dfTrain,dfTest, vY, iN, ...){
  
  # get var to permute
  vVarToPermute <- dfTrain[[sVarName]]

  # permute the var 
  vPermutedVar <- vVarToPermute[sample(iN)]
  
  # add permuted var back in the training set
  dfTrain[[sVarName]]<- vPermutedVar
  
  if(sModelType == "boosted tree"){
    
    objModel_permuted = gbm(data = dfTrain,...)
    
  }else{
    stop("Pick for sModelType one of 'boosted tree' or 'kernel ridge regression' ")
  }

  # carry out prediction with permuted var
  vPrediction <- predict(objModel_permuted, dfTest)
  
  # get mean squared error
  MSE <- mean((vY - vPrediction)^2)
  
  return(MSE)

}





create_PermutationTest_df <- function(sModelType, dfTrain, dfTest, vY, ...){
  
  # get n of observations
  iN_train <- nrow(dfTrain)
  
  lMSEresult <- mclapply(colnames(mX), calc_PermuationTest,sModelType=sModelType, dfTrain= dfTrain, dfTest = dfTest, vY = vY, iN = iN_train, ..., mc.cores = detectCores())
  
  
  dfMSEresult = data.frame(var_removed = colnames(mX), MSE_without = unlist(lMSEresult))
  
  return(dfMSEresult)
}









vLambda = 10^seq(-1, 5, length.out = 10)
vSigma = 10^seq(0.5,2, length.out = 5)
vDegree = seq(2,4,by=1)

result_krr_linear <- gridSearch_krr(vY_logged, mX_scaled, "linear", vLambda, 10)
result_krr_rbf <- gridSearch_krr(vY_logged, mX_scaled, "RBF", vLambda,10, vSigma = vSigma )
result_krr_nonhompolynom <- gridSearch_krr(vY_logged, mX_scaled, "nonhompolynom",vLambda, iK = 10, vDegree = vDegree)


###########################
#  Boosted Regression Trees
###########################


# finds best minimum of observations and interaction depth in the boosted tree
gridSearch_TreeBoost <- function(objFormula, vInteraction_depth, vMinobs, dfData, iK, n_trees =500){
  
  # create dataframe of all parameters to check
  dfParamGrid <- expand.grid(vInteraction_depth, vMinobs)
  colnames(dfParamGrid) <- c("InteractionDepth", "Minobs")
  
  # create list with all parameters
  ldfParamGrid =  split(dfParamGrid, seq(nrow(dfParamGrid)))

  # apply gbm to all combinations of the parameters
  lResult_grid <- lapply(ldfParamGrid ,gbm_inGridsearch,objFormula = objFormula, sDistribution = "gaussian", dfData = dfData, iK = iK, n_trees = n_trees)
  
  # save results of the data.frame
  lResult_grid_cleaned <- lapply(lResult_grid, create_df_TreeResult)
  
  dfResult_grid_cleaned <- do.call("rbind", lResult_grid_cleaned)
  
  return(dfResult_grid_cleaned)

}

# ensures we can apply gbm in lapply
gbm_inGridsearch <- function(lParam, objFormula, sDistribution, dfData, iK, n_trees){
  result <- gbm(formula = objFormula, distribution = sDistribution, data = dfData, cv.folds = iK, n.trees = n_trees, interaction.depth = lParam$InteractionDepth, n.minobsinnode = lParam$Minobs)
  return(result)
}

# creates df to analyse tree results
create_df_TreeResult <- function(objTreeResult){
  dfTreeResult <- data.frame(n_tree_cv_best = which.min(objTreeResult$cv.error), cv.error = objTreeResult$cv.error[which.min(objTreeResult$cv.error)], Minobs = objTreeResult$n.minobsinnode, InteractionDepth = objTreeResult$interaction.depth)
  return(dfTreeResult)
}

# clean data to use it in gbm() function
dfCrime_clean_gbm <- data.frame(cbind(mX_scaled, vY_logged))
colnames(dfCrime_clean_gbm)[ncol(dfCrime_clean_gbm)] <- "ViolentCrimesPerPop_logged"

# define interaciton depth and min. observaitonst to check
vInteraction_depth <- c(1,3,5)
vMinobs <- c(5,10,50)

# get result of the gridsearch
result_gridsearch_TreeBoost = gridSearch_TreeBoost(ViolentCrimesPerPop_logged ~., vInteraction_depth, vMinobs, dfCrime_clean_gbm, 10)

## set the seed to make your partition reproducible
set.seed(123)

# get train and test indexes in dataset
fSample_size <- floor(0.7 * nrow(dfCrime_clean_gbm))
vInd_train <- sample(seq_len(nrow(dfCrime_clean_gbm)), size = fSample_size)
# create train and test datasets
dfTrain <- dfCrime_clean_gbm[vInd_train, ]
dfTest <- dfCrime_clean_gbm[-vInd_train, ]

# train and predict test
Model_train_regressionTree <- gbm(formula = ViolentCrimesPerPop_logged ~., distribution = "gaussian", data = dfTrain, n.trees = 350, interaction.depth = 1, n.minobsinnode = 10)
prediction_test_regressionTree <- predict(Model_train_regressionTree, dfTest)

# get errors
vResidualsSquared <- (prediction_test_regressionTree - dfTest$ViolentCrimesPerPop_logged)^2
MSE_regressionTree <- mean(vResidualsSquared)
MSE_regressionTree

dfPermutationTest <- create_PermutationTest_df(sModelType = 'boosted tree', dfTrain, dfTest, dfTest$ViolentCrimesPerPop_logged,  formula = ViolentCrimesPerPop_logged ~., distribution = "gaussian", n.trees = 350, interaction.depth = 1, n.minobsinnode = 10)

dfPermutationTest$perc_worse <- (dfPermutationTest$MSE_without - MSE_regressionTree) / MSE_regressionTree




# create df to add on variables and analyse the residuals
dfResidualAnalysis <- data.frame(ViolentCrimsPerPop_logged = dfTest$ViolentCrimesPerPop_logged)


mean(vResidualsSquared)
dfResidualAnalysis$residualsSquared_RegressionTree <- vResidualsSquared

# transform back - check why loss of info
dfResidualAnalysis$ViolentCrimesPerPop_backTransformed <- exp(dfResidualAnalysis$ViolentCrimsPerPop_logged) - 1

# get in which quantile the observations are
vQuantileCrime_test <- quantcut(dfResidualAnalysis$ViolentCrimesPerPop_backTransformed ,q = 10)
dfResidualAnalysis$QuantileCrime <- vQuantileCrime_test

# show per quantile what % contributes to the errors
dfSummary_errorPerQuantile <- dfResidualAnalysis %>%
    group_by(QuantileCrime) %>%
    summarise(sumResidualsSquared = sum(residualsSquared_RegressionTree) / sum(dfTest$residualsSquared_RegressionTree), 
              sampleSize = n())


dfTrain$flooredViolentCrimesPerPop_logged<- floor(dfTrain$ViolentCrimesPerPop_logged)

dfSummary_flooredViolentCrimesPerPop_logged <- dfTrain %>% 
  group_by(flooredViolentCrimesPerPop_logged) %>%
  summarise(sampleSize = n())

dfStrat_train <- stratified(dfTrain, "flooredViolentCrimesPerPop_logged", 10)

dfStrat_train$flooredViolentCrimesPerPop_logged <- NULL


model_stratTrain_regressionTree <-  gbm(formula = ViolentCrimesPerPop_logged ~., distribution = "gaussian", data = dfStrat_train, n.trees = 85, interaction.depth = 1, n.minobsinnode = 10)
prediction_test_regressionTree_strat <- predict(model_stratTrain_regressionTree, dfTest)

vResidualsSquared_stratTrain <- dfTest$ViolentCrimesPerPop_logged - prediction_test_regressionTree_strat
mean(vResidualsSquared_stratTrain)


colnames(dfCrime_clean_gbm)
control = rpart.control(maxdepth = 2)
help(rpart)
quicktree = rpart(ViolentCrimesPerPop_logged ~.,data = dfCrime_clean_gbm, control = control)
rpart.plot(quicktree, type = 0)
