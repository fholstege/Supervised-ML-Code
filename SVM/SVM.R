
setwd("C:/Users/flori/OneDrive/Documents/GitHub/Supervised ML Code/SVM")


# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies,
               SVMMaj,
               dplyr,
               caret, 
               mclust, 
               VeryLargeIntegers,
               matlib,
               latex2exp,
               reshape2)

# create parameters for hinge errors
create_hinge_param <- function(mY, mQ, z, hinge = "absolute", epsilon = 1e-08, k_huber = 1){
  
  # when the hinge is absolute
  if(hinge == "absolute"){
  
    # add epsilon threshold to make sure its invertible
    m <- abs(1 - z)
    m <- m * (m > epsilon) + epsilon * (m <= epsilon)
  
    # calculate a and b using m
    a <- .25 * m^-1
    b <- mY * (a + 0.25)
    
    return(list(a = a,
                b = b))

  }else if (hinge == "quadratic"){
    
    # for quadratic hinge, no need to calculate a (always identity)
    m <- z* (z > 1) + (1 >= z)
    b <- m * mY
    
    return(list(b = b))
    
  }else if(hinge == "huber"){
    
    a <- 0.5 * (k_huber + 1)^-1
    
    b <- (mY == -1) * ((mQ <= -1) * (a * mQ) + 
            (mQ >= k_huber) * (a * mQ - 0.5) + 
            (mQ > -1 & mQ < k_huber) * (-a)) +
         (mY == 1) * ((mQ <= -k_huber) * (0.5 + a * mQ) + 
            (mQ >= 1) * (a * mQ) + 
            (mQ > -k_huber & mQ < 1) * (a))
    
    c <- (mY == -1) * ((mQ <= -1) * (a * mQ^2) + 
            (mQ >= k_huber) * (1 - (k_huber + 1)/2 + a * mQ^2) + 
            (mQ > -1 & mQ < k_huber) * (a)) +
         (mY == 1) * ((mQ <= -k_huber) * (1 - (k_huber + 1)/2  + a * mQ^2) + 
            (mQ >= 1) * (a * mQ^2) + 
            (mQ > -k_huber & mQ < 1) * (a))
    
    return(list(a = a, 
                b = b, 
                c = c))
  }
  
  
}




calc_loss <- function(mW, lambda, z, mY = NA, hinge = "absolute", mQ = NA, k_huber = NA ){
  
  mWtW = t(mW) %*% mW
  penalty = lambda * mWtW
  
  if(hinge == "absolute"){

    vloss <- (1 > z) * (1 - z)
    
  }else if(hinge == "quadratic"){
    
    vloss <- (1 > z) * (1 - z)^2
    
  }else if(hinge == "huber"){
    
    vloss <- (mY == 1) * (mQ <= -k_huber) * (1 - mQ -(k_huber + 1)/2) +
      (mY == 1) * (mQ > -k_huber) * (0.5 * (k_huber+1)^-1 * pmax(0,1 - mQ) ^2) +
      (mY == -1) * (mQ <= k_huber) * (0.5 * (k_huber+1)^-1 * pmax(0,1 + mQ) ^2) +
      (mY == -1) * (mQ > k_huber) * (mQ +1 -(k_huber + 1)/2)
    
  }else{
    stop("Not given an known hinge error")
  }

  return(sum(vloss) + penalty)
}

svm_mm <- function(mY, mX, hinge = "absolute", lambda = 10, epsilon = 1e-08, k_huber = 3){
  
  # set initial weights and constant. From these, create initial v = c + w
  vW = runif(ncol(mX)-1, -1, 1)
  fConstant = 0.0
  vV_current = as.matrix(c(fConstant, vW))
  
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
    mCurrent_q = mX %*% vV_current
    mCurrent_z = mCurrent_q * mY

    # get parameters given the z (absolute and quadratic) or q and y (huber)
    lHinge_param_current = create_hinge_param(mY = mY, 
                                              mQ = mCurrent_q,
                                              z = mCurrent_z,
                                              hinge = hinge, 
                                              epsilon = epsilon, 
                                              k_huber = k_huber)
    # all hinges need a b, so can define here
    b = lHinge_param_current$b

    # quick update if quadratic
    if(hinge == "quadratic"){
      
      vV_update = Z %*% b
    
    }else {
      
      # define A and b
      A = diag(as.vector(lHinge_param_current$a),n)

      # get updated v
      vV_update = solve(t(mX) %*% A %*% mX + lambda * P, t(mX) %*% b)
      
    }
    
    # get weights in previous and next 
    mW_current = tail(vV_current,-1)
    mW_update = tail(vV_update, -1)

    # get new prediction (q) and z
    mNew_q = mX %*% vV_update
    mNew_z = mNew_q * mY
    

    # calculate new, and previous loss
    fCurrent_loss = calc_loss(mW = mW_current,
                              mQ = mCurrent_q,
                              z = mCurrent_z,
                              mY = mY,
                              lambda = lambda, 
                              hinge = hinge, 
                              k_huber = k_huber)
    

    fNew_loss = calc_loss(mW = mW_update, 
                          mQ = mNew_q,
                          z = mNew_z,
                          mY = mY,
                          lambda = lambda, 
                          hinge = hinge, 
                          k_huber = k_huber)
    
    # calculate improvement
    stepscore <- (fCurrent_loss - fNew_loss)/fCurrent_loss

    # check: if all predicted correctly, turns to NaN since divided by zero
    if (is.na(stepscore)){
      stepscore = 0
    }
    
    # move to next iteration
    k = k + 1
    vV_current = vV_update
    
  }

  # get predicted category
  mY_hat = sign(mNew_q)
  
  
  
  # gather results 
  mConfusionTable <- table(mY, mY_hat)
  
  # create confusion matrix and calculate accuracy
  fAccuracy <- sum(mY_hat == mY)/length(mY)
  
  fARI <- adjustedRandIndex(mY, mY_hat)
  
  # return results object
  lResults = list(v = vV_update,
                  mY = mY,
                  mX = mX, 
                  loss = fCurrent_loss,
                  q = mNew_q,
                  yhat = mY_hat,
                  Accuracy = fAccuracy,
                  ConfusionTable = mConfusionTable,
                  ARI = fARI)
  
  return(lResults)

  
  
}


svm_mm_cv <- function(mX, mY, lambda, folds, hinge = "absolute", k_huber = NA, epsilon = 1e-08, metric = "ARI"){
  
  fTest_metric <- 0
  
  #Perform k fold cross validation
  for(i in 1:length(folds)){
    
    #Split the data according to the folds
    vTest_id = folds[[i]]
    vTrain_id = -folds[[i]]

    # define train and test set for y and x
    mY_train <- mY[vTrain_id]
    mX_train <- mX[vTrain_id,]
    mY_test <- mY[vTest_id]
    mX_test <- mX[vTest_id,]
    
    # get result from training set
    lResult = svm_mm(mY = mY_train, mX = mX_train, hinge = hinge, lambda = lambda, epsilon = epsilon, k_huber = k_huber )
    mV_train <- lResult$v
    
    # calculate the predicted y for this fold
    mQ_test <- mX_test %*% mV_train
    mY_hat <- sign(mQ_test)
    
    if(metric == "ARI"){
      
      # get the confusion matrix
      fAdjRand <- adjustedRandIndex(mY_test, mY_hat)
      fTest_metric_fold <- fAdjRand
      
    }else if (metric == "misclassification"){

      # calculate misclassification
      tab <- table(mY_test, mY_hat)
      fMisclassification <- 1-sum(diag(tab))/sum(tab)
      fTest_metric_fold <- fMisclassification
      
    }

    # add to calculate average
    fTest_metric <- fTest_metric + fTest_metric_fold

  }
  
  # calculate average
  avg_test_metric = fTest_metric/length(folds)
  
  lResult_cv <- list(lambda = lambda,
                     metric = avg_test_metric)
  
  return(lResult_cv)
  
  
}

svm_mm_gridsearch <- function(mX, mY, lambda, hinge = "absolute", k, k_huber = NA, epsilon = 1e-08, metric = "ARI"){
  
  
  if (hinge == "huber"){
    
    mParamgrid = expand.grid(lambda, k_huber)
  }else {
    
    mParamgrid = expand.grid(lambda)
  }

  # create k equally size folds
  folds = createFolds(mY, k = k, list = TRUE, returnTrain = FALSE)
  
  mParamgrid$metric <- NA

  # iterate over the grid
  for(i in 1:nrow(mParamgrid)){

    # select parameters from the grid
    param <- mParamgrid[i,]
    
    print(param[1,c(1,2)])
    
    # test these with k-fold
    cv_result <- svm_mm_cv(mX = mX,mY= mY, lambda =param[1,1], k_huber = param[1,2],folds = folds, hinge = hinge, epsilon = epsilon, metric = metric)

    print(cv_result$metric)
    
    #save the result
    mParamgrid$metric[i] <- cv_result$metric
  }
  
  return(mParamgrid)
}
  
  


#### section 1: data pre-processing

# ensure that we can reproduce the code
set.seed(123)

# load bank data
bank <- read.csv("bank-additional.csv", sep = ";")
#load("bank.Rdata")

# get dependent variable and transform to 1, -1
mY <- bank$y
mY_num <- as.matrix(ifelse(mY == "yes", 1, -1))

# get all possible independent variables
mX = bank[,-ncol(bank)]

# create dummy variables 
df_toDummy = as.matrix(data.frame(
  mX$job,
  mX$marital, 
  mX$education, 
  mX$default, 
  mX$housing, 
  mX$loan, 
  mX$contact,
  mX$month, 
  mX$day_of_week,
  mX$poutcome,
  mX$pdays))

dummy_vars <-  dummy_columns(df_toDummy)
dummy_vars <- dummy_vars[,(1+ncol(df_toDummy)):ncol(dummy_vars)]

# select numeric variables, scale these. Leave out duration for realistic predictions
numeric_vars <- as.matrix(data.frame(
  mX$age, 
  mX$campaign, 
  mX$cons.price.idx, 
  mX$cons.conf.idx, 
  mX$nr.employed, 
  mX$emp.var.rate,
  mX$euribor3m))

# combine numeric and dummy variables
mX_var <- as.matrix(cbind(1, scale(numeric_vars),dummy_vars))

# pick a random sample of 1000
sample_id <- sample(4000, 1000)
sample_x = mX_var[sample_id,]
sample_y = mY_num[sample_id]


### section 2: compare a single svm_mm with svm_maj

# We get the same predictions (see confusion matrix), and exactly the same loss, for each type of error
# but still different beta's? most likely something to do with the transformation of the data...

# our implementation
result = svm_mm(sample_y, sample_x, lambda = 0.1, hinge = "huber", k_huber = 2)
result$v
result$loss
result$ConfusionTable
result$Accuracy
result$ARI

#svmmaj implementation
result_svmmaj <- svmmaj(sample_x,sample_y, lambda = 10, scale = "none",hinge = "quadratic")
result_svmmaj$beta 
result_svmmaj$loss
result_svmmaj$q




### section 3: compare cross validation results, and create cv plots

# check out the following parameters
vLambda =10^seq(5, -3, length.out= 10) # 10
vk_huber = seq(0,3, by = 1)

# k-fold of 20 to have enough info, but not make it computationally too intensive
k = 20

# cv comparison
result_svm_mm_cv <- svm_mm_gridsearch(sample_x, sample_y, k = k, lambda = vLambda, k_huber = vk_huber, hinge = "huber", metric = "misclassification")
result_svmmaj_cv <- svmmajcrossval(sample_x,sample_y, search.grid = list(lambda = vLambda) , k = k,scale = "none",hinge = "huber")

# data cleaning of our cv results
colnames(result_svm_mm_cv) <- c("lambda", "misclassification")

# finds same lambda for quadratic
ggplot(data = result_svm_mm_cv, aes(x = log(lambda), y = misclassification, col = "red"))+
  geom_line()+
  geom_line(data = result_svmmaj_cv$param.grid, aes(y = loss, col = "blue")) +
  theme_minimal()


## now create plots per type of error for optimal lambda
result_cv_quadratic <- svm_mm_gridsearch(sample_x, sample_y, k = k, lambda = vLambda, hinge = "quadratic", metric = "ARI")
result_cv_huber <- svm_mm_gridsearch(sample_x, sample_y, k = k, lambda = vLambda, k_huber = vk_huber, hinge = "huber", metric = "ARI")
result_cv_absolute <- svm_mm_gridsearch(sample_x, sample_y, k = k, lambda = vLambda, hinge = "absolute", metric = "ARI")

# change column names
colnames(result_cv_quadratic) <- c("lambda", "ARI")
colnames(result_cv_absolute) <- c("lambda", "ARI")
colnames(result_cv_huber) <- c("lambda", "k_huber", "ARI")


df_cv_compare <- data.frame(lambda = result_cv_absolute$lambda, abs_ARI = result_cv_absolute$ARI, quad_ARI = result_cv_quadratic$ARI)
df_cv_compare <- melt(df_cv_compare, id.vars = "lambda")
df_cv_compare

# finds same lambda for quadratic
ggplot(data = df_cv_compare, aes(x = log(lambda), y = value, col = variable))+
  geom_line() + 
  theme_bw() + 
  labs(y = "Average ARI"
    
  ) + 
  theme(plot.title =element_text(size=20, face = "plain",hjust = 0.5),
        axis.title=element_text(size=15, face = "plain")) +
  scale_color_discrete(name = "Type of error", 
                      labels = c("Absolute", "Quadratic")
  ) 



optimal_lambda_quadratic <- result_cv_quadratic[which.max(result_cv_quadratic$ARI),]$lambda
optimal_lambda_absolute <- result_cv_absolute[which.max(result_cv_absolute$ARI),]$lambda


### section 4: test the optimal parameters from cv on a train and test set

## 70% of the sample size
fTrain_size <- floor(0.7 * nrow(mX_var))

# get train indexes for the dataset
vTrain_id <- sample(nrow(mX_var), size = fTrain_size)

# train and test split
mX_train <- mX_var[vTrain_id, ]
mX_test <- mX_var[-vTrain_id, ]

mY_train <- mY_num[vTrain_id]
mY_test <- mY_num[-vTrain_id]


# train for the beta's
result_quadratic <- svm_mm(mY = mY_train,mX= mX_train, lambda = optimal_lambda_quadratic, hinge = "quadratic")
result_absolute <- svm_mm(mY = mY_train,mX= mX_train, lambda = optimal_lambda_absolute, hinge = "absolute")



analyse_svm_result <- function(mY, mQ, plot_title = NaN){
  
  # get predicted category
  mY_hat <- sign(mQ)
  
  # get several key statistics
  resultOverview_quadratic <- confusionMatrix(as.factor(mY_test), as.factor(mY_hat))
  
  # get confusion matrix
  mConfusionMatrix <- resultOverview_quadratic$table
  
  # get ARI
  adjRand <- adjustedRandIndex(mY_test, mY_hat)
  
  dfComparePlot  <- data.frame(q = mQ, mY = mY_test) 
  ComparePlot <- ggplot(data = dfComparePlot, aes(x = q, fill = as.factor(mY))) +
    geom_histogram(bins = 50,alpha = 0.7) +
    labs(
      title = plot_title,
      y = "Count",
      x = TeX("$\\hat{q}$")
    ) + 
    scale_fill_discrete(name = "Result", 
                        labels = c("Did not take subscription (-1)", "Took subscription (1)")
                        )+
    scale_fill_manual(values=c("red", "blue")) + 
    theme_bw() +
    theme(plot.title =element_text(size=20, face = "plain",hjust = 0.5),
          axis.title=element_text(size=15, face = "plain"))
    
  
  
  return(list(mY_hat = mY_hat,
              ConfusionMatrix <- mConfusionMatrix,
              ARI = adjRand,
              ComparePlot = ComparePlot
  ))
  
  
}

# calculate the predicted y for this fold

### first, quadratic
mQ_test_quadratic <- mX_test %*% result_quadratic$v

analysis_quadratic <- analyse_svm_result(mY_test, mQ_test_quadratic, plot_title = "Results for quadratic error" )
analysis_quadratic$ComparePlot


 


# get said beta's, and check prediction



### section 5: test the different kernels 

result_svmmaj_cv_linKernel <- svmmajcrossval(sample_x, sample_y, 
                                             search.grid = list(lambda = 1.0e+15),
                                             k = k, scale = "none", 
                                             hinge = "quadratic", 
                                             kernel = polydot,
                                             degree = 1)
result_svmmaj_cv_linKernel






### Section 6: miscellaneous 


create_show_df <- function(vk_huber_show){
  
  df_show <- data.frame( q =seq(-3,3,by=0.2) )
  
  
  df_show$q <- seq(-3,3,by=0.2)
  col_count = 1
  
  for(ik_huber_show in vk_huber_show){
    
    vloss_plusOne <-(mQ_show <= -ik_huber_show) * (1 - mQ_show -(ik_huber_show + 1)/2) +
      (mQ_show > -ik_huber_show) * (0.5 * (ik_huber_show+1)^-1 * pmax(0,1 - mQ_show) ^2) 
    
    vloss_minusOne <-  (mQ_show <= ik_huber_show) * (0.5 * (ik_huber_show+1)^-1 * pmax(0,1 + mQ_show) ^2) +
      (mQ_show > ik_huber_show) * (mQ_show +1 -(ik_huber_show + 1)/2)
    
    col_names_huber = c(paste0("plusOne_", ik_huber_show),paste0("minusOne_", ik_huber_show) )
    
    df_show$plusOne <- vloss_plusOne
    df_show$minusOne <- vloss_minusOne
    
    colnames(df_show)[-(1:col_count)] <- col_names_huber
    col_count = col_count + 2
    
    
  }
  

  return(df_show)

}

vk_huber_show = c(0)
df_show <- create_show_df(vk_huber_show)

vAbsError_plusOne <- (1 > mQ_show) * (1 - mQ_show)
vAbsError_minusOne <- (1 > -mQ_show) * (1 - -mQ_show)

vQuadError_plusOne <- (1 > mQ_show) * (1 - mQ_show)^2
vQuadError_minusOne <- (1 > -mQ_show) * (1 - -mQ_show)^2


df_show$plusOne_abserror <- vAbsError_plusOne
df_show$minusOne_abserror <- vAbsError_minusOne
df_show$plusOne_quaderror <- vQuadError_plusOne
df_show$minusOne_quaderror <- vQuadError_minusOne



df_show_plot <- melt(df_show, measure.vars = colnames(df_show)[-1] )

ggplot(data = df_show_plot, aes(x = q, y = value, col = variable)) +
  geom_line(size = 2) + 
  theme_bw() + 
  labs(y = "Loss") + 
  lims(y = c(0,5)) + 
    scale_color_discrete(name = "Function",
                        labels = c("+1 huber loss, k = 0",
                                   "-1 huber loss, k = 0",
                                   "+1 absolute loss",
                                   "-1 absolute loss",
                                   "+1 quadratic loss",
                                   "-1 quadratic loss")
    )


