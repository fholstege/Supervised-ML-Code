

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

# define independent variables
x_smoke_train = smoke_train[,-5]
x_smoke_test = smoke_test[,-5]


##################################
# Own implementation
##################################


# define parameters:
mX <- x_smoke_train
vY <- smoke_train$intention_to_smoke



BoostLogit <- function(mX, vY, iN_classifiers){
  
  # define the initial weights
  iN = nrow(mX)
  vW <- rep(1/iN,iN)
  
  # define the initial probabilities 
  vP <- rep(0.5, iN)

  # define the probabilities as chosen by the ensemble
  vEnsembleProbabilities = rep(0, iN)
  
  # update integer with each weak classifier added
  iClassifier = 0
  
  # max value for z
  iZmax = 4
  
  # Continue updating until as many weak classifiers as defined in iN_classifiers
  while(iClassifier < iN_classifiers){
    iClassifier = iClassifier+ 1

    # Get sample weights based on probabilities - define minimum value to prevent bugs
    # Will weigh observations closer to 0.5 (e.g. uncertain ones) more
    vWeights = pmax(vP * (1-vP),1e-24)
    vWeights = vWeights/sum(vWeights)
    
    # define the 'response' variable vZ
    vZ = ifelse(vY == 1, 1/vP, ifelse(vY == 0, -1/(1-vP), (vY - vP)/vWeights))
    
    # define min and max of it
    vZ <- pmin(vZ, iZmax)
    vZ <- pmax(vZ, -iZmax)


    # get probabilites from a weak classifier
    vWeakClassifierProbabilities <- lm(vZ ~ as.matrix(mX), weights = vWeights)
    
    # update the probabilities of the ensemble with this weak classifier
    vEnsembleProbabilities = vEnsembleProbabilities +  (0.5 * vWeakClassifierProbabilities$fitted.values)

    # update current probabilities for the next set of weights
    vP <- exp(vEnsembleProbabilities)/(exp(vEnsembleProbabilities) + exp(-vEnsembleProbabilities))
    
  }
  
  resultObj = list(predicted_label = sign(vEnsembleProbabilities), 
                   predicted_raw = vEnsembleProbabilities)

  
  return(resultObj)
  
}

Result = BoostLogit(mX,vY, 10)

ctrl <- boost_control(mstop = 10)
Mboost_result = gamboost(as.factor(intention_to_smoke) ~ ., data = smoke_train, family = Binomial(), control = ctrl, baselearner = "bols")


dfResult = data.frame(our_label = Result$predicted_label, mboost_label = sign(Mboost_result$predict()), our_p = Result$predicted_raw, mboost_p = Mboost_result$predict())
dfResult$same <- dfResult$our_label == dfResult$mboost_label
                      