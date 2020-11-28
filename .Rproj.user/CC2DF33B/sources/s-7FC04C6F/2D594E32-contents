
# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sampleSelection,
               dplyr,
               rpart,
               partykit)


splitFunction <- function(vY, vX, minObs, type = "continuous"){
  
  vPossible_splits <- unique(vX)
  vOptimal_split <- c(NA,Inf)

  for(fSplitValue in vPossible_splits){
    
    if(type == "continuous"){
      
      vR1 <- vX[vX<=fSplitValue]
      vR2 <- vX[vX>fSplitValue]
      
    }else if (type == "categorical"){
      
      vR1 <- vX[vX==fSplitValue]
      vR2 <- vX[vX==fSplitValue]
      
    }
    
    fMin_obs_split <- min(length(vR1), length(vR2))
    
    if(fMin_obs_split < minObs){
      
      next
      
    }
    
    vY1 <- vY[match(vR1, vX)]
    vY2 <- vY[match(vR2, vX)]
    
    Loss <- sum((vY1 - mean(vY1))^2) + sum((vY2 - mean(vY2))^2)
    
    if(vOptimal_split[2]>Loss){
      vOptimal_split <- c(fSplitValue, Loss)
    }
    
  }

    
  return(vOptimal_split)
  
}

TreeBuilder <- function(mX, vY, maxDepth, minObs, k, node_id, bLeaf, dfSplits = NA){
  
  fBenchmark_loss <- sum((vY - mean(vY))^2)
  

  while(k < maxDepth & bLeaf == FALSE ){
    
    k = k + 1
    
    vOptimal_split_branch <- c(NA, Inf, 0)
    
    for(iBranch_index in 1:NCOL(mX)){

      vBranch = mX[,iBranch_index]

      vOptimal_split <- splitFunction(vY, vBranch, minObs) 
      
      if(is.na(vOptimal_split[2])){
        next 
      }

      if(vOptimal_split[2]<vOptimal_split_branch[2]){
        vOptimal_split_branch <- c(vOptimal_split,iBranch_index)
      }
      
    }


    # make the split
    vIndex_R1 <- mX[,vOptimal_split_branch[3]]<= vOptimal_split_branch[1]
    vIndex_R2 <- mX[,vOptimal_split_branch[3]] > vOptimal_split_branch[1]

    if(NCOL(mX)==1){
      bLeaf <- TRUE
    }else{
    
    mX_R1 <- as.matrix(mX[vIndex_R1,])
    mX_R2 <- as.matrix(mX[vIndex_R2,])
    
    vY_R1 <- vY[vIndex_R1]
    vY_R2 <- vY[vIndex_R2]
    
    fNew_loss <- vOptimal_split_branch[2]
    fImprovement <- (fNew_loss - fBenchmark_loss)/fBenchmark_loss
    
    print("Level")
    print(k)
    print("Node id")
    print(node_id)
    print("Split")
    print(vOptimal_split_branch)

    dfSplits[iNodecounter,] <- c(colnames(mX)[vOptimal_split_branch[3]], k, node_id, fImprovement, vOptimal_split_branch[1] )
    iNodecounter <<- iNodecounter + 1
    #<- rbind(dfSplits, c(colnames(mX)[vOptimal_split_branch[3]], k, node_id, fImprovement, vOptimal_split_branch[1] ))
    dfSplits <- dfSplits

    # set column names
    colnames(dfSplits) <- c("variable", "level", "id", "improvement", "index")
    
    R1_Tree <- TreeBuilder(mX <- mX_R1, vY_R1, maxDepth = maxDepth - 1, minObs, k, node_id = 1, bLeaf, dfSplits = dfSplits)
    R2_Tree <- TreeBuilder(mX <- mX_R2, vY_R2, maxDepth = maxDepth - 1, minObs, k, node_id = 2, bLeaf, dfSplits = dfSplits)
    
    }
    
    
  }
}

data(Smoke)

mSmoke <- Smoke %>%
  select(educ, age, cigs) %>%
  as.matrix()


mX_smoke <- mSmoke[,-3]
vY_smoke <- mSmoke[,3]

dfSplits <<- data.frame(variable = character(), level = numeric(), id = numeric(), improvement = numeric(), index = numeric())
iNodecounter <<- 1


TreeBuilder(mX_smoke, vY_smoke, maxDepth = 3, minObs = 5, k=0, node_id = 1, bLeaf = FALSE, dfSplits)


smokers <- rpart(cigs ~  age + educ, data = Smoke,
                      control = rpart.control(minsplit = 5, maxdepth = 3))

smokers$splits
plot(as.party(smokers))
