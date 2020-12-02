
# function to split continuous variables
split_continuous <- function(x,y, metric, minObs){
  
  # get all possible splits
  splits <- sort(unique(x))
  #final_possible_splits <- sort(unique(x))
  
  # save here the SSE
  results <- c()
  
  # go through all possible splits
  for (i in seq_along(splits)) {
    
    # get a split
    split <- splits[i]
    
    # get the dependent variable in in the possible split regions
    vY1 <- y[x < split]
    vY2 <- y[x >= split]
    
    # test if minimum of observations fulfilled
    if(min(length(vY1), length(vY2)) < minObs){
      results[i] <- NA

    # based on metric, make split
    }else if(metric == "sse"){
      
      results[i] <-calc_sse(vY1, vY2)
      
    }else if (metric =="Gini"){
      
      results[i]<- calc_Gini(vY1, vY2)
    }
  }
  

  # if na, remove 
  if(all(is.na(results))){
    return(c(metric = NA, split = NA))
  }else{
    

    # get where to split and return
    split_at <- splits[which.min(results)]
    results <- na.omit(results)
    
    return(c(metric = min(results), split = split_at))
    
  }
  
  
  
}

# split if categorical variable
split_categorical <- function(x, y, metric){
  
  # get all categories
  categories <- unique(x)

  # create all the possible combis
  nVar <- length(categories)
  
  # go over all the potential combinations
  nCombi <- floor(nVar/2)
  tot_combis <-  list()
  
  # create list of combinations
  for(i in c(1:nCombi)){
    
    set_combis <- as.matrix(combn(1:nVar, i))
    tot_combis[[i]]<- set_combis
    
  }
  
  
  # save gini scores here
  results <- c()
  counter = 1
  
  # check every combi
  for (set_combi in seq_along(tot_combis)){
    for(combi in set_combi){
      
      vY1 <- y[x %in% combi]
      vY2 <- y[! x %in% combi]
      
      # test if minimum of observations fulfilled
      if(min(length(vY1), length(vY2)) < minObs){
        next
      }
      
      # based on metric, make the split
      if (metric == "sse"){
        
        results[counter] <- calc_sse(vY1, vY2)
        
      }else if (metric == "Gini"){
        
        results[counter] <- calc_Gini(vY1, vY2)
      }
      
      counter = counter + 1
      
    }
    
  }
  
  # make split and return
  split_at <- categories[which.min(results)]
  
  return(c(metric = min(results),split = split_at ))
                     
}

# function to calculate gini index
calc_Gini <- function(vY1, vY2){
  
  Gini_R1 <- 1 - sum(prop.table(table(vY1))^2)
  Gini_R2 <- 1 - sum(prop.table(table(vY2))^2)
  
  prop_R1 <- length(vY1)/(length(vY1) + length(vY2))
  prop_R2 <- 1- prop_R1
  
  Gini <- sum(c(prop_R1, prop_R2) * c(Gini_R1, Gini_R2))
  
  return(Gini)
}

# function to calculate SSE
calc_sse <- function(vY1, vY2){
  
  sse <- sum((vY1 - mean(vY1))^2) + sum((vY2  - mean(vY2))^2) 
  return(sse)
  
}

# general function - recognizes what split to make based on dependent and independent variable
splitFunction <- function(x,y, minObs){

  if(is.numeric(x)){
    
    if(is.numeric(y)){
      
      result <- split_continuous(x,y,"sse", minObs)
    }else{
      result <- split_continuous(x,y,"Gini", minObs)
    }
    
  }else{
    
    if(is.numeric(y)){
      result <- split_categorical(x,y,"sse", minObs)
    }else{
      result <- split_categorical(x,y,"Gini", minObs)
    }
  }
  
  return(result)
}

# function to build the regression tree
RegressionTree <- function(formula, df, minObs, maxDepth) {
  
  # coerce to data.frame
  df <- as.data.frame(df)
  
  # handle formula
  formula <- terms.formula(formula)
  
  # get the independent variable matrix
  X <- model.matrix(formula, df)
  
  # get the dependent variable
  y <- df[, as.character(formula)[2]]
  
  # While TRUE, while loop continues 
  do_splits <- TRUE
  
  # Create output data.frame with splitting rules and observations
  # This dataframe also works to guide which next need to be taken
  tree_info <- data.frame(Node_id = 1,    
                          Variable = NA,  
                          Index = NA,     
                          Obs = nrow(df), 
                          Rule = NA,     
                          Status = "SPLIT", 
                          stringsAsFactors = FALSE)

  
  # keep splitting until there are only leafs left
  while(do_splits) {
    
    # Get indexes of which nodes need have to be splitted
    vIndex_split <- which(tree_info$Status == "SPLIT")

    # go through all to calculate
    for (j in vIndex_split) {
      
        # when not at the root node, select the subset based on the previous splits
        if (!is.na(tree_info[j, "Rule"])) {
          
          # subset data according to the rule
          this_data <- subset(df, eval(parse(text = tree_info[j, "Rule"])))

          # get the design matrix
          X <- model.matrix(formula, this_data)
        
        # in the first case, the split is over all the data
        } else {
          this_data <- df
        }
        
        # estimate splitting criteria
        splitting <- apply(X, MARGIN = 2, FUN = splitFunction, y = y, minObs = minObs)
        
        if(all(is.na(splitting))){
          do_splits <- FALSE
        }else{

        # get the min metric (SSE/GINI)
        tmp_splitter <- which.min(splitting[1,])

        # define maxnode
        mn <- max(tree_info$Node_id)
        
        # paste rules for the upcoming split
        tmp_rule <- c(paste(names(tmp_splitter), ">=", splitting[2,tmp_splitter]),paste(names(tmp_splitter), "<", splitting[2,tmp_splitter]))
        

        # Check if the splitting rule has already been invoked to prevent appying the same rule twice
        split_here  <- !sapply(tmp_rule,
                               FUN = function(x,y) any(grepl(x, x = y)),
                               y = tree_info$Rule)
        
        # If empty, add the splitting rules
        if (!is.na(tree_info[j, "Rule"])) {
          tmp_rule  <- paste(tree_info[j, "Rule"], 
                               tmp_rule, sep = " & ")
        } 
        
        # get the number of observations in current node
        tmp_nobs <- sapply(tmp_rule,
                           FUN = function(i, x) {
                             nrow(subset(x = x, subset = eval(parse(text = i))))
                           },
                           x = this_data)  
        
        # insufficient minObs for split
        if (any(tmp_nobs <= minObs)) {
          split_here <- rep(FALSE, 2)
        }

        # create children data frame
        children <- data.frame(Node_id = c(mn+1, mn+2),
                               Variable =names(tmp_splitter),
                               Index = splitting[2,tmp_splitter],
                               Obs = tmp_nobs,
                               Rule = tmp_rule,
                               Status = rep("SPLIT", 2),
                               row.names = NULL)[split_here,]
        

        # overwrite state of current node
        tree_info[j, "Status"] <- ifelse(all(!split_here), "LEAF", "PARENT")
        
        # bind everything
        tree_info <- rbind(tree_info, children)
        
        # check if there are any open splits left
        do_splits <- !all(tree_info$Status != "SPLIT")
        
        # check if at max depth
        k <- length(unique(tree_info$Index))-1
        
        # if no splits left, stop
        if(all(tree_info$Status != "SPLIT")){
          
          do_splits <- FALSE
        
        # also stop if max depth reached
        }else if(k>=maxDepth){
          
          do_splits <- FALSE
          
        }
      }
    }
  }
  
  # if any splits left, set to leaf since done
  if(any(tree_info$Status == "SPLIT")){
    
    tree_info[which(tree_info$Status == "SPLIT"),]$Status <- "LEAF"
  }
  
  # get the leafs to get fitted values
  leafs <- tree_info[tree_info$Status == "LEAF", ]
  fitted <- c()

  for (i in seq_len(nrow(leafs))) {
    
    # extract index
    ind <- as.numeric(rownames(subset(df, eval(parse(text = leafs[i, "Rule"])))))
    
    # if numeric dependent, estimator is the mean y value of the leaf
    if(is.numeric(y)){
      
      fitted[ind] <- mean(y[ind])
    
      # otherwise predict the category 
    }else{
      fitted[ind] <- names(sort(table(y[ind]), decreasing = TRUE)[1])
    }


  }
  
  # return everything
  return(list(tree = tree_info, fit = fitted, formula = formula, data = data))
}

# load data
load("smoking.Rdata")


# our tree
Our_smoke_tree <- RegressionTree(intention_to_smoke ~  age + alcohol_per_month, df = smoking, minObs = 10, maxDepth = 5)
Our_smoke_tree$tree

# rpart tree
Rpart_smoke_tree <- rpart(intention_to_smoke ~  age + alcohol_per_month, data = smoking,
                 control = rpart.control(minsplit = 10, maxdepth = 5))
plot(as.party(Rpart_smoke_tree))

# difference created by the fact that Rpart makes split at half points (59.5) instead of 60
# it also uses ANOVA instead of SSE, likely to lead to small difference 

