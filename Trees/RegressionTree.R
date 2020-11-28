
# Calculates which split is best according to the SSE criterion
sse_var <- function(x, y) {
  
  # get all possible splits
  splits <- sort(unique(x))
  
  # save here the SSE
  sse <- c()
  
  # go through all possible splits
  for (i in seq_along(splits)) {
    
    # get a split
    split <- splits[i]
    
    # get the dependent variable in in the possible split regions
    vY1 <- y[x < split]
    vY2 <- y[x >= split]
    
    # calculate sse 
    sse[i] <- sum((vY1 - mean(y[x < split]))^2) + sum((vY2  - mean(y[x >= split]))^2) 
  }
  # get where to split
  split_at <- splits[which.min(sse)]
  return(c(sse = min(sse), split = split_at))
}


RegressionTree <- function(formula, df, minObs, maxDepth) {
  
  # coerce to data.frame
  data <- as.data.frame(df)
  
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
          
          # subset data according to the filter
          this_data <- subset(df, eval(parse(text = tree_info[j, "Rule"])))

          # get the design matrix
          X <- model.matrix(formula, this_data)
        
        # in the first case, the split is over all the data
        } else {
          this_data <- df
        }
        
        # estimate splitting criteria
        splitting <- apply(X, MARGIN = 2, FUN = sse_var, y = y)
        
        # get the min SSE
        tmp_splitter <- which.min(splitting[1,])

        # define maxnode
        mn <- max(tree_info$Node_id)
        
        # paste rules for the upcoming split
        tmp_filter <- c(paste(names(tmp_splitter), ">=", splitting[2,tmp_splitter]),paste(names(tmp_splitter), "<", splitting[2,tmp_splitter]))
        
        # Check if the splitting rule has already been invoked to prevent appying the same rule twice
        split_here  <- !sapply(tmp_filter,
                               FUN = function(x,y) any(grepl(x, x = y)),
                               y = tree_info$Rule)
        
        # If empty, add the splitting rules
        if (!is.na(tree_info[j, "Rule"])) {
          tmp_filter  <- paste(tree_info[j, "Rule"], 
                               tmp_filter, sep = " & ")
        } 
        
        # get the number of observations in current node
        tmp_nobs <- sapply(tmp_filter,
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
                               Rule = tmp_filter,
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
  
  # calculate fitted values
  leafs <- tree_info[tree_info$Status == "LEAF", ]
  fitted <- c()
  for (i in seq_len(nrow(leafs))) {
    # extract index
    ind <- as.numeric(rownames(subset(df, eval(parse(text = leafs[i, "Rule"])))))
    # estimator is the mean y value of the leaf
    fitted[ind] <- mean(y[ind])
  }
  
  # return everything
  return(list(tree = tree_info, fit = fitted, formula = formula, data = data))
}


test_tree <- RegressionTree(cigs ~  age + educ, df = Smoke, minObs = 10, maxDepth = 5)
test_tree$tree

k = 5
k>10
