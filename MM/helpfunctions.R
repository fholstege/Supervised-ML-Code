
# Author: Floris Holstege
# Purpose: Helper functions that aid with analysing the results of several models
# Date: 27/10/2020

# dfCompare
# Parameters:
#   result1: dataframe, results from first model
#   result2: dataframe, results from second model
#   cols: vector, columns names for dataframe


dfCompare = function(result1, result2, cols){
  
  dfResult <- data.frame(result1, result2)
  colnames(dfResult) <- cols
  
  return(dfResult)
  
}



createOverviewdf <- function(result){
  
  row_sub = apply(result$Beta, 1, function(row) all(row !=0 ))
  
  df <- data.frame(result$Beta[row_sub,], 
                   result$SignificanceResults$pval[row_sub],
                   result$SignificanceResults$t[row_sub])
  colnames(df) <- c("Coefficient", "P-value", "T-val with n-p-1 degree of freedom")
  
  return(df)
  
}

