
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

