
calc_mostCorrelated <- function(mX, iThreshold){
  
  # create correlation matrix for the independent variables
  corr <- cor(mX)
  
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  
  #drop perfect correlations
  corr[corr == 1] <- NA 
  
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  
  #remove the NA values from above 
  corr <- na.omit(corr) 
  
  # only show high correlations (both negative and positive)
  high_corr <- corr %>%
    filter(abs(Freq)>iThreshold) %>%
    distinct(Var1, .keep_all = TRUE) %>%
    arrange(-Freq)
  
  return(high_corr)
  
  
}