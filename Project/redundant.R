
# dataframe of only country codes in euruope
dfISO_continent_EUR <- data.frame(con = codelist$continent, iso3c = codelist$iso3c) %>%
  filter(con == "Europe")

df_EUR <- df %>%
  filter(loc %in% dfISO_continent_EUR$iso3c)

# Ensure reproducability
set.seed(1)

# get sample to work off (for now)
sample_id <- sample(21000, 1000)

# pick small sample to check out the data
sample_df <- df_EUR[sample_id,]

colnames(sample_df)
drop_cols <- c("acctstdq",
               "add1", 
               "add2", 
               "add3", 
               "add4",
               "addzip",
               "cik",
               "city",
               "compstq",
               "conm",
               "conml",
               "county",
               "curcdq",
               "datacqtr",
               "datafqtr",
               "ein",
               "exchg",
               "fax",
               "fdateq",
               "fic",
               "ggroup",
               "gind",
               "naics",
               "incorp",
               "isin",
               "pdateq",
               "phone",
               "prican",
               "prirow",
               "priusa",
               "rp",
               "sedol",
               "sic",
               "staltq",
               "state",
               "weburl",
               "capr2q",
               "capr3q",
               "spcsrc", 
               "ipodate",
               "busdesc")

# later look at ipodate

# sample without dropped variables
sample_df_dropped <- sample_df %>% select(-drop_cols)

# dependent and independent
dfY <- sample_df_dropped$capr1q
dfX <- sample_df_dropped %>% select(-capr1q)

# create banktype variable
dummy_banktype <- dummy_columns(dfX$gsubind)
colnames(dummy_banktype) <- c("Diversified banks", "Regional banks")

# create dummy for location
dummy_loc <- dummy_columns(dfX$loc)

colSums(is.na(dfX))/nrow(dfX)
# remove columns with only na's
dfX_naRemoved <- dfX[, colSums(is.na(dfX))/nrow(dfX) < 0.1]


