

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,
               fastDummies,
               countrycode)


# load data
df <- read.csv("Compustat_BankData.csv")

# Ensure reproducability
set.seed(1)

# get sample to work off (for now)
sample_id <- sample(65000, 1000)


# pick small sample to check out the data
sample_df <- df[sample_id,]