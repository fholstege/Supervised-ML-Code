

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies,
               readxl,
               dplyr,
               nortest)



# load in dataframe, add variable names
dfCrime <- read.table("data/crime_unnormalized.txt", sep = ",")
dfVarnames <- read_xlsx("data/Varnames_unnormalized.xlsx")
colnames(dfCrime) <- dfVarnames$Variable

# remove variables that we will not consider for the regression, such as the county number
vRemoveVar <- c("communityname", "state", "countyCode","communityCode", "fold")
dfCrime_removeIrrelevant <- dfCrime %>% 
                  select(-vRemoveVar)

# remove the other potential dependent variables
vIndex_OtherDependent <- c((ncol(dfCrime_removeIrrelevant) -  18):(ncol(dfCrime_removeIrrelevant)-2))
dfCrime_removeDependent <- dfCrime_removeIrrelevant[,-vIndex_OtherDependent]

# turn all to numeric
dfCrime_num <- data.frame(apply(dfCrime_removeDependent, 2, as.numeric))

# remove rows with na in dependent
dfCrime_cleanDependent <- dfCrime_num[!is.na(dfCrime_num$ViolentCrimesPerPop),]
dfCrime_clean <- dfCrime_cleanDependent[,colSums(is.na(dfCrime_cleanDependent)) == 0]


# remove one observation with the number of violent crimes at 0
# justification; will not distort much but makes logging impossible
dfCrime_clean <- dfCrime_clean[!dfCrime_clean$ViolentCrimesPerPop == 0,]

# get the dependent and indepdent variables
vY <- as.matrix((dfCrime_clean$ViolentCrimesPerPop))
mX <- as.matrix(dfCrime_clean %>% select(-ViolentCrimesPerPop))

# scale the independent variables
mX_scaled <- scale(mX)

# check histograms
hist(vY, breaks = 30)
hist(log(vY), breaks = 30)

