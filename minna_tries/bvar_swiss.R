
## Housekeeping ----------------------------------------------------------------

library(tidyverse)
library(BVAR)


# import & prep data
getwd()
data <- read.csv("submission/data/data_quarterly.csv")
metadata <- read.csv("submission/data/metadata_quarterly.csv")

# check metadata to get vars
metadata <- metadata[, c("ts_key","variable")]

# inspect
# head(data)
# view(metadata) 

# select vars
data <- data[, c("date","gdp", "cpi", "uroff", "srate")]
# view(data)

# transform to be stationary
data <- column_to_rownames(data, var = "date")
data <- fred_transform(data, codes = c(5, 5, 5, 1), lag = 4) # transform to be stationary

# estimate model
# bayesian VAR(1): y_t = c + phi_1* y_{t-1} + epsilon_t
# with y_t = [gdp, cpi, unempl, fedfunds]
# parameters: (1+4)*4 = 20 coefficients