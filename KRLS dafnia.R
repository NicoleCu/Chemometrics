library(KRLS)
library(readxl)
library(Metrics)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51'))
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
row.names.remove <- c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% row.names.remove),]

col <- colnames(DafniaMagna)
s <- sample(1:39, 30)
DM_train = DafniaMagna[s,col]
DM_test = DafniaMagna[-s,col]

## approximate function using KRLS. Automatically compute lambda, sigma controls the bandwidth or the
# degree of spread in the Gaussian kernel weight

out <- krls(y = DM_train[,17], x = DM_train[,1:16])
summary(out)

# get fitted values and ses
fit <- predict(out,newdata = DM_test[,1:16], se.fit=TRUE)

# results
RMSEP <- rmse(DM_test[, 'Daphnia'], fit$fit)
RMSEP
