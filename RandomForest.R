library(randomForest)
library(readxl)
library(Metrics)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51'))
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
row.names.remove <- c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% row.names.remove),]

test <- c('7s','13s','11s','6b','15b','23s','19b','18b','20s')
col<-colnames(DafniaMagna)
DM_test <- DafniaMagna[test, col]
DM_train <- DafniaMagna[!(row.names(DafniaMagna) %in% test),]

# Random Forests split the data into "training" and "test" sets for you. 
## This is because Random Forests use bootstrapped data, 
# and thus, not every sample is used to build every tree. The 
## "training" dataset is the bootstrapped data and the "test" dataset is
## the remaining samples. The remaining samples are called
## the "Out-Of-Bag" (OOB) data.

RF <- randomForest(Daphnia ~., data = DafniaMagna, mtry = 5, ntree = 500, importance = T)
show(RF)
plot(RF)

for (i in 1:500) {
  if (RF[['mse']][i] == min(RF$mse)) 
    print(i)
}


#varImpPlot(RF)
#importance(RF)
# The first is based upon the mean decrease of accuracy in predictions 
#on the out of bag samples when a given variable is permuted. 
#he second is a measure of the total decrease in node impurity that results from splits over that variable, 
#averaged over all trees. In the case of regression trees, the node impurity is 
#measured by the training RSS, and for classification trees by the deviance. 


DM_pred <- predict(RF, DM_test)

st <- lm(DM_pred ~ DM_test[,17])
plot(DM_test[, 'Daphnia'], DM_pred, 
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DM_test[,17], DM_pred,
     labels = (row.names(DM_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((DM_pred - DM_test[, 'Daphnia'])^2)
RMSEP <- rmse( DM_test[, 'Daphnia'], DM_pred)
bias <- bias(DM_test[,17], DM_pred)
slope <- st$coefficients[2]
show(rbind(RMSEP, MSE, bias, slope))

