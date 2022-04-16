# support vector regression
library(readxl)
library(e1071)
library(Metrics)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51'))
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]
col <- colnames(DafniaMagna)

# remove outliers 
row.names.remove <- c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% row.names.remove),]

test <- c('7s','13s','11s','6b','15b','23s','19b','18b','20s')
col<-colnames(DafniaMagna)
DM_test <- DafniaMagna[test, col]
DM_train <- DafniaMagna[!(row.names(DafniaMagna) %in% test),]


## Tuning SVR model by varying values of maximum allowable error and cost parameter
#Tune the SVM model
#The legend on the right displays the value of Mean Square Error (MSE). 
#MSE is defined as (RMSE)2 and is also a performance indicator.

OptModelsvm <- tune(svm, Daphnia ~ ., data = DM_train, kernel = 'radial',
                    ranges=list(cost = c(1,10,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
                    tunecontrol = tune.control(cross = 30), scale = F)
print(OptModelsvm)
#Plot the performance of SVM Regression model
#plot(OptModelsvm)
#The best model is the one with lowest MSE. 

## Select the best model out of trained models and compute RMSE
#Find out the best model
BstModel <- OptModelsvm$best.model
show(BstModel)

#Predict Y using best model
pred <- predict(BstModel, DM_test)
st <- lm(pred ~ DM_test[,17])
plot(DM_test[,17], pred, 
     main = 'SVM',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DM_test[,17], pred,
     labels = (row.names(DM_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmsep <- rmse(DM_test[,17], pred)
bias <- bias(DM_test[,17], pred)
slope <- st$coefficients[2]

#R2 <- 1 - (sum((DM_test[,17] - pred)^2))/(sum((DM_test[,17] - mean(DM_test[,17]))^2))
show(rbind(rmsep, bias, slope))
