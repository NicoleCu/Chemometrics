# support vector regression
library(readxl)
library(e1071)
library(Metrics)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51'))
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]
col<-colnames(DafniaMagna)


## Tuning SVR model by varying values of maximum allowable error and cost parameter
#Tune the SVM model
#The legend on the right displays the value of Mean Square Error (MSE). 
#MSE is defined as (RMSE)2 and is also a performance indicator.

tModelsvm <- tune(svm, Daphnia ~ ., data = DaphniaMagna, kernel = 'radial',
                  ranges=list(cost = c(1:10,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
                  tunecontrol = tune.control(cross = 40))  
#Print optimum value of parameters
print(OptModelsvm)
#Plot the performance of SVM Regression model
plot(OptModelsvm)
#The best model is the one with lowest MSE. 
summary(OptModelsvm)

## Select the best model out of 1100 trained models and compute RMSE

#Find out the best model
BstModel <- OptModelsvm$best.model
show(BstModel)

BstModel$index          #numbers of support vectors
#Predict Y using best model


st <- lm(pred ~ DafniaMagna[,17])
plot(DafniaMagna[,17], pred, 
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DafniaMagna[,17], pred,
     labels = (row.names(DafniaMagna)),
     cex = 0.6, pos = 4, col = "gray")
grid()
#Calculate RMSE and other parameters of the best model 
RMSE <- rmse(pred, DafniaMagna[,17])
RMSECV <- sqrt*
bias <- bias(DafniaMagna[,17], pred)
slope <- st$coefficients[2]
show(rbind(RMSE, bias, slope))
