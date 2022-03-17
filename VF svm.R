# support vector regression
library(readxl)
library(e1071)
library(Metrics)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_fisheri",
                                       range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]


col <- colnames(VibrioFisheri)
s <- sample(1:26, 20)
VF_train = VibrioFisheri[s,col]
VF_test = VibrioFisheri[-s,col]

## Tuning SVR model by varying values of maximum allowable error and cost parameter
#Tune the SVM model
#The legend on the right displays the value of Mean Square Error (MSE). 
#MSE is defined as (RMSE)2 and is also a performance indicator.

OptModelsvm <- tune(svm, EC50 ~., data=VF_train, kernel = 'radial',
                    ranges=list(cost = c(0.01, 0.1, 1:10, 100), gamma = c(0.001, 0.01, 0.1:0.5, 1, 2)))  
#Print optimum value of parameters
print(OptModelsvm)
#Plot the performance of SVM Regression model
plot(OptModelsvm)
#The best model is the one with lowest MSE. 


## Select the best model out of 1100 trained models and compute RMSE

#Find out the best model
BstModel <- OptModelsvm$best.model
show(BstModel)

#Predict Y using best model
pred <- predict(BstModel, VF_test)
st <- lm(pred ~ VF_test[,'EC50'])
plot(VF_test[,'EC50'], pred, 
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(VF_test[,'EC50'], pred,
     labels = (row.names(VF_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmse <- rmse(VF_test[,'EC50'], pred)
bias <- bias(VF_test[,'EC50'], pred)
slope <- st$coefficients[2]

show(rbind(rmse, bias, slope))
