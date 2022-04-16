library(readxl)
library(mdatools)
library(Metrics)
library(randomForest)
library(e1074)
library(KRLS)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx", 
                                  sheet = "analysis_fisheri", range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]

VibrioFisheri <-  VibrioFisheri[!(row.names(VibrioFisheri) %in% 
                                    c('w24', 'w6')),]

test <- c('w12','w16','w22','w9','w8','w21')
xy_test <- VibrioFisheri[test,]
xy_train <- VibrioFisheri[!(row.names(VibrioFisheri) %in% test),]


#new model
m <- pls(xy_train[,1:23], xy_train[,24], 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T, cex.lab = 3)
abline(a = 0, b = 1)

p <- predict(m, xy_test[,1:23], xy_test[,24])
plotPredictions(p, show.labels = T)
abline(a = 0, b = 1)
summary(p)


################################################################################
#### SVM
OptModelsvm <- tune(svm, EC50 ~., data=xy_train, kernel = 'radial',
                    ranges=list(cost = c(1,10,20, 50,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
                    unecontrol = tune.control(cross = 18), scale = F)
#Print optimum value of parameters
print(OptModelsvm)

## Select the best model out of trained models and compute RMSE

BstModel <- OptModelsvm$best.model
show(BstModel)

#Predict Y using best model
pred <- predict(BstModel, xy_test)
st <- lm(pred ~ xy_test[,'EC50'])
plot(xy_test[,'EC50'], pred, main = 'SVM', 
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(xy_test[,'EC50'], pred,
     labels = (row.names(xy_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmse <- rmse(xy_test[,'EC50'], pred)
bias <- bias(xy_test[,'EC50'], pred)
slope <- st$coefficients[2]

show(rbind(rmse, bias, slope))

################################################################################
## RF
set.seed(23)
RF <- randomForest(EC50 ~., data = xy_train, xtest = xy_test[,1:23], ytest = xy_test[,24],  
                   mtry = 8, ntree = 500, nodesize = 5,
                   importance = T)
show(RF)
plot(RF)

#for (i in 1:1300) {
#  if (RF[['mse']][i] == min(RF$mse)) 
#    print(i)
#}

st <- lm(RF[["test"]][["predicted"]] ~ xy_test[,24])
plot(xy_test[, 'EC50'], RF[["test"]][["predicted"]],
     main = 'Random Forest',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(xy_test[,24], RF[["test"]][["predicted"]],
     labels = (row.names(xy_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((RF[["test"]][["predicted"]] - xy_test[, 'EC50'])^2)

RMSEP <- rmse( xy_test[, 'EC50'], RF[["test"]][["predicted"]])
#RMSECV <- rmse(xy_train[, 'EC50'], RF$predicted)

bias <- bias(xy_test[,24], RF[["test"]][["predicted"]])
slope <- st$coefficients[2]
show(rbind(RMSEP, MSE, bias, slope))

################################################################################
## KRLS



#KRLS <- tune(krls, train.x = xy_train[,1:23], train.y = as.vector(xy_train[,24]),
#             ranges=list(sigma = c(1:10,100)), 
#             tunecontrol = tune.control(cross = 18))
  
KRLS <- krls(xy_train[,1:23], as.vector(xy_train[,24]), vcov = T, derivative = T, sigma = 97)
summary(KRLS)

# get fitted values and ses
KRLS_pred <- predict(KRLS,newdata = xy_test[,1:23], se.fit=TRUE)


# results
st <- lm(KRLS_pred$fit ~ xy_test[,24])
#RMSECV = rmse(KRLS[["fitted"]], xy_train[,24])    #там fitted value а не предсказанные
RMSEP <- rmse(xy_test[, 'EC50'], KRLS_pred$fit)
bias <- bias(xy_test[,24], KRLS_pred$fit)
slope <- st$coefficients[2]
show(rbind(RMSEP, bias, slope))


plot(xy_test[, 'EC50'], KRLS_pred$fit,
     main = 'KRLS',
     ylim = c(0, 110), xlim = c(0, 110), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(xy_test[,24], KRLS_pred$fit,
     labels = (row.names(xy_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()


