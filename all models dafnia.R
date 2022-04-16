library(readxl)
library(mdatools)
library(Metrics)
library(writexl)

# data import 

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",
                                  sheet = "analysis_dafnia", range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% 
             c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')),]

# split on test and train set (original data)
test <- c('7s','13s','11s','6b','15b','23s','19b','18b','20s')
DM_test <- DafniaMagna[test,]
DM_train <- DafniaMagna[!(row.names(DafniaMagna) %in% test),]

#shapiro.test(as.integer(DM_test[,17]))

# data modification
output <- as.data.frame(DafniaMagna[,17])
rownames(output) <- rownames(DafniaMagna)
predictors <- DafniaMagna[,1:16]

#squared columns
squared <- as.data.frame((predictors^2))
colnames(squared) <- paste(colnames(predictors), '^2')

#multiplied columns
multiplied <- matrix(1, nrow = nrow(predictors), ncol = 1)
for (i in 1:(ncol(predictors)-1)) {
  for (j in (i+1):ncol(predictors))  {
    n <- as.data.frame(predictors[,i]*predictors[,j])
    colnames(n) <- paste(as.character(i), as.character(j), sep = "*")
    multiplied <- cbind.data.frame(multiplied, n)
  }
}
rm(i, j)
multiplied <- multiplied[,-1]


#divided columns
divided <- matrix(1, nrow = nrow(predictors), ncol = 1)
for (i in 1:(ncol(predictors)-1)) {
  for (j in (i+1):ncol(predictors)) {
    n <- as.data.frame(predictors[,i]/predictors[,j])
    colnames(n) = paste(as.character(i), as.character(j), sep = "/")
    divided <- cbind.data.frame(divided, n)
  }
}
rm(i, j)
divided <- divided[,-1]

# new dataset of predictors
predictors <- data.frame(cbind(scale(predictors), scale(squared), scale(multiplied), scale(divided), center = F))


# данные для выгрузки

predictors_names <- c(colnames(DafniaMagna[,1:16]), colnames(squared), colnames(multiplied), colnames(divided), 'Dafnia')
dataset <- cbind(predictors, DafniaMagna[,17])
colnames(dataset) <- predictors_names

write_xlsx(data.frame(rownames(dataset)), "~/R/Chemometrics/processed data Dafnia.xlsx")
write_xlsx(dataset, "~/R/Chemometrics/processed data Dafnia_data.xlsx")


# formating for building models
x_test <- predictors[test,]
y_test <- DM_test[,17]
names(y_test) <- row.names(x_test)
x_train <- predictors[!(row.names(DafniaMagna) %in% test),]
y_train <- as.data.frame(output[!(row.names(DafniaMagna) %in% test),])
rownames(y_train) <- rownames(x_train)
colnames(y_train) <-  'Dafnia'
xy_train <- cbind(x_train, y_train)


###############################################################################
## PLS на изначальных данных

# model
m <- pls(DM_train[,1:16], DM_train[,17], 10, cv = 1)
plot(m)
show(summary(m))

#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

r <- predict(m, DM_test[,1:16], DM_test[,17])
plot(r, ncomp = 4)
plotPredictions(r, show.label =T, ncomp=4)
abline(0,1)
summary(r)
show(rmse(r$y.pred[, r[["ncomp.selected"]],1], r$y.ref[,1]))     # по 4 компонентам ошибка самая малая,
# но для CV и Cal самая малая на 6 
plotRegcoeffs(m)

#remove noise

#DM_train <- DM_train[,-c(1,15)]
#DM_test <- DM_test[,-c(1,15)]

# model
#m_n <- pls(DM_train[,1:14], DM_train[,'Daphnia'], 10, cv = 1)
#plot(m_n)
#show(summary(m_n))

#look for the outliers
#plotPredictions(m_n, show.labels = T)
#abline(a = 0, b = 1)

#r_n <- predict(m_n, DM_test[,1:14], DM_test[,'Daphnia'])
#plot(r_n)
#plotPredictions(r_n, show.label =T)
#abline(0,1)
#summary(r_n)
#show(rmse(r_n$y.pred[, r_n[["ncomp.selected"]],1], r_n$y.ref[,1]))     # по 4 компонентам ошибка самая малая,
# но для CV и Cal самая малая на 6 
#plotRegcoeffs(m_n)




##############################################################################
## PLS на данных с преобразованием

#new model
new_m <- pls(x_train, y_train, 10, cv = 1, scale = F)
plot(new_m)
show(summary(new_m))    
show(summary(new_m))
#look for the outliers
plotPredictions(new_m, show.labels = T)
abline(a = 0, b = 1)

#plotYVariance(new_m)	

r_new <- predict(new_m, x_test, y_test)
plot(r_new)
plotPredictions(r_new, show.label =T)
abline(0,1)
summary(r_new)
show(rmse(r_new$y.pred[, r_new[["ncomp.selected"]],1], r_new$y.ref[,1]))
plotRegcoeffs(new_m)

coef <- new_m[["coeffs"]][["values"]]
coef <- coef[,new_m[["ncomp.selected"]],1]
vec <- vector()

for (i in 1:length(coef)){
  if (abs(coef[i]) < 0.3) {
    vec <-  append(vec, i)
  }
}

#new model
new_m_n <- pls(x_train[,-vec], y_train, 10, cv = 1, scale = F)
plot(new_m_n)
show(summary(new_m_n, ncomp = 3))
#look for the outliers
plotPredictions(new_m_n, show.labels = T)
abline(a = 0, b = 1)

#plotYVariance(new_m)	
r_newn <- predict(new_m_n, x_test[,-vec], y_test)
plot(r_newn)
plotPredictions(r_newn, show.label =T)
abline(0,1)
summary(r_newn)
show(rmse(r_newn$y.pred[, r_newn[["ncomp.selected"]],1], r_newn$y.ref[,1]))
plotRegcoeffs(new_m_n)



################################################################################
#### svm для первоначальных данных

library(e1071)

## Tuning SVM-R model by varying values of maximum allowable error and cost parameter
#The legend on the right displays the value of Mean Square Error (MSE). 
#MSE is defined as (RMSE)2 and is also a performance indicator.

OptModelsvm <- tune(svm, Daphnia ~ ., data = DM_train, kernel = 'radial',
                    ranges=list(cost = c(1:10,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
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
pred <- predict(BstModel, DM_test[,1:16])
st <- lm(pred ~ DM_test[,17])
plot(DM_test[,17], pred, main = 'SVM', ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DM_test[,17], pred,
     labels = (row.names(DM_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmsep <- rmse(DM_test[,17], pred)
bias <- bias(DM_test[,17], pred)
slope <- st$coefficients[2]

R2 <- 1 - (sum((DM_test[,17] - pred)^2))/(sum((DM_test[,17] - mean(DM_test[,17]))^2))
show(rbind(rmsep, bias, slope, R2))




################################################################################
## SVM для преобразованных данных


OptModelsvm_new <- tune(svm, Dafnia ~., data=xy_train, kernel = 'radial',
                    ranges=list(cost = c(1:10,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
                    tunecontrol = tune.control(cross = 30), scale = F)  
print(OptModelsvm_new)

BstModel_new <- OptModelsvm_new$best.model
show(BstModel_new)

#Predict Y using best model
pred_new <- predict(BstModel_new, x_test) 
st <- lm(pred_new ~ y_test)
plot(y_test, pred_new, 
     main = 'SVM for preprocessed data',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)

abline(0, 1)
abline(st, col = "Blue")
text(y_test, pred_new,
     labels = (names(y_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmsep <- rmse(y_test, pred_new)
bias <- bias(y_test, pred_new)
slope <- st$coefficients[2]
R2 <- 1 - (sum((y_test - pred_new)^2))/(sum((y_test - mean(y_test))^2))
show(rbind(rmsep, bias, slope, R2))


################################################################################
#### Random Forest for original data

# Random Forests split the data into "training" and "test" sets for you because 
# Random Forests use bootstrapped data, and thus, not every sample is used 
# to build every tree. The "training" dataset is the bootstrapped data and the 
# "test" dataset is the remaining samples. The remaining samples are called
# the "Out-Of-Bag" (OOB) data.
library(randomForest)

# Number randomly selected variables is mtry
# rf_default <- train(Daphnia~., data=DM_train, method='rf', importance = T,
#                  trControl=trainControl(method='cv', number=10))
#print(rf_default)

#tunedRF <- tuneRF(DM_train[,1:16], DM_train[,17], mtryStart = 2 , ntreeTry= 500, stepFactor=2, improve=0.05, plot=TRUE, doBest=FALSE)

set.seed(123)
RF <- randomForest(Daphnia ~., data = DM_train, xtest = DM_test[,1:16], ytest = DM_test[,17],  
                   mtry = 5, ntree = 300, nodesize = 5,
                   importance = T)
show(RF)
plot(RF)

#for (i in 1:1300) {
#  if (RF[['mse']][i] == min(RF$mse)) 
#    print(i)
#}

RMSEP <- rmse(DM_test[, 'Daphnia'], RF[["test"]][["predicted"]])
st <- lm(RF[["test"]][["predicted"]] ~ DM_test[,17])
plot(DM_test[, 'Daphnia'], RF[["test"]][["predicted"]],
     main = 'Random Forest',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DM_test[,17], RF[["test"]][["predicted"]],
     labels = (row.names(DM_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((RF[["test"]][["predicted"]] - DM_test[, 'Daphnia'])^2)

RMSEP <- rmse( DM_test[, 'Daphnia'], RF[["test"]][["predicted"]])
#RMSECV <- rmse(DM_train[, 'Daphnia'], RF$predicted)

bias <- bias(DM_test[,17], RF[["test"]][["predicted"]])
slope <- st$coefficients[2]
R2 <- 1 - (sum((DM_test[,17] - RF[["test"]][["predicted"]])))/(sum((DM_test[,17] - mean(DM_test[,17]))^2))
show(rbind(RMSEP, MSE, bias, slope))




################################################################################
##### RF for prepared data
set.seed(123)
RF_new <- randomForest(Dafnia ~., data = xy_train, xtest = x_test, ytest = y_test,  
                   mtry = 91, ntree = 300, nodesize = 5,
                   importance = T)
show(RF_new)
plot(RF_new)

st <- lm(RF_new[["test"]][["predicted"]] ~ y_test)

plot(y_test, RF_new[["test"]][["predicted"]], 
     main = 'Random Forest for preprocessed data',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_test, RF_new[["test"]][["predicted"]],
     labels = (names(y_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((RF_new[["test"]][["predicted"]] - y_test)^2)
RMSEP <- rmse(y_test, RF_new[["test"]][["predicted"]])
bias <- bias(y_test, RF_new[["test"]][["predicted"]])
slope <- st$coefficients[2]
R2 <- 1 - (sum((y_test - RF_new[["test"]][["predicted"]])))/(sum((y_test - mean(y_test))^2))
show(rbind(RMSEP, MSE, bias, slope))


################################################################################
## KRLS for original data
library(KRLS)

out <- krls(DM_train[,1:16], DM_train[,17], vcov = T, derivative = T, sigma = 34)
summary(out)


# get fitted values and ses
KRLS <- predict(out,newdata = DM_test[,1:16], se.fit=TRUE)


# results
st <- lm(KRLS$fit ~ DM_test[,17])
#RMSECV = rmse(out[["fitted"]], DM_train[,17])    #там fitted value а не предсказанные
RMSEP <- rmse(DM_test[, 'Daphnia'], KRLS$fit)
bias <- bias(DM_test[,17], KRLS$fit)
slope <- st$coefficients[2]
show(rbind(RMSEP, bias, slope))


plot(DM_test[, 'Daphnia'], KRLS$fit, main = 'KRLS',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(DM_test[,17], KRLS$fit, labels = (row.names(DM_test)), 
     cex = 0.6, pos = 4, col = "gray")
grid()


################################################################################
## KRLS for prepared


out_new <- krls(x_train, y_train, vcov = T, derivative = T)
#summary(out_new)


# get fitted values and ses
KRLS_new <- predict(out_new, newdata = x_test, se.fit=TRUE)


# results
st <- lm(KRLS_new$fit ~ y_test)
#RMSECV = rmse(out[["fitted"]], DM_train[,17])    #там fitted value а не предсказанные
RMSEP <- rmse(y_test, KRLS_new$fit)
bias <- bias(y_test, KRLS_new$fit)
slope <- st$coefficients[2]
show(rbind(RMSEP, bias, slope))


plot(y_test, KRLS_new$fit, main = 'KRLS for preprocessed data',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_test, KRLS_new$fit,
     labels = (names(y_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

