library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(randomForest)

x <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                              range ='B2:L23'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                          range ='N2:N23'))
rownames(y) <- seq(0,20, 1)
colnames(x) <- seq(1,11, 1)

set.seed(41)
s <- sample(1:21, 6)
x_test <-  x[s,]
y_test <- y[s,]
names(y_test) <- s
x_train = x[-s,]
y_train <- y[-s,]
xy_train <- cbind(x_train, y_train)

###############################################################################
###############################################################################
# PLS for original data
m_daf <- pls(x, y, 10, cv = 1)
plot(m_daf)
show(summary(m_daf))
plotYVariance(m_daf, type  = 'h', show.labels = TRUE)

# look for the outliers
plotPredictions(m_daf, show.labels = T)
abline(a = 0, b = 1)
plotRegcoeffs(m_daf, show.labels = T)

################################################################################
### SVM
library(e1071)

OptModelsvm_new <- tune(svm, y_train ~., data=xy_train, kernel = 'radial',
                        ranges=list(cost = c(1:10,100),  
                        epsilon = seq(0.01, 0.15, 0.01),
                        gamma = c(0.0001, 0.0005, 0.0001)), 
                        tunecontrol = tune.control(cross = 15), scale = F)  
print(OptModelsvm_new)

BstModel_new <- OptModelsvm_new$best.model
show(BstModel_new)

#Predict Y using best model
pred_new <- predict(BstModel_new, x_test) 
st <- lm(pred_new ~ y_test)
plot(y_test, pred_new, 
     main = 'SVM for original data',
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


###############################################################################
## Random forest for original data
library(randomForest)

set.seed(123)
RF <- randomForest(y_train ~., data = xy_train, xtest = x_test, ytest = y_test,  
                       mtry=5, ntree = 500, nodesize = 3,
                       importance = T)
show(RF)
plot(RF)

st <- lm(RF[["test"]][["predicted"]] ~ y_test)

plot(y_test, RF[["test"]][["predicted"]], 
     main = 'Random Forest for original data',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_test, RF[["test"]][["predicted"]],
     labels = (names(y_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((RF[["test"]][["predicted"]] - y_test)^2)
RMSEP <- rmse(y_test, RF[["test"]][["predicted"]])
bias <- bias(y_test, RF[["test"]][["predicted"]])
slope <- st$coefficients[2]
R2 <- 1 - (sum((y_test - RF[["test"]][["predicted"]])))/(sum((y_test - mean(y_test))^2))
show(rbind(RMSEP, MSE, bias, slope))


################################################################################
## KRLS for the original data
library(KRLS)

out <- krls(x_train, y_train, vcov = T, derivative = T, sigma = 14)
summary(out)


# get fitted values and ses
KRLS <- predict(out,newdata = x_test, se.fit=TRUE)


# results
st <- lm(KRLS$fit ~ y_test)
#RMSECV = rmse(out[["fitted"]], DM_train[,17])    #там fitted value а не предсказанные
RMSEP <- rmse(y_test, KRLS$fit)
bias <- bias(y_test, KRLS$fit)
slope <- st$coefficients[2]
show(rbind(RMSEP, bias, slope))


plot(y_test, KRLS$fit, main = 'KRLS for original data',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_test, KRLS$fit, labels = (names(y_test)), 
     cex = 0.6, pos = 4, col = "gray")
grid()






###############################################################################
###############################################################################
##PLS after removing outliers

x_n <- x[!(rownames(x) %in% c('11','19', '9', '18')), ]
y_n <- y[!(rownames(y) %in% c("11","19", "9", '18')), ]

#set.seed(41)
#c <- sample(1:17, 4)
#x_ntest <-  x[c,]
#y_ntest <- y[c,]
#names(y_ntest) <- c
#x_ntrain = x[-c,]
#y_ntrain <- y[-c,]
#xy_ntrain <- cbind(x_ntrain, y_ntrain)

m_daf_n <- pls(x_n, y_n, 10, cv = 1)
plot(m_daf_n)
show(summary(m_daf_n))
plotYVariance(m_daf_n, type  = 'h', show.labels = TRUE)
plotPredictions(m_daf_n, show.labels = T)
abline(a = 0, b = 1)
plotRegcoeffs(m_daf_n, show.labels = T)

#rem <- c(3, 10,11, 8, 6)
m_daf_n <- pls(x_n[,c(2,5,7)], y_n, 10, cv = 1)
plot(m_daf_n)
show(summary(m_daf_n))
plotYVariance(m_daf_n, type  = 'h', show.labels = TRUE)
plotPredictions(m_daf_n, show.labels = T)
abline(a = 0, b = 1)
plotRegcoeffs(m_daf_n, show.labels = T)


################################################################################
### SVM after removing outliers

OptModelsvm_new <- tune(svm, y_ntrain ~., data=xy_ntrain, kernel = 'radial',
                        ranges=list(cost = c(1:10,100),  
                                    epsilon = seq(0.01, 0.15, 0.01),
                                    gamma = c(0.0001, 0.0005, 0.0001)), 
                        tunecontrol = tune.control(cross = 15), scale = F)  
print(OptModelsvm_new)

BstModel_new <- OptModelsvm_new$best.model
show(BstModel_new)

#Predict Y using best model
pred_new <- predict(BstModel_new, x_ntest) 
st <- lm(pred_new ~ y_ntest)
plot(y_ntest, pred_new, 
     main = 'SVM after removing outliers',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_ntest, pred_new,
     labels = (names(y_ntest)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmsep <- rmse(y_ntest, pred_new)
bias <- bias(y_ntest, pred_new)
slope <- st$coefficients[2]
R2 <- 1 - (sum((y_ntest - pred_new)^2))/(sum((y_ntest - mean(y_ntest))^2))
show(rbind(rmsep, bias, slope, R2))

################################################################################
#### RF after removing outliers


set.seed(123)
RF_new <- randomForest(y_ntrain ~., data = xy_ntrain, xtest = x_ntest, ytest = y_ntest,  
                   mtry=5, ntree = 500, nodesize = 2,
                   importance = T)
show(RF_new)
plot(RF_new)

st <- lm(RF_new[["test"]][["predicted"]] ~ y_ntest)

plot(y_ntest, RF_new[["test"]][["predicted"]], 
     main = 'Random Forest for original data',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_ntest, RF_new[["test"]][["predicted"]],
     labels = (names(y_ntest)),
     cex = 0.6, pos = 4, col = "gray")
grid()

MSE <- mean((RF_new[["test"]][["predicted"]] - y_ntest)^2)
RMSEP <- rmse(y_ntest, RF_new[["test"]][["predicted"]])
bias <- bias(y_ntest, RF_new[["test"]][["predicted"]])
slope <- st$coefficients[2]
R2 <- 1 - (sum((y_ntest - RF_new[["test"]][["predicted"]])))/(sum((y_ntest - mean(y_ntest))^2))
show(rbind(RMSEP, MSE, bias, slope))

###############################################################################
### KRLS after removing outliers


out_new <- krls(x_ntrain, y_ntrain, vcov = T, derivative = T, sigma = 20)
summary(out_new)


# get fitted values and ses
KRLS_new <- predict(out_new,newdata = x_ntest, se.fit=TRUE)


# results
st <- lm(KRLS_new$fit ~ y_ntest)
#RMSECV = rmse(out[["fitted"]], DM_train[,17])    #там fitted value а не предсказанные
RMSEP <- rmse(y_ntest, KRLS_new$fit)
bias <- bias(y_ntest, KRLS_new$fit)
slope <- st$coefficients[2]
show(rbind(RMSEP, bias, slope))


plot(y_ntest, KRLS_new$fit, main = 'KRLS after removing outliers',
     ylim = c(0, 100), xlim = c(0, 100), pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(y_ntest, KRLS_new$fit, labels = (names(y_ntest)), 
     cex = 0.6, pos = 4, col = "gray")
grid()



