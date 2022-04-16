library(readxl)
library(mdatools)

removeY <- function(y, x, row.names.remove, scale = FALSE){
  y <- data.frame(y[!(row.names(y) %in% row.names.remove), ])
  row.names(y) <- row.names(x)
  colnames(y) <- "y"
  return(y)
}



VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx", 
                                 sheet = "analysis_fisheri", range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]

x <- VibrioFisheri[,1:23]
y <- data.frame(VibrioFisheri[,24])
row.names(y) <- row.names(x)

#do a regression model, cv = 1 full cross-validation

m <- pls(x, y, 5, cv = 1)
plot(m)
summary(m)

#look for the outliers
plotPredictions(m, show.label = T)
abline(a = 0, b = 1)

# remove outliers x1
row.names.remove <- c('w24')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 5, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

# remove outliers x1
row.names.remove <- c('w6')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

# remove outliers x3
#row.names.remove <- c('w17')
#x <- x[!(row.names(x) %in% row.names.remove),]
#y <- removeY(y, x, row.names.remove)

#new model
#m <- pls(x, y, 10, cv = 1)
#plot(m)
#show(summary(m))
#look for the outliers
#plotPredictions(m, show.labels = T)
#abline(a = 0, b = 1)

#r <- predict(m, xt, yt)
#plot(r)
#plotPredictions(r, show.labels = T)
#abline(a = 0, b = 1)
#par(mfrow = c(1,1))
#plotWeights(m, 2)
#plotRegcoeffs(m,2)
#summary(r)
#round(r$y.pred[, 2,1],2)  
#show(VibrioFisheri[21:26, 24])
