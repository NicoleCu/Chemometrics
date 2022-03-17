library(readxl)
library(mdatools)


removeY <- function(y, x, row.names.remove, scale = FALSE){
  y <- data.frame(y[!(row.names(y) %in% row.names.remove), ])
  row.names(y) <- row.names(x)
  colnames(y) <- "y"
  return(y)
}

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]


x <- DafniaMagna[,1:16]
y <- data.frame(DafniaMagna[,17]) 
row.names(y) <- row.names(x)
#xt <- DafniaMagna[41:50,1:16]
#yt <- data.frame(DafniaMagna[41:50,17])
#row.names(yt) <- row.names(xt)

#do a regression model, cv = 1 full cross-validation

m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
plotYVariance(m, type  = 'h', show.labels = TRUE)

#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

# remove outliers x1
row.names.remove <- c("29s", "29b")
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)
# remove outliers x2
row.names.remove <- c('28s','28b')
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
row.names.remove <- c('4s','4b')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)
# remove outliers x4
row.names.remove <- c('8b','8s')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

# remove outliers x4
row.names.remove <- c('30b','30s')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

# remove outliers x4
row.names.remove <- c('24b')
x <- x[!(row.names(x) %in% row.names.remove),]
y <- removeY(y, x, row.names.remove)

#new model
m <- pls(x, y, 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)



#plotXResiduals(m, 1, show.labels = T, norm = F )
#plotXYResiduals(m, 1, show.labels = T, norm = F )
#plotSelectivityRatio(m, 1, show.labels = T)
#plotVIPScores(m, 1, show.labels = T)
#plotXYScores(m, 1, show.labels = T)

#r <- predict(m, xt, yt)
#plot(r)
#plotPredictions(r, show.label = T, xlab = 'yt, reference', ylab = 'yt, predicted')
#abline(a = 0, b = 1)
#show(summary(r))
#round(r$y.pred[, 4,1],2)  
#show(yt)

##Preprocessing. Алгоритм же делает препроцессинг сам?
#Y <- DafniaMagna[,17]
#X <- prep.autoscale(DafniaMagna[, 1:16], scale = T, center = T)   #standardization
#X <- prep.snv(DafniaMagna[,1:16])             #normalization
#DafniaMagna <- cbind(X, Y)

##Выбор переменных
#sr <- data.frame(selratio(m))
#VIP <- data.frame(vipscores(m))
#VIPSR <- cbind(VIP,sr)
#plot(VIPSR, show.labels = T)

#select only important variables
#p <-  m$coeffs$p.values[, ncomp = 2, 1]
#xnew <- mda.subset(x, select = p < 0.05)
#plotVIPScores(m)

#new model
#model <- pls(xnew, y, 5, cv = 1, center = T, scale = T)

#plot(model)
#plotYVariance(model, type = "h")

#xt <- xt[,colnames(xnew)]