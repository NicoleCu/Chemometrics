library(readxl)
library(mdatools)
library(Metrics)


x <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                           range ='B2:L23'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                           range ='N2:Q23'))
y <- y[,-2]
rownames(y) <- seq(0,20, 1)


x <-  x[!(row.names(x) %in% c('0', '1', '5', '16')),]
y <-  y[!(row.names(y) %in% c('0','1', '5', '16')),]

# PLS for original data
m_cd <- pls(x[,], y[,2], 10, cv = 1, scale = F)
plot(m_cd)
show(summary(m_cd, ncomp = 3))
plotYCumVariance(m_cd, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_cd, type  = 'h', show.labels = TRUE)

# look for the outliers
plotPredictions(m_cd, show.labels = T, ncomp = 3)
abline(a = 0, b = 1)
plotRegcoeffs(m_cd, show.labels = T, ncomp =3)


###############################################################################
### fl sensor

fl_1 <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "FL-Aver",
                              range ='A160:Z181')) 
row.names(fl_1) <- c(fl_1[,1])
fl_1 <-  fl_1[,-c(1,3)] #remove names and Pb

fl_1 <-  fl_1[!(row.names(fl_1) %in% c('0','11', '5', '12', '16')),]


m_1 <- pls(fl_1[,2:24],fl_1[,1], 10, cv = 1, scale = F, center = T)
#plot(m_1)
show(summary(m_1))
plotRMSE(m_1)
plotYCumVariance(m_1)
#look for the outliers
plotPredictions(m_1, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_1, show.labels = T)


################################################################################
########### data fusion

fused <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "cd data fusion",
                               range ='A3:AJ24'))
row.names(fused) <- c(fused[,1])

fused <-  fused[,-c(1)] #remove names
fused <-  fused[!(row.names(fused) %in% c('0', '11', '5', '12', '16')),]


#new model
m <- pls(fused[,2:35], fused[,1], 5, cv = 1, scale = T, center = T)
plotRMSE(m)
plotYCumVariance(m)
plotXCumVariance(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h', lab.cex = 0.8)