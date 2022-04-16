library(readxl)
library(mdatools)
library(Metrics)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx", 
                                       sheet = "analysis_fisheri", range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]

VibrioFisheri <-  VibrioFisheri[!(row.names(VibrioFisheri) %in% 
                                    c('w24', 'w6')),]

colnames(VibrioFisheri) <- c(seq(1,23, 1), 'EC50')

#new model
m <- pls(VibrioFisheri[,1:23], VibrioFisheri[,24], 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T, cex.lab = 3)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T)


# remove predictors = noise

x <- VibrioFisheri[,1:23]
y <- VibrioFisheri[,24]

rem <- c(11,20)

m_new <- pls(x[,-rem], y, 10, cv = 1)
plot(m_new)
show(summary(m_new))
#look for the outliers
plotPredictions(m_new, show.labels = T)
abline(a = 0, b = 1)



plotRegcoeffs(m, show.labels = T)

