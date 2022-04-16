library(readxl)
library(mdatools)
library(Metrics)

# data import 

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/processed data VibrioFisheri.xlsx",
                                     sheet = "Sheet1", range ='A1:UH27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]
colnames(VibrioFisheri) <- read_excel("~/R/Chemometrics/processed data VibrioFisheri.xlsx",
                                    sheet = "Sheet1", range ='B1:UH1', col_names = F)

# remove outliers 
VibrioFisheri <-  VibrioFisheri[!(row.names(VibrioFisheri) %in% c('w24', 'w6')),]

VF_x <- VibrioFisheri[,1:552]
VF_y <- VibrioFisheri[,'EC50']



##############################################################################
## PLS на данных с преобразованием


#new model
new_m <- pls(VF_x, VF_y, 10, cv = 1, scale = T, center =T)
plot(new_m)
show(summary(new_m))    
#look for the outliers
plotPredictions(new_m, show.labels = T)
abline(a = 0, b = 1)

#plotYVariance(new_m)	

plotRegcoeffs(new_m, type = 'h', labels = T)

coef <- new_m[["coeffs"]][["values"]]
coef <- coef[,new_m[["ncomp.selected"]],1]
vec <- vector()

for (i in 1:length(coef)){
  if (abs(coef[i]) > 0.0050) {
    vec <-  append(vec, i)
  }
}

#new model
new_m_n <- pls(VF_x[,vec], VF_y, 10, cv = 1, scale = T)
plot(new_m_n)
show(summary(new_m_n))
#look for the outliers
plotPredictions(new_m_n, show.labels = T)
abline(a = 0, b = 1)

show(colnames(VibrioFisheri[,vec]))
plotRegcoeffs(new_m_n, type = 'h')
