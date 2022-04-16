library(readxl)
library(mdatools)
library(Metrics)

# data import 

ModelWaters <- data.frame(read_excel("~/R/Chemometrics/processed data ModelWaters.xlsx",
                                       sheet = "Sheet1", range ='A1:ED22')) 
row.names(ModelWaters) <- c(ModelWaters[,1])
ModelWaters <-  ModelWaters[,-1]
colnames(ModelWaters) <- read_excel("~/R/Chemometrics/processed data ModelWaters.xlsx",
                                      sheet = "Sheet1", range ='B1:ED1', col_names = F)
# remove outliers 
ModelWaters <-  ModelWaters[!(row.names(ModelWaters) %in% c('11','19', '9', '18')),]

MW_x <- ModelWaters[,1:132]
MW_y <- ModelWaters[,'Dafnia surv rate']



##############################################################################
## PLS на данных с преобразованием


#new model
new_m <- pls(MW_x, MW_y, 7, cv = 1, scale = T, center =T)
plot(new_m)
show(summary(new_m))    
#look for the outliers
plotPredictions(new_m, show.labels = T)
abline(a = 0, b = 1)
plotRegcoeffs(new_m, type = 'h', show.labels = T)

coef <- new_m[["coeffs"]][["values"]]
coef <- coef[,new_m[["ncomp.selected"]],1]
vec <- vector()

for (i in 1:length(coef)){
  if (abs(coef[i]) > 0.031) {
    vec <-  append(vec, i)
  }
}

#new model
new_m_n <- pls(MW_x[,vec], MW_y, 10, cv = 1, scale = T)
plot(new_m_n)
show(summary(new_m_n))
#look for the outliers
plotPredictions(new_m_n, show.labels = T)
abline(a = 0, b = 1)

show(colnames(ModelWaters[,vec]))
plotRegcoeffs(new_m_n, type = 'h', show.labels = T)
