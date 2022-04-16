library(readxl)
library(mdatools)
library(Metrics)
library(writexl)



################################################################################
################################################################################

# data import 

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",
                                     sheet = "analysis_dafnia", range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

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
predictors <- data.frame(cbind(predictors, squared, multiplied, divided))


# данные для выгрузки

predictors_names <- c(colnames(DafniaMagna[,1:16]), colnames(squared), colnames(multiplied), colnames(divided), 'Dafnia')
dataset <- cbind(predictors, DafniaMagna[,17])
colnames(dataset) <- predictors_names

write_xlsx(data.frame(rownames(dataset)), "~/R/Chemometrics/processed data Dafnia_ind.xlsx")
write_xlsx(dataset, "~/R/Chemometrics/processed data DafniaMagna.xlsx")




################################################################################
################################################################################

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx", 
                                       sheet = "analysis_fisheri", range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]
colnames(VibrioFisheri) <- c(seq(1,23, 1), 'EC50')


# data modification
output <- as.data.frame(VibrioFisheri[,'EC50'])
rownames(output) <- rownames(VibrioFisheri)
predictors <- VibrioFisheri[,1:23]

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
predictors <- data.frame(cbind(predictors, squared, multiplied, divided))


# данные для выгрузки

predictors_names <- c(colnames(VibrioFisheri[,1:23]), colnames(squared), colnames(multiplied), colnames(divided), 'EC50')
dataset <- cbind(predictors, VibrioFisheri[,24])
colnames(dataset) <- predictors_names

write_xlsx(data.frame(rownames(dataset)), "~/R/Chemometrics/processed data VF_ind.xlsx")
write_xlsx(dataset, "~/R/Chemometrics/processed data VibrioFisheri.xlsx")


###############################################################################
###############################################################################


ModelWaters <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                                     range ='B2:N23')) 
rownames(ModelWaters) <- seq(0,20, 1)
ModelWaters <- ModelWaters[,-12]


# data modification
output <- as.data.frame(ModelWaters[,12])
rownames(output) <- rownames(ModelWaters)
predictors <- ModelWaters[,1:11]

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
predictors <- data.frame(cbind(predictors, squared, multiplied, divided))


# данные для выгрузки

predictors_names <- c(colnames(ModelWaters[,1:11]), colnames(squared), colnames(multiplied), colnames(divided), 'Dafnia surv rate')
dataset <- cbind(predictors, ModelWaters[,12])
colnames(dataset) <- predictors_names

write_xlsx(data.frame(rownames(dataset)), "~/R/Chemometrics/processed data MV_ind.xlsx")
write_xlsx(dataset, "~/R/Chemometrics/processed data ModelWaters.xlsx")



