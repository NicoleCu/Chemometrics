#� PCA ������ ��� ������, ���������� � ������� �������� � ������������
library(readxl)
library(mdatools)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_fisheri",
                                     range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]


m3 <- pca(VibrioFisheri, ncomp =8, scale = T)

par(mfrow = c(1,2)) 

#������, ����� ����� ��������� ��������� ������ ����������
#����� ����� ���������
plotVariance(m3, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m3, type ='b', show.labels = TRUE)  # acum.var.

# ����� ����� ������ � ��������
plotLoadings(m3, c(1,5), show.labels = T) # � ������� �� ������� ������� ����������� ��������
plotScores(m3, c(1,5), show.labels = T) # ���������� �� �����

# ������, ����� ���������� � ������� ������� �������� �� ����������
m3$loadings

# change the number of component, ��� ������������
# ������� ��������� �����
plotResiduals(m3, show.labels = T,  norm = F, log = T, cgroup = "categories")

# !!������-�� �� ��� ��������� ������� �����
c <- categorize(m3, m3$res$cal)
names(c) <- rownames(VibrioFisheri)




