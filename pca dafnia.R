#� PCA ������ ��� ������, ���������� � ������� �������� � ������������
library(readxl)
library(mdatools)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

m <- pca(DafniaMagna, ncomp =6, scale = T)

par(mfrow = c(1,2)) 

#������, ����� ����� ��������� ��������� ������ ����������
#����� ����� ���������
plotVariance(m, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

# ����� ����� ������ � ��������
plotLoadings(m, c(1,2), show.labels = T) # � ������� �� ������� ������� ����������� ��������
plotScores(m, c(1,2), show.labels = T) # ���������� �� �����

# ������, ����� ���������� � ������� ������� �������� �� ����������
m$loadings

# change the number of component, ��� ������������
# ������� ��������� �����
plotResiduals(m, show.labels = T,  norm = F, log = T, cgroup = "categories")

# !!������-�� �� ��� ��������� ������� �����
c <- categorize(m, m$res$cal)
names(c) <- rownames(DafniaMagna)




## PCA ������ ��� ������, ���������� � ������� ������������ ��������
DafniaMagna_bioind <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='R1:R51')) 
row.names(DafniaMagna_bioind) <- row.names(DafniaMagna)

m2 <- pca(DafniaMagna, ncomp = 6, center =T, scale = T)
par(mfrow = c(1,2)) 

#������, ����� ����� ��������� ��������� ������ ����������
#����� ����� ���������
plotVariance(m2, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m2, type ='b', show.labels = TRUE)  # acum.var.

# ����� ����� ������ � ��������
plotLoadings(m2, c(1,6), show.labels = T) # � ������� �� ������� ������� ����������� ��������
plotScores(m2, c(1,6), show.labels = T) # ���������� �� �����

# ������, ����� ���������� � ������� ������� �������� �� ����������
m$loadings

# change the number of component, ��� ������������
# ������� ��������� �����
plotResiduals(m2, show.labels = T, norm = F, log = T, cgroup = "categories")
