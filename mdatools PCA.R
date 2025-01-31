library(mdatools)
data(people)

X <- people

m <- pca(X, ncomp = 4, scale = TRUE)

m$loadings

par(mfrow = c(1,2)) 
plotLoadings(m, c(1,4), show.labels = T) # � ������� �� ������� ������� ����������� ��������
# 
plotScores(m, c(1,4), show.labels = T) # ���������� �� �����

# barplo�, ����� ���������� ��� ����� ��������
plotLoadings(m, show.labels = T, type = 'h')

summary(m)
plotVariance(m, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

m$res
cumresidual <- 100 - m$res$cal$cumexpvar
cumresidual
barplot(as.numeric(cumresidual))

## ����������
# distance plot, ������������� ���������� ������ ���������, ������� ������ ����������
plotResiduals(m)
# change the number of component, ��� ������������
plotResiduals(m, c(1,2), show.labels =T, norm = F, log = T)
# ������� ��������� �����
plotResiduals(m, c(1,2), show.labels = T, cgroup = "categories")

# compute the categories
c <- categorize(m, m$res$cal)
names(c) <- rownames(people)

# �������� 2 ����� �����
P1 <- c(159, 50, +1, 36, 20, 15000, 10, 10, +1, 95, +1, 110)
P2 <- c(159, 50, +1, 36, 20, 30000, 10, 10, +1, 95, +1, 120 )
new <- rbind(P1, P2)

r <- predict(m, new) # use ready model for a new data
plotScores(m, res=list("cal" = m$res$cal, "new" = r), show.labels = TRUE)

plotResiduals(m, res=list("cal" = m$res$cal, "new" = r), show.labels = TRUE)

# �������, ��� ��� � residuals, ������ ����� ������ ������ � �������

E <- r$residuals
colnames(E) <- colnames(people)
mdaplotg(E, type ='h', show.labels = T)

# ������� �� ������ 2�� ����� ��������
cnew <- categorize(m, r)
cnew

# ��������� ������ ��� ������
m <- pca(people, ncomp = 6, scale = TRUE, lim.type = 'ddmoments', alpha = 0.05, gamma = 0.01)
plotResiduals(m,4)
