library(mdatools)
data(people)

Height <- people[,'Height']
Weight <-  people[,'Weight']

#��� ����������. �� ��������� ��� ��������� �����������
cov(Height, Weight)
cov(Weight, Height)

#��������� ����� ����� ����� �����
cov(Height, Height)

#covariance matrix ����� ����� ����� �����
C <- cov(people[,1:3])
C

#��� �����: crossprod ~ X'X (������ ���������������)
#tcrossprod ~ XX' (������ ���������������)
C2 <- crossprod(Xm)/(nrow(Xm)-1) 
C2

#mean centering
Xm <- scale(people[,1:3], center = TRUE, scale = FALSE)
head(Xm)

# �������� �������, 1 - �� ��������, 2 - �� �������
apply(people[,1:3], 2, mean)

# ����� ������ �������������� ������
plot(Xm[,1], Xm[,2],
     main = "Center data",
     xlab = 'Height',
     ylab = 'Weight')
#�������� ����� �� �������, ������ � ������ �������� intercept and slope �������������
# h =... ����� ����������� ��� �������, v =... ����������� �������
abline(h=0, col = 'Gray', lty = 2)
abline(v=0, col = 'Gray', lty = 2)
#����� ����� ������ (�������� ���������), � lm ���������� slope and intercept
abline(lm(Xm[,2]~Xm[,1]))


