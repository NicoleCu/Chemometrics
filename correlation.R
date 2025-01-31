library(mdatools)
library(pheatmap)

#������� ���������� �������

data(people)

#������� ���������� 

Age <- people[,"Age"]
Income <- people[, 'Income']

# ����� ����������� ������ �� ��������
plot(Age, Income, col ='Orange')

# ������ ������ �������� ��������� (�� �������)

abline(lm(Income ~ Age), col = 'Gray', lty =2)

# �������� ���������� (����������� ���������� �������)
cor(Age, Income)

# ������������� ��������� ��� ���������
cor.test(Age, Income)

# ������� ����������� ���������� ��� ����� 2�� ����������
R <- cor(people)
R

#����� heatmap
#image(R)  # ��-��������

#��� ����� �����
breaks <- seq(-1, 1, length.out = 100)
pheatmap(R, breaks = breaks, 
         treeheight_row = 0, 
         treeheight_col = 0)
