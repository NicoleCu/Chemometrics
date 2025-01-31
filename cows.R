library(readxl)
library(mdatools)
library(Metrics)
library(nortest)
library(nortest)


# data import 

cows <- data.frame(read_excel("~/R/1.xlsx", range ='B1:F255', col_names = T)) 



hist(cows$����, main = "����������� �������� ����", col = "tomato", freq = F)
curve(dnorm(x, mean = mean(cows$����, na.rm = TRUE), 
            sd = sd(cows$����, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)


hist(cows$���.��.���������.����..��, main = "����������� ���������� ��������� ����", col = "tomato", freq = F)
curve(dnorm(x, mean = mean(cows$���.��.���������.����..��, na.rm = TRUE), 
            sd = sd(cows$���.��.���������.����..��, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)



hist(cows$���.��.���������.�����.��, main = "����������� ���������� ��������� �����", col = "tomato", freq = FALSE)
curve(dnorm(x, mean = mean(cows$���.��.���������.�����.��, na.rm = TRUE), 
            sd = sd(cows$���.��.���������.�����.��, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)






show(lillie.test(cows$����))
show(lillie.test(cows$���.��.���������.����..��))
show(lillie.test(cows$���.��.���������.�����.��))


pearson.test(cows$����)
pearson.test(cows$���.��.���������.����..��)
pearson.test(cows$���.��.���������.�����.��)

t.test(cows$����)
t.test(cows$���.��.���������.����..��)
t.test(cows$���.��.���������.�����.��)
