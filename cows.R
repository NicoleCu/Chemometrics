library(readxl)
library(mdatools)
library(Metrics)
library(nortest)
library(nortest)


# data import 

cows <- data.frame(read_excel("~/R/1.xlsx", range ='B1:F255', col_names = T)) 



hist(cows$удой, main = "Гистограмма величины удоя", col = "tomato", freq = F)
curve(dnorm(x, mean = mean(cows$удой, na.rm = TRUE), 
            sd = sd(cows$удой, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)


hist(cows$кол.во.молочного.жира..кг, main = "Гистограмма содержания молочного жира", col = "tomato", freq = F)
curve(dnorm(x, mean = mean(cows$кол.во.молочного.жира..кг, na.rm = TRUE), 
            sd = sd(cows$кол.во.молочного.жира..кг, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)



hist(cows$кол.во.молочного.белка.кг, main = "Гистограмма содержания молочного белка", col = "tomato", freq = FALSE)
curve(dnorm(x, mean = mean(cows$кол.во.молочного.белка.кг, na.rm = TRUE), 
            sd = sd(cows$кол.во.молочного.белка.кг, na.rm = TRUE)), 
      col = "blue", lwd = 2, add = TRUE)






show(lillie.test(cows$удой))
show(lillie.test(cows$кол.во.молочного.жира..кг))
show(lillie.test(cows$кол.во.молочного.белка.кг))


pearson.test(cows$удой)
pearson.test(cows$кол.во.молочного.жира..кг)
pearson.test(cows$кол.во.молочного.белка.кг)

t.test(cows$удой)
t.test(cows$кол.во.молочного.жира..кг)
t.test(cows$кол.во.молочного.белка.кг)
