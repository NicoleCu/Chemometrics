library(mdatools)
library(pheatmap)

#вызываю встроенный датасет

data(people)

#выбираю переменные 

Age <- people[,"Age"]
Income <- people[, 'Income']

# строю зависимость дохода от возраста
plot(Age, Income, col ='Orange')

# создаю модель линейной регрессии (на графике)

abline(lm(Income ~ Age), col = 'Gray', lty =2)

# вычисл€ю коррел€цию (коэффициент коррел€ции ѕирсона)
cor(Age, Income)

# доверительные интервалы дл€ корел€ции
cor.test(Age, Income)

# матрица коэффиентов коррел€ции дл€ более 2ух переменных
R <- cor(people)
R

#строю heatmap
#image(R)  # по-простому

#или через пакет
breaks <- seq(-1, 1, length.out = 100)
pheatmap(R, breaks = breaks, 
         treeheight_row = 0, 
         treeheight_col = 0)
