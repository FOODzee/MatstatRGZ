a <- 1
disp <- 1.1

Y <- c( 1.221,  1.273, 1.158, 2.154,  0.088, -0.583, 1.418,  0.102, 0.380,  2.619,
        1.286,  2.200, 1.187, 1.324, -0.600,  0.727, 2.511,  1.943, 0.301,  2.040)

Z <- c(-0.190,  2.607, 0.810, 0.691,  1.712, -1.209, 0.472, -0.142, 0.497, -0.421,
        0.651, -0.458, 0.707, 2.245, -0.017,  0.649, 1.433,  0.886, 0.286,  0.911,
       -0.050, -0.411, 1.106, 3.506,  2.617,  0.171, 2.441,  2.983, 1.259,  1.723)

X <- c(Y,Z)

# Эти числа задают разбиение выборки для критерия хи-квадрат.
# !!!
# !!! Их следует изменить так, чтобы все значения вектора npi лежали от 8 до 10
# !!! и
# !!! число точек в каждой области было ~ npi.
# !!!
# min и max менять нет смысла, от них ничего не зависит.
zones <- c(min(X), 0.115, 1-0.27,  1+0.26, 1.87, max(X))


## Всё, что ниже, почти наверняка не надо исправлять.
## (Но желательно прочесть, чтобы знать, что происходит).

n <- length(X)
m <- length(Y)
k <- length(Z)

F <- function(x) { pnorm(x, mean=a, sd=sqrt(disp)) }
Fn <- stepfun(X[order(X)], c(0:n)/n, right = TRUE)



### Calc stats.
X_  <- mean(X)
S12 <- 1/n * sum((X - a)^2) 
S02 <- 1/(n-1) * sum((X - X_)^2)

Y_   <- mean(Y)
S02y <- 1/(m-1)*sum((Y - Y_)^2)

Z_   <- mean(Z)
S02z <- 1/(k-1)*sum((Z - Z_)^2)



### Kolmogorov criteria.
diff <- abs(Fn(X) - F(X))
Dk <- max(diff)             # Расстояние Колмогорова.
point <- X[which.max(diff)] # Точка, в которой достигается максимум расстояния.



### Draw empirical distribution function.
plot(Fn, main="Эмпирическая функция распределения", xlab="t", ylab="Fn(t)",
     cex=.3, xaxt="n", yaxt="n", verticals=FALSE)

# Отметим точку, в которой дистигается максимум расстояния Колмогорова.
abline(v=point, h=Fn(point), col="grey")

# Проведём пунктиром график истинного распределения.
points <- seq(min(X)-0.33, max(X)+0.33, length.out = 200)
lines(points, F(points), lty=3)

# Поставим на оси Х метки на точках из выборки,
axis(side=1, at=X)
# а на оси Y -- метки на скачках э.ф.р.
axis(side=2, at=c(0:n)/n)



### hi^2 criteria.
# Построим гистограмму. 
Xhist <- hist(X, breaks=zones, plot=FALSE) 
plot(Xhist, xaxt="n", yaxt="n", freq=TRUE)   # Нарисуем её,
axis(side=1, at=Xhist$breaks)                # пометив на оси X границы областей,
axis(side=2, at=Xhist$counts)                # а на оси Y - высоты столбцов.

# Функция для вычисления n*(вероятность попасть в [a, b])
np_ab <- function(a,b) {  n*(F(b) - F(a)) }

# helper-функция, чтобы интеграл от крайних областей считался до -+inf
np_vector <- function(breaks) {
  lbreaks <- breaks[2:(length(breaks)-2)] # выбираем только левые границы внутренних интервалов
  rbreaks <- breaks[3:(length(breaks)-1)] # и только правые
  c(n*F(breaks[2]), np_ab(lbreaks, rbreaks), n*(1 - F(breaks[length(breaks)-1])))
}

# Считаем теоретическую вероятность попадания в каждую из областей
npi <- np_vector(Xhist$breaks)

# Вычисляем n*(расстояние хи-квадрат)
nDhi <- sum((Xhist$counts - npi)^2/npi)

# Вычисляем реально достигнутый уровень значимости
RDUZhi <- 1 - pchisq(nDhi, length(Xhist$counts))



# Part 3
tau <- (Y_ - Z_)/sqrt((m-1)*S02y + (k-1)*S02z) * sqrt(m+k-2)/sqrt(1/m+1/k)