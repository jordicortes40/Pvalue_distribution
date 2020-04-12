##-- Jacknife: varianza dividida por n
var2 <- function(x){n <- length(x); var(x) * (n-1)/n}
n <- 4
x <- rnorm(n)
v <- var2(x)
vjack <- c()
for(i in 1:n){
  vjack[i] <- var2(x[-i]) #* (n-1)/n
}
mean(vjack)
v

##-- Jacknife: media
n <- 4
x <- rnorm(n)
v <- mean(x)
vjack <- c()
for(i in 1:n){
  vjack[i] <- mean(x[-i])
}
mean(vjack)
v






library(bootstrap)
xdata <- rnorm(10)
var2 <- function(x){n <- length(x); var(x) * (n-1)/n}
results <- jackknife(xdata,var2)
results

jackknife
