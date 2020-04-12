n <- 10
x <- rnorm(n)
s <- var(x)
aux <- 0
for(i in 1:n){
  for(j in 1:n){
    aux <- aux + (x[i] - x[j])^2
  }  
}
s2 <- aux/(2*n*(n-1))
s;s2


## New formula
# https://www.jepusto.com/distribution-of-sample-variances/
library(MASS)

covxy <- 0.7   # Cov(X,Y)
Vx <- 3        # V(X)
Vy <- .2       # V(Y)
nsim <- 3000
n <- 100
mu <- c(0,0)
Sigma <- matrix(c(Vx,covxy,covxy,Vy),nrow=2)
S2_x <- S2_y <- c()

set.seed(12345)
for (k in 1:nsim){
  mostra <- mvrnorm(n = 100, mu=mu, Sigma=Sigma)
  S2_x[k] <- var(mostra[,1])
  S2_y[k] <- var(mostra[,2])
}

# Simulated result
(res1 <- cor(S2_x,S2_y))

# Expected result
(res3 <- covxy^2/(Vx*Vy))   # Bueno!



  
  
    