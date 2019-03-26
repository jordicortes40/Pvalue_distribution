##############################################################
# List
# 1. Uniform + 1 beta
# 2. Uniform + 2 beta
# 3. Uniform + 2 triangular
# 4. Uniform + 2 exponential
##############################################################

##############################################################
# Uniform + 1 beta
# x[1]:   Proportion coming from the U(0,1)
# x[2:3]: Beta parameters
##############################################################
fn1 <- function(x, var, ...) {
  datos = var$datos
  y = log(x[1] + (1-x[1])*dbeta(datos, x[2], x[3]))
  -sum(y)
}

hin1 <- function(x, var, ...) {
  h <- c()
  h[1:3] = x-1.e-5      # x > 0
  h[4] = 1-x[1]         # p < 1
  h
}

hin.jac1 <- function(x, var, ...) {
  j <- matrix(0, 4, length(x))
  j[1:3, 1:3] = diag(c(1,1,1))
  j[4,1] = -1
  j
}

heq1 <- function(x, var, ...) {
}


##############################################################
# Uniform + 2 beta
# x[1]:   Proportion coming from the U(0,1)
# x[2]:   Proportion coming from the first beta distribution
# x[3:4]: Beta parameters for the first beta distribution
# x[5:6]: Beta parameters for the 2nd beta distribution
##############################################################

fn2 <- function(x, var, ...) {
  datos = var$datos
  y = log(x[1] + x[2]*dbeta(datos, x[3], x[4]) + (1-x[1]-x[2])*dbeta(datos, x[5], x[6]))
  -sum(y)
}

hin2 <- function(x, var, ...) {
  h <- c()
  h[1:6] = x-1.e-5    # x > 0
  h[7:8] = 1-x[1:2]   # p < 1
  h[9] = x[4]-x[3]    # a1 < b1
  h[10] = x[5]-x[6]   # a2 > b2
  h
}

hin.jac2 <- function(x, var, ...) {
  j <- matrix(0, 10, length(x))
  j[1:6, 1:6] = diag(rep(1,6))
  j[7:8, 1:2] = diag(c(-1,-1))
  j[9, c(3,4)] = c(-1, 1)
  j[10, c(5,6)] = c(1, -1)
  j
}

heq2 <- function(x, var, ...) {
} 

##############################################################
# Uniform + 2 triangular
# x[1]: Proportion coming from the U(0,1)
# x[2]: Proportion coming from the first (right) triangular distribution
# x[3]: Parameter of the first (right) triangle distribution
# x[4]: Parameter of the second (right) triangle distribution
##############################################################

dtriangle1 <- function (x,base) pmax(2/base - 2/(base^2) * x, 0.001 , na.rm = FALSE)
dtriangle2 <- function (x,base) pmax((2/(1-base)^2) * (x - base), 0.001 , na.rm = FALSE)

fn3 <- function(x, var, ...) {
  datos = var$datos
  y = log(x[1] + x[2]*dtriangle1(datos,x[3]) + (1-x[1]-x[2])*dtriangle2(datos,x[4]))
  -sum(y)
}

hin3 <- function(x, var, ...) {
  h <- c()
  h[1:4] = x-1.e-5          # x > 0
  h[5:8] = 1-x[1:4]-1.e-5   # x < 1
  h
}

hin.jac3 <- function(x, var, ...) {
  j <- matrix(0, 8, length(x))
  j[1:4, 1:4] = diag(rep(1,4))
  j[5:8, 1:4] = diag(rep(-1,4))  
  j
}

heq3 <- function(x, var, ...) {
} 

##############################################################
# Uniform + 2 exponential
# x[1]: Proportion coming from the U(0,1)
# x[2]: Proportion coming from the first exponential distribution
# x[3]: Parameter of the first exponential distribution
# x[4]: Parameter of the second translated distribution
##############################################################
##-- Exponential 1
#|\
#| \
#|  \
#|   \
#|    \
#|      \
#|         \
#|              \
#|                    \
#|                          \
#-----------------------------
#0                           1

##-- Exponential 2
#                           /|
#                         /  |
#                       /   |
#                      /     |
#                    /       |
#                 /          |
#              /             |
#         /                  |              \
#    /                       |                   
# /                          |
#-----------------------------
#0                           1


dexp1 <- function (x,lambda) dexp(x,lambda)/(1-exp(-lambda))            # Standardized to interval 0-1
dexp2 <- function (x,lambda) lambda*exp(lambda*(x-1))/(1-exp(-lambda))  # Translated

fn4 <- function(x, var, ...) {
  datos = var$datos
  y = log(x[1] + x[2]*dexp1(datos,x[3]) + (1-x[1]-x[2])*dexp2(datos,x[4]))
  -sum(y)
}

hin4 <- function(x, var, ...) {
  h <- c()
  h[1:4] = x-1.e-5    # x > 0
  h[5:6] = 1-x[1:2]   # p < 1
  h
}

hin.jac4 <- function(x, var, ...) {
  j <- matrix(0, 6, length(x))
  j[1:4, 1:4] = diag(rep(1,4))
  j[5:6, 3:4] = diag(rep(-1,2))
  j
}

heq4 <- function(x, var, ...) {
} 




##############################################################
# Function results (AIC, Kolmogorov,...)
##############################################################
results <- function(datos,ans,nparam,F,f,tit,hist=TRUE){
  
  param <- ans$par
  
  ##-- AIC
  LL <- (-ans$value)       # Log-likelihood
  AIC <- 2*nparam - 2*LL   # AIC
  
  ##-- KS --> Method 1
  x <- seq(0.0001,0.9999,0.001)
  y <- f(x,param)
  Fteorica <- F(x,param)
  Fd <- ecdf(datos)
  Fdatos <- Fd(x)
  KS1 <- max(abs(Fdatos-Fteorica))
  
  ##-- KS --> Method 2
  nd <- 0
  ysim <- c()
  ntry <- 0
  ymax <- max(y)
  set.seed(12345)
  while (nd<5000){
    u1 <- runif(1,0,ymax)
    x1 <- runif(1)
    IN <- f(x1,param)>u1
    if(IN) ysim <- c(ysim,x1)
    ntry <- ntry + 1
    nd <- length(ysim)
  }
  kst <- ks.test(datos,ysim)
  KS2 <- as.numeric(kst$sta)
  
  ## IC95%
  SE <- sqrt(diag(solve(ans$hessian)))
  LI <- ans$par - 1.96*SE
  LS <- ans$par + 1.96*SE
  
  ##-- Uniform
  pu <- ans$par[1]
  
  ##-- Plot
  if(hist){
    hist(datos,col=4,border='grey80',freq=FALSE,main=tit)
  }else{
    plot(NA,xlim=c(0,1),ylim=c(0,2.14),main='Model',bty='n',xlab = 'p-value',ylab='Density')
  }
  x0 <- c(0,x,1,0)
  y0 <- c(0,y,0,0)
  polygon(x0,y0,col=rgb(0,0,1,0.5),border=NA,ylim=c(0,5))
  # plot(ecdf(datos),verticals=F,pch=0)
  # plot(ecdf(y),verticals=F,pch=0,add=TRUE,col=2)
  
  return(c(nparam,LL,AIC,KS1,KS2,pu,LI[1],LS[1],max(y),ntry))
}

##############################################################
# Density functions
##############################################################
f1 <- function(datos,x) x[1] + (1-x[1])*dbeta(datos, x[2], x[3])
f2 <- function(datos,x) x[1] + x[2]*dbeta(datos, x[3], x[4]) + (1-x[1]-x[2])*dbeta(datos, x[5], x[6])
f3 <- function(datos,x) x[1] + x[2]*dtriangle1(datos,x[3]) + (1-x[1]-x[2])*dtriangle2(datos,x[4])
f4 <- function(datos,x) x[1] + x[2]*dexp1(datos,x[3]) + (1-x[1]-x[2])*dexp2(datos,x[4])

##############################################################
# Distribution functions
##############################################################
F1 <- function(d,x) x[1]*d + (1-x[1])*pbeta(d, x[2], x[3])
F2 <- function(d,x) x[1]*d + x[2]*pbeta(d, x[3], x[4]) + (1-x[1]-x[2])*pbeta(d, x[5], x[6])
F3 <- function(d,x){
  n <- length(d)
  r <- c()
  for(i in 1:n){
    if(d[i] < x[3])               r[i] <- x[1]*d[i] + x[2]*d[i]/x[3]*(2-d[i]/x[3]) + (1-x[1]-x[2])*0
    if(x[3] <= d[i] & d[i] <x[4]) r[i] <- x[1]*d[i] + x[2]*1 + (1-x[1]-x[2])*0
    if(d[i] >= x[4])              r[i] <- x[1]*d[i] + x[2]*1 + (1-x[1]-x[2])/(1-exp(-x[4]))*(exp(x[4]*(d[i]-1)) - exp(-x[4]))
  }
  
  return(r)
}
F4 <- function(d,x) x[1]*d + x[2]*pexp(d,x[3])/(1-exp(-x[3])) + (1-x[1]-x[2])/(1-exp(-x[4]))*(exp(x[4]*(d-1)) - exp(-x[4]))
