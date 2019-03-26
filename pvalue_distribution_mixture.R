rm(list=ls())

##-- Load packages
library(alabama)

##-- Read data and functions
path <- '...'                                          # Path for the data 
setwd(path)                                            # Set the path for the data
source('pvalue_distribution_mixture_functions.R')      # Functions for the likelihood
datos <- read.table('pvalues.txt',header=TRUE)[,1]     # Read data

##############################################################
# Models estimation
##############################################################

var = list(datos=datos)     # If some pvalues are almost 0 or 1, some bounds could be added --> pmin(pmax(datos,10^-5),1-10^-5)

##-- Uniform + Beta 
P0 = c(0.4, 1, 1)
ans1 = auglag(par=P0, fn=fn1, heq=heq1,hin=hin1, hin.jac=hin.jac1, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Betas
P0 = c(0.4, 0.3, 1, 1, 1, 1)
ans2 = auglag(par=P0, fn=fn2, heq=heq2,hin=hin2, hin.jac=hin.jac2, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Triangulars
P0 = c(0.4, 0.3, 0.05, 0.95)
ans3 = auglag(par=P0, fn=fn3, heq=heq3,hin=hin3, hin.jac=hin.jac3, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Exponential
P0 = c(0.8, 0.1, 0.1, 0.1)
ans4 = auglag(par=P0, fn=fn4, heq=heq4,hin=hin4, hin.jac=hin.jac4, control.outer=list(trace=FALSE), var=var)

##############################################################
# Check distributions using a QQplot
##############################################################
graphics.off()
windows()
par(las=1,mfrow=c(2,2))
d <- seq(0,1,0.001)
plot(ecdf(datos),cex=0.2,main='Uniform + Beta')
lines(d,F1(d,ans1$par),col=2)
plot(ecdf(datos),cex=0.2,main='Uniform +  2 Betas')
lines(d,F2(d,ans2$par),col=2)
plot(ecdf(datos),cex=0.2,main='Uniform +  2 Triangular')
lines(d,F3(d,ans3$par),col=2)
plot(ecdf(datos),cex=0.2,main='Uniform +  2 Exponentials')
lines(d,F4(d,ans4$par),col=2)

##############################################################
# Check distributions using:
# 1. Histogram with mixture distribution density overlapped
# 2. Statistics, such as AIC and Komogorov- Smirnov
# results also provide the proportion (pu) of pvalues comming 
# from the uniform distribution
##############################################################
graphics.off()
windows()
par(las=1,mfrow=c(2,2))
RES <- rbind(results(datos,ans1,nparam=3,F1,f1,'Uniform + Beta'),
             results(datos,ans2,nparam=6,F2,f2,'Uniform + 2 Betas'),
             results(datos,ans3,nparam=4,F3,f3,'Uniform + 2 Triangular'),
             results(datos,ans4,nparam=4,F4,f4,'Uniform + 2 Exponentials'))
colnames(RES) <- c('nparam','LL','AIC','KS1','KS2','pu','95%LL(pu)','95%UL(pu)','max_f(x)','iterations')
RES
