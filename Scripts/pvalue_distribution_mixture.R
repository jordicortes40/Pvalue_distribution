rm(list=ls())

##-- Load packages
library(alabama)

##-- Read data and functions
path <- '...'                                          # Path for the data 
setwd(path)                                            # Set the path for the data
source('pvalue_distribution_mixture_functions.R')      # Functions for the likelihood
datos <- read.table('pvalues.txt',header=TRUE)[,1]     # Read data

##############################################################
# Descriptive
##############################################################
library(ggplot2)
library(gridExtra)
d <- data.frame(x=datos)
gg1 <- ggplot(d,aes(x=x,y=stat(density))) + geom_histogram(col='white',bins=10,breaks=seq(0,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous() + xlab('p-values') + ylab('n') + labs(title='Histogram') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
gg2 <- ggplot(d,aes(sample=x)) + stat_qq(distribution = qunif,size=2,alpha=0.5,color=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  stat_qq_line(distribution = qunif) +
  scale_x_continuous() + xlab('Theoretical Uniform') + ylab('Sample') + labs(title='QQplot') +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
grid.arrange(gg1,gg2,nrow=1)


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
## -- Base
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

x <- seq(0,1,0.001)

##-- Uniform + Beta
d1 <- data.frame(x=x,y=pmin(2,f1(x,ans1$par)))
gg4 <- ggplot(d,aes(x=x)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + Beta') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d1,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Betas
d2 <- data.frame(x=x,y=pmin(2,f2(x,ans2$par)))
gg3 <- ggplot(d,aes(x=x)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Betas') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d2,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Triangular
d3 <- data.frame(x=x,y=pmin(2,f3(x,ans3$par)))
gg1 <- ggplot(d,aes(x=x)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Triangulars') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d3,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Exponentials
d4 <- data.frame(x=x,y=pmin(2,f4(x,ans4$par)))
gg2 <- ggplot(d,aes(x=x)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Exponentials') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d4,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

grid.arrange(gg1,gg2,gg3,gg4,nrow=2)
