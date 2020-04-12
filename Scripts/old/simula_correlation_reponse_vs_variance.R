rm(list=ls())
###############################
# PARAMETERS
###############################
n <- 1000        # Tamaño muestral estudio
NSIM <- 5000     # Numero de simulaciones para distintas correlación
nsim <- 100      # Numero de simulaciones para misma correlación
ruido.min <- 0   # Mínima desviación del ruido
ruido.max <- 20  # Máxima desviación del ruido
cor.outcome <- cor.var <- c()
cov.outcome <- cov.var <- c()
sta <- c()
EFFECT <- c(1,10,100)

###############################
# SIMULATION
###############################
set.seed(12345)
graphic.off()
windows(10,10)
par(mfrow=c(1,2))

for(k in 3:3){
  for(i in 1:NSIM){
    cat('iteration:',i,'\n')
    vb <- c()                  # Varianza estimada basal
    vf <- c()                  # Varianza estimada final
    cor.out <- cov.out <- c()  # Correlación y covarianza basal-final respuesta (se podia estimar analíticamente)
    coeff <- runif(1,-EFFECT[k],EFFECT[k])
    #ruido <- runif(1,ruido.min,ruido.max)
    ruido <- max(0,rnorm(1,2,1))
    for(j in 1:nsim){
      z <- rnorm(n)             # Valores basales
      y <- z + rnorm(n,0,ruido) # Valores finales
      
      vb[j] <- var(z)
      vf[j] <- var(y)
      cor.out[j] <- cor(z,y)
      cov.out[j] <- cov(z,y)
    }
    cor.outcome[i] <- mean(cor.out)
    cov.outcome[i] <- mean(cov.out)
    sta[i] <- log(1+2*cor.outcome[i]^2/(n-2))
    cor.var[i] <- cor(vb,vf)
    cov.var[i] <- cov(vb,vf)
  }
  
  ##-- Correlation plot
  plot(cor.outcome,cor.var,pch=19,col=rgb(0,0,1,0.4))
  lines(lowess(cor.outcome,cor.var,f=1/8),col=rgb(0,1,0,0.6),lwd=2)
  curve(x^2,col=rgb(1,0,0,0.6),lwd=2,add=TRUE)
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=0,v=0,lty=2)
  
  ##-- Covariance plot
  # plot(cov.outcome,cov.var,pch=19,col=rgb(0,0,1,0.4))
  plot(sta,cov.var,pch=19,col=rgb(0,0,1,0.4))
  lines(lowess(sta,cov.var),col=rgb(0,1,0,0.6),lwd=2)
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=0,v=0,lty=2)
  
}

###############################
# PLOT
###############################
library(ggplot2)
library(data.table)

fun.x2 <- function(x) x^2
d <- data.table(x=cor.outcome,y=cor.var)
ggplot(d,aes(x=x,y=y)) +
  geom_point(colour='blue',alpha=0.4) +
  stat_function(fun = fun.x2,colour='green',size=1) +
  geom_smooth(colour='red',linetype=2) +
  xlab(expression(bold(Corr(Y[BT],Y[OT])))) + ylab(expression(bold(Corr(S[BT]^2,S[OT]^2))))
  
