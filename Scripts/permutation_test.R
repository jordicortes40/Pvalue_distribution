rm(list=ls())
set.seed(1234)
NSIM <- 100
nsim <- 500
n <- 200
minE <- 0
maxE <- 10
Th <- 0:70
LS <- c()
fn <- function(x,d) formatC(x,digits=d,format='f')
for (sim in 1:NSIM){
  ##-- Generate potential outcomes
  Y1 <- round(rnorm(n,100,10))
  effect <- sample(minE:maxE,n,rep=TRUE)
  Y2 <- Y1 + effect
  # plot(Y1,Y2)
  # abline(0,1)
  
  ##-- Allocation
  treated <- sample(c(FALSE,TRUE),n,rep=TRUE)
  y1 <- Y1[treated]
  y2 <- Y2[!treated]
  d <- data.frame(y=c(y1,y2),x=c(rep(0,length(y1)),rep(1,length(y2))))
  
  ##-- Observed statistic
  #t_obs <- t.test(y~x,d)$st
  t_obs <- with(d,diff(tapply(y,x,mean)))
  
  
  ##-- Permutations
  j <- 1
  p <- 0
  th <- 0
  while(th<max(Th) & p<0.05){
    cat('=')
    th <- Th[j]
    y1_imp <- c(y2 - th, y1)
    y2_imp <- c(y2,      y1 + th)
    y_imp <- ifelse(sample(0:1,length(y1_imp),rep=TRUE),y1_imp,y2_imp)
    sta <- c()
    for (k in 1:nsim){
      daux <- data.frame(y=y_imp,x=sample(0:1,length(y_imp),rep=TRUE))  
      sta[k] <- with(daux,diff(tapply(y,x,mean))) #t.test(y~x,daux)$st
    }
    p <- sum(sta>t_obs)/nsim
    j <- j+1  
  }
  cat('\n')
  #summary(mod <- lm(p~Th))
  #plot(p~Th)
  #abline(h=0.05,lty=2)
  #abline(h=0.025,lty=2)
  #abline(v=maxE,lty=2)
  #abline(mod,lty=3)
  LS[sim] <- th
  bn <- binom.test(sum(LS<maxE),length(LS))
  cat('sim:',sim,'LS:',LS[sim],'t_obs:',formatC(t_obs,digits=2,format='f'),
      'p:',p,'confidence:',
      paste0(round(bn$estimate,3),' [',bn$con[1],',',bn$con[2],']','\n')
}



