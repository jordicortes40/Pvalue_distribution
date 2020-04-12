rm(list=ls())
set.seed(12345)
NSIM <- 1000   # Number of simulations
nsim <- 500    # number of permutation tests
n <- 200       # sample size
minE <- 0      # minimum effect
maxE <- 10     # maximum effect  
Th <- 0:100    # 70 # (maxE+1)
LS <- c()
fn <- function(x,d=2) formatC(x,digits=d,format='f')
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
    sta <- c()
    w <- sample(0:1,length(y1_imp),rep=TRUE)
    y_imp <- ifelse(w,y1_imp,y2_imp)
    
    for (k in 1:nsim){
      #w <- sample(0:1,length(y1_imp),rep=TRUE)
      #y_imp <- ifelse(w,y1_imp,y2_imp)
      daux <- data.frame(y=y_imp,x=sample(0:1,length(y1_imp),rep=TRUE))  
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
  cat('sim:',sim,'LS:',LS[sim],'t_obs:',fn(t_obs,d=2),
      'p:',fn(p,3),'confidence:',
      paste0(fn(bn$estimate,3),' [',fn(bn$con[1],3),',',fn(bn$con[2],3),']'),'\n')
}
write.table(data.frame(LS=LS),
           file='C:/Users/jcortes/Google Drive/Tesis/Memoria/Scripts/LS_Confidence_intervals.txt',
           append=FALSE,col.names = FALSE,row.names = FALSE, quote = FALSE, sep='')
hist(LS)
summary(LS)

# Alternativa 1 --> Proporciona un alpha mejor: 15/300 pero todos los LS son o max(Th) o 0
# for (k in 1:nsim){
#   w <- sample(0:1,length(y1_imp),rep=TRUE)
#   y_imp <- ifelse(w,y1_imp,y2_imp)
#   daux <- data.frame(y=y_imp,x=w)  
#   sta[k] <- with(daux,diff(tapply(y,x,mean))) #t.test(y~x,daux)$st
# }

# Alternativa 2 --> Proporciona un alpha peor: 23/300 pero los LS varian
#  w <- sample(0:1,length(y1_imp),rep=TRUE)
#  y_imp <- ifelse(w,y1_imp,y2_imp)
# for (k in 1:nsim){
#   daux <- data.frame(y=y_imp,x=sample(0:1,length(y1_imp),rep=TRUE))  
#   sta[k] <- with(daux,diff(tapply(y,x,mean))) #t.test(y~x,daux)$st
# }
