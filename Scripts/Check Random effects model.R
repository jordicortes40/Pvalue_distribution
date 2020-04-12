rm(list=ls())

####################################################################
# Load libraries
####################################################################
library(metafor)
library(ggplot2)
library(gridExtra)
library(ggpubr)

####################################################################
# Read data
####################################################################
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
data <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")

nsim <- 100                                  # Number of simulations
N1 <- data$final_cases_T1                    # Sample size in one group
N2 <- data$final_cases_T2                    # Sample size in one group
V1 <- data$final_sd_T1^2
V2 <- data$final_sd_T2^2                     # Actual variances

#-------------------------------------------------------------------
#
#
#                          Graphic
#
#
#-------------------------------------------------------------------
####################################################################
# Parameters of the simulation
####################################################################
EFFECT <- seq(0.1,2,0.1)
TAU <- seq(0,1,0.1)
RES <- 2/N1+2/N2
V1 <- data$final_sd_T2^2
M <- matrix(ncol=4,nrow=length(EFFECT)*length(TAU))
set.seed(12345)
j=1
for(e in EFFECT){
  for(t in TAU){
    V2 <- e*V1
    mu <- tau <- c()
    for(i in 1:nsim){
      y <- log(e) + rnorm(208,0,t) + rnorm(208,0,sqrt(RES))
      rma.model <- rma(yi=y,vi=RES)  
      mu <- c(mu,coef(rma.model)[1])
      tau <- c(tau,sqrt(rma.model$tau2))
    }
    M[j,] <- c(e,t,mean(mu),mean(tau))
    j <- j+1
  }
  cat('Effect:',e,'\n')
}

##-- Plot
graphics.off()
windows(8,5)
par(las=1,font.axis=4,font.lab=2,mfrow=c(1,2))
plot(M[,1],exp(M[,3]),col='darkblue',lwd=1,xlab='Real variance ratio',ylab='Estimated variance ratio',log='xy',main=expression(bold(mu)))
rect(0.001,0.001,10000,10000,col='grey90')
abline(v=c(0.1,0.2,0.5,1,2),col='white')
abline(h=c(0.1,0.2,0.5,1,2),col='white')
abline(0,1,lty=2,lwd=2)
points(M[,1],exp(M[,3]),col='darkblue',lwd=1)
plot(M[,2],M[,4],col='darkblue',lwd=1,xlab='Real Heterogeneity',ylab='Estimated Heterogeneity',main=expression(bold(tau)))
rect(-2,-2,10000,10000,col='grey90')
abline(0,1,lty=2,lwd=2)
abline(v=seq(0,2,0.2),col='white')
abline(h=seq(0,2,0.2),col='white')
points(M[,2],M[,4],col='darkblue',lwd=1)

#-------------------------------------------------------------------
#
#
# Simulation under equal variances                          
#
# 
#-------------------------------------------------------------------
####################################################################
# Read data
####################################################################
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
data <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")

####################################################################
# Parameters of the simulation
####################################################################
nsim <- 10000                    # Number of simulations
N1 <- data$final_cases_T1        # Sample size in one group
N2 <- data$final_cases_T2        # Sample size in other group
V1 <- data$final_sd_T1^2         # Variances in treated group
V2 <- V1                         # Variances in control group=treated group


####################################################################
# Simulation
####################################################################
set.seed(12345)

M <- matrix(ncol=3,nrow=nsim)
for (i in 1:nsim){
  est <- c()
  cat('i:',i,'/',nsim,'\n')                     # Print iteration
  for(j in 1:208){
    # n1 <- N1[j]                                 # Sample size in experimental arm
    # n2 <- N2[j]                                 # Sample size in control arm
    # v1 <- V1[j]                                 # Variance in experimental arm
    # v2 <- V1[j]                                 # Variance in control arm = Variance in experimental arm
    y01 <- rnorm(N1[j],0,sqrt(V1[j]))           # Generate normal data for one arm
    y02 <- rnorm(N2[j],0,sqrt(V1[j]))           # Generate normal data for the other arm
    est[j] <- log(var(y01)/var(y02))            # Response
  }
  try(rma.model <- rma(yi=est,vi=2/(N1-2)+2/(N2-2)))
  M[i,] <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2)
}

colnames(M) <- c('mu','tau','I2')
labs <- c(expression(bold(mu^'')),expression(bold(tau^'')),expression(I^2))
M
##-- Normal plot
graphics.off()
windows()
par(mfrow=c(1,3),las=1)
REAL.BASELINE <- c(0.02,0.31,100*0.583)
REAL.FINAL <- c(-0.12,0.59,100*0.833)
ymin <- c(-0.15,0,0)
ymax <- c(0.15,0.65,90)
for (i in 1:3){
  boxplot(M[,i],main=labs[i],ylim=c(ymin[i],ymax[i]),cex.main=1.8)
  rect(-100,-100,100,100,border=NA, col='grey90')
  abline(h=pretty(c(ymin[i],ymax[i])),col='white')
  points(1,REAL.BASELINE[i],col=4,lwd=2,pch=3,cex=2)
  points(1,REAL.FINAL[i],col=2,lwd=2,pch=4,cex=2)
  boxplot(M[,i],add=TRUE,col='white')
  legend('topleft',c('Reference model','Outcome model'),bg = 'white',
         pch=3:4,col=c(4,2),cex=1,pt.cex=1.1,pt.lwd = 2)
}
summary(M)
with(data,rma(2*log(final_sd_T1/final_sd_T2),vi=2/(N1-2) + 2/(N2-2)))                      # Final model
with(data,rma(2*log(base_sd_T1/base_sd_T2),vi=2/(base_cases_T1-2) + 2/(base_cases_T2-2)))  # Baseline model

##-- Plot with ggplot2
dd <- as.data.frame(M)

common.theme <- theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title = element_text(size=15,face='bold'),
                      legend.text = element_text(face='bold'))

cols <- c("BASELINE"="blue","OUTCOME"="red")
sha <- c("BASELINE"=3,"OUTCOME"=4)

gg1 <- ggplot(dd,aes(y=mu)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[1],colour="BASELINE",shape="BASELINE"),size=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[1],colour="OUTCOME",shape="OUTCOME"), size=3,stroke =2) +
  xlab('') + ylab(expression(bold(hat(mu)))) + common.theme +
  scale_colour_manual(name="",values=cols) + scale_shape_manual(name="",values=sha)
gg2 <- ggplot(dd,aes(y=tau)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[2]), colour="blue",size=3,shape=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[2]), colour="red",size=3,shape=4,stroke =2) +
  xlab('') + ylab(expression(bold(hat(tau)))) + common.theme     # "\u03c4"
gg3 <- ggplot(dd,aes(y=I2)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[3]), colour="blue",shape=3,size=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[3]), colour="red",shape=4,size=3,stroke =2) +
  xlab('') + ylab(expression(bold(I^2))) + common.theme


ggarrange(gg1,gg2,gg3,nrow=1, common.legend = TRUE,legend = 'bottom')

#-------------------------------------------------------------------
#
#
# Jacknife                        
#
#
#-------------------------------------------------------------------
####################################################################
# Read data
####################################################################
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
data <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")


M <- matrix(ncol=3,nrow=208)
for (j in 1:208){
  print(j)
  n1 <- data$final_cases_T1[-j]
  n2 <- data$final_cases_T2[-j]
  v1 <- (data$final_sd_T1[-j])^2
  v2 <- (data$final_sd_T2[-j])^2
  yi <- log(v1/v2)
  vi <- 2/(n1-2)+2/(n2-2)  
  rma.model <- rma(yi=yi,vi=vi)
  M[j,] <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2)
}

colnames(M) <- c('mu','tau','I2')
M
graphics.off()
windows()
par(mfrow=c(1,3),las=1)
REAL.BASELINE <- c(0.02,0.31,100*0.583)
n1 <- data$final_cases_T1
n2 <- data$final_cases_T2
v1 <- (data$final_sd_T1)^2
v2 <- (data$final_sd_T2)^2
yi <- log(v1/v2)
vi <- 2/(n1-2)+2/(n2-2)  
rma.model <- rma(yi=yi,vi=vi)
REAL.FINAL <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2)
labs <- c(expression(bold(mu^'')),expression(bold(tau^'')),expression(I^2))
ymin <- c(-0.12,0.54,81)
ymax <- c(-0.08,0.66,84)
for (i in 1:3){
  boxplot(M[,i],main=labs[i],ylim=c(ymin[i],ymax[i]),cex.main=1.8)
  rect(-100,-100,100,100,border=NA, col='grey90')
  abline(h=pretty(c(ymin[i],ymax[i])),col='white')
  boxplot(M[,i],add=TRUE,col='white')
  points(1,REAL.FINAL[i],col=2,lwd=2,pch=4,cex=2)
  abline(h=mean(M[,i]),lty=2)
  
  legend('topleft',c(expression(bold(hat(theta[' ']))),NA,expression(bold(hat(theta)[('.')]))),
         bg = 'white',pch=c(4,NA,NA),col=c(2,NA,1),lty=c(NA,NA,2),cex=1.05,pt.cex=1.1,
         ncol=3,pt.lwd = 2)

}

summary(M)

#-------------------------------------------------------------------
#
#
# Normality                        
#
#
#-------------------------------------------------------------------

##-- Normality
par(mfrow=c(1,1))
with(data,qqnorm(2*log(final_sd_T1/final_sd_T2)))
with(data,qqline(2*log(final_sd_T1/final_sd_T2),col=2))

#-------------------------------------------------------------------
#
#
# Forest-plot                        
#
#
#-------------------------------------------------------------------
##-- Simulated forest under NULL hypotheses
# Data
set.seed(12345)
est0 <- c()
for(j in 1:208){
  n1 <- data$final_cases_T1[j]
  n2 <- data$final_cases_T2[j]
  v1 <- data$final_sd_T1[j]^2
  v2 <- data$final_sd_T1[j]^2
  y01 <- rnorm(n1,0,sqrt(v1))                 # Generate normal data for one arm
  y02 <- rnorm(n2,0,sqrt(v2))                 # Generate normal data for the other arm
  est0[j] <- log(var(y01)/var(y02))           # Response
}

# Graphic
graphics.off()
windows(15,10)
par(mfrow=c(1,3),las=1,cex.main=1.8)

##-- Forestplot simulated data  --------------------
# Data
ord <- order(est0)
est <- est0[ord]
se <- with(data,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- est-1.96*se
LS <- est+1.96*se
rma.mod1 <- rma(est,sei=se)

# Graphic
# plot(NA,xlim=c(-4.5,4.5),ylim=c(1,208),yaxt='n',ylab='',xlab='log(Vot/Voc)',
#      main=expression(bold('Simulated data under'~H[0])),sub=paste0('Tau=',round(sqrt(rma.mod1$tau2),2)))
# rect(-100,-20,100,1000,col='grey90',border=NA)
# abline(v=seq(-4,4,2),col='white',lty=2)
# abline(v=0,col='white',lwd=2)
# co <- (LI>0 | LS<0) + 1
# points(est,1:208,pch=15,cex=0.5,col=co)
# segments(LI,1:208,LS,1:208,col=co)

# GGplot2
dd <- data.frame(y=1:208,est=exp(est),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

gg1 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[BT]^2~"/"~S[BC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle(expression(bold('Simulated data under'~H[0]))) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))


  # stat_summary(geom = "pointrange",fun.ymax = LI)

##-- Forestplot baseline  --------------------
# Data
p0 <- with(data,2*log(base_sd_T1/base_sd_T2))
ord <- order(p0)
p <- p0[ord]
se <- with(data,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- p-1.96*se
LS <- p+1.96*se
rma.mod2 <- rma(p,sei=se)

# Graphic
# plot(NA,xlim=c(-4.5,4.5),ylim=c(1,208),yaxt='n',ylab='',xlab='log(Vot/Voc)',
#      main='Baseline data',sub=paste0('Tau=',round(sqrt(rma.mod2$tau2),2)))
# rect(-100,-20,100,1000,col='grey90',border=NA)
# abline(v=seq(-4,4,2),col='white',lty=2)
# abline(v=0,col='white',lwd=2)
# co <- (LI>0 | LS<0) + 1
# points(p,1:208,pch=15,cex=0.5,col=co)
# segments(LI,1:208,LS,1:208,col=co)

# GGplot2
dd <- data.frame(y=1:208,est=exp(p),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

gg2 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[BT]^2~"/"~S[BC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle('Baseline data') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))


##-- Forestplot outcome --------------------
# Data
p0 <- with(data,2*log(final_sd_T1/final_sd_T2))
ord <- order(p0)
p <- p0[ord]
se <- with(data,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- p-1.96*se
LS <- p+1.96*se
rma.mod2 <- rma(p,sei=se)

# Graphic
# plot(NA,xlim=c(-4.5,4.5),ylim=c(1,208),yaxt='n',ylab='',xlab='log(Vot/Voc)',
#      main='Outcome Data',sub=paste0('Tau=',round(sqrt(rma.mod2$tau2),2)))
# rect(-100,-20,100,1000,col='grey90',border=NA)
# abline(v=seq(-4,4,2),col='white',lty=2)
# abline(v=0,col='white',lwd=2)
# co <- (LI>0 | LS<0) + 1
# points(p,1:208,pch=15,cex=0.5,col=co)
# segments(LI,1:208,LS,1:208,col=co)

# GGplot2
dd <- data.frame(y=1:208,est=exp(p),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

gg3 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[OT]^2~"/"~S[OC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle('Outcome data') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))

ggarrange(gg1,gg2,gg3,nrow=1,common.legend = TRUE,legend = "bottom")