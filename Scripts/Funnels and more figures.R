#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#
# All scripts are in 'http://www-eio.upc.es/teaching/best/variability_data/'
# You can see them at:
# http://www-eio.upc.es/teaching/best/variability_data/Main.R
# http://www-eio.upc.es/teaching/best/variability_data/functions.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models_reduced_data.R
# http://www-eio.upc.es/teaching/best/variability_data/subgroups.R
#
#-----------------------------------------------------------------

##-- Remove objects in memory
rm(list=ls())

#-----------------------------------------------------------------
#
# Install and load packages
#
#-----------------------------------------------------------------
##-- Install packages
list.of.packages <- c('ggplot2','weights','catspec','alabama','metafor')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##-- Load packages
library(ggplot2)
library(weights)
library(catspec)
library(alabama)
library(metafor)

##-- Penalize scientific notation
options(scipen=3)

#-----------------------------------------------------------------
#
# Load specific functions
#
#-----------------------------------------------------------------
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/'
source(paste0(URL,'functions.R'))
source('C:/Users/jcortes/Google Drive/Tesis/1stPaper/Science/parameters.R')

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
datos1 <- read.table(url(paste0(URL,'data.csv')),header=TRUE,sep=";",stringsAsFactors = TRUE,quote = "")
closeAllConnections()

#-----------------------------------------------------------------
#
# Create new variables: outcomes and their standard errors
#
#-----------------------------------------------------------------
##-- Between arms --> Baseline log variance ratio
datos1$yBaselineRatio <- with(datos1,2*log(base_sd_T1/base_sd_T2))                               # Outcome
datos1$seBaselineRatio <- with(datos1,sqrt(2*(1/(base_cases_T1-2)+1/(base_cases_T2-2))))         # Standard error (n-2 = df-1 --> best approximation)


##-- Between arms  --> Outcome log variance ratio
datos1$yBetweenArmsRatio <- with(datos1,2*log(final_sd_T1/final_sd_T2))                          # Outcome
datos1$seBetweenArmsRatio <- with(datos1,sqrt(2*(1/(final_cases_T1-2)+1/(final_cases_T2-2))))    # Standard error (n-2 = df-1 --> best approximation)


##-- Over-time --> Experimental log variance ratio
datos1$yOverTimeRatioT <- with(datos1,2*log(final_sd_T1/base_sd_T1))                             # Outcome Treateds

cor.var <- datos1$rho^2                                                                          # Correlation between variances are aprox. the squared correlation between measures (97 available)
cov.log <- with(datos1,log(1 + 2*cor.var/(final_cases_T1-2)))                                    # Covariance between log variances
datos1$seOverTimeRatioT <- with(datos1,sqrt(2*(1/(final_cases_T1-2) +                            # Standard error (adding covariance) for over-time ratio
                                                 1/(final_cases_T1-2) -
                                                 cov.log)))       


##-- Over-time --> Reference log variance ratio
datos1$yOverTimeRatioC <- with(datos1,2*log(final_sd_T2/base_sd_T2))                             # Outcome Controls
datos1$seOverTimeRatioC <- with(datos1,sqrt(2*(1/(final_cases_T2-2) +                            # Standard error in control group
                                                 1/(final_cases_T2-2) -
                                                 cov.log)))

## Models
source(paste0(URL,'rma_models.R'))                # Models shown in table S1 with all data
source(paste0(URL,'rma_models_reduced_data.R'))   # Models shown in table S1 with reduced data   

#-----------------------------------------------------------------
#
# Relationship between significance and the outcome Between arms
#
#-----------------------------------------------------------------
d <- datos1[datos1$pvalue_exact=='Exact',]
yl1 <- expression(log~bgroup("(",frac(S[OT]^2,S[OC]^2),")"))
yl2 <- expression(bgroup('|',log~bgroup("(",frac(S[OT]^2,S[OC]^2),")"),'|'))

color <- rgb(0,0,1,0.5)

graphics.off()
windows(10,8)
par(mfrow=c(2,2),las=1,mar=c(5,7,5,1),font.lab=2,font.axis=4)
plot(2*log(final_sd_T1/final_sd_T2)~pvalue,d,log='x',pch=19,cex=1.2,
     col=color,
     ylab=yl1,xlab='P-value (log-scale)',
     main='Log variance ratio as function of p-value')
mod1 <- lm(2*log(final_sd_T1/final_sd_T2)~pvalue,d)
abline(mod1,col=2,lty=2)
abline(h=0,col=1,lty=2)
plot(abs(2*log(final_sd_T1/final_sd_T2))~pvalue,d,log='x',pch=19,cex=1.2,
     col=color,
     ylab=yl2,xlab='P-value (log-scale)',
     main='Absolute log variance ratio as function of p-value')
mod2 <- lm(abs(2*log(final_sd_T1/final_sd_T2))~pvalue,d)
abline(mod2,col=2,lty=2)
summary(mod1)
summary(mod2)

##-- Boxplots
boxplot(2*log(final_sd_T1/final_sd_T2)~significant,datos1,
        col=color,xlab='Statistically significant',ylab=yl1,
        main='Log variance ratio versus significance')
abline(h=0,lty=2)
boxplot(abs(2*log(final_sd_T1/final_sd_T2))~significant,datos1,
        col=color,xlab='Statistically significant',ylab=yl2,
        main='Absolute log variance ratio versus significance')

t.test(yBetweenArmsRatio~significant,datos1,var.equal=TRUE)
t.test(abs(yBetweenArmsRatio)~significant,datos1,var.equal=TRUE)


#-----------------------------------------------------------------
#
# Two funnel plots
#
#-----------------------------------------------------------------

##-- Labels
xl <- bquote(bold(frac(S[OT]^2,S[OC]^2)))
yl <- 'Standard error'


##-- Plot region
zoom <- TRUE
graphics.off()
windows(7,12)
if(zoom){xmax <- max(abs(exp(datos1$yBetweenArmsRatio)));ymax <- 1.02;ymin <- 0.05}else{
         xmax <- 100;ymax <- 1.1;ymin <- 0}
par(las=1,mfrow=c(2,1))
plot(seBetweenArmsRatio~yBetweenArmsRatio,datos1,pch=19,xlab='',ylab=yl,col=0,xaxt='n',
     main='Between arms (Non-significant studies)',xlim=log(c(1/xmax,xmax)),ylim=c(ymax,ymin))
mtext(xl,1,at=0,line=4,adj=0.5)
rect(-10,-10,300,100,col='grey85')
ticksy <- seq(0,1,0.2) # c(0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50)
abline(h=ticksy,lwd=2,col='white')
x0 <- c(4,0,-4,4) 
y0 <- c(2,0,2,2)
polygon(x0,y0,col='white',lty=3)
abline(v=0)
abline(v=coef(rma.unadj)[1],lwd=2,lty=2,col=4)
ticksx <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticksx),lab=ticksx)

##-- Points
d1 <- datos1[datos1$significant=='No',]
co <- with(d1,ifelse(pvalue>=0.05 | significant=='No',rgb(0,0,0,0.8),
                     ifelse(pvalue>0.001,rgb(227,66,52,255/2,maxColorValue = 255),
                            rgb(128,0,0,255/2,maxColorValue = 255))))#'#E34234','#800000'))))
pc <- with(d1,ifelse(significant=='No',1,19))
co[is.na(co)] <- ifelse(d1$significant[is.na(co)]=='No',rgb(0,0,0,0.8),rgb(128,0,0,255/2,maxColorValue = 255))


points(d1$yBetweenArmsRatio,d1$seBetweenArmsRatio,pch=pc,col=co,lwd=2,cex=1.1)

##-- More labels
if(zoom){
  mtext("Greater Treated",1,adj=1,at=log(xmax),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=1,at=log(xmax),line=3.5,cex=.8,font=2,las=0)
  mtext("Greater Control",1,adj=0,at=log(1/xmax),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=0,at=log(1/xmax),line=3.5,cex=.8,font=2,las=0) 
}else{
  mtext("Greater Treated",1,adj=1,at=log(100),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=1,at=log(100),line=3.5,cex=.8,font=2,las=0)
  mtext("Greater Control",1,adj=0,at=log(0.01),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=0,at=log(0.01),line=3.5,cex=.8,font=2,las=0)
}

##-- 2nd plot
plot(seOverTimeRatioT~yOverTimeRatioT,datos1,pch=19,xlab='',ylab=yl,col=0,xaxt='n',
     main='Between arms (significant studies)',xlim=log(c(1/xmax,xmax)),ylim=c(ymax,ymin))
mtext(xl,1,at=0,line=4,adj=0.5)
rect(-10,-10,300,100,col='grey85')
ticksy <- seq(0,1,0.2) # c(0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50)
abline(h=ticksy,lwd=2,col='white')
x0 <- c(4,0,-4,4) 
y0 <- c(2,0,2,2)
polygon(x0,y0,col='white',lty=3)
abline(v=0)
abline(v=coef(rma.unadj)[1],lwd=2,lty=2,col=4)
ticksx <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticksx),lab=ticksx)

##-- Points
d2 <- datos1[datos1$significant=='Yes',]
co <- with(d2,ifelse(pvalue>0.05 | significant=='No',rgb(0,0,0,0.8),
                     ifelse(pvalue>0.001,rgb(227,66,52,255/2,maxColorValue = 255),
                            rgb(128,0,0,255/2,maxColorValue = 255))))#'#E34234','#800000'))))
pc <- with(d2,ifelse(significant=='No',1,19))
co[is.na(co)] <- ifelse(d2$significant[is.na(co)]=='No',rgb(0,0,0,0.8),rgb(128,0,0,255/2,maxColorValue = 255))
points(d2$yBetweenArmsRatio,d2$seBetweenArmsRatio,pch=pc,col=co,lwd=2,cex=1.1)
legend('topright',c('0.001 < p < 0.05','p < 0.001'),
       pch=c(19,19),pt.lwd = 2,text.font = 2,pt.cex=1.3,
       co=c(rgb(227,66,52,255/2,maxColorValue = 255),rgb(128,0,0,255/2,maxColorValue = 255)))

##-- More labels
if(zoom){
  mtext("Greater Treated",1,adj=1,at=log(xmax),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=1,at=log(xmax),line=3.5,cex=.8,font=2,las=0)
  mtext("Greater Control",1,adj=0,at=log(1/xmax),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=0,at=log(1/xmax),line=3.5,cex=.8,font=2,las=0) 
}else{
  mtext("Greater Treated",1,adj=1,at=log(100),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=1,at=log(100),line=3.5,cex=.8,font=2,las=0)
  mtext("Greater Control",1,adj=0,at=log(0.01),line=2.7,cex=.8,font=2,las=0)
  mtext("Arm Variability",1,adj=0,at=log(0.01),line=3.5,cex=.8,font=2,las=0)
}


#------------------------------------------------------------------------------------------------------------------------

##-- Funnel (non-centered)
par(mfrow=c(1,1))
funnel(rma.unadj)

#-- Other funnels
par(mfrow=c(1,3))
funnel(rma.unadj)
funnel(rma.unadj,addtau2=TRUE)
funnel(rma.adj,addtau2=TRUE)



##-- Funnel (non-centered)
par(mfrow=c(1,1))
funnel(rma.unadj2)

############################################################
# Two funnels
############################################################
source('C:/Users/jcortes/Google Drive/Tesis/1stPaper/Science/Random effects model/funnel.R')







##-- Plot: Outcome versus baseline
windows(10,10)
par(las=1,mfrow=c(1,1),mar=c(6,8,2,1),font.axis=4,font.lab=2,mgp=c(4.5,1,0))
m <- sqrt(apply(datos1[,c('final_cases_T1','final_cases_T2',
                         'base_cases_T1','base_cases_T2')],1,mean))   # Point size (m) based on the sample size average
plot(yBetweenArmsRatio~yBaselineRatio,datos1,pch=1,
     xlab=expression(log~bgroup("(",frac(S[BT]^2,S[BC]^2),")")),
     ylab=expression(log~bgroup("(",frac(S[OT]^2,S[OC]^2),")")),
     cex=0.15*m,xlim=c(-1.6,1.6),ylim=c(-1.6,1.6))                    # Plot
abline(h=0,v=0,lty=2)
abline(0,1,lty=1)

##-- Plot: Outcome versus baseline (with ggplot2)
datos1$size <- m
ggplot(datos1,aes(x=yBaselineRatio,y=yBetweenArmsRatio,size=size)) + 
  geom_point(color='darkblue',alpha=0.5) +
  geom_vline(xintercept = 0,linetype=2) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_abline(intercept = 0,slope=1,linetype=2) +
  xlab(expression(bold(log~bgroup("(",frac(S[BT]^2,S[BC]^2),")")))) +
  ylab(expression(bold(log~bgroup("(",frac(S[OT]^2,S[OC]^2),")")))) +
  theme(legend.position = 'none',
        axis.title = element_text(face = 'bold',size=13),
        axis.text = element_text(face = 'bold',size=13))


##-- Correlation between outcome and baseline
with(datos1,cor(yBaselineRatio,yBetweenArmsRatio))

############################################################
# Previous calculations for funnels
############################################################
##-- Axis
y <- with(datos1,log(final_sd_T1)-log(final_sd_T2))
w <- with(datos1,sqrt((final_cases_T1+final_cases_T2)/2))
w <- with(datos1,sqrt((final_cases_T1+final_cases_T2)))

##-- Weighted t-test
wtd.t.test(y,weight=w)

##-- F statistic
Fest <- with(datos1,(final_sd_T1/final_sd_T2)^2)
p <- with(datos1,pf(Fest,final_cases_T1,final_cases_T2))
write.table(data.frame(pvalue=p),'pvalor.txt',row.names = FALSE,quote=FALSE)
LL <- with(datos1,qf(0.025,final_cases_T1,final_cases_T2))
UL <- with(datos1,qf(0.975,final_cases_T1,final_cases_T2))
x <- log(sqrt(UL))
sign <- Fest<LL | Fest>UL
signL <- Fest<LL
signU <- Fest>UL
cat('Studies with lower variability in Treated arm:',sum(signL),'\n')
cat('Studies with greater variability in Treated arm:',sum(signU),'\n')
cat('Studies with different variability between arms:',sum(sign),'\n')

############################################################
# Different funnel-plots to choose for print the triangle
############################################################
# Option 1: Statistic
# Option 2: Regression
# Option 3: Lowess
source('C:/Users/jcortes/Google Drive/Tesis/1stPaper/Science/choose_funnel.R')

############################################################
# Funnel plot between arms (horizontal)
############################################################
graphics.off()
windows(7,5)
par(las=1,mfrow=c(1,1),mar=c(5,6,2,1),font.axis=4,font.lab=2)
xl <- 'Uncertainty'
yl <- bquote(bold(frac(SD[O],SD[C])))
plot(1/w,y,pch=19,xlab=xl,ylab=yl,col=0,ylim=c(log(0.1),log(10)),yaxt='n') # xlim=c(0.03,0.42),
rect(-10,-10,300,100,col=colbg)
ticks <- c(0.1,0.2,0.5,1,2,5,10)
axis(2,at=log(ticks),lab=ticks)
abline(h=log(ticks),lwd=2,col='white')

xL <- xU <- 1/w
yL <- log(sqrt(LL))
yU <- log(sqrt(UL))
abline(h=0,col=4,lwd=2,lty=2)
abline(lm(yL~xL),col=4,lwd=2,lty=2)
abline(lm(yU~xU),col=4,lwd=2,lty=2)
points(1/w,y,pch=19,col=co1[sign+1])
mtext("Greater Treated",2,adj=1,at=log(11),line=3.3,cex=.8,font=2,las=0)
mtext("Arm Variability",2,adj=1,at=log(11),line=2.5,cex=.8,font=2,las=0)
mtext("Greater Control",2,adj=0,at=log(0.09),line=3.3,cex=.8,font=2,las=0)
mtext("Arm Variability",2,adj=0,at=log(0.09),line=2.5,cex=.8,font=2,las=0)

############################################################
# Funnel plot between arms (vertical)
############################################################
graphics.off()
windows(9,7)
par(las=1,mfrow=c(1,1),mar=c(6,7,2,1),font.axis=4,font.lab=2,mgp=c(4.2,1,0))
yl <- bquote(bold(Uncertainty~~~bgroup("(",bold(frac(1,sqrt(n[OT]+n[CT]))),")"))) # 'Uncertainty'
xl <- bquote(bold(frac(S[OT]^2,S[CT]^2)))
y2 <- log(exp(y)^2)
plot(y2,1/w,pch=19,xlab=xl,ylab=yl,col=0,xlim=c(log(0.01),log(100)),ylim=c(0.42,0.03),xaxt='n') 
rect(-10,-10,100,300,col=colbg)
ticks <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticks),lab=ticks)
ticks2 <- c(0.1,0.2,0.3,0.4)
abline(h=ticks2,lwd=2,col='white')

xL <- xU <- 1/w
yL <- log(LL)
yU <- log(UL)
modL <- lm(xL~yL)  # Model lower
modU <- lm(xU~yU)  # Model upper
xpoly <- c(log(0.01),log(1),log(100))
ypoly <- as.numeric(c(predict(modL,data.frame(yL=xpoly[1:2])),predict(modU,data.frame(yU=xpoly[3]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.42+ 0.04*(0.42-0.03)
segments(log(0.01),xaxis,log(100),xaxis,col=1,lwd=1,lty=1)
abline(v=0,col=1,lwd=1,lty=1)
points(y2,1/w,pch=19,col=co1[sign+1])
mtext("Greater Treated",1,adj=1,at=log(100),line=2.5,cex=.8,font=2,las=1)
mtext("Arm Variability",1,adj=1,at=log(100),line=3.3,cex=.8,font=2,las=1)
mtext("Greater Control",1,adj=0,at=log(0.01),line=2.5,cex=.8,font=2,las=1)
mtext("Arm Variability",1,adj=0,at=log(0.01),line=3.3,cex=.8,font=2,las=1)


##########################################################################
##-- Over-time (previous calculations)
##########################################################################
##-- Paired variance test (pag 331 de Lothar Sachs)  --> Only 95 papers with covariance

sum(!is.na(datos1$rho))
datos2 <- datos1[!is.na(datos1$rho) & datos1$rho<1 & datos1$rho>-1,]
sdx <- datos2$base_sd_T1
sdy <- datos2$final_sd_T1
n <- datos2$final_cases_T1
rho <- datos2$rho
Qest <- var.paired.test(sdx,sdy,n,rho)[[1]]
num <- var.paired.test(sdx,sdy,n,rho)[[2]]
den <- var.paired.test(sdx,sdy,n,rho)[[3]]
LI2 <- qt(0.025,n-2)
LS2 <- qt(0.975,n-2)
sign2 <- Qest<LI2 | Qest>LS2

cat('Studies with lower variability at the end of the sudy:',sum(Qest>LS2),'\n')
cat('Studies with greater variability at the end of the sudy:',sum(Qest<LI2),'\n')
cat('Studies with different variability over time:',sum(sign2),'\n')

############################################################
# Funnel plot over-time (horizontal) --> Function of the statistic
############################################################
##-- Funnel con Q
x2 <- 1/sqrt(n)
y2 <- -Qest
plot(x2,y2,col=co1[sign2+1],pch=19,ylim=c(-10,10),xlab=bquote(bold(frac(1,sqrt(n)))),ylab='Q')
rect(0,-15,1,15,col=colbg)
abline(h=seq(-10,10,2),col='white',lwd=2)
abline(lm(LI2~x2),lty=2,col=4,lwd=2)
abline(lm(LS2~x2),lty=2,col=4,lwd=2)
points(x2,y2,col=co1[sign2+1],pch=19)
#points(x2,LI2,pch='-')
#points(x2,LS2,pch='-')

##-- Funnel con Q (option 2 --> MAL)
# x2 <- den
# y2 <- with(datos2,log(final_sd_T1/base_sd_T1))
# plot(x2,y2,col=co1[sign2+1],pch=19,xlab=bquote(bold(frac(1,sqrt(n)))),ylab='Q',log='x') # ylim=c(-10,10),
# rect(0,-15,1,15,col=colbg)
# abline(h=seq(-10,10,2),col='white',lwd=2)
# abline(lm(LI2~x2),lty=2,col=4,lwd=2)
# abline(lm(LS2~x2),lty=2,col=4,lwd=2)
# points(x2,y2,col=co1[sign2+1],pch=19)
# #points(x2,LI2,pch='-')
# #points(x2,LS2,pch='-')

############################################################
# Funnel plot over-time (vertical) --> Function of the statistic
############################################################
##-- Funnel con Q
y2 <- 1/sqrt(n)
x2 <- -Qest
par(mgp=c(3,1,0))
plot(x2,y2,col=co1[sign2+1],pch=19,xlim=c(-10,10),ylim=c(0.35,0.05),ylab=bquote(bold("Uncertainty"~~~bgroup("(",frac(1,sqrt(n[BT])),")"))),xlab='Q')
rect(-15,0,15,1,col=colbg)
abline(h=seq(0.05,0.35,0.05),col='white',lwd=2)
ypoly <- c(0,.4,.4,0)
xpoly <- as.numeric(c(predict(lm(LI2~x2),data.frame(x2=ypoly[1:2])),predict(lm(LS2~x2),data.frame(x2=ypoly[3:4]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.35 + 0.04*(0.35-0.05)
segments(-10,xaxis,10,xaxis,col=1,lwd=1,lty=1)
abline(v=0,lty=1,col=1,lwd=1)
points(x2,y2,col=co1[sign2+1],pch=19)
mtext("Greater Outcome",1,adj=1,at=10,line=2.5,cex=.8,font=2,las=0)
mtext("Variability",1,adj=1,at=10,line=3.3,cex=.8,font=2,las=0)
mtext("Greater Baseline",1,adj=0,at=-10,line=2.5,cex=.8,font=2,las=0)
mtext("Variability",1,adj=0,at=-10,line=3.3,cex=.8,font=2,las=0)


############################################################
# Funnel plot over-time (horizontal) --> Function of variance ratio
############################################################
##-- Funnel con cociente de desviaciones
QLIM <- function(Qx,Qy,Qxy,n){
  num <- (Qx-Qy)*sqrt(n-2)
  aux <- Qx*Qy-Qxy^2
  den <- 2*sqrt(max(-aux,aux))
  Qest <- num/den
  return(Qest - qt(c(0.975,0.025),n-2))
}
QLIM1 <- function(Qx,Qy,Qxy,n) QLIM(Qx,Qy,Qxy,n)[1]
QLIM2 <- function(Qx,Qy,Qxy,n) QLIM(Qx,Qy,Qxy,n)[2]
n.datos2 <- nrow(datos2)
LS2 <- LI2 <- c()
for (i in 1:n.datos2){
  n <- datos2$final_cases_T1[i]
  Qy <- datos2$base_sd_T1[i]^2*(n-1)
  Qxy <- with(datos2,rho[i]*base_sd_T1[i]*final_sd_T1[i]*(n-1))
  
  roo <- uniroot(QLIM1,c(Qy/100,100*Qy),Qy=Qy,Qxy=Qxy,n=n)$root 
  LS2[i] <- sqrt(roo/(n-1))
  
  roo <- try(uniroot(QLIM2,c(Qy/1000,100*Qy),Qy=Qy,Qxy=Qxy,n=n)$root)
  if(class(roo)!='try-error') LI2[i] <- sqrt(roo/(n-1))
  
  print(i)
}
y2 <- with(datos2,log(final_sd_T1/base_sd_T1))
y_pos <- y2>=0
y_neg <- y2<0
x2_1 <- log(LS2/datos2$base_sd_T1)[y_pos] #1/sqrt(nest) # LS2/datos2$base_sd_T1
x2_2 <- -log(LI2/datos2$base_sd_T1)[y_neg]
x2<- log(LS2/datos2$base_sd_T1)

# x2 <- pmax(abs(log(LS2/datos2$base_sd_T1)),abs(log(LI2/datos2$base_sd_T1)))
# x2 <- log(LS2/datos2$base_sd_T1)-log(LI2/datos2$base_sd_T1)
# x2 <- 1/datos2$rho
x2 <- 1/sqrt(datos2$base_cases_T1)
par(las=1,mfrow=c(1,1),mar=c(5,6,2,1))
plot(x2,y2,pch=19,col=sign2+1,ylim=c(log(0.1),log(10)),yaxt='n',
     ylab=bquote(bold(frac(SD[O],SD[B]))),xlab='Uncertainty') # bquote(bold(log~bgroup("(",frac(SD[O],SD[B]),")")))
ticks <- c(0.1,0.2,0.5,1,2,5,10)
axis(2,at=log(ticks),lab=ticks)
rect(-10,-10,300,100,col=colbg)
abline(h=log(ticks),lwd=2,col='white')
#points(x2,log(LS2/datos2$base_sd_T1),col=4,pch='-')
#points(x2,log(LI2/datos2$base_sd_T1),col=4,pch='-')
# segments(x2,log(LS2/datos2$base_sd_T1),x2,log(LI2/datos2$base_sd_T1),lty=2,col='grey')
# a1 <- log(LS2/datos2$base_sd_T1)
# a2 <- log(LI2/datos2$base_sd_T1)
a1 <- log(LS2/datos2$base_sd_T1)[y_pos]
a2 <- log(LI2/datos2$base_sd_T1)[y_neg]
abline(lm(a1~x2_1),lty=2,lwd=2,col=4)
abline(lm(a2~x2_2),lty=2,lwd=2,col=4)
abline(h=0,lty=2,lwd=2,col=4)
points(x2_1,y2[y_pos],pch=19,col=co1[sign2+1][y_pos])
points(x2_2,y2[y_neg],pch=19,col=co1[sign2+1][y_neg])
mtext("Greater Outcome",2,adj=1,at=log(11),line=3.3,cex=.8,font=2,las=0)
mtext("Variability",2,adj=1,at=log(11),line=2.5,cex=.8,font=2,las=0)
mtext("Greater Baseline",2,adj=0,at=log(0.09),line=3.3,cex=.8,font=2,las=0)
mtext("Variability",2,adj=0,at=log(0.09),line=2.5,cex=.8,font=2,las=0)

############################################################
# Funnel plot over-time (vertical) --> Function of variance ratio
############################################################
##-- Funnel con cociente de desviaciones
yl <- bquote(bold(Uncertainty~~~bgroup("(",bold(frac(1,sqrt(n[BT]))),")"))) # 'Uncertainty'
xl <- bquote(bold(frac(S[OT]^2,S[BT]^2)))
y2 <- with(datos2,log(final_sd_T1^2/base_sd_T1^2))
y_pos <- y2>=0
y_neg <- y2<0
x2_1 <- log(LS2/datos2$base_sd_T1)[y_pos] #1/sqrt(nest) # LS2/datos2$base_sd_T1
x2_2 <- -log(LI2/datos2$base_sd_T1)[y_neg]
# x2<- log(LS2/datos2$base_sd_T1)
x2 <- 1/sqrt(datos2$base_cases_T1)

par(las=1,mfrow=c(1,1),mar=c(6,7,2,1),font.axis=4,font.lab=2,mgp=c(4.2,1,0))
plot(y2,x2,pch=19,col=sign2+1,xlim=c(log(0.01),log(100)),ylim=c(0.42,0.03),xaxt='n',
     xlab=xl,ylab=yl) # bquote(bold(log~bgroup("(",frac(SD[O],SD[B]),")")))
ticks <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticks),lab=ticks)
rect(-10,-10,100,300,col=colbg)
abline(h=seq(0.1,0.4,0.1),lwd=2,col='white')
a1 <- log(LS2/datos2$base_sd_T1)[y_pos]*2
a2 <- log(LI2/datos2$base_sd_T1)[y_neg]*2
xpoly <- c(log(0.01),log(1),log(100))
ypoly <- -as.numeric(c(predict(lm(x2_1~a1),data.frame(a1=xpoly[1:2])),predict(lm(x2_2~a2),data.frame(a2=xpoly[3]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.42+ 0.04*(0.42-0.03)
segments(log(0.01),xaxis,log(100),xaxis,col=1,lwd=1,lty=1)
abline(lm(x2_1~a1),lty=3,lwd=1,col=1)
abline(lm(x2_2~a2),lty=3,lwd=1,col=1)
abline(v=0,lty=1,lwd=1,col=1)
# points(y2[y_pos],x2_1,pch=19,col=co1[sign2+1][y_pos])
# points(y2[y_neg],x2_2,pch=19,col=co1[sign2+1][y_neg])
points(y2[y_pos],x2[y_pos],pch=19,col=co1[sign2+1][y_pos])
points(y2[y_neg],x2[y_neg],pch=19,col=co1[sign2+1][y_neg])
mtext("Greater Outcome",1,adj=1,at=log(100),line=2.5,cex=.8,font=2,las=1)
mtext("Variability",1,adj=1,at=log(100),line=3.3,cex=.8,font=2,las=1)
mtext("Greater Baseline",1,adj=0,at=log(0.01),line=2.5,cex=.8,font=2,las=1)
mtext("Variability",1,adj=0,at=log(0.01),line=3.3,cex=.8,font=2,las=1)

## Non-significative study with less variance ratio
which.min(y2[!sign2])
exp(y2[!sign2][11])
datos2[!sign2,][11,] # rho 0.306

## Significative with variance ratio closest to 1. It's by chance to obtain the same number (11)
which.max(y2[sign2 & y2<0])
exp(y2[sign2 & y2<0][11])
datos2[sign2 & y2<0,][11,] # rho 0.997
max(datos2$rho)

##########################################################################
##-- Over-time comparing between arms
##########################################################################
##-- Axis
y <- with(datos1,log(change_sd_T1)-log(change_sd_T2))
w <- with(datos1,sqrt((final_cases_T1+final_cases_T2)/2))

##-- Weighted t-test
wtd.t.test(y,weight=w)

##-- F statistic
Fest <- with(datos1,(change_sd_T1/change_sd_T2)^2)
p <- with(datos1,pf(Fest,final_cases_T1,final_cases_T2))
write.table(data.frame(pvalue=p),'pvalor.txt',row.names = FALSE,quote=FALSE)
LL <- with(datos1,qf(0.025,final_cases_T1,final_cases_T2))
UL <- with(datos1,qf(0.975,final_cases_T1,final_cases_T2))
x <- log(sqrt(UL))
sign <- Fest<LL | Fest>UL
signL <- Fest<LL
signU <- Fest>UL

##-- Elegir funnel plot
source('choose_funnel.R')

##-- Funnel plot between arms
graphics.off()
windows(7,5)
par(las=1,mfrow=c(1,1),mar=c(5,6,2,1),font.axis=4,font.lab=2)
xl <- bquote(bold(frac(1,sqrt(n))))
yl <- bquote(bold(log~bgroup("(",frac(SD[T],SD[C]),")")))
plot(1/w,y,pch=19,xlab=xl,ylab=yl,xlim=c(0.03,0.42),col=0)
rect(-10,-10,300,100,col='grey85')
abline(h=seq(-2,2,0.5),lwd=2,col='white')
xL <- xU <- 1/w
yL <- log(sqrt(LL))
yU <- log(sqrt(UL))
abline(h=0,col=4,lwd=2,lty=2)
abline(lm(yL~xL),col=4,lwd=2,lty=2)
abline(lm(yU~xU),col=4,lwd=2,lty=2)
points(1/w,y,pch=19,col=co1[sign+1])

#-----------------------------------------------------------------
#
# Mosaic plot
#
#-----------------------------------------------------------------
par(mfrow=c(1,2))
sign11 <- ifelse(Fest<LL | Fest>UL,'Significant','Non-significant')
sign22 <- ifelse(Fest<LL,'Control greater',ifelse(Fest>UL,'Treated greater','Non evidence'))
mosaicplot(datos1$Area_WOS~sign11,col=1:2,main='',las=2,xlab='',ylab='')
mosaicplot(datos1$Area_WOS~sign22,col=c(3,1:2),main='',las=2,xlab='',ylab='')


#-----------------------------------------------------------------
#
# Write data
#
#-----------------------------------------------------------------
datos.low <- datos1[Fest<LL,]
datos.up <- datos1[Fest>UL,]
# write.table(datos.low,'datoslow.txt',sep='\t',quote=FALSE,row.names=FALSE)
# write.table(datos.up,'datosup.txt',sep='\t',quote=FALSE,row.names=FALSE)

#-----------------------------------------------------------------
#
# P-values
#
#-----------------------------------------------------------------
##########################################################################
##-- Histogram with different breaks
##########################################################################
graphics.off()
windows(10,10)
par(mfrow=c(2,2),las=1,font.lab=2,font.axis=4,cex.lab=1,cex.axis=1)
hist(p,freq=F,br=5,main='Histogram of p-value (br=5)',col=4,border='grey50')
abline(h=1,lwd=2,col=2)
hist(p,freq=F,br=10,main='Histogram of p-value (br=10)',col=4,border='grey50')
abline(h=1,lwd=2,col=2)
hist(p,freq=F,br=15,main='Histogram of p-value (br=15)',col=4,border='grey50')
abline(h=1,lwd=2,col=2)
ru <- runif(1000)
qqplot(ru,p,xlab='Theorethical quantiles',pch=1,col=4,
       ylab='Sample quantiles',main='QQplot for Uniform distribution')
qqline(p,distribution = function(p) qunif(p,0,1),col=2,lty=2)

##########################################################################
##-- Histogram with probability bars
##########################################################################
graphics.off()
windows(8,5)
par(mfrow=c(1,2),las=1,font.lab=2,font.axis=4,cex.lab=1,cex.axis=1,mar=c(5,5,5,1))
hist(p,freq=F,br=10,main='Histogram',col=4,border='grey50',xlab='p-value')
abline(h=1,lwd=2,col=2)
## CI95%
n <- 208;p1 <- 1/10;q1 <- 1 - p1
V <- n*p1*q1
s <- sqrt(V)
CI95 <- (20.8 + c(-1,1)*2*s)/20.8
abline(h=CI95,lwd=2,col=2,lty=2)

##-- Modelo escogido
var = list(datos=p)   ######### pon aqu? los p valores
P0 = c(0.5, 1, 1)
ans1 = auglag(par=P0, fn=fn1, heq=heq1,hin=hin1, hin.jac=hin.jac1, control.outer=list(trace=FALSE), var=var)
ans1
results(p,ans1,nparam=3,F1,f1,tit='Uniform and Beta distributions',hist=FALSE)

#-----------------------------------------------------------------
#
# Four standard deviations at time
#
#-----------------------------------------------------------------
OT <- datos1$final_sd_T1 
OC <- datos1$final_sd_T2
BT <- datos1$base_sd_T1
BC <- datos1$base_sd_T2
x <- log(OT/BT)
y <- log(OT/OC)

datos1$sign2 <- NA
sel <- with(datos1,!is.na(rho) & abs(rho)<=1)
datos1$sign2[sel] <- sign2
sign2bis <- datos1$sign2
co4 <- ifelse(sign & sign2bis & !is.na(sign2bis),co3[1],
              ifelse(sign & (!sign2bis | is.na(sign2bis)),co3[2],
                     ifelse(!sign & sign2bis & !is.na(sign2bis),co3[3],co3[4])))
pc1 <- ifelse(is.na(datos1$rho) & sign,2,
              ifelse(is.na(datos1$rho) & !sign,1,
                     ifelse(sign & sign2bis & !is.na(sign2bis),15,
                            ifelse(sign & (!sign2bis | is.na(sign2bis)),17,
                                   ifelse(!sign & sign2bis & !is.na(sign2bis),18,19)))))
pc2 <- ifelse(is.na(datos1$rho),NA,
              ifelse(sign & sign2bis & !is.na(sign2bis),15,
                     ifelse(sign & (!sign2bis | is.na(sign2bis)),17,
                            ifelse(!sign & sign2bis & !is.na(sign2bis),18,19))))

#-----------------------------------------------------------------
#
# Normal plot
#
#-----------------------------------------------------------------
graphics.off()
windows(8,8)
par(mfrow=c(1,1),las=1,mar=c(6,6,5,1))
plot(x,y,pch=pc1,col=co4,main='All studies',cex=1.1,
     xlab=bquote(bold(log~bgroup("(",frac(SD[TO],SD[TB]),")"))),
     ylab=bquote(bold(log~bgroup("(",frac(SD[TO],SD[CO]),")"))))
mtext('NS:Non-Significant',1,at=-2.2,line=2,adj=0,cex=0.7)
mtext('S:Significant',1,at=-2.2,line=2.5,adj=0,cex=0.7)
mtext('NI:No Information',1,at=-2.2,line=3,adj=0,cex=0.7)
legend('topleft',c('NS between arms / NI over-time',
                   'S between arms / NI over-time',
                   'NS between arms / NS over-time',
                   'NS between arms / S over-time',
                   'S between arms / NS over-time',
                   'S between arms / S over-time'),
       col=c(co3[4],co3[2],co3[4],co3[3],co3[2],co3[1]),
       pch=c(1,2,19,18,17,15),pt.cex=1.3)
abline(v=0,lty=2)
abline(h=0,lty=2)

windows(8,8)
par(mfrow=c(1,1),las=1,mar=c(6,6,5,1))
plot(x,y,pch=pc2,col=co2,main='Studies with covariance',cex=1.2,
     xlab=bquote(bold(log~bgroup("(",frac(SD[TO],SD[TB]),")"))),
     ylab=bquote(bold(log~bgroup("(",frac(SD[TO],SD[CO]),")"))))
abline(v=0,lty=2)
abline(h=0,lty=2)
mtext('NS:Non-Significant',1,at=-2.2,line=2,adj=0,cex=0.7)
mtext('S:Significant',1,at=-2.2,line=2.5,adj=0,cex=0.7)
legend('topleft',c('NS between arms / NS over-time',
                   'NS between arms / S over-time',
                   'S between arms / NS over-time',
                   'S between arms / S over-time'),
       col=c(co3[4],co3[3],co3[2],co3[1]),
       pch=c(19,18,17,15),pt.cex=1.3)
#-----------------------------------------------------------------
#
# ggplot
#
#-----------------------------------------------------------------
dgg <- data.frame(variance_ratio_time=exp(x)^2,variance_ratio_arms=exp(y)^2,
          precision=1/sqrt(1/datos1$final_cases_T1 + 1/datos1$final_cases_T2),        
          significance_time=as.character(with(datos1,yOverTimeRatioT > (2*seOverTimeRatioT) | yOverTimeRatioT < (-2*seOverTimeRatioT))),
          significance_arms=as.character(with(datos1,yBetweenArmsRatio < (-2*seBetweenArmsRatio) | yBetweenArmsRatio> (2*seBetweenArmsRatio))),
          stringsAsFactors = FALSE)
dgg$significance_time[is.na(dgg$significance_time)] <- "NA"
dgg$significance_time <- factor(dgg$significance_time,levels=c('NA','FALSE','TRUE'))
dgg$Significance <- factor(with(dgg,ifelse(significance_time=='NA' & significance_arms=='FALSE','NS between arms/NI over-time',
                                    ifelse(significance_time=='NA' & significance_arms=='TRUE','S between arms/NI over-time',
                                           ifelse(significance_time=='FALSE' & significance_arms=='FALSE','NS between arms/NS over-time',
                                                  ifelse(significance_time=='FALSE' & significance_arms=='TRUE','S between arms/NS over-time',
                                                         ifelse(significance_time=='TRUE' & significance_arms=='FALSE','NS between arms/S over-time',
                                                                                                                      'S between arms/S over-time')))))),
                          levels=c('NS between arms/NI over-time','S between arms/NI over-time','NS between arms/NS over-time',
                                   'S between arms/NS over-time','NS between arms/S over-time','S between arms/S over-time'))


ggplot(dgg,aes(x=variance_ratio_time,y=variance_ratio_arms)) + 
  geom_point(alpha=0.8,mapping=aes(color=Significance,size=precision)) + 
  geom_smooth(method='lm',se=FALSE,linetype=2) +
  scale_x_log10(limits=c(0.01,100)) + scale_y_log10(limits=c(0.01,100)) +
  xlab(expression(bold(frac(S[OT]^2,S[BT]^2)))) + ylab(expression(bold(frac(S[OT]^2,S[OC]^2)))) +
  geom_vline(xintercept = 1,linetype=2) + geom_hline(yintercept = 1,linetype=2) +
  theme(legend.position='bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.7,'cm'),
        legend.text = element_text(face='bold')) +
  scale_color_manual(values=c("white","pink","yellow","red","purple", "black")) +
  guides(size=FALSE,color = guide_legend(override.aes = list(size = 3)))
with(dgg,cor(variance_ratio_time,variance_ratio_arms))
with(dgg,cor(variance_ratio_time,variance_ratio_arms,method = 'spe'))

table(dgg$Significance)
prop.table(table(dgg$Significance))

#-----------------------------------------------------------------
#
# Bland-Altman
#
#-----------------------------------------------------------------

##########################################################################
##-- Bland-Altman (Baseline - Outcome) --> Treated
##########################################################################
datos2 <- datos1

qplot(log(base_sd_T1)+log(final_sd_T1),
      log(final_sd_T1)-log(base_sd_T1),data=datos1,log='')

qplot(sqrt(final_sd_T1*base_sd_T1),
      final_sd_T1/base_sd_T1,data=datos1,log='xy')


sg <- datos1[,c('Intervention_type','Outcome_type','Condition_type','Measurement_type')]
sg$All <- rep(NA,nrow(sg))
p <- length(sg)

ylim <- log(8)

##-- Tractats
graphics.off()
pdf('Comparing Outcome and Baseline in Experimental arm.pdf',10,7,paper='USr')
MT <- matrix(NA,nrow=0,ncol=3)
for (i in 1:p){
  #windows(11,7)
  m <- with(datos1,bland.altman(v1=base_sd_T1,v2=final_sd_T1,z=factor(sg[,i]),
                                s1=base_cases_T1,s2=final_cases_T1,
                                l1='TO',l2='TB',sub='Final vs. Baseline (Experimental arms)',
                                log=TRUE,dup.axis=FALSE,ylim=ylim,rm=0,main=sub('_',' ',names(sg)[i])))
  MT <- rbind(MT,m)
}
sum(datos1$final_sd_T1>datos1$base_sd_T1)/208
prop.test(96,208)

#windows(11,7)
forestM(M=MT,xl=bquote(bold(SD[TO]/SD[TB])),lab1,lab2)
dev.off()

##-- Controls
graphics.off()
pdf('Comparing Outcome and Baseline in Reference arm.pdf',10,7,paper='USr')
MC <- matrix(NA,nrow=0,ncol=3)
for (i in 1:p){
  #windows(10,7)
  m <- with(datos1,bland.altman(v1=base_sd_T2,v2=base_sd_T2,z=factor(sg[,i]),
                                s1=base_cases_T2,s2=final_cases_T2,
                                l1='CO',l2='CB',sub='Final vs. Baseline (Reference arms)',
                                log=TRUE,dup.axis=FALSE,ylim=ylim,rm=0,main=sub('_',' ',names(sg)[i])))
  MC <- rbind(MC,m)
}

#windows(10,7)
forestM(M=MC,xl=bquote(bold(SD[CO]/SD[CB])),lab1,lab2)
dev.off()


##-- Placebos
datos.placebo <- subset(datos1,grepl('Placebo',Treatment_2) | grepl('placebo',Treatment_2))
graphics.off()

MP <- matrix(NA,nrow=0,ncol=3)
for (i in 1:p){
  windows(10,7)
  m <- with(datos.placebo,bland.altman(v1=base_sd_T2,v2=final_sd_T2,z=factor(sg[,i]),
                                       s1=base_cases_T2,s2=final_cases_T2,
                                       l1='PO',l2='PB',sub='Final vs. Baseline (Placebo arms)',
                                       log=TRUE,dup.axis=FALSE,ylim=ylim,rm=0,main=sub('_',' ',names(sg)[i])))
  MP <- rbind(MP,m)
}

##-- Tractats i controls (forest)
graphics.off()
pdf('Comparing Outcome and Baseline in both groups.pdf',10,7,paper='USr')
forestM(M=MT,xl=bquote(bold(SD[O]/SD[B])),lab1,lab2,add=FALSE,sep=0.15,col=c(rep(c(NA,pal[1:2]),4),NA,1),lty=1)
forestM(M=MC,xl=bquote(bold(SD[O]/SD[B])),lab1,lab2,add=TRUE,sep=-0.15,col=c(rep(c(NA,pal[1:2]),4),NA,1),lty=2)
legend('topright',c('Treated','Control'),col=1,lwd=2,lty=1:2,text.font = 2)
dev.off()


##########################################################################
##-- Bland-Altman (Treatment - Control)
##########################################################################
##-- Tractats
graphics.off()
pdf('Comparing Experimental and Reference outcome_1.pdf',10,7,paper='USr')
M <- matrix(NA,nrow=0,ncol=3)
for (i in 1:p){
  #windows(10,7)
  m <- with(datos1,bland.altman(v1=final_sd_T2,v2=final_sd_T1,z=factor(sg[,i]),
                                s1=final_cases_T2,s2=final_cases_T1,
                                l1='TO',l2='CO',sub='Experimental vs. control outcomes',
                                log=TRUE,dup.axis=FALSE,ylim=ylim,rm=0,main=sub('_',' ',names(sg)[i])))
  M <- rbind(M,m)
}
sum(datos1$final_sd_T2>datos1$final_sd_T1)/208
prop.test(117,208)

windows(10,7)
lab1[5:6] <- c('Measured','Scored') 
forestM(M=M,xl=bquote(bold(SD[TO]/SD[CO])),lab1,lab2)
dev.off()

##-- F-test
pvalue <- c()
n <- nrow(datos1)
for(i in 1:n){
  Fsta <- with(datos1,final_sd_T1[i]^2/final_sd_T2[i]^2)
  pr <- with(datos1,pf(Fsta,base_cases_T1[i]-1,base_cases_T2[i]-1))
  pvalue[i] <- 2*min(pr,1-pr)
}
sum(pvalue<0.05)/n
with(datos1,sum(final_sd_T1[pvalue<0.05]>final_sd_T2[pvalue<0.05]))
with(datos1,sum(final_sd_T1[pvalue<0.05]<final_sd_T2[pvalue<0.05]))
##########################################################################
##-- Canvi Final-Basal Treated vs Controls
##########################################################################


graphics.off()
pdf('Change Final-Baseline comparing Experimental and Reference_sin_un_punto.pdf',10,7,paper='USr')
M <- matrix(NA,nrow=0,ncol=3)
for (i in 1:p){
  #windows(10,7)
  m <- with(datos1[-154,],bland.altman(v1=final_sd_T2/base_sd_T2,v2=final_sd_T1/base_sd_T1,z=factor(sg[,i]),
                                       s1=final_cases_T2,s2=final_cases_T1,
                                       l1=c('TO','TB'),l2=c('CO','CB'),sub='Final/Baseline comparing Experimental vs. Reference',
                                       log=TRUE,dup.axis=FALSE,ylim=ylim,rm=0,main=sub('_',' ',names(sg)[i])))
  M <- rbind(M,m)
}

windows(10,7)
par(mgp=c(3,1,0))
forestM(M=M,xl='',lab1,lab2)
mtext(bquote(bold(frac(SD[TO]/SD[TB],SD[CO]/SD[CB]))),1,at=1,line=4,adj=0.5)
dev.off()
