######################################
# Triangle
######################################

##-- Triangle density (general)
x0 <- c(0,6,10,0)
y0 <- c(0,6,0,0)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,10),ylim=c(0,10),xaxt='n',xlab='',ylab='f(x)')
polygon(x0,y0,border=NA,col='grey90')
segments(0,0,6,6,lwd=2)
segments(6,6,10,0,lwd=2)
segments(6,0,6,6,lwd=1,lty=2)
axis(1,at=c(0,6,10),lab=c('a','c','b'))



##-- Triangle density (for our study)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,10),ylim=c(0,10),xaxt='n',xlab='',ylab='f(x)')
x0 <- c(0.05,0.05,2,0.05)
y0 <- c(0,6,0,0)
polygon(x0,y0,border=1,col='grey90',lty=2,lwd=2)

x0 <- c(9.95,9.95,7,9.95)
y0 <- c(0,4,0,0)
polygon(x0,y0,border=1,col='grey80',lty=3,lwd=2)

axis(1,at=c(0,2,7,10),lab=c('0','b1','a2','1'))

text(1.2,3,expression(bold(f[T]^1~(x))),adj=0)
text(7.6,3,expression(bold(f[T]^2~(x))),adj=0)

######################################
# Exponential
######################################

##-- Exponential (general)
par(mar=c(2.5,2.5,1,1),mgp=c(1,1,0))
xx <- seq(0,10,0.01)
x0 <- c(0,xx,10)
y0 <- c(0,dexp(xx,1/2),0)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,10),ylim=c(0,0.5),xaxt='n',xlab='',ylab='f(x)')
polygon(x0,y0,border=1,col='grey90',lwd=1.5)
axis(1,at=c(0,10),lab=c('0','1'))



##-- Exponential (inverted)
dexp2 <- function(x,mu,rate) rate*exp(-rate*(mu-x)) 
xx <- seq(0,10,0.01)
x1 <- c(0,xx,10)
y1 <- c(0,dexp2(xx,10,1/2),0)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,10),ylim=c(0,0.5),xaxt='n',xlab='',ylab='f(x)')
polygon(x1,y1,border=1,col='grey90',lwd=1.5)
axis(1,at=c(0,10),lab=c('0','1'))

##-- Exponential (both)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,10),ylim=c(0,0.5),xaxt='n',xlab='',ylab='f(x)')
polygon(x0,y0,border=1,lty=2,lwd=1.5,col=rgb(0,0,0,0.1))
polygon(x1,y1,border=1,lty=3,lwd=1.5,col=rgb(0,0,0,0.1))
axis(1,at=c(0,10),lab=c('0','1'))
text(1.5,0.25,expression(bold(f[E]^1~(x))),adj=0)
text(8.5,0.25,expression(bold(f[E]^2~(x))),adj=1)

######################################
# Beta
######################################

##-- Beta (one beta)
par(mar=c(2.5,2.5,1,1),mgp=c(1,1,0))
xx <- seq(0.001,0.999,0.001)
x0 <- c(0,xx,1)
y0 <- c(0,dbeta(xx,1/2,1/2),0)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,1),ylim=c(0,2),xaxt='n',xlab='',ylab='f(x)',xaxs="i")
polygon(x0,y0,border=1,col='grey90',lwd=1.5)
axis(1,at=c(0,1),lab=c('0','1'))



##-- Beta (two betas)
par(mar=c(2.5,2.5,1,1),mgp=c(1,1,0))
xx <- seq(0.001,0.999,0.001)
x0 <- c(0,xx,1)
y0 <- c(0,dbeta(xx,1,4),0)
x1 <- c(0,xx,1)
y1 <- c(0,dbeta(xx,4,1),0)
par(las=1,cex.lab=1,cex.axis=1,font.lab=2,font.axis=4,yaxs='i',yaxt='n')
plot(NA,xlim=c(0,1),ylim=c(0,4),xaxt='n',xlab='',ylab='f(x)',xaxs="i")
polygon(x0,y0,border=1,col=rgb(0,0,0,0.1),lwd=1.5,lty=2)
polygon(x1,y1,border=1,col=rgb(0,0,0,0.1),lwd=1.5,lty=3)
axis(1,at=c(0,1),lab=c('0','1'))
text(.25,2,expression(bold(f[beta]^1~(x))),adj=0)
text(.75,2,expression(bold(f[beta]^2~(x))),adj=1)
