#-------------------------------------------------------------------
#
#
# Jacknife                        
#
#
#-------------------------------------------------------------------

rm(list=ls())
library(metafor)
library(ggplot2)
library(data.table)

####################################################################
# Read data
####################################################################
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
data <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")

####################################################################
# Jacknife
####################################################################
M <- matrix(ncol=3,nrow=208)
for (j in 1:208){
  print(j)
  n1 <- data$final_cases_T1[-j]
  n2 <- data$final_cases_T2[-j]
  v1 <- (data$final_sd_T1[-j])^2
  v2 <- (data$final_sd_T2[-j])^2
  yi <- log(v1/v2)
  vi <- 2/(n1-2) + 2/(n2-2)  
  rma.model <- rma(yi=yi,vi=vi)
  M[j,] <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2)
}

####################################################################
# Plot
####################################################################
REAL.FINAL <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2)
colnames(M) <- c('mu','tau','I2')
M
d <- melt(as.data.table(M),measure.vars = 1:3)
d$real_value <- ifelse(d$variable=='mu',REAL.FINAL[1],ifelse(d$variable=='tau',REAL.FINAL[2],REAL.FINAL[3]))
to_expression <- as_labeller(c('mu' = expression(mu), 'tau' = expression(tau), 'I2'=expression(I^2)), label_parsed)
ggplot(d,aes(y=value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = real_value),linetype=2,size=1) +
  stat_summary(aes(x=0,y=value,group=variable),fun.y = "mean", geom='point', colour = "red", shape=4, size=2, stroke = 2) +
  facet_wrap(~variable,scales = 'free',labeller = to_expression) +
  xlab('') + ylab('') +
  theme(axis.title = element_text(face='bold'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size=15))
ggsave('C:/Users/jcortes/Google Drive/Tesis/Memoria/Figures/jacknife_good.png')





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

for (i in 1:3){
  boxplot(M[,i],main=colnames(M)[i])#,ylim=c(ymin[i],ymax[i]))
  abline(h=mean(M[,i]),lty=2)
  points(1,REAL.FINAL[i],col=2,lwd=2,pch=4,cex=2)
  legend('bottomleft',c('Our final mod'),pch=4,col=2,cex=0.9)
}

summary(M)
