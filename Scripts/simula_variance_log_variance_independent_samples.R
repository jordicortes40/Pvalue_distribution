##-- Remove all objects
rm(list=ls())

#-------------------------------------------------------------------
#
#
#  Simulation A                         
#
#
#-------------------------------------------------------------------

####################################################################
# Parameters of the simulation
####################################################################

nsim <- 10000                     # Number of simulations
N1 <- c(30,300,3000)              # Sample size in one group
FRAC <- c(1/3,1/2,1)              # Ratio between arms
V1 <- V2 <- 1                     # Actual variances

####################################################################
# Store information
####################################################################
M <- matrix(ncol=5,nrow=length(N1)*length(FRAC)) # Matrix to store infiormation
j <- 1                                           # Counter

####################################################################
# Simulation
####################################################################
set.seed(12345)
for(n1 in N1){
  for (frac in FRAC){
    cat('n1:',n1,'frac:',frac,'\n')              # Print the process
    n2 <- n1*frac                                # Sample size in second group
    
    v1 <- v2 <- est1 <- est2 <- c()              # Arrays for storing information
    
    for (i in 1:nsim){
      
      ##-- Method 1-----------------------------------
      y1 <- rnorm(n1,0,sqrt(V1))                 # Generate normal data for one arm
      y2 <- rnorm(n2,0,sqrt(V2))                 # Generate normal data for the other arm
      v1[i] <- var(y1)                           # Sample variance for one arm
      v2[i] <- var(y2)                           # Sample variance for the other arm
      est1[i] <- log(v1[i]/v2[i])                # First estimation
      
      ##-- Method 2-----------------------------------
      v1 <- rchisq(1,df=n1-1)*V1/(n1-1)          # Sample variance for one arm          
      v2 <- rchisq(1,df=n2-1)*V2/(n2-1)          # Sample variance for the other arm
      est2[i] <- log(v1/v2)                      # Second estimation
      
    }
    
    E1 <- sd(est1)                               # First estimation: E1
    E2 <- sd(est2)                               # Second estimation: E2
    E3 <- sqrt(2*(1/(n1-1)+1/(n2-1)))            # Estimation from delta method: DM
    
    ##-- Store information and update the counter
    M[j,] <- c(E3,E1,E1/E3,E2,E2/E3)
    j <- j+1
  }

}
colnames(M) <- c('SE Real','SE_from_var','E1divE3','SE_from_norm','E2divE3')
rownames(M) <- paste0(rep(N1,each=3),'_',round(FRAC,2))
M

#-------------------------------------------------------------------
#
#
#                          Simulation B
#
#
#-------------------------------------------------------------------
####################################################################
# Parameters of the simulation
####################################################################

nsim <- 10000                                # Number of simulations
N <- exp(seq(log(12),log(6144/4),length=8))  # Total sample size
FRAC <- c(1/2,1/3,1/4)                       # Fraction in one group over the total sample size
V1 <- V2 <- 1                                # Actual variances

####################################################################
# Store information
####################################################################
M1 <- M2 <- M3 <- matrix(ncol=7,nrow=length(N))
M <- list(M1,M2,M3) 

####################################################################
# Simulation
####################################################################
set.seed(12345)
for(n in N){
  for (frac in FRAC){
    cat('n:',n,'frac:',round(frac,2),'\n')
    n2 <- n*frac                                 # Sample size in second group
    n1 <- n-n2                                   # Sample size in first group   
    
    v1 <- v2 <- est1 <- est2 <- c()              # Arrays for storing information
    
    for (i in 1:nsim){
      
      ##-- Method 1-----------------------------------
      y1 <- rnorm(n1,0,sqrt(V1))                 # Generate normal data for one arm
      y2 <- rnorm(n2,0,sqrt(V2))                 # Generate normal data for the other arm
      v1[i] <- var(y1)                           # Sample variance for one arm
      v2[i] <- var(y2)                           # Sample variance for the other arm
      est1[i] <- log(v1[i]/v2[i])                # First estimation
      
      ##-- Method 2-----------------------------------
      v1 <- rchisq(1,df=n1-1)*V1/(n1-1)          # Sample variance for one arm          
      v2 <- rchisq(1,df=n2-1)*V2/(n2-1)          # Sample variance for the other arm
      est2[i] <- log(v1/v2)                      # Second estimation
      
    }
    
    E1 <- sd(est1)                               # First estimation: E1
    E2 <- sd(est2)                               # Second estimation: E2
    E3 <- sqrt(2*(1/(n1-1)+1/(n2-1)))            # Estimation from delta method: DM 
    E4 <- sqrt(2*(1/(n1-2)+1/(n2-2)))            # Estimation from corrected delta method: DM' (n-2) in the denominator instead of (n-1)
    
    ##-- Store information
    i <- which(N==n)
    j <- which(FRAC==frac)  
    M[[j]][i,] <- c(E3,E1,E1/E3,E1/E4,E2,E2/E3,E2/E4)
  }
}

##-- Print M
M

#-------------------------------------------------------------------
#
#
#                         Graphic                         
#
#
#-------------------------------------------------------------------
####################################################################
# Read data (for obtaining sample sizes)
####################################################################

URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
datos <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")


####################################################################
# Graphic
####################################################################
##-- Parameters for the graphic
TIT <- paste('Ratio',c('1:1','2:1','3:1'),'between arms')    # Titles
tck0 <- c(1.02,1.05,1.1,1.2,1.3)                             
tck <- c(1/tck0,1,tck0)                                      # Ticks for Y-axis

##-- Create window
graphics.off()
windows(9.5,4)
par(mfrow=c(1,3),las=1)

##-- Plot
for(i in 1:3){
  
  ##-- Background
  plot(NA,xlim=c(12,1500),ylim=c(3/4,4/3),log='xy',yaxt='n',
       xlab='Total Sample Size',ylab='Ratio of of estimations',main=TIT[i])
  rect(0.1,0.1,10^5,10,col='grey85')
  abline(h=1,lty=1,col='white',lwd=1)
  abline(h=tck,lty=2,col='white',lwd=1)
  
  ##-- Vertical axis
  axis(2,at=tck,lab=paste(c(100*(tck0-1),0,100*(tck0-1)),'%'))
  
  ##-- Marks in upper axis with actual sample sizes
  sel <- with(datos,round(pmax(final_cases_T1,final_cases_T2)/pmin(final_cases_T1,final_cases_T2)))==i
  cs <- with(datos[sel,],final_cases_T1+final_cases_T2)
  mtext('|',3,at=jitter(cs),line=0)
  
  ##-- Plot lines
  for (j in c(3,6)){
    lines(N,M[[i]][,j],col=j-1,lwd=2,lty=2)   # /DM
    lines(N,M[[i]][,j+1],col=j-1,lwd=2,lty=1) # /DM corrected
  }
  
  ##--Legend
  legend('topright',c(expression(bold(SE[1]/SE[DM])),    expression(bold(SE[2]/SE[DM])),
                      expression(bold(SE[1]/SE[DM]^"'")),expression(bold(SE[2]/SE[DM]^"'"))),lwd=2,col=c(2,5,2,5),lty=c(2,2,1,1),
         border = 'white',box.lty=0)
}

##-- Summary of the real sample sizes
summary(with(datos,final_cases_T1+final_cases_T2))
sum(with(datos,final_cases_T1<10 | final_cases_T2<10))/208

library(data.table)
library(ggplot2)
df0 <- as.data.table(rbind(M[[1]],M[[2]],M[[3]]))
colnames(df0) <- c('E3','E1','E1/E3','E1/E4','E2','E2/E3','E2/E4')
df0$SS <- rep(N,3)
df0$Ratio <- rep(paste0('Ratio ',c('1:1','1:2','1:3')),each=8)
df1 <- melt.data.table(df0,id.vars = c('Ratio','SS'),measure.vars = c(3,4,6,7) )
mm <- matrix(unlist(strsplit(as.character(df1$variable),'/',fixed = TRUE)),ncol=2,byrow=TRUE)
df1$Comparator <- mm[,1]
df1$Estimator <- mm[,2]

ggplot(df1,aes(x=SS,y=value,colour=Estimator)) + # linetype=Numerator,
  geom_hline(yintercept = 1) +
  geom_line(size=1.2) +
  scale_y_log10(limits=c(0.9,1.1),breaks=seq(0.85,1.15,0.05)) +
  xlim(c(0,769)) +
  facet_wrap(Comparator~Ratio,labeller = label_bquote(bold(Estimator: E[.(substr(Comparator,2,2))]~~~~~~~'Allocation'~.(Ratio)))) +
  xlab('Total sample size') +
  ylab('Ratio') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=13,face='bold'),
        strip.text = element_text(size=10,face='bold'),
        axis.title = element_text(size=13,face='bold')) +
  scale_color_discrete(labels = c(expression(bold(SE[DM])), expression(bold(SE[DM]^"C"))))

##-- Remove all objects
rm(list=ls())

#-------------------------------------------------------------------
#
#
#                            Simulation C                         
#
#
#-------------------------------------------------------------------
####################################################################
# Parameters of the simulation
####################################################################
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
data <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")

nsim <- 10000                                # Number of simulations
N1 <- data$final_cases_T1                    # Sample size in one group
N2 <- data$final_cases_T2                    # Sample size in other group
V1 <- data$final_sd_T1^2
V2 <- data$final_sd_T2^2                     # Actual variances

####################################################################
# Store information
####################################################################
M <- matrix(ncol=7,nrow=208) # Matrix to store infiormation


####################################################################
# Simulation
####################################################################
set.seed(12345)
for(j in 1:208){
  n1 <- N1[j]
  n2 <- N2[j]
  v1 <- V1[j]
  v2 <- V2[j]
  cat('Iteration:',j,'of 208\n')                                  # Print the process
  y01 <- rnorm(n1,0,sqrt(v1))                       # Generate normal data for one arm
  y02 <- rnorm(n2,0,sqrt(v2))                       # Generate normal data for the other arm
  
  var1 <- var2 <- var01 <- var02 <- est1 <- est2 <- c()
  
  for (i in 1:nsim){
    
    ##-- Method 1-----------------------------------
    y1 <- rnorm(n1,0,sqrt(v1))                    # Generate normal data for one arm
    y2 <- rnorm(n2,0,sqrt(v2))                    # Generate normal data for the other arm
    var1[i] <- var(y1)                            # Sample variance for one arm
    var2[i] <- var(y2)                            # Sample variance for the other arm
    est1[i] <- log(var1[i]/var2[i])               # First estimation
    
    ##-- Method 2-----------------------------------
    y1 <- sample(y01,n1,rep=TRUE)                 # Generate normal data for one arm
    y2 <- sample(y02,n2,rep=TRUE)                 # Generate normal data for the other arm
    var01[i] <- var(y1)                           # Sample variance for one arm
    var02[i] <- var(y2)                           # Sample variance for the other arm
    est2[i] <- log(var01[i]/var02[i])             # Second estimation
    
  }
  
  E1 <- var(est1,na.rm=TRUE)                      # First estimation: E1
  E2 <- var(est2,na.rm=TRUE)                      # Second estimation: E2
  E3 <- 2/(n1-1)+2/(n2-1)                         # Estimation from delta method: DM
  E4 <- 2/(n1-2)+2/(n2-2)                         # Estimation from corrected delta method: DM'
  
  ##-- Store information and update the counter
  M[j,] <- c(E3,E1,E1/E3,E1/E4,E2,E2/E3,E2/E4)
}
M <- cbind(M,1/(M[,4]/M[,2]))
#colnames(M) <- c('SE Real','SE_from_var','E1divE3','E1divE4','SE_from_norm','E2divE3','E2divE4','E4')
colnames(M) <- c('E3','E1','E1divE3','E1divE4','E2','E2divE3','E2divE4','E4')
rownames(M) <- paste0(rep(N1,each=3),'_',round(FRAC,2))

# Plot
windows(8,6)
par(mfrow=c(1,1),font.axis=2,font.lab=4,mar=c(5,5,5,1),las=1)
plot(M[,4]~I(N1+N2),cex=sqrt(V1/V2),log='xy',xlab='Total Sample size',ylab=expression(bold(SE[1]/SE[DM]^"'")))
rect(0.0001,0.1,10000,10,border=NA,col='grey90')
abline(h=seq(0.94,1.06,0.02),col='white')
abline(v=c(20,50,100,200,500,1000),col='white')
abline(h=1,lty=2)
points(M[,4]~I(N1+N2),cex=sqrt(V1/V2),lwd=2,col='darkblue')

####################################################################
# Graphic with ggplot
####################################################################
library(data.table)
library(ggplot2)
df0 <- as.data.table(M)
df0$SS <- N1+N2
df2 <- melt.data.table(df0,id.vars = c('SS'),measure.vars = c(3,4,6,7))
mm <- matrix(unlist(strsplit(as.character(df2$variable),'div',fixed = TRUE)),ncol=2,byrow=TRUE)
df2$Comparator <- mm[,1]
df2$Estimator <- mm[,2]

ggplot(df2[Comparator=='E1'],aes(x=SS,y=value,colour=Estimator)) + 
  geom_hline(yintercept = 1) +
  geom_smooth(size=1.2,se = FALSE) +
  scale_y_log10(limits=c(0.9,1.1),breaks=seq(0.9,1.1,0.05)) +
  xlim(c(0,900)) +
  xlab('Total sample size') +
  ylab('Ratio') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=11,face='bold')) +
  scale_color_discrete(labels = c(expression(bold(SE[DM])), expression(bold(SE[DM]^C))))


# for (i in 1:5) points(M[,4]~I(N1+N2),cex=i*sqrt(V1/V2),lwd=2,col=rgb(0,0,1,0.1),pch=19)
# abline(lm(log(M[,4])~I(log(N1+N2))),col='orange',lwd=4)
# for (i in 1:100) abline(runif(1,-0.05,0.05),runif(1,-0.05,0.05),col='orange',lwd=4)
# summary(M)
# pairs(M[,c(5,2,1,8)])  # Comparing variance of the studies with the two estimations
# plot(M[,8]~N1)
# plot(M[,3]~V1,log='x')
# plot(M[,3]~M[,1])
# abline(0,1)
# plot(M[,4]~M[,1])
# abline(0,1)
# boxplot(M[,c(3)])

# library(PairedData)
# x=M[,2]
# y=1/(M[,4]/M[,2])
# p <- paired(x,y)
# plot(p,type='BA')
# plot(M[,8]~M[,2],xlim=c(0,0.2),ylim=c(0,0.2))
# 
# par(mfrow=c(2,2))
# plot(lm(M[,8]~M[,2]))
# plot(lm(log(M[,8])~log(M[,2])))





