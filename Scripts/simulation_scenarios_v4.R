rm(list=ls())

####################################################################
# Load libraries
####################################################################
library(metafor)
library(ggplot2)
library(data.table)

####################################################################
# Define colors
####################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####################################################################
# Read data
####################################################################
# Read data in order to extract the sample size (n) of each study
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/data.csv'
d <- read.table(url(URL),header=TRUE,sep=";",stringsAsFactors = FALSE,quote = "")
summary(d)

####################################################################
# 
# Simulations
#
####################################################################
####################################################################
# Information from data
####################################################################
N <- nrow(d)                   # Number of rows
n_treat <- d$final_cases_T1    # Sample size treated
n_ctrl <- d$final_cases_T2     # Sample size controls
n <- n_treat + n_ctrl          # Total sample size



#############################################################################
##-- Scenario 1: All Additive effect and fixed  
#############################################################################
##-- Parameters
SIGMA <- c(0,0.1,0.5,1,2)        # Sigma: baseline variability
DELTA_C <- 0                     # Change in controls 
DELTA_T <- 0                     # Change in treated 
EFFECT <- DELTA_T - DELTA_C      # Effect
nsim <- 10                       # Number of simulations

##-- Martix to store the results
M <- matrix(nrow=nsim*length(SIGMA),ncol=5)
colnames(M) <- c('sigma','iteration','mu','tau','I2')

##-- Simulations
set.seed(12345)

for(sigma in SIGMA){
    for(k in 1:nsim){
      cat('Sigma:',sigma,'Iteration:',k,'\n')
      y <- se <- yBaselineRatio <- c()
      
      ##-- Generate data for each single study
      for(j in 1:N){
        Y_B <- rnorm(n[j],sigma)    # Baseline values
        Y_OC <- Y_B + DELTA_C       # Potential outcome controls
        Y_OT <- Y_B + DELTA_T       # Potential outcome treateds

        ##-- Random allocation
        sel.trt <- sample(1:n[j],n_treat[j],rep=FALSE)
        sel.ctrl <- (1:n[j])[-sel.trt]
        y_bc <- Y_B[sel.ctrl]        # patients allocated to control
        y_bt <- Y_B[sel.trt]         # patients allocated to treated
        y_oc <- Y_OC[sel.ctrl]       # patients allocated to control
        y_ot <- Y_OT[sel.trt]        # patients allocated to treated
        
        y[j] <- log(var(y_ot)/var(y_oc))
        yBaselineRatio[j] <- log(var(y_bt)/var(y_bc))
        se[j] <- sqrt(2/(n_treat[j]-2) + 2/(n_ctrl[j]-2))
      }
      
      data <- data.frame(y=y,se=se)                                     # data to fit the model
      mod <- rma(y,sei=se,data=data,mods=~yBaselineRatio,method='REML') # fit the model                                     
      
      ##-- Store results
      s <- which(SIGMA==sigma)
      pos <- (s-1)*nsim
      M[pos + k,1] <- sigma                    # sigma
      M[pos + k,2] <- k                        # number of simulation
      M[pos + k,3] <- as.numeric(mod$beta[1,1])# coefficient that estimates log(mu)
      M[pos + k,4] <- sqrt(mod$tau2)           # estimate of tau 
      M[pos + k,5] <- mod$I2                   # estimate of I^2 
    }
}

##-- Summary of the results
summary(M)
df <- as.data.frame(M)
df$I2 <- df$I2/100
df2 <- melt(df,id.vars=1, measure.vars=3:5,variable.name = "stat")

##-- Plot
ggplot(df2,aes(x=sigma,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='pointrange',size=0.5,alpha=0.5,
               fun.ymin=function(x) mean(x) - 2*sd(x)/sqrt(nsim),
               fun.ymax=function(x) mean(x) + 2*sd(x)/sqrt(nsim)) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  geom_hline(yintercept = c(-0.1213,0.55,0.8147),linetype=2, color=gg_color_hue(3)) +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~~~I^2)),values=gg_color_hue(3)) +
  labs(title = "All studies with additive treatment effect",
       subtitle = expression(Y[OT[i]]==Y[CT[i]]+1+epsilon[i]~~~~~~~~~epsilon[i]~"~"~N(0,sigma)))+
  theme(legend.title = element_blank(),
        legend.position = 'bottom') +
  xlab(expression(sigma)) + ylab('Value')

#############################################################################
##-- Scenario 2: All Additive effect and some random --> % Fixed and % Random
#############################################################################
##-- Parameters
THETA_MAX <- c(1,3,5,7)          # theta=sigma_C/sigma_T; theta_max = max(theta)
PROP_RANDOM <- seq(0,0.5,0.05)   # proportion of studies with random treatment effect
DELTA_C <- 0                     # Change in controls 
DELTA_T <- 0                     # Change in treated 
EFFECT <- DELTA_T - DELTA_C      # Effect
nsim <- 100                      # number of simulations

##-- Martix to store the results
M <- matrix(nrow=nsim*length(THETA_MAX)*length(PROP_RANDOM),ncol=9)
colnames(M) <- c('theta','prop_random','iteration',
                 'mu','tau','I2',        # Between arms
                 'mu_2','tau_2','I2_2')  # Over time

##-- Simulations
set.seed(12345)
rowi <- 1                    # Row indicator
for(pa in PROP_RANDOM){
  for(theta in THETA_MAX){
    for(k in 1:nsim){
      y <- se <- yBaselineRatio <- c()
      y2 <- se2 <- yControlRatio <- c()
      
      ##-- Generate data for each single study
      for(j in 1:N){
        effect.type   <- sample(c('fixed','random'),1,prob=c(1-pa,pa))   # generate type of treatment effect: random or fixed
        sigma.control <- ifelse(effect.type=='fixed',0,runif(1,0,theta)) # generate standard deviation of the effect in control group. If fixed, sigma=0, otherwise sigma=runif(0,sigma_max)
        sigma.treated <- as.numeric(effect.type=='random')               # generate standard deviation of the effect in treated group. If fixed, sigma=0, otherwise sigma=1
        
        # Does not work
        # sigma.control <- 7*as.numeric(effect.type=='random')             # generate standard deviation of the effect in control group. If fixed, sigma=0, otherwise sigma=runif(0,sigma_max)
        # sigma.treated <- ifelse(effect.type=='fixed',0,runif(1,0,theta)) # generate standard deviation of the effect in treated group. If fixed, sigma=0, otherwise sigma=1
          
        Y_B <- rnorm(n[j],0,1)                            # generate baseline values
        Y_OC <- Y_B + rnorm(n[j],DELTA_C,sigma.control)   # generate potential outcomes in control arm
        Y_OT <- Y_B + rnorm(n[j],DELTA_T,sigma.treated)   # generate potential outcomes in treatment arm
        
        sel.trt <- sample(1:n[j],n_treat[j],rep=FALSE)                    # indicator of treated patients 
        sel.ctrl <- (1:n[j])[-sel.trt]                                    # indicator of reference patients
        y_bc <- Y_B[sel.ctrl]                                                # observed outcome of treated patients
        y_bt <- Y_B[sel.trt]                                                 # observed outcome of treated patients
        y_oc <- Y_OC[sel.ctrl]                                               # observed outcome of treated patients
        y_ot <- Y_OT[sel.trt]                                                # observed outcome of treated patients
        
        # Between arms comparison
        y[j] <- log(var(y_ot)/var(y_oc))                                     # response of the model
        yBaselineRatio[j] <- log(var(y_bt)/var(y_bc))                        # response of the model at baseline
        se[j] <- sqrt(2/(n_treat[j]-2) + 2/(n_ctrl[j]-2))                    # within standard error of the response
        
        # Over time comparison
        y2[j] <- log(var(y_ot)/var(y_bt))                                    # response of the model
        yControlRatio[j] <- log(var(y_oc)/var(y_bc))                         # response of the model at baseline
        se2[j] <- sqrt(2/(n_treat[j]-2) + 
                         2/(n_treat[j]-2) - 
                         2*log(1 + 2*cor(y_ot,y_bt)^2/(n_treat[j]-1)))       # within standard error of the response
        
        
      }
      
      # Between arms
      data <- data.frame(y=y,se=se,yBaselineRatio=yBaselineRatio)                        # data to fit the model
      mod <- try(rma(y,sei=se,data=data,mods=~yBaselineRatio,method='REML'),silent=TRUE) # fit the model                                     
      
      # Over time
      data2 <- data.frame(y=y2,se=se2,yBaselineRatio=yControlRatio)                        # data to fit the model
      mod2 <- try(rma(y2,sei=se2,data=data2,mods=~yControlRatio,method='REML'),silent=TRUE) # fit the model                                     
      
      
      M[rowi,1] <- theta                     # THETA_MAX                            
      M[rowi,2] <- pa                        # proportion of random studies
      M[rowi,3] <- k                         # Iteration
      
      # Between arms
      if(class(mod)[1]!="try-error"){
        ##-- Store results
        M[rowi,4] <- as.numeric(mod$beta[1,1]) # coefficient that estimates log(mu)
        M[rowi,5] <- sqrt(mod$tau2)            # estimate of tau^2 
        M[rowi,6] <- mod$I2                    # estimate of I^2 
      }else{
        cat('Fisher scoring algorithm did not converge in between arms comparison.\n')
      }
      
      # Over time
      if(class(mod2)[1]!="try-error"){
        ##-- Store results
        M[rowi,7] <- as.numeric(mod2$beta[1,1]) # coefficient that estimates log(mu)
        M[rowi,8] <- sqrt(mod2$tau2)            # estimate of tau^2 
        M[rowi,9] <- mod2$I2                    # estimate of I^2 
      }else{
        cat('Fisher scoring algorithm did not converge in over time comparison.\n')
      }
      
      ##-- Update indicator
      rowi <- rowi+1
    }
    cat('Theta_max:',theta,'Prop. Random:',pa,'\n')
  }
}

##-- Summary of the results
summary(M)

##-- Store the data
df <- as.data.frame(M)
save.image(file='C:/Users/jordi/Google Drive/Tesis/Memoria/Scripts/simulation_REM.RData')
# load('C:/Users/jordi/Google Drive/Tesis/Memoria/Scripts/simulation_REM.RData')

##-- Arrange the data
df$I2 <- df$I2/100
df$I2_2 <- df$I2_2/100
df2 <- as.data.table(df)[,.(mu=mean(mu,na.rm=TRUE),sd.mu=sd(mu,na.rm=TRUE),
                            tau=mean(tau,na.rm=TRUE),sd.tau=sd(tau,na.rm=TRUE),
                            I2=mean(I2,na.rm=TRUE),sd.I2=sd(I2,na.rm=TRUE),
                            mu2=mean(mu_2,na.rm=TRUE),sd.mu2=sd(mu_2,na.rm=TRUE),
                            tau2=mean(tau_2,na.rm=TRUE),sd.tau2=sd(tau_2,na.rm=TRUE),
                            I22=mean(I2_2,na.rm=TRUE),sd.I22=sd(I2_2,na.rm=TRUE)),
                         by=.(theta,prop_random)]
df3 <- melt(df2,id.vars=1:2, measure.vars=c(3,5,7),variable.name = "stat")
df3[,theta:=factor(theta)]
df3[,theta_string:=paste0('theta[M]==',theta)]
df3$stat <- factor(df3$stat,levels=c('mu','tau','I2'))

df4 <- melt(df2,id.vars=1:2, measure.vars=c(9,11,13),variable.name = "stat")
df4[,theta:=factor(theta)]
df4[,theta_string:=paste0('theta[M]==',theta)]

##-- Between arms
ggplot(df3,aes(x=prop_random,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='point',size=1.5) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  facet_wrap(~theta_string,nrow = 2,labeller = label_parsed) +
  scale_x_continuous(limits=c(0,0.5)) +
  geom_hline(yintercept = c(-0.1213,0.5456,0.8085),linetype=2,color=rep(gg_color_hue(3),length(THETA_MAX))) + 
  xlab(expression('Proportion ('~pi[R]~') of studies with random effect')) +
  labs(title = expression('All studies with additive treatment effect:'~pi[R]~'of them with random effect')) + # subtitle = expression("Fixed effect:"~Y[OT[i]]==Y[OC[i]]~~~~~"Random effect:"~Y[OT[i]]==Y[OC[i]] + T~~~~~~~~T~"~"~N(0,tau^2)~~~~~~~~tau^2~"~ ["~0~","~theta[M]~"]")
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  ylab('Value') +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~I^2)),values=gg_color_hue(3))

##-- Over time --> Does not work
ggplot(df4,aes(x=prop_random,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='point',size=1.5) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  facet_wrap(~theta_string,nrow = 2,labeller = label_parsed) +
  scale_x_continuous(limits=c(0,0.5)) +
  geom_hline(yintercept = c(-0.1510,0.5897,0.9365),linetype=2,color=rep(gg_color_hue(3),length(THETA_MAX))) + 
  xlab(expression('Proportion ('~pi[R]~') of studies with random effect')) +
  labs(title = expression('All studies with additive treatment effect:'~pi[R]~'of them with random effect')) + # subtitle = expression("Fixed effect:"~Y[OT[i]]==Y[OC[i]]~~~~~"Random effect:"~Y[OT[i]]==Y[OC[i]] + T~~~~~~~~T~"~"~N(0,tau^2)~~~~~~~~tau^2~"~ ["~0~","~theta[M]~"]")
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  ylab('Value') +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~I^2)),values=gg_color_hue(3))


########################### Revisado hasta aquí #################3


#############################################################################
##-- Scenario 3: All Multiplicative effect and fixed
#############################################################################
SIGMA <- seq(0,0.9,.1)
M <- matrix(nrow=nsim*length(SIGMA),ncol=5)
colnames(M) <- c('sigma','ite','mu','tau','I2')

for(sigma in SIGMA){
  for(k in 1:nsim){
    cat('Sigma:',sigma,'Iteration:',k,'\n')
    y <- c()
    se <- c()
    for(j in 1:N){
      OC <- rnorm(n[j])
      T <- sqrt(exp(-0.1213)-sigma^2)
      #OT <- OC * sqrt(exp(-0.1213)-sigma^2) + rnorm(n[j],0,sigma)
      OT <- OC * T + rnorm(n[j],0,sigma)
      
      sel.trt <- sample(1:n[j],n_treat[j],rep=FALSE)
      sel.ctrl <- (1:n[j])[-sel.trt]
      oc <- OC[sel.ctrl]
      ot <- OT[sel.trt]
      
      y[j] <- log(var(ot)/var(oc))
      se[j] <- sqrt(2/(n_treat[j]-2) + 2/(n_ctrl[j]-2))
    }
    
    data <- data.frame(y=y,se=se)
    mod <- rma(y,sei=se,data=data)
    
    s <- which(SIGMA==sigma)
    pos <- (s-1)*nsim
    M[pos + k,1] <- sigma
    M[pos + k,2] <- k
    M[pos + k,3] <- as.numeric(mod$beta)
    M[pos + k,4] <- sqrt(mod$tau2)
    M[pos + k,5] <- mod$I2
  }
}


summary(M)
df <- as.data.frame(M)
df$I2 <- df$I2/100
df2 <- melt(df,id.vars=1, measure.vars=3:5,variable.name = "stat")
ggplot(df2,aes(x=sigma,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='point',size=2,
               fun.ymin=function(x) mean(x) - 2*sd(x)/sqrt(nsim),
               fun.ymax=function(x) mean(x) + 2*sd(x)/sqrt(nsim)) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~~~I^2)),values=gg_color_hue(3)) +
  geom_hline(yintercept = c(-0.1213,0.55,0.8147),linetype=2, color=gg_color_hue(3)) +
  # ggtitle('Multiplicative: Yt = Yc * sqrt(exp(mu)-sigma^2) + N(0,sigma)') +
  labs(title = "All studies with multiplicative treatment effect",
       subtitle = expression(Y[OT[i]]==T%.%Y[CT[i]] + epsilon[i]~~~~~~~~T^2==e^hat(mu)-sigma^2~~~~~~~~~epsilon[i]~"~"~N(0,sigma))) +
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  xlab(expression(sigma)) + ylab('Value')
# "Yt = sqrt(exp(mu)-sigma^2) * Yc + N(0,sigma)")  --> lo que ponia antes


#############################################################################
##-- Scenario 4: All multiplicative effect --> % Fixed and % Random
#############################################################################
PER <- seq(0,1,.2)
THETA_MAX <- c(0.1,0.5,1,2)
M <- matrix(nrow=nsim*length(PER)*length(THETA_MAX),ncol=6)
colnames(M) <- c('per','sigma','ite','mu','tau','I2')
i <- 1
for(THETA_MAX in THETA_MAX){
  for(per in PER){
    for(k in 1:nsim){
      cat('per:',per,'THETA_MAX',THETA_MAX,'Iteration:',k,'\n')
      y <- c()
      se <- c()
      for(j in 1:N){
        OC <- rnorm(n[j])
        effect.type <- sample(c('fixed','random'),1,prob=c(1-per,per))
        sigma <- ifelse(effect.type=='fixed',0,runif(1,0,THETA_MAX))
        effect <- ifelse(effect.type=='fixed',sqrt(0.8),1)
        # OT <- effect*OC + rnorm(n[j],0,sigma)
        OT <- 0.8*OC + rnorm(n[j],0,sigma)
        
        sel.trt <- sample(1:n[j],n_treat[j],rep=FALSE)
        sel.ctrl <- (1:n[j])[-sel.trt]
        oc <- OC[sel.ctrl]
        ot <- OT[sel.trt]
        
        y[j] <- log(var(ot)/var(oc))
        se[j] <- sqrt(2/(n_treat[j]-2) + 2/(n_ctrl[j]-2))
      }
      
      data <- data.frame(y=y,se=se)
      mod <- rma(y,sei=se,data=data)
      
      M[i,1] <- per
      M[i,2] <- THETA_MAX
      M[i,3] <- k
      M[i,4] <- as.numeric(mod$beta)
      M[i,5] <- sqrt(mod$tau2)
      M[i,6] <- mod$I2
      i <- i+1
    }
  }
}

summary(M)
df <- as.data.frame(M)
df$I2 <- df$I2/100
df2 <- melt(df,id.vars=1:2, measure.vars=4:6,variable.name = "stat")
df2$sigma <- factor(df2$sigma)
theta_names <- list(
  '0.1'=expression(theta[M]==0.1),
  '0.5'=expression(theta[M]==0.5),
  '1'=expression(theta[M]==1),
  '2'=expression(theta[M]==2)
)
theta_labeller <- function(variable,value) return(theta_names[value])


ggplot(df2,aes(x=per,y=value,color=stat)) + 
  geom_hline(yintercept = c(-0.1213,0.55,0.8147),linetype=2,color=rep(gg_color_hue(3),4)) +
  # ggtitle('Multiplicative effect: per% Random and (1-per)% Fixed') + 
  xlab('Percentatge (per) of studies with random effect') +
  facet_wrap(~sigma,nrow = 2,labeller = theta_labeller) +
  stat_summary(fun.y=mean,geom='point',size=2) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~~~I^2)),values=gg_color_hue(3)) +
  labs(title = "All studies with multiplicative treatment effect",
       subtitle = expression('Fixed effect: '~Y[OT[i]]==0.8%.%Y[OC[i]]~~~~~'Random effect:'~Y[OT[i]]==0.8%.%Y[OC[i]]+epsilon[i]~~~~~~~~epsilon[i]~'~'~N(0,sigma)~~~~~~~~sigma~'~'~U~"["~0~','~theta[M]~"]")) +
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  xlab(expression('Proportion ('~pi[R]~') of studies with random effect')) + ylab('Value')


