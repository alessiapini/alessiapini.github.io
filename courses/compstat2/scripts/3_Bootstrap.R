library(boot)
help(boot)

# Bootstrap on test score data -------------------------------------------------
library(MVT)
data(examScor)
View(examScor)

theta <- function(data,indexes){
  Shat = cov(data[indexes,])
  lambdas = eigen(Shat)$values
  return(max(lambdas)/sum(lambdas))
}

set.seed(1)
boot.score <- boot(data=examScor,statistic=theta,R=500)      
boot.score
plot(boot.score)

pdf('Boot_score.pdf',height=5)
hist(boot.score$t,col=rgb(231,243,243,maxColorValue = 255),probability=TRUE,main='Histogram of Bootstrap replications')
dev.off()

# Exploring how the estimate change with B --------------------------------
B = c(10,100,500,1000,2000,5000)
boot.score = list()
boot.sd = numeric(length(B))
for(i in 1:length(B)){
  #set.seed(1)
  boot.score[[i]] <- boot(data=examScor,statistic=theta,R=B[i])      
  boot.sd[i] = sd(boot.score[[i]]$t)
  #plot(boot.score[[i]])
}
boot.sd


# Choosing B - very small simulation ---------------------------------------------------
MC = 10
B.explore = c(10, 100, 500, 1000, 2000, 5000)
se.estimates = matrix(nrow=MC,ncol=length(B.explore))
set.seed(1)
seeds = round(runif(MC,0,1000000))

for(sim in 1:MC){
  set.seed(seeds[sim])
  for(i in 1:length(B.explore)){
    boot.score <- boot(data=examScor,statistic=theta,R=B.explore[i])   
    se.estimates[sim,i] = sd(boot.score$t)
  }
}

matplot(B.explore,t(se.estimates),type='l',lty=1)
se.estimates[,6]


# Parametric bootstrap ----------------------------------------------------
# We use the parametric Bootstrap for the test score data
# We assume that the distribution of the scores is multivariate Normal
set.seed(1)
B = 1000
boot.replications = numeric(B)
muhat = colMeans(examScor)
Shat = cov(examScor)
n = dim(examScor)[1]

library(mvtnorm)
   
for(b in 1:B){
  xstar = rmvnorm(n,muhat,Shat)
  boot.replications[b] = theta(xstar,1:n)
}
sd(boot.replications)


pdf('ParBoot_Testscore.pdf',height=5)
hist(boot.replications,main='Parametric Bootstrap Replications',col=rgb(214,234,244,maxColorValue = 255))
dev.off()

sd(boot.replications)
mean(boot.replications) - theta(examScor,1:n)

# Parametric Bootstrap Failure ----------------------------------------------------
# We simulate a sample of size 10 from an Exponential distribution with rate lambda=1
# We use Parametric and nonparametric Bootstrap to estimate the first quartile of the distribution
# In the case of parametric Bootstrap, we simulate data from a Normal distribution.
n = 10
set.seed(3)
lambda = 0.5
data = rexp(n,lambda)

thetahat = quantile(data,probs=0.25)
thetahat
B = 1000

# Parametric Bootstrap
mu = mean(data)
sigma = sd(data)
param.bootstrap = numeric(B)

shapiro.test(data)

for(i in 1:B){
  xstar = rnorm(n,mu,sigma)
  param.bootstrap[i] = quantile(xstar,probs=0.25)
}


# Nonparametric Bootstrap
Q1 = function(data,indexes){
  return(quantile(data[indexes],probs=0.25))
}
nonparam.bootstrap <- boot(data=data,statistic=Q1,R=B)   

# MC Simulation
MC.sim = numeric(B)
for(i in 1:B){
  xstar = rexp(n,lambda)
  MC.sim[i] = quantile(xstar,probs=0.25)
}

pdf('Paramfailure.pdf',width=8,height=5)
layout(cbind(1,2,3))
hist(param.bootstrap,probability=TRUE,main='Parametric Bootstrap',xlim=c(-2,3),ylim=c(0,1.7),col=rgb(214,234,244,maxColorValue = 255))
hist(nonparam.bootstrap$t,probability=TRUE,main='Nonparametric Bootstrap',xlim=c(-2,3),ylim=c(0,1.7),col=rgb(162,217,214,maxColorValue = 255))    
hist(MC.sim,probability=TRUE,main='Monte Carlo',xlim=c(-2,3),ylim=c(0,1.7),col=rgb(221,166,13,maxColorValue = 255))
dev.off()

layout(cbind(1))
plot(ecdf(data))
xaxis = seq(0,6,len=100)
lines(xaxis,pexp(xaxis,lambda),type='l',col=2,lwd=3)

sd(param.bootstrap)
sd(nonparam.bootstrap$t)
sd(MC.sim)

# Nonparametric Bootstrap Failure ----------------------------------------------------
# data generated from uniform distribution on (0,theta) with theta = 1
# we want to estimate the upper bound of the domain theta

theta = 1
n = 500
set.seed(1)
x = runif(n,0,theta)
thetahat = max(x)
thetahat

B = 1000
# nonparametric:
thetahat.fn = function(data,indexes){
  return(max(data[indexes]))
}

# parametric:
set.seed(1)
thetahat = thetahat.fn(x,1:n)
boot.parametric = numeric(B)

for(b in 1:B){
  xstar = runif(n,0,thetahat)
  boot.parametric[b] = thetahat.fn(xstar,1:n)
}
layout(cbind(1,2,3))
hist(boot.parametric,xlim=c(0.97,1))




set.seed(1)
boot.nonparametric <- boot(data=x,statistic=thetahat.fn,R=B)      
hist(boot.nonparametric$t,xlim=c(0.97,1))


# MC simulation
MC.simulation = numeric(B)
for(b in 1:B){
  xstar = runif(n,0,1)
  MC.simulation[b] = max(xstar)
}
hist(MC.simulation,xlim=c(0.97,1))


pdf('NonParamfailure.pdf',width=8,height=5)
layout(cbind(1,2,3))
hist(boot.parametric,probability=TRUE,main='Parametric Bootstrap',xlim=c(0.8,1),ylim=c(0,40),col=rgb(214,234,244,maxColorValue = 255))
hist(boot.nonparametric$t,probability=TRUE,main='Nonparametric Bootstrap',xlim=c(0.8,1),ylim=c(0,40),col=rgb(162,217,214,maxColorValue = 255))    
hist(MC.simulation,probability=TRUE,main='Monte Carlo',xlim=c(0.8,1),ylim=c(0,40),col=rgb(221,166,13,maxColorValue = 255))
dev.off()

sd(boot.parametric)
sd(boot.nonparametric$t)
sd(MC.simulation)

# Bootstrap is consistent for the median ------------------------------------
# data simulated from a Binomial distribution with p=0.5 n = 10
p.binom = 0.5
n.binom = 10
median.true = n.binom*p.binom
n = 20

set.seed(2)
data = rbinom(n,size = n.binom, prob = p.binom)
sample.median = function(data,i){
  return(median(data[i]))
}

boot.result = boot(data,sample.median,R=1000)
boot.result

hist(boot.result$t,xlim=c(3,7),probability=TRUE,ylim=c(0,4),col=rgb(162,217,214,maxColorValue = 255))


#MC
MC.medians = numeric(B)
for(sim in 1:B){
  data = rbinom(n,size = n.binom, prob = p.binom)
  MC.medians[sim] = median(data)
}

hist(MC.medians,xlim=c(3,7),probability=TRUE,ylim=c(0,4),col=rgb(221,166,13,maxColorValue = 255))

sd(boot.result$t)
sd(MC.medians)

# bigger sample size
n = 100
set.seed(2)
data = rbinom(n,size = n.binom, prob = p.binom)

#Boot
boot.result = boot(data,sample.median,R=1000)
sd(boot.result$t)
#MC
MC.medians = numeric(B)
for(sim in 1:B){
  data = rbinom(n,size = n.binom, prob = p.binom)
  MC.medians[sim] = median(data)
}
sd(boot.result$t)
sd(MC.medians)


pdf('Median.pdf',height=5)
layout(cbind(1,2))
hist(boot.result$t,xlim=c(3,7),probability=TRUE,ylim=c(0,4),col=rgb(162,217,214,maxColorValue = 255),main='Bootstrap')
hist(MC.medians,xlim=c(3,7),probability=TRUE,ylim=c(0,4),col=rgb(221,166,13,maxColorValue = 255),main='Monte Carlo')
dev.off()

# Regression --------------------------------------------------------------

# generating data from a linear regression with normal errors
n = 50
set.seed(15072017)
x1 = runif(n,0,1)
x2 = runif(n,-2,3)

b0 = 1
b1 = 4
b2 = -2

y = b0 + x1*b1 + x2*b2 + rnorm(n,0,3)
classical.lm = lm(y~x1+x2)
summary(classical.lm)
betahat = classical.lm$coefficients

B = 1000

# Bootstrapping residuals:
# we first create a data frame containing data, residuals, and fitted values
x.frame <- data.frame(y=y,x1=x1,x2=x2,res=classical.lm$residuals,fit=classical.lm$fitted)

lm.statistic <- function(data,i){
  d <- data
  d$y <- d$fit+d$res[i] # applies bootstrap to the residuals
  lm.boot <- lm(y~x1+x2,data=d)
  return(c(coef(lm.boot)))
}
boot.regression = boot(x.frame,statistic=lm.statistic,R=B)
boot.regression

# Bootstrapping pairs:
lm.statistic2 <- function(data,i){
  d<-data[i,]
  lm.boot <- lm(y~x1+x2,data=d)
  return(c(coef(lm.boot)))
}
x.frame2 = data.frame(y=y,x1=x1,x2=x2)
boot.regression2 = boot(x.frame2,statistic=lm.statistic2,R=B)
boot.regression2

pdf('Regression0.pdf',height=6,width=8)
layout(cbind(1,2))
hist(boot.regression$t[,1],probability=1,main='Bootstrapping residuals - beta0',xlim=c(-3,4),ylim=c(0,0.5),col=rgb(221,166,13,maxColorValue = 255))
grid.x = seq(-3,7,len=100)
lines(grid.x,dnorm(grid.x,summary(classical.lm)$coeff[1,1],summary(classical.lm)$coeff[1,2]),col='darkblue',lwd=3)
hist(boot.regression2$t[,1],probability=1,main='Bootstrapping pairs - beta0',xlim=c(-3,4),ylim=c(0,0.5),col=rgb(162,217,214,maxColorValue = 255))
lines(grid.x,dnorm(grid.x,summary(classical.lm)$coeff[1,1],summary(classical.lm)$coeff[1,2]),col='darkblue',lwd=3)
dev.off()


pdf('Regression1.pdf',height=6,width=8)
layout(cbind(1,2))
hist(boot.regression$t[,2],probability=1,main='Bootstrapping residuals - beta1',xlim=c(-3,9),ylim=c(0,0.3),col=rgb(221,166,13,maxColorValue = 255))
grid.x = seq(-3,9,len=100)
lines(grid.x,dnorm(grid.x,summary(classical.lm)$coeff[2,1],summary(classical.lm)$coeff[2,2]),col='darkblue',lwd=3)
hist(boot.regression2$t[,2],probability=1,main='Bootstrapping pairs - beta1',xlim=c(-3,9),ylim=c(0,0.3),col=rgb(162,217,214,maxColorValue = 255))
lines(grid.x,dnorm(grid.x,summary(classical.lm)$coeff[2,1],summary(classical.lm)$coeff[2,2]),col='darkblue',lwd=3)
dev.off()
