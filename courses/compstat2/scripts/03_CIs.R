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

# using B= 100
set.seed(1)
boot.score <- boot(data=examScor,statistic=theta,R=100)      
boot.score

plot(boot.score)

# using B= 100000
set.seed(1)
boot.score2 <- boot(data=examScor,statistic=theta,R=10000)      
boot.score2

plot(boot.score2)


# Bootstrap confidence intervals:
# A unique function for everything.
help(boot.ci)
intervals1 = boot.ci(boot.score, conf = 0.95)
intervals1

intervals2 = boot.ci(boot.score2, conf = 0.95)
intervals2

# Asymptotic Normal  ------------------------------------------------------
intervals1$normal
#pdf('Boot_normal1.pdf',height=5,width=8)
plot(boot.score)
#dev.off()

#pdf('Boot_normal2.pdf',height=5,width=8)
plot(boot.score2)
#dev.off()


# Bootstrap-t interval ----------------------------------------------------
# in this case we need to find a formula - or compute with Bootstrap - the variance.
# Authomatically, the function boot.ci only computes the approximate interval
# that uses as Pivot the quantity \hat \theta - theta (unstudentized).
intervals1$basic
intervals2$basic

# in order to compute the Bootstrap-t interval we need to modify
# the function theta. The function needs to return two values: the estimate
# and its variance.
thetavariance = function(data,indexes){
  datamod = data[indexes,]
  Shat = cov(datamod)
  lambdas = eigen(Shat)$values
  
  sdboot = boot(data=datamod,statistic=theta,R=100)
  
  return(c(max(lambdas)/sum(lambdas), sd(sdboot$t)^2 ))
}
boot.varestimate = boot(data=examScor,statistic=thetavariance,R=100)
boot.varestimate

intervals1_student = boot.ci(boot.varestimate)
pdf('Student1.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals1_student$percent[4],col='darkblue',lwd=3)
abline(v=intervals1_student$percent[5],col='darkblue',lwd=3)
plot(ecdf(boot.score$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=intervals1_student$percent[4],col='darkblue',lwd=3)
abline(v=intervals1_student$percent[5],col='darkblue',lwd=3)
dev.off()

# very long!
boot.varestimate2 = boot(data=examScor,statistic=thetavariance,R=10000)
boot.varestimate2
intervals2_student = boot.ci(boot.varestimate2)

pdf('Student2.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score2$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals2_student$percent[4],col='darkblue',lwd=3)
abline(v=intervals2_student$percent[5],col='darkblue',lwd=3)
plot(ecdf(boot.score2$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=intervals2_student$percent[4],col='darkblue',lwd=3)
abline(v=intervals2_student$percent[5],col='darkblue',lwd=3)
dev.off()


# percentile method -------------------------------------------------------
intervals1$percent

# we can also use the R function quantile. We obtain slightly different values since the empirical quantiles are approximated
quantile(boot.score$t,0.025)
quantile(boot.score$t,0.975)

#pdf('Percentile1.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals1$percent[4],col='darkblue',lwd=3)
abline(v=intervals1$percent[5],col='darkblue',lwd=3)
plot(ecdf(boot.score$t),main='ECDF of Bootstrap replications')
abline(v=intervals1$percent[4],col='darkblue',lwd=3)
abline(v=intervals1$percent[5],col='darkblue',lwd=3)
#dev.off()

# B = 10000
intervals2$percent

# now we obtain closer values with the quantile function:
quantile(boot.score2$t,0.025)
quantile(boot.score2$t,0.975)

#pdf('Percentile2.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score2$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals2$percent[4],col='darkblue',lwd=3)
abline(v=intervals2$percent[5],col='darkblue',lwd=3)
plot(ecdf(boot.score2$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=intervals2$percent[4],col='darkblue',lwd=3)
abline(v=intervals2$percent[5],col='darkblue',lwd=3)
#dev.off()



# BCa intervals -----------------------------------------------------------
intervals1$bca

#pdf('BCa1.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals1$bca[4],col='darkblue',lwd=3)
abline(v=intervals1$bca[5],col='darkblue',lwd=3)
plot(ecdf(boot.score$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=intervals1$bca[4],col='darkblue',lwd=3)
abline(v=intervals1$bca[5],col='darkblue',lwd=3)
#dev.off()

intervals2$bca

#pdf('BCa2.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score2$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=intervals2$bca[4],col='darkblue',lwd=3)
abline(v=intervals2$bca[5],col='darkblue',lwd=3)
plot(ecdf(boot.score2$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=intervals2$bca[4],col='darkblue',lwd=3)
abline(v=intervals2$bca[5],col='darkblue',lwd=3)
#dev.off()

# ABC ---------------------------------------------------------------------
# we need to modify the theta function.
# We need to introduce theta as a function of Pj
# covariance matrix computed with weights:
# Cov(X) = E(XX') - mux mux'
thetamod = function(data,w){ # the function needs to take as input the weights
  n = dim(data)[1]
  p = dim(data)[2]
  datamatrix = as.matrix(data)
  # mean vector of X:
  mux = t(datamatrix) %*% (w)
  
  computecov = function(xvector){
    xvector = matrix(data=xvector,ncol=1)
    return(xvector %*% t(xvector))
  }
  XtX = apply(datamatrix,1,computecov)
  meanXtX = matrix(XtX %*% w,nrow=p,ncol=p)
  
  Shat = meanXtX - mux %*% t(mux)
  
  lambdas = eigen(Shat)$values
  return(max(lambdas)/sum(lambdas))
}
abc1 = abc.ci(examScor,thetamod)
# does not depend on B


pdf('ABC.pdf',width=8,height=5)
layout(cbind(1,2))
hist(boot.score2$t,probability=TRUE,col=rgb(231,243,243,maxColorValue = 255))
abline(v=abc1[2],col='darkblue',lwd=3)
abline(v=abc1[3],col='darkblue',lwd=3)
plot(ecdf(boot.score2$t),main='ECDF of Bootstrap replications',lwd=2)
abline(v=abc1[2],col='darkblue',lwd=3)
abline(v=abc1[3],col='darkblue',lwd=3)
dev.off()


# Comparison --------------------------------------------------------------
# Let us consider a simpler problem: finding a confidence interval for the mean of a population
# We will simulate a sample of size 10 from:
# a) a Normal distribution of mean 0 and sd 1
# b) a translated Exponential distribution of mean 0 and sd 1 (lambda = 1)
# we repeat the computation of confidence intervals 1000 times and compute the coverage probability of a 95% CI


M = 1000 # Simulation to estimate coverage probability
B = 1000 # MC simulation for Bootstrap
n=10 # sample size 

# scenario a) Normal distribution

# matrices where we store all confidence intervals:
ICs_t = matrix(nrow=M,ncol=2) # classical CI based on t distribution
ICs_AsN = matrix(nrow=M,ncol=2) # asymptotic Normal
ICs_Boott = matrix(nrow=M,ncol=2) # Bootstrap t
ICs_Perc = matrix(nrow=M,ncol=2) # Percentile
ICs_BCa = matrix(nrow=M,ncol=2) # BCa
ICs_ABC = matrix(nrow=M,ncol=2) # ABC

# now for returning the variance (needed for the Studentized interval) we can use the closed formula
computemean = function(data,indexes){
  xbar = mean(data[indexes])
  variance = (n-1)*var(data[indexes])/(length(data)^2)
  return(c(xbar,variance))
}

computemeanw = function(data,w){
  xbar = sum(data * w)
  return(xbar)
}

set.seed(15072017)
for(sim in 1:M){
  data = rnorm(n)
  
  # classical:
  ICs_t[sim,] = t.test(data)$conf.int[1:2]
  
  # running bootstrap 
  boot.result = boot(data,computemean,R=B)
  intervals = boot.ci(boot.result)
  
  # Asymptotic Normal
  ICs_AsN[sim,] = intervals$normal[2:3] # the first argument is the conf level
  
  # Bootstrap-t
  ICs_Boott[sim,] = intervals$student[4:5] # the forst three arguments are conf level and positions
  
  # Percentile
  ICs_Perc[sim,] = intervals$percent[4:5]
  
  # BCa
  ICs_BCa[sim,] = intervals$bca[4:5]
  
  # ABC
  ICs_ABC[sim,] = abc.ci(data,computemeanw)[2:3] 
}

# how many times 0 is included?
coverage.t = mean( ICs_t[,1]<0 & ICs_t[,2]>0 )
coverage.AsN = mean( ICs_AsN[,1]<0 & ICs_AsN[,2]>0 )
coverage.Boott = mean( ICs_Boott[,1]<0 & ICs_Boott[,2]>0 )
coverage.Perc = mean( ICs_Perc[,1]<0 & ICs_Perc[,2]>0 )
coverage.BCa = mean( ICs_BCa[,1]<0 & ICs_BCa[,2]>0 )
coverage.ABC = mean( ICs_ABC[,1]<0 & ICs_ABC[,2]>0 )

coverage.t
coverage.AsN 
coverage.Boott 
coverage.Perc
coverage.BCa 
coverage.ABC 

coverage.ABC - 1.96*sqrt(coverage.ABC*(1-coverage.ABC)/M)
coverage.ABC + 1.96*sqrt(coverage.ABC*(1-coverage.ABC)/M)

# scenario b) Exponential distribution

set.seed(15072017)
for(sim in 1:M){
  data = rexp(n) - 1
  
  # classical:
  ICs_t[sim,] = t.test(data)$conf.int[1:2]
  
  # running bootstrap 
  boot.result = boot(data,computemean,R=B)
  intervals = boot.ci(boot.result)
  
  # Asymptotic Normal
  ICs_AsN[sim,] = intervals$normal[2:3] # the first argument is the conf level
  
  # Bootstrap-t
  ICs_Boott[sim,] = intervals$student[4:5] # the forst three arguments are conf level and positions
  
  # Percentile
  ICs_Perc[sim,] = intervals$percent[4:5]
  
  # BCa
  ICs_BCa[sim,] = intervals$bca[4:5]
  
  # ABC
  ICs_ABC[sim,] = abc.ci(data,computemeanw)[2:3] 
}

# how many times 0 is included?
coverage.t = mean( ICs_t[,1]<0 & ICs_t[,2]>0 )
coverage.AsN = mean( ICs_AsN[,1]<0 & ICs_AsN[,2]>0 )
coverage.Boott = mean( ICs_Boott[,1]<0 & ICs_Boott[,2]>0 )
coverage.Perc = mean( ICs_Perc[,1]<0 & ICs_Perc[,2]>0 )
coverage.BCa = mean( ICs_BCa[,1]<0 & ICs_BCa[,2]>0 )
coverage.ABC = mean( ICs_ABC[,1]<0 & ICs_ABC[,2]>0 )

coverage.t
coverage.AsN 
coverage.Boott 
coverage.Perc
coverage.BCa 
coverage.ABC 


# scenario c) Exponential distribution, big n
n = 100

set.seed(15072017)
for(sim in 1:M){
  data = rexp(n) - 1
  
  # classical:
  ICs_t[sim,] = t.test(data)$conf.int[1:2]
  
  # running bootstrap 
  boot.result = boot(data,computemean,R=B)
  intervals = boot.ci(boot.result)
  
  # Asymptotic Normal
  ICs_AsN[sim,] = intervals$normal[2:3] # the first argument is the conf level
  
  # Bootstrap-t
  ICs_Boott[sim,] = intervals$student[4:5] # the forst three arguments are conf level and positions
  
  # Percentile
  ICs_Perc[sim,] = intervals$percent[4:5]
  
  # BCa
  ICs_BCa[sim,] = intervals$bca[4:5]
  
  # ABC
  ICs_ABC[sim,] = abc.ci(data,computemeanw)[2:3] 
}

# how many times 0 is included?
coverage.t = mean( ICs_t[,1]<0 & ICs_t[,2]>0 )
coverage.AsN = mean( ICs_AsN[,1]<0 & ICs_AsN[,2]>0 )
coverage.Boott = mean( ICs_Boott[,1]<0 & ICs_Boott[,2]>0 )
coverage.Perc = mean( ICs_Perc[,1]<0 & ICs_Perc[,2]>0 )
coverage.BCa = mean( ICs_BCa[,1]<0 & ICs_BCa[,2]>0 )
coverage.ABC = mean( ICs_ABC[,1]<0 & ICs_ABC[,2]>0 )

coverage.t
coverage.AsN 
coverage.Boott 
coverage.Perc
coverage.BCa 
coverage.ABC 


