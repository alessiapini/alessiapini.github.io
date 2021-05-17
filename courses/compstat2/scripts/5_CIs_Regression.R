# Test on regression parameters --------------------------------------------------------------

# simulation

# Normal data
M = 1000
n = 50
set.seed(15072017)
x1 = runif(n,0,1)
x2 = runif(n,-2,3)
b0 = 1
b1 = 4
b2 = -2
sigma = 3
B = 1000

ICs.parametric = matrix(nrow=M,ncol=2)
ICs.boot1 = matrix(nrow=M,ncol=2)
ICs.boot2 = matrix(nrow=M,ncol=2)

lm.statistic2 <- function(data,i){
  d<-data[i,]
  lm.boot <- lm(y~x1+x2,data=d)
  return(c(coef(lm.boot),summary(lm.boot)$coefficients[,2]^2))
}

lm.statistic <- function(data,i){
  d <- data
  d$y <- d$fit+d$res[i] # applies bootstrap to the residuals
  lm.boot <- lm(y~x1+x2,data=d)
  return(c(coef(lm.boot),summary(lm.boot)$coefficients[,2]^2))
}

for(sim in 1:M){
  # generating data from a linear regression with normal errors
  y = b0 + x1*b1 + x2*b2 + rnorm(n,0,sigma)
  classical.lm = lm(y~x1+x2)
  
  ICs.parametric[sim,1] = summary(classical.lm)$coefficients[2,1]-summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  ICs.parametric[sim,2] = summary(classical.lm)$coefficients[2,1]+summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  
  
  # Bootstrapping residuals:
  x.frame <- data.frame(y=y,x1=x1,x2=x2,res=classical.lm$residuals,fit=classical.lm$fitted)
  boot.regression = boot(x.frame,statistic=lm.statistic,R=B)
  ICs.boot1[sim,] = boot.ci(boot.regression,index=c(2,5))$student[4:5]
  
  # Bootstrapping pairs:
  x.frame2 = data.frame(y=y,x1=x1,x2=x2)
  boot.regression2 = boot(x.frame2,statistic=lm.statistic2,R=B)
  ICs.boot2[sim,] = boot.ci(boot.regression2,index=c(2,5))$student[4:5]
}

load('simNorm.RData')
mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )
mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )

# all those values are estimated probabilities over M replications of the experiment
# they are also affected by an error:
p.param = mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
p.boot1 = mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )
p.boot2 = mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )

c(p.param - qnorm(0.975)*sqrt(p.param*(1-p.param)/M),p.param + qnorm(0.975)*sqrt(p.param*(1-p.param)/M))
c(p.boot1 - qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M),p.param + qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M))
c(p.boot2 - qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M),p.boot2 + qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M))


#save.image('simNorm.RData')


# Exponential data
for(sim in 1:M){
  # generating data from a linear regression with normal errors
  y = b0 + x1*b1 + x2*b2 + rexp(n,1/sigma) - 1/sigma
  classical.lm = lm(y~x1+x2)
  
  ICs.parametric[sim,1] = summary(classical.lm)$coefficients[2,1]-summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  ICs.parametric[sim,2] = summary(classical.lm)$coefficients[2,1]+summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  
  
  # Bootstrapping residuals:
  x.frame <- data.frame(y=y,x1=x1,x2=x2,res=classical.lm$residuals,fit=classical.lm$fitted)
  boot.regression = boot(x.frame,statistic=lm.statistic,R=B)
  ICs.boot1[sim,] = boot.ci(boot.regression,index=c(2,5))$student[4:5]
  
  # Bootstrapping pairs:
  x.frame2 = data.frame(y=y,x1=x1,x2=x2)
  boot.regression2 = boot(x.frame2,statistic=lm.statistic2,R=B)
  ICs.boot2[sim,] = boot.ci(boot.regression2,index=c(2,5))$student[4:5]
}

load('simExp.RData')
mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )
mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )
#save.image('simExp.RData')

# all those values are estimated probabilities over M replications of the experiment
# they are also affected by an error:
p.param = mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
p.boot1 = mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )
p.boot2 = mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )

c(p.param - qnorm(0.975)*sqrt(p.param*(1-p.param)/M),p.param + qnorm(0.975)*sqrt(p.param*(1-p.param)/M))
c(p.boot1 - qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M),p.param + qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M))
c(p.boot2 - qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M),p.boot2 + qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M))


# Exponential data - non linear model
for(sim in 1:M){
  # generating data from a linear regression with normal errors
  y = exp(b0 + x1*b1 + x2*b2 + rexp(n,1/sigma) - 1/sigma)
  classical.lm = lm(y~x1+x2)
  
  ICs.parametric[sim,1] = summary(classical.lm)$coefficients[2,1]-summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  ICs.parametric[sim,2] = summary(classical.lm)$coefficients[2,1]+summary(classical.lm)$coefficients[2,2]*qt(0.975,n-3)
  
  
  # Bootstrapping residuals:
  x.frame <- data.frame(y=y,x1=x1,x2=x2,res=classical.lm$residuals,fit=classical.lm$fitted)
  boot.regression = boot(x.frame,statistic=lm.statistic,R=B)
  ICs.boot1[sim,] = boot.ci(boot.regression,index=c(2,5))$student[4:5]
  
  # Bootstrapping pairs:
  x.frame2 = data.frame(y=y,x1=x1,x2=x2)
  boot.regression2 = boot(x.frame2,statistic=lm.statistic2,R=B)
  ICs.boot2[sim,] = boot.ci(boot.regression2,index=c(2,5))$student[4:5]
}

load('simExp2.RData')
mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )
mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )

#save.image('simExp2.RData')

# all those values are estimated probabilities over M replications of the experiment
# they are also affected by an error:
p.param = mean( ICs.parametric[,1]<b1 & ICs.parametric[,2]>b1 )
p.boot1 = mean( ICs.boot1[,1]<b1 & ICs.boot1[,2]>b1 )
p.boot2 = mean( ICs.boot2[,1]<b1 & ICs.boot2[,2]>b1 )

c(p.param - qnorm(0.975)*sqrt(p.param*(1-p.param)/M),p.param + qnorm(0.975)*sqrt(p.param*(1-p.param)/M))
c(p.boot1 - qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M),p.param + qnorm(0.975)*sqrt(p.boot1*(1-p.boot1)/M))
c(p.boot2 - qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M),p.boot2 + qnorm(0.975)*sqrt(p.boot2*(1-p.param)/M))

