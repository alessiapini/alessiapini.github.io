# Monte Carlo Simulation --------------------------------------------------

# We observe a sample of size 10 drawn from a bivariate Normal distribution with Mean (0,0) and covariance matrix (2,1 ; 1,2). 
# We are interested in estimating the eigenvalues of the covariance matrix, using the plugin estimator.
# How variable is the estimator? 

library(mvtnorm)
help(rmvnorm)

mu = c(0,0)
Sigma = rbind(c(2,1),c(1,2))
n = 10

# one sample from the data:
set.seed(1)
data = rmvnorm(n,mu,Sigma)
plot(data)

# true eigenvalues
theta = eigen(Sigma)$values
theta

# plug-in estimator
thetahat = eigen(cov(data))$values
thetahat

B = 10000 # number of MC iterations
thetastar = matrix(nrow=B,ncol=2) # simulated values of theta
for(sim in 1:B){
  data.sim = rmvnorm(n,mu,Sigma)
  thetastar[sim,] = eigen(cov(data.sim))$values
}

# Estimated values
mean(thetastar[,1]) 
mean(thetastar[,2]) 

# standard deviation:
sd(thetastar[,1])
sd(thetastar[,2])

# distribution:
layout(cbind(1,2))
hist(thetastar[,1],main='First Eigenvalue',xlim=c(0,10),probability=TRUE,ylim=c(0,1))
hist(thetastar[,2],main='Second Eigenvalue',xlim=c(0,10),probability=TRUE,ylim=c(0,1))

