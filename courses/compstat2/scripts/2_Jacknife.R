# We can directly use the jacknife function of the boostrap package
library(bootstrap)
help(jackknife)

x <- rnorm(20)               
results <- jackknife(x,mean)      
results

theta = function(x){
  return(mean(x))
}
results2 <- jackknife(x,theta)      
results2

# in theory:
sd(x)/sqrt(length(x))


# Test score data ---------------------------------------------------------
library(MVT)
data(examScor)
View(examScor)

Shat = cov(examScor)
eigen(Shat)
lambdas = eigen(Shat)$values
max(lambdas)/sum(lambdas)
n = dim(examScor)[1]

# parameter that we want to compute:
theta <- function(x,datamatrix){
  Shat = cov(datamatrix[x,])
  lambdas = eigen(Shat)$values
  return(max(lambdas)/sum(lambdas))
}

# Estimated value on data:
theta(1:n,examScor)

# jacknife estimate
results.score <- jackknife(1:n,theta,datamatrix=examScor)      
results.score

layout(1)
pdf('Jacknife_Scoredata.pdf',height=5)
hist(results.score$jack.values,main='Jacknife replications for test score data',col='darkslategray2')
dev.off()



# Jacknife is inconsistent for the median ------------------------------------
# data simulated from a Binomial distribution with p=0.5 n = 10
p.binom = 0.5
n.binom = 10
median.true = n.binom*p.binom
n = 20

set.seed(2)
data = rbinom(n,size = n.binom, prob = p.binom)
jack.result = jackknife(data,median)
jack.result

# is the Jacknife estimate consistent? 
# We can repeat this example increasing the sample size n

# MC Simulation:
n.sim = seq(10,200,by=20)
estimate.sd = numeric(length(n.sim))
estimate.bias = numeric(length(n.sim))

M = 1000 # number of simulated data sets
for(scenario in 1:length(n.sim)){
  n = n.sim[scenario]
  estimated.medians = numeric(M)
  for(sim in 1:M){
    data = rbinom(n,size = n.binom, prob = p.binom)
    estimated.medians[sim] = median(data)
  }
  estimate.sd[scenario] = sd(estimated.medians)
  estimate.bias[scenario] = (mean(estimated.medians) - median.true)^2
}
estimate.bias
estimate.sd

# Jacknife:
jacknife.sd = numeric(length(n.sim))
jacknife.bias = numeric(length(n.sim))

for(scenario in 1:length(n.sim)){
  n = n.sim[scenario]
  data = rbinom(n,size = n.binom, prob = p.binom)
  jack.result = jackknife(data,median)
  jacknife.sd[scenario] = jack.result$jack.se
  jacknife.bias[scenario] = jack.result$jack.bias
}

jacknife.sd
jacknife.bias

