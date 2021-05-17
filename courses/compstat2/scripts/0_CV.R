n = 30
set.seed(15072017)
epsilon = rnorm(n)
X = runif(n,0,5)

Y = 2 + 3*X - X^2 + epsilon
plot(X,Y,pch=16)
regr = lm(Y ~ X)
abline(regr$coefficients,col='blue',lwd=3)

regr2 = lm(Y~poly(X,20,raw=TRUE))
xnew = data.frame(X=seq(0,5,len=500))
ynew = predict(regr2,newdata = xnew)
lines(xnew$X,ynew,lwd=3,col=2)
points(X,Y,pch=16)



# validation set
n = 30
set.seed(15072017)
epsilon = rnorm(n)
X = runif(n,0,5)

Y = 2 + 3*X - X^2 + epsilon
plot(X,Y,pch=16)


set.seed(200)
training = sample(1:n,size = 10,replace = FALSE)
test = (1:30)[-training]

degmax = 6
Xpoly = poly(X,degmax)
MSEtest = numeric(degmax)

for(degree in 1:degmax){
  Ytr = Y[training]
  Xtr = Xpoly[training,1:degree]
  regr2 = lm(Ytr~Xtr)
  Xnew = Xpoly[test,1:degree]
  ypred = cbind(1,(Xnew)) %*% regr2$coefficients
  MSEtest[degree] = mean((Y[test]-ypred)^2)
}
plot(1:degmax,MSEtest,type='b',lwd=2,xlab='degree of polynomial',ylab='Test set MSE',ylim=c(0,10))


# retry several times
set.seed(200)
training = sample(1:n,size = 10,replace = FALSE)
test = (1:30)[-training]

degmax = 6
Xpoly = poly(X,degmax)
MSEtest = numeric(degmax)

for(degree in 1:degmax){
  Ytr = Y[training]
  Xtr = Xpoly[training,1:degree]
  regr2 = lm(Ytr~Xtr)
  Xnew = Xpoly[test,1:degree]
  ypred = cbind(1,(Xnew)) %*% regr2$coefficients
  MSEtest[degree] = mean((Y[test]-ypred)^2)
}
plot(1:degmax,MSEtest,type='b',lwd=2,xlab='degree of polynomial',ylab='Test set MSE',ylim=c(0,10))

for(ii in 1:10){
  set.seed(ii)
  training = sample(1:n,size = 10,replace = FALSE)
  test = (1:30)[-training]
  
  degmax = 6
  Xpoly = poly(X,degmax)
  MSEtest = numeric(degmax)
  
  for(degree in 1:degmax){
    Ytr = Y[training]
    Xtr = Xpoly[training,1:degree]
    regr2 = lm(Ytr~Xtr)
    Xnew = Xpoly[test,1:degree]
    ypred = cbind(1,(Xnew)) %*% regr2$coefficients
    MSEtest[degree] = mean((Y[test]-ypred)^2)
  }
  
  lines(1:degmax,MSEtest,type='b',lwd=2,col=ii)
  
}


# CV

set.seed(200)
K = 5
folds = sample(1:K,size = n,replace = TRUE)


degmax = 6
Xpoly = poly(X,degmax)
MSEtest = matrix(ncol=degmax,nrow=K)

for(k in 1:K){
  for(degree in 1:degmax){
    test = which(folds==k)
    training = (1:n)[-test]
    Ytr = Y[training]
    Xtr = Xpoly[training,1:degree]
    regr2 = lm(Ytr~Xtr)
    Xnew = Xpoly[test,1:degree]
    ypred = cbind(1,(Xnew)) %*% regr2$coefficients
    MSEtest[k,degree] = mean((Y[test]-ypred)^2)
  }
  
}
MSEtest = colMeans(MSEtest)

plot(1:degmax,MSEtest,type='b',lwd=2,xlab='degree of polynomial',ylab='Test set MSE',ylim=c(0,10))


set.seed(200)
K = 5
folds = sample(1:K,size = n,replace = TRUE)
degmax = 6
Xpoly = poly(X,degmax)
MSEtest = matrix(ncol=degmax,nrow=K)

for(k in 1:K){
  for(degree in 1:degmax){
    test = which(folds==k)
    training = (1:n)[-test]
    Ytr = Y[training]
    Xtr = Xpoly[training,1:degree]
    regr2 = lm(Ytr~Xtr)
    Xnew = Xpoly[test,1:degree]
    ypred = cbind(1,(Xnew)) %*% regr2$coefficients
    MSEtest[k,degree] = mean((Y[test]-ypred)^2)
  }
  
}
MSEtest = colMeans(MSEtest)
plot(1:degmax,MSEtest,type='b',lwd=2,xlab='degree of polynomial',ylab='Test set MSE',ylim=c(0,10))

for(ii in 1:10){
  set.seed(ii)
  folds = sample(1:K,size = n,replace = TRUE)
  
  
  degmax = 6
  Xpoly = poly(X,degmax)
  MSEtest = matrix(ncol=degmax,nrow=K)
  
  for(k in 1:K){
    for(degree in 1:degmax){
      test = which(folds==k)
      training = (1:n)[-test]
      Ytr = Y[training,drop=FALSE]
      Xtr = Xpoly[training,1:degree,drop=FALSE]
      regr2 = lm(Ytr~Xtr)
      Xnew = Xpoly[test,1:degree,drop=FALSE]
      ypred = cbind(1,(Xnew)) %*% regr2$coefficients
      MSEtest[k,degree] = mean((Y[test]-ypred)^2)
    }
    
  }
  MSEtest = colMeans(MSEtest)
  
  lines(1:degmax,MSEtest,type='b',lwd=2,col=ii)
  
}

# LOOCV
set.seed(200)
K = 30
folds = 1:30
degmax = 6
Xpoly = poly(X,degmax)
MSEtest = matrix(ncol=degmax,nrow=K)

for(k in 1:K){
  for(degree in 1:degmax){
    test = which(folds==k)
    training = (1:n)[-test,drop=FALSE]
    Ytr = Y[training,drop=FALSE]
    Xtr = Xpoly[training,1:degree,drop=FALSE]
    regr2 = lm(Ytr~Xtr)
    Xnew = Xpoly[test,1:degree,drop=FALSE]
    ypred = cbind(1,(Xnew)) %*% regr2$coefficients
    MSEtest[k,degree] = mean((Y[test]-ypred)^2)
  }
  
}
MSEtest = colMeans(MSEtest)
plot(1:degmax,MSEtest,type='b',lwd=2,xlab='degree of polynomial',ylab='Test set MSE',ylim=c(0,10))

