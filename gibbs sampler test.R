
library(mvtnorm) #for multivariate normal distribution
library(MCMCpack) #for inverse gamma distribution
library(coda) #for trace/density plots

#setup parameters/data
b0.real <- 4.5
b1.real <- 0.25
sigma <- 0.3
n <- 100
x <- runif(n,0,10)
X <- cbind(1,x)

y <- rnorm(n,b0.real+b1.real*x,sigma)

plot(x,y)

#define priors
b.prior.mean <- c(0,0)
b.prior.covar <- matrix(c(100,0,0,100),nrow=2)
b.prior.prec <- solve(b.prior.covar)

sigma2.prior.s1 <- 1
sigma2.prior.s2 <- 1

#initialization
b.curr <- c(0,0)
sigma2.curr <- 1

iter <- 5000

b <- matrix(nrow=iter,ncol=2)
sigma2 <- numeric(iter)

#gibbs sampling
for (i in 1:5000) {
  
  #conditional posterior for b
  v <- sigma2.curr^-1 * t(X) %*% y + b.prior.prec %*% b.prior.mean
  V.inv <- sigma2.curr^-1 * t(X) %*% X + b.prior.prec
  V <- solve(V.inv)
  
  #sample b
  b.curr <- t(rmvnorm(1,mean=V%*%v,sigma=V))
  
  pred.curr <- X %*% b.curr
  resid.curr <- y - pred.curr
  
  #conditional posterior for sigma2
  mu1 <- sigma2.prior.s1 + n/2
  mu2 <- sigma2.prior.s2 + 0.5 * t(resid.curr) %*% resid.curr
  
  #sample sigma2
  sigma2.curr <- rinvgamma(1,mu1,mu2)
  
  #save sample
  b[i,] <- b.curr
  sigma2[i] <- sigma2.curr
    
}

#plot results
colnames(b) <- c("b0","b1")
mcmc.out <- as.mcmc(cbind(b,sigma2))
plot(mcmc.out)


