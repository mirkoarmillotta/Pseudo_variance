

###############################################################################
## This file contains the code used in the paper "Pseudo-variance 
## quasi-maximum likelihood estimation of semi-parametric time series models" 
## by Mirko Armillotta and Paolo Gorgi.
##
## This file contains the code for Wald test simulation of Section 5
## regarding simulation with INAR(1).
##
## Note: R version 4.3.2 used throughout.
##
###############################################################################


## Poisson thinning
thinning_pois <- function(alfa,ni)
{
  y<-rpois(ni,alfa)
  thin<-sum(y)
  return(thin)
}

## Negative Binomial (NB) thinning
thinning_nb<-function(alfa, ni, ri)
{
  y <- rnbinom(ni, size=ri, mu=alfa)
  thin<-sum(y)
  return(thin)
}

## Binomial-Negative Binomial (BiNB) thinning
thinning_binb <- function(mu, ni, pi)
{
  y_1 <- rbinom(ni, size=1, prob=pi)
  y_2 <- rnbinom(ni, size=1, mu=mu)
  y = y_1 + y_2
  thin<-sum(y)
  return(thin)
}

# Generate data from INAR with Poisson thinning and errors
rpinar<-function(n,alfa, lamda)
{
  x<-matrix(0,n)
  epsilon<-rpois(n,lamda)
  x0<-rpois(1,lamda/(1-alfa))
  x[1]=thinning_pois(alfa,x0)+epsilon[1]
  
  for (i in 2:n){
    x[i]=thinning_pois(alfa,x[i-1])+epsilon[i]
  }
  
  return(x)
}

# Generate data from INAR with NB thinning and Poisson errors
rp_inar<-function(n,alfa, lamda, rr)
{
  x<-matrix(0,n)
  epsilon<-rpois(n,lamda)
  x0<-rpois(1,lamda/(1-alfa))
  x[1]=thinning_nb(alfa, x0, rr)+epsilon[1]
  
  for (i in 2:n){
    x[i]=thinning_nb(alfa, x[i-1], rr)+epsilon[i]
  }
  
  return(x)
}


# Generate data from INAR with BiNB thinning and Poisson errors
rbinb_inar<-function(n, mu, lamda, pp)
{
  x<-matrix(0,n)
  epsilon<-rpois(n,lamda)
  x0<-rpois(1,lamda/(1-mu-pp))
  x[1]=thinning_binb(mu, x0, pp)+epsilon[1]
  
  for (i in 2:n){
    x[i]=thinning_binb(mu, x[i-1], pp)+epsilon[i]
  }
  
  return(x)
}


### Unrestricted PQMLE
QML2 = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  v=exp(theta[3])
  b=exp(theta[4])
  
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=b*y1+v
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

### Function for computation of Hessian matrix of unrestricted PQMLE
QML2_hess = function(y,theta){
  
  a=theta[1]
  la=theta[2]
  v=theta[3]
  b=theta[4]
  
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=b*y1+v
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

### Function for computation of the score of unrestricted PQMLE
sc_QML2 = function(theta, y){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  v=exp(theta[3])
  b=exp(theta[4])
  
  n=length(y)
  ss = vector()
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=b*y1+v
  dmt = cbind(y1,1,0,0)
  dvt = cbind(0,0,1,y1)
  for(j in 1:4) ss[j] = sum( -1/(2*vt)*dvt[,j]+(y0-mt)/vt*dmt[,j]+(y0-mt)^2/(2*vt^2)*dvt[,j] )
  -ss
}

### Function for computation of Fisher information matrix of unrestricted PQMLE
out_QML2 = function(theta, y){
  
  a=theta[1]
  la=theta[2]
  v=theta[3]
  b=theta[4]
  
  n=length(y)
  out = matrix(0,4,4)
  s = matrix(0,n-1,4)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=b*y1+v
  dmt = cbind(y1,1,0,0)
  dvt = cbind(0,0,1,y1)
  for(j in 1:4){
    s[,j] = -1/(2*vt)*dvt[,j]+(y0-mt)/vt*dmt[,j]+(y0-mt)^2/(2*vt^2)*dvt[,j]
  }
  for(t in 1:(n-1)){
    out = out + s[t,]%*%t(s[t,])
  }
  out
}


### Function for computation of empirical size and power

freq <- function(x, alpha){
  z <- vector()
  s <- sum(!is.na(x))
  z[1]<-sum(x>qchisq(alpha[1],1), na.rm = T)/s
  z[2]<-sum(x>qchisq(alpha[2],1), na.rm = T)/s
  z[3]<-sum(x>qchisq(alpha[3],1), na.rm = T)/s
  return(z)
}



# select sample size

n=100
n=250
n=500
n=1000
n=2000

# True value of coefficients
a=0.75
#a=0.99   # close to non-stationary region
la=1

#r=seq(0.5, 3, by=0.5)  # for NB  thinning
r=seq(0, 0.15,length.out=6) # for BiNB thinning


### Simulations for empirical size of Wald specification test

B=5000  # number of simulations
QML2_est <- list()
QML2_se <- list()
QML2_t <- list()
W=matrix(NA,B,length(r))
score <- list()
I <- list()

j=1

set.seed(1234)

ri <- 0.5
score[[j]] <- matrix(NA, B, 4)
QML2_est[[j]] <- matrix(NA, B, 4)
I[[j]] <- list()
QML2_se[[j]] <- matrix(NA, B, 4)
QML2_t[[j]] <- matrix(NA, B, 4)
  
for(i in 1:B){
    
    y=rpinar(n,a,la)
    
    # run to generate data with outlier
    #y[50] = 3*sd(y)+mean(y)
    
    theta=c(a,la)
    ini2=c(log(theta[1]/(1-theta[1])),log(theta[2]), log(la),log(a+a^2/ri))
    
    est2 <- optim(par = ini2, fn = QML2, method = "BFGS", hessian=F, y=y,
                  control = list(ndeps=rep(1e-5,4), reltol=1e-12, maxit=1000))
    
    score[[j]][i,] = sc_QML2(y=y, theta=est2$par)
    
    alpha = exp(est2$par[1])/(1+exp(est2$par[1]))
    omega1 = exp(est2$par[2])
    omega2 = exp(est2$par[3])
    beta = exp(est2$par[4])
    
    QML2_est[[j]][i,] = c(alpha, omega1, omega2, beta)
    
    try(hess <- optimHess(par = QML2_est[[j]][i,], fn = QML2_hess, y=y))
    hes = solve(hess)
    
    out = out_QML2(y=y, theta = QML2_est[[j]][i,])
    
    I[[j]][[i]] = hes%*%out%*%hes
    
    QML2_se[[j]][i,] = sqrt(diag(I[[j]][[i]]))
    QML2_t[[j]][i,] = QML2_est[[j]][i,]/QML2_se[[j]][i,]
    
    
    # Wald test for Poisson thinning
    
    g = beta - alpha
    G = c(-1,0,0,1)
    W[i,j] = g^2/(t(G)%*%I[[j]][[i]]%*%G)
    
    print(i)
    
}


# summary estimation

est = se = t = matrix(NA,1,4) 

est <- apply(QML2_est[[1]], 2, mean)
se <- apply(score[[1]], 2, mean)
t <- apply(QML2_t[[1]], 2, mean)

#est100 = est
#est250 = est
#est500 = est
#est1000 = est
est2000 = est

estt <- cbind(est100, est250, est500, est1000, est2000)
estt = t(estt)

library(xtable)
print(xtable(estt, type = "latex", digits=c(4,4,4,4,4)),
      include.rownames=FALSE, file = "Tab_est.tex")



## Results empirical size (10%, 5%, 1%)
## (run the chit corresponding to the chosen value of n)

#chit100 <- apply(W, 2, freq, alpha=c(0.9,0.95,0.99))[,1]
#chit250 <- apply(W, 2, freq, alpha=c(0.9,0.95,0.99))[,1]
#chit500 <- apply(W, 2, freq, alpha=c(0.9,0.95,0.99))[,1]
#chit1000 <- apply(W, 2, freq, alpha=c(0.9,0.95,0.99))[,1]
#chit2000 <- apply(W, 2, freq, alpha=c(0.9,0.95,0.99))[,1]


chit <- cbind(chit100, chit250, chit500, chit1000, chit2000)

library(xtable)
print(xtable(chit, type = "latex", digits=c(4,4,4,4,4,4)),
      include.rownames=FALSE, file = "Tab_size.tex")



### Simulations for empirical power of Wald specification test

# select sample size

n=100
n=250
n=500
n=1000
n=2000


B=5000  # number of simulations
QML2_est <- list()
QML2_se <- list()
QML2_t <- list()
W=matrix(NA,B,length(r))
score <- list()
I <- list()



set.seed(1234)

for(j in 1:length(r)){
  
  ri <- r[j]
  
  print("--------------------")
  
  score[[j]] <- matrix(NA, B, 4)
  QML2_est[[j]] <- matrix(NA, B, 4)
  I[[j]] <- list()
  QML2_se[[j]] <- matrix(NA, B, 4)
  QML2_t[[j]] <- matrix(NA, B, 4)
  
  for(i in 1:B){
    
    #y=rp_inar(n,a,la,ri) # for NB thinning
    y=rbinb_inar(n,a,la,ri) # for BiNB thinning
    
    #run to generate data with outlier
    #y[50] = 3*sd(y)+mean(y)
    
    theta=c(a,la)
    
    # for NB thinning
    #ini2=c(log(theta[1]/(1-theta[1])),log(theta[2]), log(la),log(a+a^2/ri))
    
    # for BiNB thinning
    ini2=c(log(theta[1]/(1-theta[1])),log(theta[2]), log(la),log(ri*(1-ri)+a*(1+a)))
    
    est2 <- optim(par = ini2, fn = QML2, method = "BFGS", hessian=F, y=y,
                  control = list(ndeps=rep(1e-5,4), reltol=1e-12, maxit=1000))
    
    score[[j]][i,] = sc_QML2(y=y, theta=est2$par)
    
    alpha = exp(est2$par[1])/(1+exp(est2$par[1]))
    omega1 = exp(est2$par[2])
    omega2 = exp(est2$par[3])
    beta = exp(est2$par[4])
    
    QML2_est[[j]][i,] = c(alpha, omega1, omega2, beta)
    
    try(hess <- optimHess(par = QML2_est[[j]][i,], fn = QML2_hess, y=y))
    hes = solve(hess)
    
    out = out_QML2(y=y, theta = QML2_est[[j]][i,])
    
    I[[j]][[i]] = hes%*%out%*%hes
    
    QML2_se[[j]][i,] = sqrt(diag(I[[j]][[i]]))
    QML2_t[[j]][i,] = QML2_est[[j]][i,]/QML2_se[[j]][i,]
    
    
    # Wald test for Poisson thinning
    
    g = beta - alpha
    G = c(-1,0,0,1)
    W[i,j] = g^2/(t(G)%*%I[[j]][[i]]%*%G)
    
    print(i)
    
  }
  
}


## (run the pow corresponding to the chosen value of n)

#pow100 = cbind(apply(W, 2, freq, alpha=c(0.9,0.95,0.99)),chit100)
#pow250 = cbind(apply(W, 2, freq, alpha=c(0.9,0.95,0.99)),chit250)
#pow500 = cbind(apply(W, 2, freq, alpha=c(0.9,0.95,0.99)),chit500)
#pow1000 = cbind(apply(W, 2, freq, alpha=c(0.9,0.95,0.99)),chit1000)
#pow2000 = cbind(apply(W, 2, freq, alpha=c(0.9,0.95,0.99)),chit2000)



## plots


#ov = 1-a/(a+a^2/r) # For NB thinning
ov = 1-1/(1+a-r)   # For BiNB thinning
ov = c(ov, 0)



windows(width=8, height=3)
par(mfrow=c(1,3))
par(mar=c(4,4.1,2.1,1.1))

plot(pow100[1,]~ov, type="b", pch=17, col="purple", ylim=c(0,1), ylab="Power",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.1")))
lines(pow250[1,]~ov, type="b", pch=19, col="green", lty=2, ylim=c(0,1), ylab="Power",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.1")))
lines(pow500[1,]~ov, type="b", pch=20, col="red", lty=3, ylim=c(0,1), ylab="Power",
      main=expression(paste(alpha, "=0.1")))
lines(pow1000[1,]~ov, type="b", pch=21, col="orange", lty=4, ylim=c(0,1), ylab="Power",
      main=expression(paste(alpha, "=0.1")))
lines(pow2000[1,]~ov, type="b", pch=18, col="blue", lty=5, ylim=c(0,1), 
      ylab="Power", main=expression(paste(alpha, "=0.1")))
#legend(0.35, 0.35, legend=c("T=100", "T=250", "T=500", "T=1000", "T=2000"),
#       col=c("purple", "green", "red", "orange", "blue"), lty=1:5, text.font=80, cex = 0.9, box.lty=0)
#legend(0.46, 0.25, legend=c("T=100", "T=250", "T=500", "T=1000", "T=2000"),
#       col=c("purple", "green", "red", "orange", "blue"), lty=1:5, text.font=80, cex = 0.7, box.lty=0)
legend(0.25, 0.35, legend=c("T=100", "T=250", "T=500", "T=1000", "T=2000"),
       col=c("purple", "green", "red", "orange", "blue"), lty=1:5, text.font=80, cex = 0.9, box.lty=0)


plot(pow100[2,]~ov, type="b", pch=17, col="purple", ylim=c(0,1), ylab="Power",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.05")))
lines(pow250[2,]~ov, type="b", pch=19, col="green", lty=2, ylim=c(0,1), ylab="",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.05")))
lines(pow500[2,]~ov, type="b", pch=20, col="red", lty=3, ylim=c(0,1), ylab="",
      main=expression(paste(alpha, "=0.05")))
lines(pow1000[2,]~ov, type="b", pch=21, col="orange", lty=4, ylim=c(0,1), ylab="",
      main=expression(paste(alpha, "=0.05")))
lines(pow2000[2,]~ov, type="b", pch=18, col="blue", lty=5, ylim=c(0,1), 
      ylab="", main=expression(paste(alpha, "=0.05")))

plot(pow100[3,]~ov, type="b", pch=17, col="purple", ylim=c(0,1), ylab="Power",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.01")))
lines(pow250[3,]~ov, type="b", pch=19, col="green", lty=2, ylim=c(0,1), ylab="",
     xlab="% overdispersion", main=expression(paste(alpha, "=0.01")))
lines(pow500[3,]~ov, type="b", pch=20, col="red", lty=3, ylim=c(0,1), ylab="",
      main=expression(paste(alpha, "=0.01")))
lines(pow1000[3,]~ov, type="b", pch=21, col="orange", lty=4, ylim=c(0,1), ylab="",
      main=expression(paste(alpha, "=0.01")))
lines(pow2000[3,]~ov, type="b", pch=18, col="blue", lty=5, ylim=c(0,1), 
      ylab="", main=expression(paste(alpha, "=0.01")))
