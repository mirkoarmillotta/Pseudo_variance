
###############################################################################
## This file contains the code used in the paper "Pseudo-variance 
## quasi-maximum likelihood estimation of semi-parametric time series models" 
## by Mirko Armillotta and Paolo Gorgi.
##
## This file contains the code to replicate Table 1.
##
## Note: R version 4.3.2 used throughout.
##
###############################################################################

# Install spINAR package for MLE estimation of INAR model
install.packages("spINAR")


# Binomial thinning function
thinning<-function(alfa,ni)
{
  y<-rbinom(ni,1,alfa)
  thin<-sum(y)
  return(thin)
}

# Generate data from Poisson INAR(1) model 
rpinar<-function(n, alfa, lamda)
{
  x<-matrix(0,n)
  epsilon<-rpois(n,lamda)
  x0<-rpois(1,lamda/(1-alfa))
  x[1]=thinning(alfa,x0)+epsilon[1]
  
  for (i in 2:n){
    x[i]=thinning(alfa,x[i-1])+epsilon[i]
  }
  
  return(x)
}

### Fully restricted PVQMLE
QML = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1-a)*y1+la
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

### Partially restricted PVQMLE (intercept unrestricted)
QML1 = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  v=exp(theta[3])
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1-a)*y1+v
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

### Unrestricted PVQMLE
QML2 = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  v=exp(theta[3])
  a1=exp(theta[4])
  
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a1*y1+v
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

### Partially restricted PVQMLE (a unrestricted)
QML3 = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  a1=exp(theta[3])
  
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a1*y1+la
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}


### Poisson QMLE
pQML = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  n=length(y)
  
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  llik=sum(dpois(y0,mt,log=T))
  -llik
}


B=1000
QML_est=matrix(NA,B,2)
QML1_est=matrix(NA,B,3)
QML3_est=matrix(NA,B,3)
QML2_est=matrix(NA,B,4)
pQML_est=matrix(NA,B,2)
OLS_est=matrix(NA,B,2)
WLS2_est=matrix(NA,B,2)
WLS2s_est=matrix(NA,B,2)
WLS_un=matrix(NA,B,2)
MLE_est=matrix(NA,B,2)

n=100
n=500
n=2000

a=0.85
la=3



# Generate simulation results of Table 2 
# with specified sample dimension n

set.seed(1234)

for(i in 1:B){
  
  y=rpinar(n,a,la)
  theta=c(a,la)
  ini=c(log(theta[1]/(1-theta[1])),log(theta[2]))
  ini1=c(log(theta[1]/(1-theta[1])),log(theta[2]),log(la))
  ini2=c(log(theta[1]/(1-theta[1])),log(theta[2]),log(la),log(a*(1-a)))
  ini3=c(log(theta[1]/(1-theta[1])),log(theta[2]),log(a*(1-a)))
  
  
  est_pQML=optim(ini, function(x) pQML(y,x))
  pQML_est[i,] = c(exp(est_pQML$par[1])/(1+exp(est_pQML$par[1])),exp(est_pQML$par[2]))
  
  est=optim(ini, function(x) QML(y,x))
  QML_est[i,] = c(exp(est$par[1])/(1+exp(est$par[1])),exp(est$par[2]))
  
  est1=optim(ini1, function(x) QML1(y,x))
  QML1_est[i,] = c(exp(est1$par[1])/(1+exp(est1$par[1])),exp(est1$par[2]),exp(est1$par[3]))
  
  est2=optim(ini2, function(x) QML2(y,x))
  QML2_est[i,] = c(exp(est2$par[1])/(1+exp(est2$par[1])),exp(est2$par[2]),exp(est2$par[3]),exp(est2$par[4]))
  
  est3=optim(ini3, function(x) QML3(y,x))
  QML3_est[i,] = c(exp(est3$par[1])/(1+exp(est3$par[1])),exp(est3$par[2]),exp(est3$par[3]))
  
  
  y1=y[-n]
  y0=y[-1]
  X=cbind(rep(1,(n-1)),y1)
  OLS_est[i,] =solve(t(X)%*%X)%*%t(X)%*%y0
  
  
  W=diag(1/(OLS_est[i,2]*(1-OLS_est[i,2])*y1+OLS_est[i,1]))
  WLS2_est[i,] =solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y0
  
  s2 = mean((y0-OLS_est[i,2]*y1-OLS_est[i,1])^2-OLS_est[i,2]*(1-OLS_est[i,2])*y1)
  W=diag(1/(OLS_est[i,2]*(1-OLS_est[i,2])*y1+s2))
  WLS2s_est[i,] =solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y0
  
  Wun=diag(1/(a*(1-a)*y1+la))
  WLS_un[i,]=solve(t(X)%*%Wun%*%X)%*%t(X)%*%Wun%*%y0
  
  MLE_est[i,]=spINAR::spinar_est_param(y, p=1, "ml", "poi")
  
  print(i)
}

# RMSE
sqrt(mean((MLE_est[,1]-a)^2))    # MLE
sqrt(mean((QML_est[,1]-a)^2))    # fully constrained
sqrt(mean((QML1_est[,1]-a)^2))   # partially constrained (free intercept)
sqrt(mean((QML3_est[,1]-a)^2))   # partially constrained (free alpha)
sqrt(mean((QML2_est[,1]-a)^2))   # unconstrained
sqrt(mean((pQML_est[,1]-a)^2))   # Poisson QMLE
sqrt(mean((OLS_est[,2]-a)^2))    # OLS
sqrt(mean((WLS2_est[,2]-a)^2))   # two step WLS Francq
sqrt(mean((WLS2s_est[,2]-a)^2))  # two step WLS Francq (est. variance)
sqrt(mean((WLS_un[,2]-a)^2))     # unfeasible WLS

# Bias
mean(MLE_est[,1])-a    # MLE
mean(QML_est[,1])-a    # fully constrained
mean(QML1_est[,1])-a   # partially constrained (free intercept)
mean(QML3_est[,1])-a   # partially constrained (free alpha)
mean(QML2_est[,1])-a   # unconstrained
mean(pQML_est[,1])-a   # Poisson QMLE
mean(OLS_est[,2])-a    # OLS
mean(WLS2_est[,2])-a   # two step WLS Francq
mean(WLS2s_est[,2])-a  # two step WLS Francq (est. variance)
mean(WLS_un[,2])-a     # unfeasible WLS


# RMSE
sqrt(mean((MLE_est[,2]-la)^2))   # MLE
sqrt(mean((QML_est[,2]-la)^2))   # fully constrained
sqrt(mean((QML1_est[,2]-la)^2))  # partially constrained (free intercept)
sqrt(mean((QML3_est[,2]-la)^2))  # partially constrained (free alpha)
sqrt(mean((QML2_est[,2]-la)^2))  # unconstrained
sqrt(mean((pQML_est[,2]-la)^2))  # Poisson QMLE
sqrt(mean((OLS_est[,1]-la)^2))   # OLS
sqrt(mean((WLS2_est[,1]-la)^2))  # two step WLS Francq
sqrt(mean((WLS2s_est[,1]-la)^2))  # two step WLS Francq (est. variance)
sqrt(mean((WLS_un[,1]-la)^2))    # unfeasible WLS

# Bias
mean(MLE_est[,2])-la    # MLE
mean(QML_est[,2])-la    # fully constrained
mean(QML1_est[,2])-la   # partially constrained (free intercept)
mean(QML3_est[,2])-la   # partially constrained (free alpha)
mean(QML2_est[,2])-la   # unconstrained
mean(pQML_est[,2])-la   # Poisson QMLE
mean(OLS_est[,1])-la    # OLS
mean(WLS2_est[,1])-la   # two step WLS Francq
mean(WLS2s_est[,1])-la  # two step WLS Francq (est. variance)
mean(WLS_un[,1])-la     # unfeasible WLS



library(xtable)


bias = matrix(0, 10, 2)

bias[1,1] = mean(pQML_est[,2])-la   # Poisson QMLE
bias[2,1] = mean(OLS_est[,1])-la    # OLS
bias[3,1] = mean(WLS2_est[,1])-la   # two step WLS Francq
bias[4,1] = mean(WLS2s_est[,1])-la  # two step WLS Francq (est. variance)
bias[5,1] = mean(WLS_un[,1])-la     # unfeasible WLS
bias[6,1] = mean(QML2_est[,2])-la   # unconstrained
bias[7,1] = mean(QML1_est[,2])-la   # partially constrained (free intercept)
bias[8,1] = mean(QML3_est[,2])-la   # partially constrained (free alpha)
bias[9,1] = mean(QML_est[,2])-la    # fully constrained
bias[10,1] = mean(MLE_est[,2])-la   # MLE

bias[1,2] = mean(pQML_est[,1])-a   # Poisson QMLE
bias[2,2] = mean(OLS_est[,2])-a    # OLS
bias[3,2] = mean(WLS2_est[,2])-a   # two step WLS Francq
bias[4,2] = mean(WLS2s_est[,2])-a  # two step WLS Francq (est. variance)
bias[5,2] = mean(WLS_un[,2])-a     # unfeasible WLS
bias[6,2] = mean(QML2_est[,1])-a   # unconstrained
bias[7,2] = mean(QML1_est[,1])-a   # partially constrained (free intercept)
bias[8,2] = mean(QML3_est[,1])-a   # partially constrained (free alpha)
bias[9,2] = mean(QML_est[,1])-a    # fully constrained
bias[10,2] = mean(MLE_est[,1])-a   # MLE

rmse = matrix(0, 10, 2)

rmse[1,1] = sqrt(mean((pQML_est[,2]-la)^2))  # Poisson QMLE
rmse[2,1] = sqrt(mean((OLS_est[,1]-la)^2))   # OLS
rmse[3,1] = sqrt(mean((WLS2_est[,1]-la)^2))  # two step WLS Francq
rmse[4,1] = sqrt(mean((WLS2s_est[,1]-la)^2))  # two step WLS Francq (est. variance)
rmse[5,1] = sqrt(mean((WLS_un[,1]-la)^2))    # unfeasible WLS
rmse[6,1] = sqrt(mean((QML2_est[,2]-la)^2))  # unconstrained
rmse[7,1] = sqrt(mean((QML1_est[,2]-la)^2))  # partially constrained (free intercept)
rmse[8,1] = sqrt(mean((QML3_est[,2]-la)^2))  # partially constrained (free alpha)
rmse[9,1] = sqrt(mean((QML_est[,2]-la)^2))   # fully constrained
rmse[10,1] = sqrt(mean((MLE_est[,2]-la)^2))  # MLE

rmse[1,2] = sqrt(mean((pQML_est[,1]-a)^2))   # Poisson QMLE
rmse[2,2] = sqrt(mean((OLS_est[,2]-a)^2))    # OLS
rmse[3,2] = sqrt(mean((WLS2_est[,2]-a)^2))   # two step WLS Francq
rmse[4,2] = sqrt(mean((WLS2s_est[,2]-a)^2))  # two step WLS Francq (est. variance)
rmse[5,2] = sqrt(mean((WLS_un[,2]-a)^2))     # unfeasible WLS
rmse[6,2] = sqrt(mean((QML2_est[,1]-a)^2))   # unconstrained
rmse[7,2] = sqrt(mean((QML1_est[,1]-a)^2))   # partially constrained (free intercept)
rmse[8,2] = sqrt(mean((QML3_est[,1]-a)^2))   # partially constrained (free alpha)
rmse[9,2] = sqrt(mean((QML_est[,1]-a)^2))    # fully constrained
rmse[10,2] = sqrt(mean((MLE_est[,1]-a)^2))   # MLE


tab <- cbind(bias[,1], rmse[,1], bias[,2], rmse[,2])
tab2 <- cbind(bias[,1], rmse[,1], bias[,2], rmse[,2])
tab3 <- cbind(bias[,1], rmse[,1], bias[,2], rmse[,2])
tabb <- cbind(tab, tab2, tab3)
print(xtable(tabb, type = "latex", digits=4),
      include.rownames=FALSE, file = "Tab_INAR_eff.tex")
