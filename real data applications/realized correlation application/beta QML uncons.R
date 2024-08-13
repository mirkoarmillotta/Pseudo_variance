
### Functions for unrestricted PVQMLE of interval data
### with mean and variance of beta distribution


### unrestricted PVQMLE
QMLu = function(y,theta){
  d=exp(theta[1])/(1+exp(theta[1]))
  a=exp(theta[2])/(1+exp(theta[2]))
  b=exp(theta[3])/(1+exp(theta[3]))
  phi=exp(theta[4])
  dd=exp(theta[5])/(1+exp(theta[5]))
  aa=exp(theta[6])/(1+exp(theta[6]))
  bb=exp(theta[7])/(1+exp(theta[7]))
  n=length(y)
  lambda = mu = vector()
  lambda[1] = mu[1] = mean(y)
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    mu[t] = dd + aa*y[t-1] + bb*mu[t-1]
  }
  
  if(all(lambda<1 & mu<1)==FALSE){
    llik=10e6
  }else{
    y0=y[-1]
    mt=lambda[-1]
    mt2=mu[-1]
    vt=mt2*(1-mt2)/(1+phi)
    llik=-sum(dnorm(y0,mt,sqrt(vt),log=T))
  }
  
  llik
}

# Score unrestricted PVQMLE
QMLu_sc = function(theta, y){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  dd=theta[5]
  aa=theta[6]
  bb=theta[7]
  n=length(y)
  lambda = mu = vector()
  lambda[1] = mu[1] = mean(y)
  dla =  dmu = matrix(NA,3,n)
  dla[,1] = dmu[,1] = rep(0,3)
  ss = vector()
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    mu[t] = dd + aa*y[t-1] + bb*mu[t-1]
    la = c(1, y[t-1], lambda[t-1])
    dla[,t] = la + b*dla[,t-1]
    muu = c(1, y[t-1], mu[t-1])
    dmu[,t] = muu + bb*dmu[,t-1]
  }
  
  y0=y[-1]
  mt=lambda[-1]
  mt2=mu[-1]
  vt=mt2*(1-mt2)/(1+phi)
  dmt = cbind(t(dla[,-1]), 0, 0,0,0)
  dmt2 = cbind(0,0,0, 0, t(dmu[,-1]))
  dpsi = (1-2*mt2)/(1+phi)*t(dmu[,-1])
  dvt = cbind(0,0,0, -mt2*(1-mt2)/(1+phi)^2, dpsi)
  ss = -apply( -1/(2*vt)*dvt+(y0-mt)/vt*dmt+(y0-mt)^2/(2*vt^2)*dvt, 2, sum )
  ss
}

# Fisher matrix unrestricted PVQMLE
QMLu_out = function(theta, y){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  dd=theta[5]
  aa=theta[6]
  bb=theta[7]
  n=length(y)
  lambda = mu = vector()
  lambda[1] = mu[1] = mean(y)
  dla =  dmu = matrix(NA,3,n)
  dla[,1] = dmu[,1] = rep(0,3)
  s = matrix(0,n-1,7)
  out = matrix(0,7,7)
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    mu[t] = dd + aa*y[t-1] + bb*mu[t-1]
    la = c(1, y[t-1], lambda[t-1])
    dla[,t] = la + b*dla[,t-1]
    muu = c(1, y[t-1], mu[t-1])
    dmu[,t] = muu + bb*dmu[,t-1]
  }
  
  y0=y[-1]
  mt=lambda[-1]
  mt2=mu[-1]
  vt=mt2*(1-mt2)/(1+phi)
  dmt = cbind(t(dla[,-1]), 0, 0,0,0)
  dmt2 = cbind(0,0,0, 0, t(dmu[,-1]))
  dpsi = (1-2*mt2)/(1+phi)*t(dmu[,-1])
  dvt = cbind(0,0,0, -mt2*(1-mt2)/(1+phi)^2, dpsi)
  s = -1/(2*vt)*dvt+(y0-mt)/vt*dmt+(y0-mt)^2/(2*vt^2)*dvt
  
  for(t in 1:(n-1)){
    out = out + s[t,]%*%t(s[t,])
  }
  out
}

# Function to compute Hessian of unrestricted PVQMLE
QMLu_hess = function(y,theta){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  dd=theta[5]
  aa=theta[6]
  bb=theta[7]
  n=length(y)
  lambda = mu = vector()
  lambda[1] = mu[1] = mean(y)
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    mu[t] = dd + aa*y[t-1] + bb*mu[t-1]
  }
  
    y0=y[-1]
    mt=lambda[-1]
    mt2=mu[-1]
    vt=mt2*(1-mt2)/(1+phi)
    llik=-sum(dnorm(y0,mt,sqrt(vt),log=T))
  
  llik
}
