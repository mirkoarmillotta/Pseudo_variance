
### Functions for fully restricted PVQMLE of interval data with
### mean and variance of beta distribution


### fully restricted PVQMLE
QML = function(y,theta){
  d=exp(theta[1])/(1+exp(theta[1]))
  a=exp(theta[2])/(1+exp(theta[2]))
  b=exp(theta[3])/(1+exp(theta[3]))
  phi=exp(theta[4])
  n=length(y)
  lambda = vector()
  lambda[1] = mean(y)
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
  }
  
  if(all(lambda<1)==FALSE){
    llik=10e6
  }else{
    y0=y[-1]
    mt=lambda[-1]
    vt=mt*(1-mt)/(1+phi)
    llik=-sum(dnorm(y0,mt,sqrt(vt),log=T))
  }
  
  llik
}

# Score fully restricted PVQMLE
QML_sc = function(theta, y){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  n=length(y)
  lambda = vector()
  lambda[1] = mean(y)
  dla = matrix(NA,3,n)
  dla[,1] = rep(0,3)
  ss = vector()
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    la = c(1, y[t-1], lambda[t-1])
    dla[,t] = la + b*dla[,t-1]
  }
  
  y1=y[-n]
  la1 = lambda[-n]
  y0=y[-1]
  mt=lambda[-1]
  vt=mt*(1-mt)/(1+phi)
  dmt = cbind(t(dla[,-1]), 0)
  dpsi = (1-2*mt)/(1+phi)*t(dla[,-1])
  dvt = cbind(dpsi,-mt*(1-mt)/(1+phi)^2)
  ss = -apply( -1/(2*vt)*dvt+(y0-mt)/vt*dmt+(y0-mt)^2/(2*vt^2)*dvt, 2, sum )
  ss
}

# Fisher matrix fully restricted PVQMLE
QML_out = function(theta, y){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  n=length(y)
  lambda = vector()
  lambda[1] = mean(y)
  dla = matrix(NA,3,n)
  dla[,1] = rep(0,3)
  s = matrix(0,n-1,4)
  out = matrix(0,4,4)
  
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
    la = c(1, y[t-1], lambda[t-1])
    dla[,t] = la + b*dla[,t-1]
  }
  
  y1=y[-n]
  la1 = lambda[-n]
  y0=y[-1]
  mt=lambda[-1]
  vt=mt*(1-mt)/(1+phi)
  dmt = cbind(t(dla[,-1]), 0)
  dpsi = (1-2*mt)/(1+phi)*t(dla[,-1])
  dvt = cbind(dpsi,-mt*(1-mt)/(1+phi)^2)
  
  s = -1/(2*vt)*dvt+(y0-mt)/vt*dmt+(y0-mt)^2/(2*vt^2)*dvt
  for(t in 1:(n-1)){
    out = out + s[t,]%*%t(s[t,])
  }
  out
}

# Function to compute Hessian of  fully restricted PVQMLE
QML_hess = function(y,theta){
  d=theta[1]
  a=theta[2]
  b=theta[3]
  phi=theta[4]
  n=length(y)
  lambda = vector()
  lambda[1] = mean(y)
  
  for (t in 2:n){
    lambda[t] = d + a*y[t-1] + b*lambda[t-1]
  }
  
  y0=y[-1]
  mt=lambda[-1]
  vt=mt*(1-mt)/(1+phi)
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}
