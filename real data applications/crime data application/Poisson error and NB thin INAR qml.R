
### Functions for the estimation of Poisson error and NB thinning INAR PVQMLE

QMLPNB = function(y,theta,r){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  n=length(y)
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1+a/r)*y1+la
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}

QMLPNB_hess = function(y,theta,r){
  a=theta[1]
  la=theta[2]
  n=length(y)
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1+a/r)*y1+la
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}


sc_QMLPNB = function(theta, y, r){
  a=exp(theta[1])/(1+exp(theta[1]))
  la=exp(theta[2])
  n=length(y)
  ss = vector()
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1+a/r)*y1+la
  dmt = cbind(y1,1)
  dvt = cbind(y1+2*a/r*y1,1)
  for(j in 1:2) ss[j] = sum( -1/(2*vt)*dvt[,j]+(y0-mt)/vt*dmt[,j]+(y0-mt)^2/(2*vt^2)*dvt[,j] )
  -ss
}

out_QMLPNB = function(theta, y, r){
  a=theta[1]
  la=theta[2]
  n=length(y)
  out = matrix(0,2,2)
  s = matrix(0,n-1,2)
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1+a/r)*y1+la
  dmt = cbind(y1,1)
  dvt = cbind(y1+2*a/r*y1,1)
  for(j in 1:2){
    s[,j] = -1/(2*vt)*dvt[,j]+(y0-mt)/vt*dmt[,j]+(y0-mt)^2/(2*vt^2)*dvt[,j]
  }
  for(t in 1:(n-1)){
    out = out + s[t,]%*%t(s[t,])
  }
  out
}
