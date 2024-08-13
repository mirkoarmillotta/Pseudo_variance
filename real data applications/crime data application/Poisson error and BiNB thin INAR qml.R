
### Functions for the estimation of Poisson error and BiNB thinning INAR PVQMLE

QMLPBiNB = function(y,theta){
  a=exp(theta[1])/(1+exp(theta[1]))
  p=exp(theta[2])/(1+exp(theta[2]))
  
  if(a+p<1){
    
    la=exp(theta[3])
    n=length(y)
    y1=y[-n]
    y0=y[-1]
    mt=(a+p)*y1+la
    vt=(p*(1-p)+a*(1+a))*y1+la
    llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
    -llik
    
   } else{
     10e+7
   }
  
}

QMLPBiNB_hess = function(y,theta){
  a=theta[1]
  p=theta[2]
  la=theta[3]
  n=length(y)
  y1=y[-n]
  y0=y[-1]
  mt=(a+p)*y1+la
  vt=(p*(1-p)+a*(1+a))*y1+la
  llik=sum(dnorm(y0,mt,sqrt(vt),log=T))
  -llik
}



out_QMLPBiNB = function(theta, y, r){
  a=theta[1]
  p=theta[2]
  la=theta[3]
  n=length(y)
  out = matrix(0,3,3)
  s = matrix(0,n-1,3)
  y1=y[-n]
  y0=y[-1]
  mt=(a+p)*y1+la
  vt=(p*(1-p)+a*(1+a))*y1+la
  dmt = cbind(y1,y1,1)
  dvt = cbind(y1+2*a*y1, y1-2*p*y1, 1)
  for(j in 1:2){
    s[,j] = -1/(2*vt)*dvt[,j]+(y0-mt)/vt*dmt[,j]+(y0-mt)^2/(2*vt^2)*dvt[,j]
  }
  for(t in 1:(n-1)){
    out = out + s[t,]%*%t(s[t,])
  }
  out
}
