
### Functions for estimation of unrestricted model

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