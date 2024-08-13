
###############################################################################
## This file contains the code used in the paper "Pseudo-variance 
## quasi-maximum likelihood estimation of semi-parametric time series models" 
## by Mirko Armillotta and Paolo Gorgi.
##
## This file contains the code to replicate Figure 1-2.
##
## Note: R version 4.3.2 used throughout.
##
###############################################################################



# Binomial thinning function
thinning<-function(alfa,ni)
{
  y<-rbinom(ni,1,alfa)
  thin<-sum(y)
  return(thin)
}

# Generate data from Poisson INAR(1) model 
rpinar<-function(n,alfa, lamda)
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

### Function for fully restricted PVQMLE
sQML = function(y,theta){
  la=theta[1]
  a=theta[2]
  n=length(y)
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a*(1-a)*y1+la
  one = rep(1,n-1)
  dm = y1
  Dm = rbind(one,t(dm))
  dv = y1-2*a*y1
  Dv = rbind(one,t(dv))
  e = y0 - mt
  h =  Dm %*% diag(1/vt) %*% t(Dm) + Dv %*% diag(1/(2*vt^2)) %*% t(Dv)
  h = 1/(n-1) * h
  i = Dm %*% diag(1/vt) %*% t(Dm) + 
    (Dm + Dv) %*% diag(e^3/(2*vt^3)) %*% t(Dv + Dm) +
    Dv %*% diag(e^4/vt^2-1) %*% diag(1/(4*vt^2)) %*% t(Dv)
  i = 1/(n-1) * i
  h1 = solve(h)
  sig = h1 %*% i %*% h1
  return(list(h = h, i = i, h1 = h1, sig = sig))
}

### Function for unrestricted PVQMLE
sQML2 = function(y,theta){
  la=theta[1]
  v=theta[2]
  a= theta[3]
  a1=theta[4]
  n=length(y)
  y1=y[-n]
  y0=y[-1]
  mt=a*y1+la
  vt=a1*y1+v
  one = rep(1,n-1)
  dm = y1
  Dm = rbind(one,t(dm))
  h =  Dm %*% diag(1/vt) %*% t(Dm)
  h = 1/(n-1) * h
  sig = solve(h)
  return(list(h = h, sig = sig))
}


### Simulation for a

n = 10000


la = v = 3

av = seq(0.1, 0.95, by=0.05) 

cQML = list()
cQML2 = list()
vqml = vector()
vqml2 = vector()
eig = matrix(NA, length(av), 2)


for(j in 1:length(av)){
  
  a = av[j]
  a1 = a*(1-a)
  
  set.seed(1234)
  y = rpinar(n,a,la)
  
  theta = c(la, a)
  thetau = c(la, v, a, a1)
  
  c = sQML(y, theta)
  unc = sQML2(y, thetau)
  
  cQML = c$sig
  cQML2 = unc$sig
  
  eig[j,] = eigen(unc$sig - c$sig)$value
  
  vqml[j] = c$sig[2,2]
  vqml2[j] = unc$sig[2,2]
  
  print(j)
}


### Simulations for omega

a2 = 0.85
a12 = a2*(1-a2)
lav = seq(0.5, 15, length.out = 30) 
wcQML = list()
wcQML2 = list()
wvqml = vector()
wvqml2 = vector()
weig = matrix(NA, length(lav), 2)

for(j in 1:length(lav)){
  
  la2 = lav[j]
  v2 = la2
  
  set.seed(1234)
  y = rpinar(n,a2,la2)
  
  theta2 = c(la2, a2)
  thetau2 = c(la2, v2, a2, a12)
  
  wc = sQML(y, theta2)
  wunc = sQML2(y, thetau2)
  
  wcQML = wc$sig
  wcQML2 = wunc$sig
  
  weig[j,] = eigen(wunc$sig - wc$sig)$value
  
  wvqml[j] = wc$sig[1,1]
  wvqml2[j] = wunc$sig[1,1]
  
  print(j)
}

## Figure 2

par(mfrow = c(1,2))
par(mar=c(4,4.1,2.1,1.1))
plot(av, log10(vqml2/vqml), pch = 1, ylim = c(-0.1,1.2), xlab = "a", ylab="log10 Variance ratio") #, xlim = c(0.2,0.9), ylim=c(0,2))
abline(h=0, type="l", lty=2, col=2)
plot(lav, log10(wvqml2/wvqml), pch = 1, ylim = c(-0.1,1.2), xlab=expression(paste(omega)), ylab="")
abline(h=0, type="l", lty=2, col=2)


## Contour plots Figure 1

avd = seq(0.1, 0.95, length.out = 10) 
lavd = seq(0.5, 15, length.out = 10) 

zz=length(avd)
EE = matrix(NA, zz,zz)


z_fun <- function(avd, lavd){
  
  for(i in 1:zz){
    
    a = avd[i]
    a1 = a*(1-a)
    
    for(j in 1:zz){
      
      la = lavd[j]
      v = la
      
      set.seed(1234)
      y = rpinar(n,a,la)
      
      theta = c(la, a)
      thetau = c(la, v, a, a1)
      
      c = sQML(y, theta)
      unc = sQML2(y, thetau)
      cQML = c$sig
      cQML2 = unc$sig
      Vqml = c$sig[2,2]
      Vqml2 = unc$sig[2,2]
      EE[i,j] = Vqml2/Vqml
    }
  }
  return(EE)
}

z_values = z_fun(avd, lavd)



EE2 = matrix(NA, zz,zz)

z_fun_la <- function(avd, lavd){
  
  for(i in 1:zz){
    
    a = avd[i]
    a1 = a*(1-a)
    
    for(j in 1:zz){
      
      la = lavd[j]
      v = la
      
      set.seed(1234)
      y = rpinar(n,a,la)
      
      theta = c(la, a)
      thetau = c(la, v, a, a1)
      
      c = sQML(y, theta)
      unc = sQML2(y, thetau)
      cQML = c$sig
      cQML2 = unc$sig
      Wqml = c$sig[1,1]
      Wqml2 = unc$sig[1,1]
      EE2[i,j] = Wqml2/Wqml
    }
  }
  return(EE2)
}

z_values_la = z_fun_la(avd, lavd)


lz = log10(t(z_values))
lz_la = log10(t(z_values_la))

cols <- hcl.colors(20,"Blue-Red")
par(mfrow=c(1,2))
par(mar=c(4,4.1,2.1,1.1))
par(mar=c(4,3.7,1.1,1.1))
filled.contour(lavd, avd, lz,
               nlevels = 7,  color.palette = function(n) rainbow(n,rev = T,start=0.85,end=0.25),
               xlab=expression(paste(omega)), ylab="a",
               main=expression(paste("log10"~"Variance"~"ratio"~"for"~a)),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(lavd, avd, lz, nlevels = 40, add = TRUE, lwd = 1)
               }
)
filled.contour(lavd, avd, lz_la,
               nlevels = 7,  color.palette = function(n) rainbow(n,rev = T,start=0.85,end=0.25),
               xlab=expression(paste(omega)),
               main=expression(paste("log10"~"Variance"~"ratio"~"for"~omega)),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(lavd, avd, lz_la, nlevels = 40, add = TRUE, lwd = 1)
               }
)

