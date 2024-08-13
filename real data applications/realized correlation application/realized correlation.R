

###############################################################################
## This file contains the code used in the paper "Pseudo-variance 
## quasi-maximum likelihood estimation of semi-parametric time series models" 
## by Mirko Armillotta and Paolo Gorgi.
##
## This file contains the code for the application on 
## realized correlation data (Section 6).
##
## Note: R version 4.3.2 used throughout.
##
###############################################################################

# packages
library(tseries)

source("beta QML.R")
source("beta QML uncons.R")

#data
current_dir <- getwd()
setwd(current_dir)
load("rk.dta") 

## Select data from the lists using the index "m"
# The index m can be selected from 1 to 21. From 1 to 10 there are the bivariate series and realized covariances as in the paper Gorgi, P., Hansen, P. R., Janus, P., and Koopman, S. J. (2019). Realized Wishart-GARCH: A Score-driven Multi-Asset Volatility Model. Journal of Financial Econometrics, 17, 1-32. (same order). 
# From 11 to 20 there are the 5x5 series (same order as in the paper). For m=21 there is the 15x15 series.
m <- 4

rv <- rk[[m]]         
rc <- rv[1,2,]/(sqrt(rv[1,1,])*sqrt(rv[2,2,]))
data <- (rc+1)/2


y = data


par(mfrow=c(1,1),mar=c(4.1,4.1,1.1,2.1))
time <- seq(from=2001,to=2010,length.out = length(data))
plot(time, 2*y-1, type="l", ylab = "Realized correlation", xlab="Time")

min(y)
max(y)
my = mean(y)
my
vy = var(y)
vy
n=length(y)


# compute starting values by ARMA
arma = arma(y-my, order=c(1,1), include.intercept = F)
b=-arma$coef[2]
a=arma$coef[1]-b
d=my*(1-a-b)
psi = c(d,a,b)
psi
ini0=log(psi/(1-psi))

mu  = vector()
mu[1] = my
for (t in 2:n){
  mu[t] = d + a*y[t-1] + b*mu[t-1]
}

#moment estimator for phi
phi=(mean((y-mu)^2/(mu-mu^2)))^(-1)-1
ini = c(ini0, log(phi))



## Estimate fully restricted PVQMLE

test2 = optim(par=ini, fn=QML, method = "BFGS", y=y,
              control = list(ndeps=rep(1e-5,4), reltol=1e-30, maxit=1000))
QML_est = c(exp(test2$par[1:3])/(1+exp(test2$par[1:3])), exp(test2$par[4]))

hessq <- optimHess(par = QML_est, fn = QML_hess, y=y)
hesq = solve(hessq)
out = QML_out(theta = QML_est, y=y)
sigma = hesq%*%out%*%hesq
QML_se=sqrt(diag(sigma))
QML_t=QML_est/QML_se


## Estimate unrestricted PVQMLE

ini2=c(ini,ini0)

test5 = optim(par=ini2, fn=QMLu, method = "BFGS", y=y,
              control = list(ndeps=rep(1e-5,7), reltol=1e-30, maxit=1000))

QMLu_est = c(exp(test5$par[1:3])/(1+exp(test5$par[1:3])), exp(test5$par[4]),
             exp(test5$par[5:7])/(1+exp(test5$par[5:7])))

hessqu <- optimHess(par = QMLu_est, fn = QMLu_hess, y=y)
hesqu = solve(hessqu)
outu = QMLu_out(theta = QMLu_est, y=y)
sigmau = hesqu%*%outu%*%hesqu

QMLu_se=sqrt(diag(sigmau))
QMLu_t=QMLu_est/QMLu_se


## results


round(QML_est,3)
round(QML_se,3)

round(QMLu_est,3)
round(QMLu_se,3)



## Wald tests


# test for equality of intercepts: omega_1=omega_2

gd = c(QMLu_est[1]-QMLu_est[5])
Gd = cbind(1,0,0,0,-1,0,0)
Wd = t(gd) %*% solve(Gd%*%sigmau%*%t(Gd)) %*% gd
round(1-pchisq(Wd,1),2)


# test for alpha_1=alpha_2

ga = c(QMLu_est[2]-QMLu_est[6])
Ga = cbind(0,1,0,0,0,-1,0)
Wa = t(ga) %*% solve(Ga%*%sigmau%*%t(Ga)) %*% ga
round(1-pchisq(Wa,1),2)


# test for beta_1=beta_2

gb = c(QMLu_est[3]-QMLu_est[7])
Gb = cbind(0,0,1,0,0,0,-1)
Wb = t(gb) %*% solve(Gb%*%sigmau%*%t(Gb)) %*% gb
round(1-pchisq(Wb,1),2)

# Joint test

g = c(QMLu_est[1]-QMLu_est[5] , QMLu_est[2]-QMLu_est[6],
      QMLu_est[3]-QMLu_est[7])
G = cbind(diag(3),rep(0,3),diag(-1,3))
W = t(g) %*% solve(G%*%sigmau%*%t(G)) %*% g

round(1-pchisq(W,3),2)


