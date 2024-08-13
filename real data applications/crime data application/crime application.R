


###############################################################################
## This file contains the code used in the paper "Pseudo-variance 
## quasi-maximum likelihood estimation of semi-parametric time series models" 
## by Mirko Armillotta and Paolo Gorgi.
##
## This file contains the code for the application on crime data (Section 6).
##
## Note: R version 4.3.2 used throughout.
##
###############################################################################

# data
current_dir <- getwd()
setwd(current_dir)
crime <- read.delim("Dataset.txt", header=FALSE)
y = t(crime)


par(mfrow=c(3,1),mar=c(4.1,4.1,1.1,2.1))
time <- seq(from=1995,to=2015,length.out = length(y))
plot(time, y, type="l", ylab="Crime counts", xlab="Time")
acf(y)
pacf(y)

my = mean(y)
vy = var(y)
# overdispersion
my
vy
n=length(y)
# zero-inflation
sum(y==0)/n
exp(-my)

source("uncons model.R")
source("Poisson INAR qml.R")
source("Poisson error and thin INAR qml.R")
source("Poisson error and NB thin INAR qml.R")
source("Poisson error and BiNB thin INAR qml.R")

### Estimation of unrestricted model

# starting values (method of moments)
a=cov(y[-1],y[-n])/vy
b=a*(1-a)
w1 = (1-a)*my
w2 = (1-a^2)*vy-b*my 
theta=c(a, w1, w2, b)
ini2=c(log(theta[1]/(1-theta[1])),log(theta[2]), log(theta[3]),log(theta[4]))

# estimation
est2 <- optim(par = ini2, fn = QML2, method = "BFGS", hessian=F, y=y,
              control = list(ndeps=rep(1e-5,4), reltol=1e-12, maxit=1000))

alpha = exp(est2$par[1])/(1+exp(est2$par[1]))
omega1 = exp(est2$par[2])
omega2 = exp(est2$par[3])
beta = exp(est2$par[4])
QML2_est = c(alpha, omega1, omega2, beta)

hess <- optimHess(par = QML2_est, fn = QML2_hess, y=y)
hes = solve(hess)
out = out_QML2(y=y, theta = QML2_est)
I = hes%*%out%*%hes
QML2_se = sqrt(diag(I))
QML2_t = QML2_est/QML2_se

# results estimation unrestricted PVQMLE

round(QML2_est,3)
round(QML2_se,3)
round(QML2_t,3)

# Wald test for Poisson thinning

g = beta - alpha
G = c(-1,0,0,1)
Wp2 = g^2/(t(G)%*%I%*%G)
1-pchisq(Wp2,1)


# Wald test for Binomial thinning

g = beta - alpha*(1-alpha)
G = c(-1+2*alpha,0,0,1)
Wb2 = g^2/(t(G)%*%I%*%G)
1-pchisq(Wb2,1)


# Wald test for Geometric thinning

g = beta - alpha*(1+alpha)
G = c(-1-2*alpha,0,0,1)
Wg2 = g^2/(t(G)%*%I%*%G)
1-pchisq(Wg2,1)



# Wald test for equidispersed error

g = omega2 - omega1
G = c(0,-1,1,0)
Wep = g^2/(t(G)%*%I%*%G)
1-pchisq(Wep,1)


################ models estimation

### Poisson INAR PVQMLE

# starting values (method of moments)
a=cov(y[-1],y[-n])/vy
w = (1-a)*my
theta=c(a, w)
inip=c(log(theta[1]/(1-theta[1])),log(theta[2]))

# estimation
estp <- optim(par = inip, fn = QMLP, method = "BFGS", hessian=F, y=y,
              control = list(ndeps=rep(1e-5,2), reltol=1e-12, maxit=1000))

alpha = exp(estp$par[1])/(1+exp(estp$par[1]))
omega = exp(estp$par[2])
QMLP_est = c(alpha, omega)

hess <- optimHess(par = QMLP_est, fn = QMLP_hess, y=y)
hes = solve(hess)
out = out_QMLP(y=y, theta = QMLP_est)
I = hes%*%out%*%hes
QMLP_se = sqrt(diag(I))
QMLP_t = QMLP_est/QMLP_se

# results 
round(QMLP_est,3)
round(QMLP_se,3)
round(QMLP_t,3)

### Poisson error and thinning INAR PVQMLE

# starting values (method of moments)
a=cov(y[-1],y[-n])/vy
w = (1-a)*my
theta=c(a, w)
inip=c(log(theta[1]/(1-theta[1])),log(theta[2]))

# estimation
estpp <- optim(par = inip, fn = QMLPP, method = "BFGS", hessian=F, y=y,
               control = list(ndeps=rep(1e-5,2), reltol=1e-12, maxit=1000))

alpha = exp(estpp$par[1])/(1+exp(estpp$par[1]))
omega = exp(estpp$par[2])
QMLPP_est = c(alpha, omega)

hess <- optimHess(par = QMLPP_est, fn = QMLPP_hess, y=y)
hes = solve(hess)
out = out_QMLPP(y=y, theta = QMLPP_est)
I = hes%*%out%*%hes
QMLPP_se = sqrt(diag(I))
QMLPP_t = QMLPP_est/QMLPP_se

# results 
round(QMLPP_est,3)
round(QMLPP_se,3)
round(QMLPP_t,3)


### Poisson error Geometric thinning INAR PVQMLE

# starting values (method of moments)
rr=1
a=cov(y[-1],y[-n])/vy
w = (1-a)*my
theta=c(a, w)
inip=c(log(theta[1]/(1-theta[1])),log(theta[2]))

# estimation 
estpg <- optim(par = inip, fn = QMLPNB, method = "BFGS", hessian=F, y=y, r=rr,
               control = list(ndeps=rep(1e-5,2), reltol=1e-12, maxit=1000))

alpha = exp(estpg$par[1])/(1+exp(estpg$par[1]))
omega = exp(estpg$par[2])
QMLPG_est = c(alpha, omega)
QMLPG_est

hess <- optimHess(par = QMLPG_est, fn = QMLPNB_hess, y=y, r=rr)
hes = solve(hess)
out = out_QMLPNB(y=y, theta = QMLPG_est, r=rr)
I = hes%*%out%*%hes
QMLPG_se = sqrt(diag(I))
QMLPG_t = QMLPG_est/QMLPG_se

# results
round(QMLPG_est,3)
round(QMLPG_se,3)
round(QMLPG_t,3)


### Poisson error Binomial-Negative Binomial (BiNB) thinning INAR PVQMLE

# starting values (method of moments)
a=cov(y[-1],y[-n])/vy
w = (1-a)*my
theta=c(a/2, a/2, w)
inipbinb=c(log(theta[1]/(1-theta[1])), log(theta[2]/(1-theta[2])),log(theta[3]))

# estimation 
estpbinb <- optim(par = inipbinb, fn = QMLPBiNB, method = "BFGS", hessian=F, y=y,
               control = list(ndeps=rep(1e-5,3), reltol=1e-12, maxit=1000))

mu = exp(estpbinb$par[1])/(1+exp(estpbinb$par[1]))
pi = exp(estpbinb$par[2])/(1+exp(estpbinb$par[2]))
omega = exp(estpbinb$par[3])
QMLPBiNB_est = c(mu, pi, omega)

hess <- optimHess(par = QMLPBiNB_est, fn = QMLPBiNB_hess, y=y)
hes = solve(hess)
out = out_QMLPBiNB(y=y, theta = QMLPBiNB_est)
I = hes%*%out%*%hes
QMLPBiNB_se = sqrt(diag(I))
QMLPBiNB_t = QMLPBiNB_est/QMLPBiNB_se

# results
round(QMLPBiNB_est,3)
round(QMLPBiNB_se,3)
round(QMLPBiNB_t,3)


# residuals unrestricted model

mt_QML2 = alpha*y[-n] + omega1
vt_QML2 = beta*y[-n] + omega2
et_QML2 = (y[-1] - mt_QML2)/sqrt(vt_QML2)


par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
time <- seq(from=1995,to=2015,length.out = length(y))
plot(time, y, type="l", ylab="Crime counts", xlab="Time")
acf(et_QML2)

