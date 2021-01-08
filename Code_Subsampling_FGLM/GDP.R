# parameters setting
domain = c(0,1)
T = 101
d = 3
sigma = sqrt(0.1)
n = 10^5
r = 5000
r0 = r/2
# K = 15
K = ceiling(1.25*n^(0.25))
lambda = seq(0.1*n^(-3/8),1.5*n^(-3/8),length.out = 10)
lambda = lambda[10]

# Scenario I
aind = 3

# generate x & y0
data = FLMR.data.generator.bsplines(n=n,nknots=66,norder=4,T = T, domain=c(0,1), aind=aind)
X   = data$X

y0   = data$y0
psi = function(x) exp(x)/(1+exp(x))
prob = psi(y0)
Y = rbinom(n = n, prob = prob, size = 1)
# Y[Y==0] = -1

tobs   = seq(domain[1],domain[2],length.out = T)

# new design matrix and smoothness matrix
# NV = compute.NV(xc,K,d,domain)
NV = compute.NV(X,K,d,domain)
N = NV$N
V = NV$V
N_norm = NV$N_norm
basismat = NV$basismat

#==== define PQL ====
PQL = function(y,N,V,lambda){
  q = dim(N)[2]
  
  # quasi-likelihood equation
  f = function(c) colSums(as.vector(y-psi(N%*%c))*N)-lambda*V%*%c
  
  # Jacobian of f
  J = function(c) -(t(N)%*%diag(as.vector(dpsi(N%*%c)))%*%N+lambda*V)
  
  # initial value
  cstart = seq(1,1,length=q)/sqrt(q)
  
  # solve the nonlinear equation f 
  chat = nleqslv(cstart,f,J, method="Newton",control=list(trace=1,allowSingular=TRUE))
  
  # obtain the result c 
  cpql = chat$x
  
  return(cpql)
}
library(nleqslv)
result = PQL(Y,N,V,lambda)

# beta_true = exp(-32*(tobs-0.5)^2)+2*tobs-1
beta_true = sin(pi*tobs)+1
beta_full = basismat%*%result

par(mfrow = c(1,2))
plot(tobs,beta_full)
# lines(tobs,beta_true)
plot(tobs,beta_true)

#=== package-rms ====
library(rms)
# ?lrm
lrm1 = lrm(Y~N,x=TRUE, y=TRUE)
lrm1$coefficients
beta_full1 = basismat%*%lrm1$coefficients[-1]

beta_true = tobs^2

IMSE_full = sqrt(mean((beta_full1-beta_true)^2))
par(mfrow = c(1,2))
plot(tobs,beta_full1)
# lines(tobs,beta_true)
plot(tobs,beta_true)
par(mfrow = c(1,1))


# choosing lambda
pen = pentrace(lrm1, 2*lambda,V,
         # method=c('grid','optimize'),
         which='bic')
pen = pentrace(lrm1, 2*lambda, V,
               method ='optimize')

lrm1n = update(lrm1, penalty=pen$penalty)
lrm1nn = lrm(Y~N,penalty = pen$penalty, penalty.matrix = V,x=TRUE, y=TRUE)

beta_full = basismat%*%lrm1nn$coefficients[-1]
IMSE_full1 = sqrt(mean((beta_full-beta_true)^2))
plot(tobs,beta_full)
lines(tobs,beta_true)
lines(tobs,beta_full1)