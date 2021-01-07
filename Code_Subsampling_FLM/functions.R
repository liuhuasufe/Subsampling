#==== FLMR in Big Data ====
library(fda)
library(MASS)
library(mnormt) ## rmnorm: multinorm & rmt : multivariate t distribution
library(mvtnorm) ## Mvnorm: Multivariate Normal Density & Mvt: The Multivariate t Distribution
library(psych) ## the trace of matrix
library(wordspace) ## the rowNorms function
#=== functions ========================================

# function to generate a(the coefficient of basis)
# a_fun = function(n,nbasis, ii)
# {
#   sigma = 2*0.5^(toeplitz(0:(n-1)))
#   if(ii == 1) {a = t(rmnorm(nbasis,mean = rep(0,n),sigma))}
#   else if(ii == 2){a = t(rmt(nbasis,mean = rep(0,n),sigma,df=1))}
#   else if(ii == 3){a = t(rmt(nbasis,mean = rep(0,n),sigma,df=3))}
#   else{print("model does not exit")}
#   return(a)
# }

a_fun = function(n,nbasis, ii)
{
  #sigma = 2*0.5^(toeplitz(0:(n-1)))
  if(ii == 1) {a = matrix(rnorm(nbasis*n,0,1),n,nbasis)}
  # else if(ii == 2){a = matrix(rt(nbasis*n,df=1),n,nbasis)}
  # else if(ii == 3){a = matrix(rt(nbasis*n,df=3),n,nbasis)}
  else if(ii == 2) {a = matrix(rmt(nbasis*n,mean = rep(0,1),1,df=3),n,nbasis)}
  else if(ii == 3) {a = matrix(rmt(nbasis*n,mean = rep(0,1),1,df=2),n,nbasis)}
  # else if(ii == 2) {a = matrix(rmvt(nbasis*n,sigma = diag(1),df=3),n,nbasis)}
  # else if(ii == 3) {a = matrix(rmvt(nbasis*n,sigma = diag(1),df=1),n,nbasis)}
  else{print("model does not exit")}
  return(a)
}

# function of generate data
FLMR.data.generator.bsplines = function(n,nknots,norder,T,domain=c(0,1),aind)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)

  tobs = seq(domain[1],domain[2],length.out = T)
  basismat = eval.basis(tobs, basis)

  x = a_fun(n,nbasis,aind) %*% t(basismat)

  # betaeval = sqrt(2)*(sin(2*pi*tobs)+cos(2*pi*tobs))
  betaeval = exp(-32*(tobs-0.5)^2)+2*tobs-1

  # y0 the signals
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = h/3*x%*%diag(cef)%*%betaeval
  # eps0= sd(y0)
  # y   = y0 + rnorm(n,mean = 0, sd = sigma)

  return(list(X=x,y0=y0))
}

# The function to get the design matrix N
compute.NV = function(X, K, d, domain)
{
  n = dim(X)[1]
  norder   = d+1 # order = degree +1
  nknots   = K+2 # the number of all knots = the number of inner knots + 2
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2 # the number of basis = number of inner knots + order
  basis    = create.bspline.basis(knots,nbasis,norder)
  T   = dim(X)[2]
  tobs = seq(domain[1], domain[2], length.out = T)
  basismat = eval.basis(tobs, basis) # T x (K+d+1)
  V = eval.penalty(basis,int2Lfd(2))

  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  u   = h/3*X%*%diag(cef)%*%basismat # n x nbasis
  u_norm = rowNorms(u, "euclidean")
  return(list(N=u,V = V,basismat = basismat,N_norm = u_norm))
}

# BIC
BIC.Lopt = function(Y, N, V, n, lambda,W)
{
  A = t(N)%*%W%*%N+lambda*V
  B = t(N)%*%W%*%Y
  c = ginv(A)%*%B
  Yhat = N%*%c
  hat = N%*%ginv(A)%*%t(N)
  df  = tr(hat)
  SSE = t(Y - Yhat)%*%(Y - Yhat)
  BIC.Lopt = n*log(SSE/n) + log(n)*df
  return(BIC.Lopt = BIC.Lopt)
}

# The estimator of c using the subsample data selected by the L_opt
Lopt = function(N, N_norm, yc, r, r0, lambda, V )
{
  n = dim(N)[1]

  ###############################################
  #    Step 1: subsampel using the uniform      #
  #            sampling probabilities 1/n       #
  ###############################################
  index_uni0= sample(1:n,r0, replace = FALSE)
  N_uni0 = N[index_uni0,]
  y_uni0 = yc[index_uni0]

  # calculate c0
  A0 = t(N_uni0)%*%N_uni0
  B0 = t(N_uni0)%*%y_uni0
  # c0 = ginv(t(N_uni0)%*%N_uni0+lambda*V)%*%t(N_uni0)%*%y_uni0
  c0 = ginv(A0)%*%B0

  ###############################################
  #    Step 2: calculate the sampling           #
  #            probabilities p_Lopt             #
  ###############################################
  # N_norm = apply(N,1,norm,"2")
  # N_norm = rowNorms(N, "euclidean")
  res = yc - N%*%c0
  p_Lopt = abs(res)*N_norm/sum(abs(res)*N_norm)

  # subsample using p_Lopt
  index_Lopt= sample(1:n,r, prob =  p_Lopt,replace = TRUE)
  N_Lopt = N[index_Lopt,]
  y_Lopt = yc[index_Lopt]
  p_s_Lopt = p_Lopt[index_Lopt]
  W_Lopt = diag(1/(r*p_s_Lopt))

  ###############################################
  #    Step 3: estimate using the               #
  #            subsample data                   #
  ###############################################

  # estimate using the subsample data
  # c_Lopt = ginv(t(N_Lopt)%*%W_Lopt%*%N_Lopt+lambda*V)%*%t(N_Lopt)%*%W_Lopt%*%y_Lopt
  A_Lopt = t(N_Lopt)%*%W_Lopt%*%N_Lopt+lambda*V
  B_Lopt = t(N_Lopt)%*%W_Lopt%*%y_Lopt
  c_Lopt = ginv(A_Lopt)%*%B_Lopt
  return(c = c_Lopt)
}

Lopt_BIC = function(N, N_norm, yc, r, r0, lambda, V )
{
  n = dim(N)[1]

  ###############################################
  #    Step 1: subsampel using the uniform      #
  #            sampling probabilities 1/n       #
  ###############################################
  index_uni0= sample(1:n,r0, replace = FALSE)
  N_uni0 = N[index_uni0,]
  y_uni0 = yc[index_uni0]

  # calculate c0
  A0 = t(N_uni0)%*%N_uni0
  B0 = t(N_uni0)%*%y_uni0
  # c0 = ginv(t(N_uni0)%*%N_uni0+lambda*V)%*%t(N_uni0)%*%y_uni0
  c0 = ginv(A0)%*%B0

  ###############################################
  #    Step 2: calculate the sampling           #
  #            probabilities p_Lopt             #
  ###############################################
  # N_norm = apply(N,1,norm,"2")
  # N_norm = rowNorms(N, "euclidean")
  res = yc - N%*%c0
  p_Lopt = abs(res)*N_norm/sum(abs(res)*N_norm)

  # subsample using p_Lopt
  index_Lopt= sample(1:n,r, prob =  p_Lopt,replace = TRUE)
  N_Lopt = N[index_Lopt,]
  y_Lopt = yc[index_Lopt]
  p_s_Lopt = p_Lopt[index_Lopt]
  W_Lopt = diag(1/(r*p_s_Lopt))


  ###############################################
  #    Step 3:  choosing the optimal lambda     #
  #        based on the optimal subsample data  #
  ###############################################
  nlambda   = length(lambda)
  bic.lopt = array(NA,nlambda)
  for (i in 1:nlambda)
  {
    bic.lopt[i] = BIC.Lopt(y_Lopt, N_Lopt, V, r, lambda[i], W_Lopt)
  }

  idx = which(bic.lopt == min(bic.lopt), arr.ind = TRUE)
  lambda.lopt = lambda[idx]

  ###############################################
  #    Step 4: estimate using the               #
  #            subsample data                   #
  ###############################################
  # estimate using the subsample data
  A_Lopt = t(N_Lopt)%*%W_Lopt%*%N_Lopt+lambda.lopt*V
  B_Lopt = t(N_Lopt)%*%W_Lopt%*%y_Lopt
  c_Lopt = ginv(A_Lopt)%*%B_Lopt

  return(list(c_Lopt = c_Lopt,lambda.lopt=lambda.lopt))
}

# The estimator of c using the subsample data selected by the leverage score
Lscore = function(N, yc, r, lambda, V)
{
  n = dim(N)[1]

  ###############################################
  #    Step 1: calculate the leverage score     #
  ###############################################
  lscore = apply(svd(N)$u,1,norm,"2")^2/norm(svd(N)$u,"F")^2

  # subsample using the sampling probabilities leverage score
  index_LEV= sample(1:n,r,prob = lscore, replace = TRUE)
  N_LEV = N[index_LEV,]
  y_LEV = yc[index_LEV]
  p_s_LEV = lscore[index_LEV]

  # calculate c_LEV
  W_LEV = diag(1/(r*p_s_LEV))
  A_LEV = t(N_LEV)%*%W_LEV%*%N_LEV+lambda*V
  B_LEV = t(N_LEV)%*%W_LEV%*%y_LEV
  c_LEV = ginv(A_LEV)%*%B_LEV
  return(c = c_LEV)
}

# The estimator of c using the subsample data selected by uniform probability
Unisub = function(N, yc, r, lambda, V)
{
  n = dim(N)[1]
  index_uni0= sample(1:n,r, replace = TRUE)
  N_uni0 = N[index_uni0,]
  y_uni0 = yc[index_uni0]
  W_uni = diag(n/r,r)

  # calculate c0
  A0 = t(N_uni0)%*%W_uni%*%N_uni0+lambda*V
  B0 = t(N_uni0)%*%W_uni%*%y_uni0
  c0 = ginv(A0)%*%B0
  return(c = c0)
}

# calculate IMSE
IMSE = function(x,y,domain)
{
  T = dim(x)[1]
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  Imse = h/3*cef%*%(x-y)^2
  return(Imse)
}

# sim function
sim = function(n,T,domain,lambda,r,r0,sigma,y0,NV)
{
  # generate data
  Y   = y0 + rnorm(n,mean = 0, sd = sigma)
  tobs   = seq(domain[1],domain[2],length.out = T)
  # matplot(tobs,t(X[sample(1:n,10,replace = FALSE),]),xlab="t",ylab="X(t)",type="l",main= expression(t[1](0,1)))

  # # centered the data
  # Xmean = apply(X,2,mean)
  # Xmean = colMeans(X)
  # xc = X - matrix(rep(Xmean,n),nrow=n,byrow= T)
  yc = Y - mean(Y)
  rm(Y)

  # new design matrix and smoothness matrix
  # NV = compute.NV(xc,K,d,domain)
  N = NV$N
  V = NV$V
  N_norm = NV$N_norm
  basismat = NV$basismat

  # estimate c
  # c_full = ginv(t(N)%*%N+lambda*V)%*%t(N)%*%yc
  # c_Lopt = Lopt(N,yc,r,r0,K,d,domain,lambda,V)
  # c_Lev = Lscore(N, yc, r, lambda, V)
  # c_uni = Unisub(N, yc, r, lambda, V)

  # BIC
  Lopt_result = Lopt_BIC(N,N_norm,yc,r,r0,lambda,V)
  lambda.lopt = Lopt_result$lambda.lopt
  c_Lopt = Lopt_result$c_Lopt

  # c_Lev = Lscore(N, yc, r, lambda.lopt, V)
  c_uni = Unisub(N, yc, r, lambda.lopt, V)


  # estimate beta(t)
  # beta_true = sqrt(2)*(sin(2*pi*tobs)+cos(2*pi*tobs))
  beta_true = exp(-32*(tobs-0.5)^2)+2*tobs-1
  # beta_full = basismat%*%c_full
  beta_Lopt = basismat%*%c_Lopt
  # beta_Lev = basismat%*%c_Lev
  beta_uni = basismat%*%c_uni

  # IMSE
  # IMSE_Lopt = IMSE(beta_Lopt,beta_true,domain = c(0,1))
  # IMSE_Lev = IMSE(beta_Lev,beta_true,domain = c(0,1))
  # IMSE_uni = IMSE(beta_uni,beta_true,domain = c(0,1))

  IMSE_Lopt = sqrt(mean((beta_Lopt-beta_true)^2))
  # IMSE_Lev = sqrt(mean((beta_Lev-beta_true)^2))
  IMSE_uni = sqrt(mean((beta_uni-beta_true)^2))
  # IMSE_full = sqrt(mean((beta_full-beta_true)^2))

  # result = c(IMSE_Lopt, IMSE_Lev, IMSE_uni)
  result = c(IMSE_Lopt,  IMSE_uni)
  return(result)

}
