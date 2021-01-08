#==== FGLMR in Big Data ====
library(fda)
library(MASS)
library(mnormt) ## rmnorm: multinorm & rmt : multivariate t distribution
library(mvtnorm) ## Mvnorm: Multivariate Normal Density & Mvt: The Multivariate t Distribution
library(psych) ## the trace of matrix
library(wordspace) ## the rowNorms function
library(CVXR)
#=== functions ========================================

# inverse of the link function 
psi = function(x) exp(x)/(1+exp(x))

# function to generate a(the coefficient of basis)
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
FGLMR.data.generator.bsplines = function(n,nknots,norder,T,domain=c(0,1),aind)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  
  tobs = seq(domain[1],domain[2],length.out = T)
  basismat = eval.basis(tobs, basis)
  
  x = a_fun(n,nbasis,aind) %*% t(basismat)
  
  # betaeval = sqrt(2)*(sin(2*pi*tobs)+cos(2*pi*tobs))
  # betaeval = exp(-32*(tobs-0.5)^2)+2*tobs-1
  betaeval = sin(pi*tobs)
  
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

# The estimator of c using the subsample data selected by the L_opt
Lopt_FGLM = function(N, N_norm, Y, r, r0, lambda, V)
{
  n = dim(N)[1]
  p = dim(N)[2]
  beta = CVXR::Variable(p)
  
  ###############################################
  #    Step 1: subsampel using the uniform      #
  #            sampling probabilities 1/n       #
  ###############################################
  index_uni0= sample(1:n,r0, replace = FALSE)
  N_uni0 = N[index_uni0,]
  Y_uni0 = Y[index_uni0]
  
  # calculate c0(without penalty)
  obj0 <- (-sum(CVXR::logistic(N_uni0[Y_uni0<= 0, ] %*% beta)) - 
             sum(CVXR::logistic(-N_uni0[Y_uni0>0,] %*% beta)))
  prob0 <- CVXR::Problem(CVXR::Maximize(obj0))
  result0 <- solve(prob0)
  c0 <- result0$getValue(beta)
  
  ###############################################
  #    Step 2: calculate the sampling           #
  #            probabilities p_Lopt             #
  ###############################################
  y_prob0 <- 1/(1 + exp(-N %*% c0))
  res = Y - y_prob0
  p_Lopt = abs(res)*N_norm/sum(abs(res)*N_norm)
  
  # subsample using p_Lopt
  index_Lopt= sample(1:n,r, prob =  p_Lopt,replace = TRUE)
  N_Lopt = N[index_Lopt,]
  Y_Lopt = Y[index_Lopt]
  p_s_Lopt = p_Lopt[index_Lopt]
  W_Lopt = 1/(r*p_s_Lopt)
 
  ###############################################
  #    Step 3: estimate using the               #
  #            subsample data                   #
  ###############################################
  obj_Lopt <- (-sum(multiply(W_Lopt[Y_Lopt<= 0],CVXR::logistic(N_Lopt[Y_Lopt<= 0, ] %*% beta))) - 
                 sum(multiply(W_Lopt[Y_Lopt>0],CVXR::logistic(-N_Lopt[Y_Lopt>0,] %*% beta)))-
                 lambda * CVXR::quad_form(beta,V))
  prob_Lopt <- CVXR::Problem(CVXR::Maximize(obj_Lopt))
  result_Lopt <- solve(prob_Lopt)
  c_Lopt <- result_Lopt$getValue(beta)
  
  return(c_Lopt = c_Lopt)
}

Lopt_FGLM_BIC = function(N, N_norm, Y, r, r0, lambda, V)
{
  n = dim(N)[1]
  p = dim(N)[2]
  beta = CVXR::Variable(p)
  
  ###############################################
  #    Step 1: subsampel using the uniform      #
  #            sampling probabilities 1/n       #
  ###############################################
  index_uni0= sample(1:n,r0, replace = FALSE)
  N_uni0 = N[index_uni0,]
  Y_uni0 = Y[index_uni0]
  
  # calculate c0(without penalty)
  obj0 <- (-sum(CVXR::logistic(N_uni0[Y_uni0<= 0, ] %*% beta)) - 
             sum(CVXR::logistic(-N_uni0[Y_uni0>0,] %*% beta)))
  prob0 <- CVXR::Problem(CVXR::Maximize(obj0))
  result0 <- solve(prob0)
  c0 <- result0$getValue(beta)
  
  ###############################################
  #    Step 2: calculate the sampling           #
  #            probabilities p_Lopt             #
  ###############################################
  y_prob0 <- 1/(1 + exp(-N %*% c0))
  res = Y - y_prob0
  p_Lopt = abs(res)*N_norm/sum(abs(res)*N_norm)
  
  # subsample using p_Lopt
  index_Lopt= sample(1:n,r, prob =  p_Lopt,replace = TRUE)
  N_Lopt = N[index_Lopt,]
  Y_Lopt = Y[index_Lopt]
  p_s_Lopt = p_Lopt[index_Lopt]
  W_Lopt = 1/(r*p_s_Lopt)
  
  ###############################################
  #    Step 3:  choosing the optimal lambda     #
  #        based on the optimal subsample data  #
  ###############################################
  nlambda   = length(lambda)
  bic.lopt = array(NA,nlambda)
  for (i in 1:nlambda)
  {
    A = t(N_Lopt)%*%diag(W_Lopt)%*%N_Lopt+lambda[i]*V
    hat_matrix = N_Lopt%*%ginv(A)%*%t(N_Lopt)
    df = tr(hat_matrix)
    
    # calculate beta for each lambda
    obj_bic <- (-sum(multiply(W_Lopt[Y_Lopt<= 0],CVXR::logistic(N_Lopt[Y_Lopt<= 0, ] %*% beta))) - 
                  sum(multiply(W_Lopt[Y_Lopt>0],CVXR::logistic(-N_Lopt[Y_Lopt>0,] %*% beta)))-
                  lambda[i] * CVXR::quad_form(beta,V))
    prob_bic <- CVXR::Problem(CVXR::Maximize(obj_bic))
    result_bic <- solve(prob_bic)
    c_bic <- result_bic$getValue(beta)
    
    # log-likelihood
    l = (-sum(log(1+exp(N_Lopt[Y_Lopt<= 0, ] %*% c_bic))) - 
           sum(log(1+exp(-N_Lopt[Y_Lopt>0,] %*% c_bic))))
    bic.lopt[i] = -2*l/r + log(r)*df
  }
  idx = which(bic.lopt == min(bic.lopt), arr.ind = TRUE)
  lambda.lopt = lambda[idx]
  
  ###############################################
  #    Step 4: estimate using the               #
  #            subsample data                   #
  ###############################################
  obj_Lopt <- (-sum(multiply(W_Lopt[Y_Lopt<= 0],CVXR::logistic(N_Lopt[Y_Lopt<= 0, ] %*% beta))) - 
                 sum(multiply(W_Lopt[Y_Lopt>0],CVXR::logistic(-N_Lopt[Y_Lopt>0,] %*% beta)))-
                 lambda.lopt * CVXR::quad_form(beta,V))
  prob_Lopt <- CVXR::Problem(CVXR::Maximize(obj_Lopt))
  result_Lopt <- CVXR::solve(prob_Lopt)
  c_Lopt <- result_Lopt$getValue(beta)
  
  return(list(c_Lopt = c_Lopt,lambda.lopt=lambda.lopt))
}

# The estimator of c using the subsample data selected by uniform probability
Unisub_FGLM = function(N, Y, r, lambda, V)
{
  n = dim(N)[1]
  index_uni0= sample(1:n,r, replace = TRUE)
  N_uni0 = N[index_uni0,]
  Y_uni0 = Y[index_uni0]
  # W_uni = rep(n/r,r)
  
  # calculate c0
  obj0 <- (-sum(CVXR::logistic(N_uni0[Y_uni0<= 0, ] %*% beta)) - 
             sum(CVXR::logistic(-N_uni0[Y_uni0>0,] %*% beta))-
    (r/n)*lambda * CVXR::quad_form(beta,V))
  prob0 <- CVXR::Problem(CVXR::Maximize(obj0))
  result0 <- CVXR::solve(prob0)
  c0 <- result0$getValue(beta)
  return(c = c0)
}

# sim function
sim_FGLM = function(n,T,domain,lambda,r,r0,sigma,y0,NV)
{
  # generate data
  prob = psi(y0)
  Y = rbinom(n = n, prob = prob, size = 1)
  
  tobs   = seq(domain[1],domain[2],length.out = T)
  
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
  Lopt_result = Lopt_FGLM_BIC(N,N_norm,Y,r,r0,lambda,V)
  lambda.lopt = Lopt_result$lambda.lopt
  c_Lopt = Lopt_result$c_Lopt
  
  # c_Lev = Lscore(N, yc, r, lambda.lopt, V)
  c_uni = Unisub_FGLM(N, Y, r, lambda.lopt, V)
  
  
  # estimate beta(t)
  # beta_true = sqrt(2)*(sin(2*pi*tobs)+cos(2*pi*tobs))
  beta_true = sin(pi*tobs)
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
  
  # Proportions of correct classification (PCC)
  Y_prob_Lopt <- 1/(1 + exp(-N %*% c_Lopt))
  Y_hat_Lopt <- 1*(Y_prob_Lopt>0.5)
  PCC_Lopt = sum(Y==Y_hat_Lopt )/n
  Y_prob_unif <- 1/(1 + exp(-N %*% c_uni))
  Y_hat_unif <- 1*(Y_prob_unif>0.5)
  PCC_unif = sum(Y==Y_hat_unif )/n
  
  # result = c(IMSE_Lopt, IMSE_Lev, IMSE_uni)
  result = list(IMSE = c(IMSE_Lopt,  IMSE_uni),PCC = c(PCC_Lopt,PCC_unif))
  return(result)
  
}
