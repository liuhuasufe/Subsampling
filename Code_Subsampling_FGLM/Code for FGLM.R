library(CVXR)
p <- dim(N)[2]
beta <- Variable(p)
obj2 <- (-sum(CVXR::logistic(N[Y <= 0, ] %*% beta)) - 
           sum(CVXR::logistic(-N[Y>0,] %*% beta))-
           lambda * CVXR::quad_form(beta,V))
prob <- CVXR::Problem(CVXR::Maximize(obj2))
result <- solve(prob)
beta_res <- result$getValue(beta)



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
beta0 <- result0$getValue(beta)

###############################################
#    Step 2: calculate the sampling           #
#            probabilities p_Lopt             #
###############################################
y_prob0 <- 1/(1 + exp(-N %*% beta_res))
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
  prob_bic <- CVXR::Problem(CVXR::Maximize(obj_Lopt))
  result_bic <- solve(prob_bic)
  beta_bic <- result_bic$getValue(beta)
  
  # log-likelihood
  l = (-sum(log(1+exp(N_Lopt[Y_Lopt<= 0, ] %*% beta_bic))) - 
         sum(log(1+exp(-N_Lopt[Y_Lopt>0,] %*% beta_bic))))
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
               lambda * CVXR::quad_form(beta,V))
prob_Lopt <- CVXR::Problem(CVXR::Maximize(obj_Lopt))
result_Lopt <- solve(prob_Lopt)
beta_Lopt <- result_Lopt$getValue(beta)
y_p_Lopt <- 1/(1 + exp(-N %*% beta_res))
y_hat_Lopt = 1*(y_p_Lopt>0.5)
sum(Y==y_hat_Lopt)/n


prop.table(table(Y,y_hat_Lopt))






y_prob <- 1/(1 + exp(-N %*% beta_res))

res = Y - y_prob
p_Lopt = abs(res)*N_norm/sum(abs(res)*N_norm)

y_hat = 1*(y_prob>0.5)

sum(Y==y_hat)/n

table(Y,y_hat)
pp = as.matrix(prop.table(table(Y,y_hat)))



beta_true = sin(pi*tobs)
# beta_true = exp(-32*(tobs-0.5)^2)+2*tobs-1
beta_full = basismat%*%beta_Lopt
IMSE_full1 = sqrt(mean((beta_full-beta_true)^2))
# plot(tobs,beta_true)
plot(tobs,beta_Lopt)
lines(tobs,beta_true)
lines(tobs,beta_uni)
# lines(tobs,beta_full)
