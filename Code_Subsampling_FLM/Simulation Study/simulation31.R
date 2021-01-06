setwd("C:/Users/liuhua/Desktop/Study In SFU/Big Data/code_for_subsample/simulationIIInew_2")

source("functions.R")

# parameters setting
domain = c(0,1)
T = 101
d = 3
sigma = sqrt(0.1)
n = 5*10^6
r = c(5000,6000,7000,8000,9000,10000)
r0 = r/2
K = ceiling(1.25*n^(0.25))
lambda = seq(0.1*n^(-3/8),1.5*n^(-3/8),length.out = 10)

# Scenario I
aind = 1

# generate x & y0
data = FLMR.data.generator.bsplines(n=n,nknots=66,norder=4,T = T, domain=c(0,1), aind=aind)
X   = data$X
y0   = data$y0
rm(data)
xc = scale(X,center = TRUE,scale = FALSE)
rm(X)
NV = compute.NV(xc,K,d,domain)
rm(xc)

# result = sim(n,T,domain,lambda,r,r0,sigma,y0,NV)

# IMSE result
IMSE_1 = array(0,dim = c(100,2,6))
for (i in 1:100){
  for (j in 1:6){
    IMSE_1[i,,j] = sim(n,T,domain,lambda,r[j],r0[j],sigma,y0,NV)
    print(c(i,j))
  }
}
save.image("result31_0809.RData")

mIMSE_1 = apply(IMSE_1,c(2,3),mean)
matplot(r,t(mIMSE_1),type = "l",xlab = "subsample size",ylab = "IMSE",col = 1:2,lty = 1:2,lwd = 2,main = "Scenario I")
legend("topright",legend= c("Lopt","UNIF"),lty = 1:2,col = 1:2,lwd =rep(2,2),bty = "n")

