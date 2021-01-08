
source("functions.R")

# parameters setting
domain = c(0,1)
T = 101
d = 3
sigma = sqrt(0.1)
n = 10^5
r = c(600,800,1000,1200,1400,1600)
r0 = r/2
K = ceiling(1.25*n^(0.25))
lambda = seq(0.1*n^(-3/8),1.5*n^(-3/8),length.out = 10)

# Scenario I
aind = 1 

IMSE_1 = array(0,dim = c(100,2,6))
r = c(600,800,1000,1200,1400,1600)
r0 = r/2
for (i in 1:100){
  for (j in 1:6){
    IMSE_1[i,,j] = sim(n,T,domain,K,d,aind,lambda,r[j],r0[j],sigma)
    print(c(i,j))
  }
}
save.image("result11.RData")

# Scenario II
aind = 2 
IMSE_2 = array(0,dim = c(500,2,6))
r = c(600,800,1000,1200,1400,1600)
r0 = r/2
for (i in 1:500){
  for (j in 1:6){
    IMSE_2[i,,j] = sim(n,T,domain,K,d,aind,lambda,r[j],r0[j],sigma)
    print(c(i,j))
  }
}

save.image("result12.RData")