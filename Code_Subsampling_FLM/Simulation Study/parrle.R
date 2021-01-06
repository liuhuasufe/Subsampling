x1 = matrix(rnorm((5*10^5)*500, 0,1),5*10^5,500)
a = Sys.time()
n1 = rowNorms(x1,method="euclidean")
m1 = n1*n1
a = Sys.time()-a

x2 = matrix(rnorm((5*10^5)*500, 0,1),5*10^5,500)
n2 = rowNorms(x2,method="euclidean")


x3 = matrix(rnorm((5*10^5)*500, 0,1),5*10^5,500)
n3 = rowNorms(x3,method="euclidean")

x4 = matrix(rnorm((5*10^5)*500, 0,1),5*10^5,500)
n4 = rowNorms(x4,method="euclidean")

for (i in 1:10)
{
  x = matrix(rnorm((5*10^5)*500, 0,1),5*10^5,500)
  assign(paste("n",i,sep = ""),rowNorms(x,method="euclidean"))
  
}


data = lapply(paste0("n", 1:10), function(x) eval(as.name(x)))
datanew = do.call("cbind",data)
p = datanew/sum(datanew)

time = matrix(0,6,2)
r = (5:10)*1000
for (i in 1:6)
{
  datax = matrix(rnorm(r[i]*504,0,1),r[i],504)
  datay = rnorm(r[i],0,1)
  t1 = Sys.time()
  indexnew = sample(5*10^6,r[i],replace = TRUE,prob = p)
  a = t(datax)%*%datax
  b = ginv(a)
  c = t(datax)%*%datay
  d = b%*%c
  t1 = Sys.time()-t1
  t2 = Sys.time()
  indexnew = sample(5*10^6,r[i],replace = TRUE)
  a = t(datax)%*%datax
  b = ginv(a)
  c = t(datax)%*%datay
  d = b%*%c
  t2 = Sys.time()-t2
  result = c(t1,t2)
  time[i,] = result
  print(i)
}


datax = matrix(rnorm(r*504,0,1),r,504)
datay = rnorm(r,0,1)

t1 = Sys.time()
indexnew = sample(5*10^6,5000,replace = TRUE,prob = p)
a = t(datax)%*%datax
b = ginv(a)
c = t(datax)%*%datay
d = b%*%c
t1 = Sys.time()-t1



install.packages("data.table")
library(data.table)
index <- seq(0,5*10^4, by = 15000)
