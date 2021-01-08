rm(list = ls())
library(ggplot2)
Method = c(rep("Lopt",6),rep("Unif",6))
mIMSE_1 = apply(IMSE_1,c(2,3),mean)
datag = data.frame(r = rep(r,2), Method = Method, value = as.vector(t(mIMSE_1)))
ggplot(data=datag, aes(x=r, y=value, group=Method,colour=Method)) + 
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("black", "red"))+
  # ggtitle("Scenario I")+
  labs(x="subsample size", y="IMSE")+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())


# p1 = ggplot(datag, aes(x=r)) + 
#   geom_point(aes(y=Lopt)) + 
#   geom_line(aes(y=Lopt , color="black")) +
#   geom_point(aes(y=Unif)) + 
#   geom_line(aes(y=Unif, color="red"))
# p1+scale_fill_discrete(name = "Method", labels = c("Lopt", "Unif"))+
#   theme(legend.position=c(0.9,0.5))+
#   theme(legend.key = element_blank())+
#   theme(legend.background = element_blank())
# 
#   # theme(legend.title = element_blank())