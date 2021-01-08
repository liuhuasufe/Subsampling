library(reshape2)
datapn = melt(datap, id.vars="tobs")
ggplot(datapn, aes(x=tobs, y=value)) + 
  geom_line(aes(color=variable,linetype = variable))+
  labs(x="t", y="x(t)")+
  theme(legend.position=c(2,1))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())