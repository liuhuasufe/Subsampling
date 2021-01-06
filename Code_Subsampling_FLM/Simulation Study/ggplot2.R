x = 1:10
y = 2:11
z = 3:12
datag = data.frame(cbind(x,y,z))
p1 = ggplot(datag, aes(x=x)) + 
  geom_point(aes(y=y)) + 
  geom_line(aes(y=y , color="cyan")) +
  geom_point(aes(y=z)) + 
  geom_line(aes(y=z, color="red"))
p1+theme(legend.position=c(0.9,0.5))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(legend.title = element_blank())
  # guides(fill = guide_legend(title = NULL))