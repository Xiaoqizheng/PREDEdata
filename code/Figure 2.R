
load("Figure2.RData")

p1 = ggplot(dat_A, aes(x = step, y = aic)) +
  geom_line(aes(linetype = K,color = K,size = K))+
  scale_size_manual(values=c(1,1,1))+
  scale_color_manual(values = brewer.pal(9, "Set1")[c(2,3,4)]) +
  scale_x_continuous(breaks=seq(2,15,2))+
  geom_point(x = 3,y = min(v1,na.rm=T),color="#E41A1C",size=2)+
  geom_point(x = 6,y = min(v2,na.rm=T),color="#E41A1C",size=2)+
  geom_point(x = 10,y = min(v3,na.rm=T),color="#E41A1C",size=2)+
  theme_classic()+
  xlab("Predicted number of cell types")+
  ylab("AIC")+
  labs(fill = "K")

p2 = ggplot(dat_B, aes(x = step, y = aic, group = K1,color = K1)) +
  geom_line(aes(linetype = K1,color = K1,size = K1))+
  scale_size_manual(values=c(1,1,1,1,1))+
  scale_color_manual(values = c(brewer.pal(9, "Set1")[c(2,3,4,7)],"#4D4D4D")) +
  scale_x_continuous(breaks=seq(2,15,2))+
  scale_y_continuous(breaks=seq(-60000,-28000,8000))+
  geom_vline(aes(xintercept=6), colour="#E41A1C", linetype=1,size=1.1)+
  theme_classic()+
  xlab("Predicted number of cell types")+
  ylab("AIC") + 
  labs(fill = "K1")


p3 = ggplot(dat_C, aes(x = step, y = aic))+
  geom_line(aes(linetype = K1,color = K1,size = K1))+
  scale_size_manual(values=c(1,1,1,1,1))+
  scale_color_manual(values = c(brewer.pal(9, "Set1")[c(2,3,4,7)],"#4D4D4D")) +
  scale_x_continuous(breaks=seq(2,15,2))+
  scale_y_continuous(breaks=seq(-40000,-2000,4000))+
  geom_vline(aes(xintercept=10), colour="#E41A1C", linetype=1,size=1.1)+
  theme_classic()+
  xlab("Predicted number of cell types")+
  ylab("AIC") + 
  labs(fill = "K1")


p =ggarrange(p1,p2,p3,ncol=3,nrow=1,labels=c("A","B","C"))
ggsave(plot = p,filename="Figure2.pdf",width=9, height = 2.5)

