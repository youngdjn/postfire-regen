library(ggplot2)

setwd("C:/Users/DYoung/Documents/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

s <- seq(-2.2,2.2,by=.1)


dry <- dnorm(s,mean=1)
wet <- dnorm(s,mean=-1)

d.dry <- data.frame(s,prob=dry,norm="Hot and dry")
d.wet <- data.frame(s,prob=wet,norm="Cool and wet")

d <- rbind(d.dry,d.wet)

p <-ggplot(d,aes(x=s,y=prob,color=norm)) +
  geom_line(size=1.5) +
  scale_color_manual(values=c("Cool and wet"="turquoise4","Hot and dry"="darkorange1")) +
  theme_bw(16)+
  labs(x="Post-fire weather anomaly",y="Seedling recruitment",color="Normal climate") +
  geom_vline(xintercept=0,linetype="dashed") +
  theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
  scale_x_continuous(breaks=c(-2.2,0,2.2), labels= c("unusually hot or dry","0","unusually cool or wet")) +
  scale_y_continuous(breaks=c(0,0.4),labels=c("0","high")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.8, 0.4),plot.margin = unit(c(.5,2,.5,.5), "cm"))

tiff(file=paste0("../Figures/Fig1_hypothetical_recruitment_",Sys.Date(),".tiff"),width=1500,height=1000,res=200) 
p
dev.off()
