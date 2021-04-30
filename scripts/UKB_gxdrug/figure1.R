library(ggplot2)
tmp <- data.frame(group=c(rep("Medication Users",11),rep("Non-users",11)),
                  percentile=rep(seq(0,1,by=0.1),2),
                  incidence=c((seq(0.05,0.2,length.out = 11)),rep(0.05,11))
)
g1 <- ggplot(tmp,aes(x=percentile,y=incidence,col=group)) + geom_line(aes(lty=group)) + theme_bw() + 
  theme(panel.grid=element_blank()) + 
  labs(x='Genetic risk',y='Disease incidence') +
  guides(col=F,lty=F) + 
  scale_x_continuous(breaks=c(0,1),labels=c("Low","High")) +
  scale_y_continuous(breaks=c(0,0.3),labels=c("Low","High"),limits = c(0,.3))
library(ggplot2)
tmp2 <- data.frame(group=c(rep("Medication Users",11),rep("Non-users",11)),
                  percentile=rep(seq(0,1,by=0.1),2),
                  incidence=c((seq(0.05,0.2,length.out = 11)),rep(seq(0.05,0.3,length.out = 11)))
)
g2 <- ggplot(tmp2,aes(x=percentile,y=incidence,col=group)) + geom_line(aes(lty=group)) + theme_bw() + 
  theme(panel.grid=element_blank()) + 
  labs(x='Genetic risk',y='Disease incidence') +
  guides(col=F,lty=F) + 
  scale_x_continuous(breaks=c(0,1),labels=c("Low","High")) +
  scale_y_continuous(breaks=c(0,0.3),labels=c("Low","High"),limits = c(0,.3))

library(cowplot)
plot_grid(g1,g2)