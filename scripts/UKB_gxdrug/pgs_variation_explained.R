# Creates a plot

library(ggplot2)
tmp <- data.frame(type=c("Users","Non-users"),R2=c(0.0904657,0.0317438))
g2 <- ggplot(tmp,aes(x=type,y=R2,fill=type)) + geom_bar(stat = "identity") + scale_fill_manual(values=c('darkorange','blue3')) + theme_bw() + theme(panel.grid=element_blank()) +
  labs(x='Group',y='Nagelkerke R-squared') + guides(fill=F) + coord_flip(); g2
