library(data.table)
dataf <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/sim/results3.txt",data.table = F,stringsAsFactors = F)
plot(dataf$r,dataf$int_beta)
dataf2 <- as.data.frame(aggregate(dataf$int_p,list(dataf$beta_cor),function(x) {mean(x<.05)}))
dataf2$nindiv <- 5000
dataf2$Sig <- NA; for (i in 1:nrow(dataf2)) dataf2$Sig[i] <- binom.test(dataf2$x[i]*100,100,.05)$p.value < 0.05; dataf2$Sig <- ifelse(dataf2$Sig,'*','')
library(ggplot2)

g0 <- ggplot(dataf,aes(x=as.factor(round(r,1)),y=int_beta,fill=as.factor(round(r,1)))) + geom_violin() + theme_bw() + theme(panel.grid = element_blank()) + scale_fill_brewer(palette="Greens")+
  labs(x='Covariance(Beta, Omega)',y='Observed interaction effects') + guides(fill=F) + geom_jitter(width=0.1,alpha=0.2); g0

g <- ggplot(dataf2,aes(Group.1,nindiv,fill=x)) +
  geom_tile(color = "white")+
  geom_text(aes(label=Sig),size=rel(10)) +
  scale_x_continuous(breaks = dataf2$Group.1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
                       name="Correlation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid = element_blank(),
        plot.title=element_blank()) + 
  guides(fill=F) +
  labs(y='Power');g

library(cowplot)
plot_grid(g0,g,ncol=1,rel_heights = c(0.8,0.2))
