library(data.table)
polyscore.full <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_drug.txt",data.table = F,stringsAsFactors = F)
polyscore.full <- polyscore.full[,21:ncol(polyscore.full)]
colnames(polyscore.full) <- c('Mich: P < 0.05','Mich: P < 0.005','Mich: P < 0.0005','Mich: P < 0.00001', 'Mich: P < 5e-8','Mav: Stepwise Forward','Mav: LASSO')
cor.mat <- cor(polyscore.full)
cor.mat[lower.tri(cor.mat)] <- NA
cor_df <- melt(cor.mat)

library(ggplot2)
g <-ggplot(cor_df,aes(Var1,Var2,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient(low = "white", high = "red",
                      name="Correlation") +
  theme_minimal() +#+ # minimal theme
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(angle=30,hjust = 1),
        panel.grid = element_blank()) + 
  # guides(fill=F) +
  # coord_fixed() +
  labs(title='Polygenic score correlations');g

nagelr2 <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/pgs_nagelkerke_r2.txt",data.table = F,stringsAsFactors = F)
nagelr2[,1] <- c('Mich: P < 0.05','Mich: P < 0.005','Mich: P < 0.0005','Mich: P < 0.00001', 'Mich: P < 5e-8','Mav: Stepwise Forward','Mav: LASSO')
g2 <- ggplot(nagelr2,aes(x=pgs_name,y=r2,fill=pgs_name)) + geom_bar(stat = "identity") + theme_bw() + theme(panel.grid=element_blank()) +
  labs(x='PGS type',y='Nagelkerke R-squared') + guides(fill=F) + coord_flip(); g2

library(cowplot)
png("~/Documents/Research/druggene/output/pgs_int/suppfig3.png",height=2*5970/3,width=2*5510/3,res=500)
plot_grid(g,g2,ncol=1,rel_heights = c(0.5,0.3))
dev.off()
