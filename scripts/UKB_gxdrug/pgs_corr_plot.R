library(data.table)

cor_mat <- fread('/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/cor_mat.txt',data.table = F,stringsAsFactors = F)
rownames(cor_mat) <- cor_mat[,1]
cor_mat <- as.matrix(cor_mat[,-1])
cor_df <- fread('/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/cor_df.txt',data.table = F,stringsAsFactors = F)

# library(corrplot)
# corrplot(cor_mat,type='lower',is.corr = FALSE)
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
        axis.text.x=element_text(angle=70,hjust = 1)) + 
  guides(fill=F) +
  # coord_fixed() +
  labs(title='Polygenic score correlations');g


