library(data.table)
df <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/interaction_results.txt",data.table = F,stringsAsFactors = F)
library(qqman)
qq(df$p)
