library(data.table)
library(cowplot)
library(ggfortify)
library(survival)

tmp <- fread("~/Documents/Research/druggene/output/pgs_drug_plusSNP.txt",data.table = F,stringsAsFactors = F)
tmp$snp1 <- tmp$rs62119267_chr19
tmp$P <- tmp$disease
tmp$drug <- tmp$S01BA
tmp1 <- subset(tmp,snp1>=1)
fit <- survfit(Surv(days, P) ~ drug, data = tmp1)
g1 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") + 
  lims(y=c(0.8,1))
tmp0 <- subset(tmp,snp1==0)
fit <- survfit(Surv(days, P) ~ drug, data = tmp0)
g2 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") +
  lims(y=c(0.8,1))
png("~/Documents/Research/druggene/output/snp_int/survival_plots.png",height = 3950/2,width=8500/2,res = 300*2)
plot_grid(g2,g1,ncol = 2)
dev.off()

fit <- survfit(Surv(days, P) ~ drug, data = tmp1)
g1 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank(),axis.text=element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") + 
  lims(y=c(0.8,1))
tmp0 <- subset(tmp,snp1==0)
fit <- survfit(Surv(days, P) ~ drug, data = tmp0)
g2 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank(),axis.text=element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") +
  lims(y=c(0.8,1))
png("~/Documents/Research/druggene/output/snp_int/survival_plots2.png",height = 3950/2,width=8500/2,res = 300*2)
plot_grid(g2,g1,ncol = 2)
dev.off()

g1

