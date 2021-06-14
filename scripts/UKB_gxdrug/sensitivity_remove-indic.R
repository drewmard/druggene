
df.mg.tmp <- merge(df.mg,df3[,c('eid','first_diagnosis_days')],by='eid')
df.mg.tmp <- subset(df.mg.tmp,first_diagnosis_days>0 & asthma==0 & allergy==0 & crohns_disease==0 &
                    ulcerative_colitis==0 & rheumatoid_arthritis==0 & sle==0)
table(df.mg.tmp$S01BA)
mod0 <- coxph(Surv(days,disease) ~ rs62119267_chr19*S01BA+
                 bmi+age+#CRP+townsend+
                 # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
                 menopause+number_live_birth+one_birth+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg.tmp)
mod1 <- coxph(Surv(days,disease) ~ rs62119267_chr19*S01BA+
               bmi+age+#CRP+townsend+
               # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
               menopause+number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg)
summary(mod0)$coef["rs62119267_chr19:S01BA",]
summary(mod1)$coef["rs62119267_chr19:S01BA",]

mod0 <- coxph(Surv(days,disease) ~ pgs*S01BA+
               bmi+age+#CRP+townsend+
               # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
               menopause+number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg.tmp)
mod1 <- coxph(Surv(days,disease) ~ pgs*S01BA+
               bmi+age+#CRP+townsend+
               # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
               menopause+number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg)
summary(mod0)$coef["pgs:S01BA",]
summary(mod1)$coef["pgs:S01BA",]

mod0 <- coxph(Surv(days,disease) ~ rs4784227_chr16*S01BA+
               bmi+age+#CRP+townsend+
               # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
               menopause+number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg.tmp)
mod1 <- coxph(Surv(days,disease) ~ rs4784227_chr16*S01BA+
               bmi+age+#CRP+townsend+
               # asthma + allergy + crohns_disease + ulcerative_colitis+rheumatoid_arthritis + sle+
               menopause+number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg)
summary(mod0)$coef["rs4784227_chr16:S01BA",]
summary(mod1)$coef["rs4784227_chr16:S01BA",]
