pgs_name <- type.lst[i]
pgs_name <- 'PGS000007'

dataf.mg$combined <- scale(apply(data.frame(scale(dataf.mg[,'PGS000007']),scale(dataf.mg[,'P_0.00001.R2_0.1.KB_250'])),1,mean))
ind <- sample(1:10,nrow(dataf.mg),replace=TRUE)
res <- list()
for (ind.sub in 1:10) {
  print(ind.sub)
  df.train <- dataf.mg[ind!=ind.sub,]
  df.test <- dataf.mg[ind==ind.sub,]
  mod1 <- glm(disease ~ `PGS000007`+bmi+age+menopause+
               number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.train,family = binomial(link='logit'))
  pred1 <- predict(mod1,newdata = df.test,type='response')
  mod2 <- glm(disease ~ P_0.00001.R2_0.1.KB_250+bmi+age+menopause+
               number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.train,family = binomial(link='logit'))
  pred2 <- predict(mod2,newdata = df.test,type='response')
  mod3 <- glm(disease ~ combined+bmi+age+menopause+
                number_live_birth+one_birth+
                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.train,family = binomial(link='logit'))
  pred3 <- predict(mod3,newdata = df.test,type='response')
  res[[ind.sub]] <- data.frame(r1=cor(df.test$disease,pred1),r2=cor(df.test$disease,pred2),r3=cor(df.test$disease,pred3))
}
res.df <- do.call(rbind,res)
t.test(res.df[,1],res.df[,3],paired=TRUE)
t.test(res.df[,2],res.df[,3],paired=TRUE)
t.test(res.df[,1],res.df[,2],paired=TRUE)

for (i in 1:length(type.lst)) {
  print(paste0(i,'/',length(type.lst)))
  pgs_name <- type.lst[i]
  mod <- glm(disease ~ dataf.mg[,type.lst[i]]+bmi+age+menopause+
               number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg,family = binomial(link='logit'))
  param <- nagelkerke(fit=mod,null=null.mod)$Pseudo.R.squared.for.model.vs.null[3,]
  param.vec <- c(param.vec,param)
}
