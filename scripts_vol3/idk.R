# create DGEobject 
y <- DGEList(counts=samples2, samples=si, genes=rownames(samples2))
names(y)


# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]


# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y) 

norm_samples<-function.norm(y$counts)

# batch effect corrections

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

#sampleInfo = DATA FRAME samples treatment replicate
comsva<- function(samples,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  adjustedc <- ComBat(samples, batch=sampleInfo$replicate)
  svaob<-adjustedc - min(adjustedc)
  svobj <- svaseq(svaob, mod1, mod0) 
  comsvare <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(list(comsvare,svobj))
}

comsva1<- function(samples,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  adjustedc <- ComBat(samples, batch=sampleInfo$replicate)
  svaob<-adjustedc - min(adjustedc)
  svobj <- svaseq(svaob, mod1, mod0) 
  comsvare <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(svobj)
}

svacom<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(as.matrix(svaob), mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv[c(1,2)])) 
  svacomre <- ComBat(cleaned_count, batch=sampleInfo$replicate)
  return(svacomre)
}

svAll<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(svaob, mod1, mod0) 
  svAllre <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(svAllre)
}

com<- function(samples,sampleInfo,replicate){
  comre <- ComBat(samples, batch=sampleInfo$replicate)
  return(comre)
}

############################ DIFFERENTIAL EXPRESSION ANALYSIS EdgeR
#each gene will get its own unique dispersion estimate, but the common dispersion is still used in the calculation.
si$groups<-paste(si$treatment,si$hours,sep="_")
design <- model.matrix(~0 + groups, si)
design <- cbind(design, sv$sv)
y <- DGEList(counts=norm_samples, samples=si, genes=rownames(norm_samples))
y <- estimateDisp(y, design)


# mean-variance plot. raw variances of the counts (grey dots), the variances using the tagwise

meanVarPlot <- plotMeanVar( y , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

#
fit <- glmFit(y, design)
colnames(design)<-make.names(colnames(design))
contr.matrix <- makeContrasts(
  pacdmso8 = groupspac_8-groupsdmso_8, 
  pacdmso24 = groupspac_24 - groupsdmso_24,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(y, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))


tfit <- treat(vfit, lfc=0.5)
dt <- decideTests(tfit)
summary(dt)


## 24 hours
D24vsM24 <- glmLRT(fit, contrast = c(-1,0,1,0,0,0)) # -1 is the denominator (check design matrix)
res_D24vsM24 <- topTags(D24vsM24, n=nrow(y))
resultados_D24vsM24 <- as.data.frame(res_D24vsM24)
topTags(D24vsM24) # just the top 10 by default

## 8 hours
D8vsM8 <- glmLRT(fit, contrast = c(0,-1,0,1,0,0)) # -1 is the denominator (check design matrix)
res_D8vsM8 <- topTags(D8vsM8, n=nrow(y))
resultados_D8vsM8 <- as.data.frame(res_D8vsM8)
topTags(D8vsM8) # just the top 10 by default

write.table(resultados_D24vsM24, "results_D24vsM24_ALL.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
write.table(resultados_D8vsM8, "results_D8vsM8_ALL.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
