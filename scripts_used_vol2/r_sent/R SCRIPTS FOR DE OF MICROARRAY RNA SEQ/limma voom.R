library(limma)
library(edgeR)
counts<-read.delim("ALLTOGETHER.csv")
rownames(counts)<-counts[,1]
counts<-counts[,-1]
head(counts)

exp.info <-data.frame(batch=rep(c("01","02","03"),3), condition=c(rep("WT",3),rep("HET", 3),rep("HOM",3)))
dds1<-DGEList(counts = counts, lib.size = colSums(counts),
        norm.factors = rep(1,ncol(counts)), samples = exp.info,
        group = NULL, genes = NULL, remove.zeros = FALSE)
dge <- calcNormFactors(dds1, method = "TMM")
keep <- rowSums(cpm(dge)>=1) >= 1
d1<-dge[keep,]
dim(d1)
batch<-rep(c("01","02","03"),3)
genotype<-c(rep("WT",3), rep("HET",3), rep("HOM",3))
mm<-relevel(factor(mm, ref="WT"))
#mm<-model.matrix(~0+design+batch)
#y2<-voom(d1, mm, plot=T)
#fit<-lmFit(y2, mm)
#head(coef(fit))
#contr<-makeContrasts(genotypeHET-genotypeWT,genotypeHOM-genotypeWT, levels=colnames(coef(fit)))
#contr
#tfit<-contrasts.fit(fit, contr)
#eb<-eBayes(tfit)
#hist(eb$p.value)

mm<-model.matrix(~0+genotype+batch)
y<-voom(d1, mm, plot=T)
fit1<-lmFit(y, mm)
head(coef(fit1))
contr1<-makeContrasts(genotypeHET-genotypeWT,genotypeHOM-genotypeWT, levels=colnames(coef(fit1)))
tfit1<-contrasts.fit(fit1, contr1)
eb1<-eBayes(tfit1)
hist(eb1$p.value)
a<-decideTests(eb1, p.value=0.05, lfc=0.5)
vennDiagram(a)

tT1<-topTable(eb1, coef=1, n="Inf.",adjust.="BH", p.value=0.05)
tT2<-topTable(eb1, coef=2, n="Inf.",adjust.="BH", p.value=0.05)
dim(tT1)
dim(tT2)
infominerHET<-cbind(rownames(tT1),tT1[,2],tT1[,6])
infominerHOM<-cbind(rownames(tT2),tT2[,2],tT2[,6])
write.table(infominerHET,file="infoHETlimma.csv", sep=",",quote=F,row.names=F)
write.table(tT1,file="infoHOMlimma.csv", sep=",",quote=F,row.names=F)

hist(eb1$p.value)
res<-decideTests(eb1)
vennDiagram(res, include=c("up", "down"))
volcanoplot(eb1, coef=1, highlight=tT1[,6])
volcanoplot(eb1, coef=2,highlight=tT2[,6])


probemeans <- apply(d1, 1, mean, na.rm=T)
probesd <- apply(d1, 1, sd)
plot(probemeans, probesd)
probemeans <- apply(dds1, 1, mean, na.rm=T)
probesd <- apply(dds1, 1, sd)
plot(probemeans, probesd)
plot(pca$x)
text(pca$x, labels=colnames(d), cex=0.5)
