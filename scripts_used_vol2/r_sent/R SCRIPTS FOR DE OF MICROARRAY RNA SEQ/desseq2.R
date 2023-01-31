counts<-read.delim("ALLTOGETHER.csv")
rownames(counts)<-counts[,1]
counts<-counts[,-1]
exp.info <-data.frame(batch=rep(c("01","02","03"),3),condition=c(rep("WT",3),rep("HET",3),rep("HOM",3)))
rownames(exp.info) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=exp.info,design=~batch+condition)
keep <- rowSums(counts(dds)) >= 10
dds<- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "WT")
dds$condition <- factor(dds$condition, levels = c("WT","HET","HOM"))
design(dds) <- formula(~ batch + condition)
des <- DESeq(dds)

des<-DESeq(dds)
resHOM <- results(des,alpha=0.05, contrast=c("condition", "HOM", "WT"))
res.sigHOM <- resHOM[which(resHOM$padj<0.05),]
res.sigHOM <- res.sigHOM[order(res.sigHOM$padj),]
dim(res.sigHOM)
head(res.sigHOM)
infominer<-cbind(rownames(res.sigHOM),res.sigHOM[,2],res.sigHOM[,6])
write.table(res.sigHOM,file="wt-hom.csv", sep=",",quote=F,row.names=F)
write.table(infominer,file="INFOwt-homcpm12.csv", sep=",",quote=F,row.names=F)

sigHOM <- lfcHOM[which(lfcHOM$svalue<0.01),]
sigHOM <- sigHOM[order(sigHOM$padj),]


res <- results(des,alpha=0.05, contrast=c("condition", "HET", "WT"))
res.sig <- res[which(res$padj<0.05),]
res.sig <- res.sig[order(res.sig$padj),]
dim(res.sig)
head(res.sig)
infominer1<-cbind(rownames(res.sig),res.sig[,2],res.sig[,6])
write.table(res.sig,file="wt-het.csv", sep=",",quote=F,row.names=F)
write.table(infominer1,file="INFOwt-homcpm12.csv", sep=",",quote=F,row.names=F)