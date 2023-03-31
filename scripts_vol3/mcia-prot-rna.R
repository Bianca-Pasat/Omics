##### LOAD MATRICES THAT DON'T CONTAIN 0: REMOVED LOWLY EXPRESSED VALUES
library(omicade4)
prot<-readRDS("/home/bianca/Desktop/prot.RDS")
x_b=rownames(prot)
prot1=as.data.frame(sapply(prot, function(x) as.numeric(x)))
rownames(prot)=x_b

#rna<-readRDS("/home/bianca/Desktop/rna.RDS")
rna<-readRDS("/home/bianca/Desktop/rna_filt_norm.RDS")
#rna=rna[,c(1,5,4,6,3,2,9,8,10,12,11,7)]
#rna<-read.table("/home/bianca/Downloads/rna_normalized_final.txt", sep="\t")
#rna=rna[,c(1,5,4,6,3,2,9,8,10,12,11,7)]
#saveRDS(rna,"/home/bianca/Desktop/rna_filt_norm.RDS")
#write.table(sample_info_RNAseq,"/home/bianca/Downloads/rna_normalized_final_sample_info.txt", sep="\t",quote=F,  row.names=T)

#### FIX MICROARRAY TO HAVE HGNC SYMBOL!!!!!!!!
micro<-readRDS("/home/bianca/Desktop/micro.RDS")
#colnames(micro)=rep("MKC_24",3)


pr1=prot[,c(1,2,3,7,8,9)]
r1=rna[,c(1,2,3,7,8,9)]
mda<-list(as.matrix(pr1), as.matrix(r1))
mda<-list(as.matrix(prot), as.matrix(rna))

#### all 8 de
demirna8<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
colnames(demirna8)<-c("log-ratio","Pvalue")
demirna8$Gene=rownames(demirna8)
demirna8=demirna8[,c(3,1,2)]

mirna8<-readRDS("/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_8")
mrna8<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")

demi8<-mirna8[rownames(mirna8) %in% demirna8$Gene, ]
dedemi8<-demirna8[demirna8$Gene %in% rownames(mirna8) , ]

mrna8<-mrna8[,colnames(mirna8)]

prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot8")
lipids6<-readRDS("/home/ben/Desktop/lipidsb/lipids6")
lipred<-lipids6[,c(1,3,4,5,6,8)]
colnames(lipred)=colnames(prot8)


### start with clean datasets but not de
mrna8=readRDS("/home/ben/Desktop/rna_bianca/clean_AllRna_8")
mirna8=readRDS("/home/ben/Desktop/mirna_bianca/clean_mirna_comsva")
mirna8<-mirna8[,c(1:6)]
mrna8<-mrna8[,colnames(mirna8)]
lipids6<-readRDS("/home/ben/Desktop/lipidsb/clean6")
lipred<-lipids6[,c(1,3,4,5,6,8)]
prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot8")
colnames(lipred)=colnames(prot8)


mda<-list(as.matrix(mrna8), as.matrix(demi8))
mda<-list(mRNA=as.matrix(mrna8), miRNA=as.matrix(demi8), lipids=as.matrix(lipred), proteins=as.matrix(prot8))


all(apply((x <- sapply(mda, colnames)), 2, function(y)
  + identical(y, x[,1])))


layout(matrix(1:4, 2, 2))
par(mar=c(2, 1, 0.1, 6))
for (df in mda) {
  d <- dist(t(df))
  hcl <- hclust(d)
  dend <- as.dendrogram(hcl)
  plot(dend, horiz=TRUE)
}


mcoin <- mcia(mda)
cond=c(rep("DMSO",3),rep("MKC",3))
### plot to see which genes differentiate your conditions!
plot(mcoin, axes=1:2, phenovec=cond, sample.lab=FALSE, df.color=1:2)
plot(mcoin, axes=1:2, phenovec=cond, sample.lab=FALSE, df.color=1:4)

## choose genes~proteins according to your pca
goi <- selectVar(mcoin, a1.lim=c(-Inf,-2), a2.lim=c(-Inf, Inf))
goi <- selectVar(mcoin, a1.lim=c(-Inf,), a2.lim=c(-Inf, Inf))

goi2 <- selectVar(mcoin, a1.lim=c(2,Inf), a2.lim=c(-Inf, Inf))
goi05 <- selectVar(mcoin, a1.lim=c(0.5,Inf), a2.lim=c(-Inf, Inf))
goi1 <- selectVar(mcoin, a1.lim=c(1,Inf), a2.lim=c(-Inf, Inf))
goi0 <- selectVar(mcoin, a1.lim=c(0,Inf), a2.lim=c(-Inf, Inf))

write.table(goi2$var,"/home/ben/Desktop/integrated_mcia_2_features", sep="\t", quote=F, row.names = F)
write.table(goi05$var,"/home/ben/Desktop/integrated_mcia_05_features", sep="\t", quote=F, row.names = F)
write.table(goi1$var,"/home/ben/Desktop/integrated_mcia_1_features_with_all_features", sep="\t", quote=F, row.names = F)
write.table(goi0$var,"/home/ben/Desktop/integrated_mcia_0_features", sep="\t", quote=F, row.names = F)


geneStat <- plotVar(mcoin, var=c("S100B", "S100A1"), var.lab=TRUE)
