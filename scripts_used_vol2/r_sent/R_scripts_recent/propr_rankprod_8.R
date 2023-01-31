library("sva")
library(tidyverse)
library(ggplot2)


### proteins full and info
X<-read.table("/home/bianca/Desktop/proteomics_all/proteins_inputs/proteins_final.txt",sep="\t")
X_info<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final_info_correct_order.txt",sep="\t")
colnames(X_info)=X_info[1,]
X_info=X_info[-1,]
X_info$sample=colnames(X)
x_b=X
X=as.data.frame(lapply(X, function(x) as.numeric(x)))
rownames(X)=rownames(x_b)

# ## exclude low counts
# keep=apply(X,1,function(x) sum(x>=0) >= 3)
# # 
# # keep=apply(X,1,function(x) sum(x>=15) >= 3)
# 
# X1=X[keep,]


## another exculsion method
X1=X[!(rowSums(X==0)>=2),]



## rankprod differentially expressed
de_prot_8<-read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
tr_prop_8=read.table("/home/bianca/Desktop/transcripts_8.txt",sep="\t")

a=c(rownames(de_prot_8), tr_prop_8$x)
de_8_ex=X1[a,]
de_8_ex=X1[rownames(X1) %in% a,]
de_8_ex=de_8_ex[,c(10:12,22:24)]
# 
# all_8=X[,c(10:12,22:24)]
# X_info=X_info[c(10:12,22:24),]
# 
# keep=apply(de_8_ex,1,function(x) sum(x>=15) >=2)



de_8_ex1=de_8_ex[!(rowSums(de_8_ex==0)>=2),]
de_8_ex1=de_8_ex[!(rowSums(de_8_ex==0)),]


trans_cts=as.data.frame(na.omit(de_8_ex1))
trans_cts=na.omit(trans_cts)



#8 hours
de_8_ex=de_8_ex + 1
X3=as.data.frame(lapply(de_8_ex, function(x) as.numeric(x)))
rownames(X3)=rownames(de_8_ex)

trans_cts=as.data.frame(X3)

trans_cts=de_prot_8

sample_info=X_info
sample_info$sample=sample_info$short_name


### try
X1=X1[,c(10:12,22:24)]
X_info=X_info[c(10:12,22:24),]

trans_cts=as.data.frame(all_8)

trans_cts=as.data.frame(X1)
sample_info=X_info
sample_info$sample=sample_info$short_name

################## visual
demod<-read.table("/home/bianca/Desktop/gene_prio_propr_rankprod_8.tsv",sep="\t")
demod<-read.table("/home/bianca/Desktop/prior_genes_deModules_rankprod_8.tsv",sep="\t")
de8<-read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
fplot=de8[rownames(de8) %in% demod$V2,]
fplot$genes=rownames(fplot)
colnames(fplot)=c("logFC",'pval',"genes")
fplot$direction=ifelse(fplot$logFC < 0, "downregulated", "upregulated")

ggplot(fplot, aes(y=rownames(fplot), x=logFC)) + 
  geom_point(aes(size=  pval, col=direction)) + 
  scale_color_manual(values=c("blue3", "brown2")) + 
  ggtitle("Differentially regulated proteins after propd and RankProd") +
  theme(axis.title.y =element_blank(), panel.background = element_blank(), panel.grid.major = element_line(color = "gray"))
