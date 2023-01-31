library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library(preprocessCore)

load("/home/ben/Desktop/hell.RData")

raw_counts <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/countsaitor.txt", sep="", quote = "") 
raw_counts <- read.delim("/home/bianca/Downloads/countsaitor.txt", sep="", quote = "") 
rownames(raw_counts)<-raw_counts[,1]
exp <- raw_counts[,c(7:18)] # first 6 columns is annotation
sample_info_RNAseq <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/sample_info_RNAseq.txt", sep="") 
sample_info_RNAseq <- read.delim("/home/bianca/Downloads/sample_info_RNAseq.txt", sep="", row.names = 1) 
colnames(sample_info_RNAseq)<-c("sampleID", "replicate", "times","treatment")
sample_info_RNAseq$condition<-c(rep(0,6), rep(1,6))
#colnames(exp)<-sample_info_RNAseq$sample
sample_info<-sample_info_RNAseq
sample_info$sample<-c("D8R1","D24R3","D24R2","D8R3","D8R2","D24R1","M24R3","M8R2","M8R1","M8R3","M24R2","M24R1")
colnames(exp)<-sample_info$sample

library(biomaRt)
# #annotation of genes
httr::set_config(httr::config(ssl_verifypeer = FALSE))
values=rownames(exp)
attributes = c("ensembl_gene_id", "hgnc_symbol")
filters="ensembl_gene_id"

homo.anno = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")
# #listAttributes(homo.anno)

genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
exp$ens=rownames(exp)
exp1<-merge(genes, exp, by.x=1,by.y=13)

library(plyr)
exp2<-ddply(exp1, .(hgnc_symbol), numcolwise(mean))
rownames(exp2)=exp2$hgnc_symbol
exp2<-exp2[,-1]

y <- DGEList(counts=exp, samples=sample_info_RNAseq , genes=rownames(exp))
y <- DGEList(counts=exp2, samples=sample_info_RNAseq)
y <- DGEList(counts=exp2, samples=sample_info)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- aveLogCPM(y) > 0
y <- y[keep.exprs, keep.lib.sizes=F]
write.table(y$counts,"/home/bianca/Desktop/mRNA/rna_filtered.txt",sep = "\t",quote = F, row.names = F)

# keep.exprs <- filterByExpr(y)
# y <- y[keep.exprs, keep.lib.sizes=T]
# 
# keep.exprs <- rowSums(y$counts==0) < 12
# y <- y[keep.exprs,, keep.lib.sizes=FALSE]


#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  
  return(norm)
}


norm_samples = function.norm(log2(y$counts))

norm_samples = function.norm(log2(exp2))
write.table(norm_samples,"/home/bianca/Downloads/rna_normalized_final.txt", sep="\t",quote=F,  row.names=T)
write.table(sample_info_RNAseq,"/home/bianca/Downloads/rna_normalized_final_sample_info.txt", sep="\t",quote=F,  row.names=T)


# 8 hours
norm_samples<-norm_samples[,c(1,4,5,8,9,10)]
y$samples<-y$samples[c(1,4,5,8,9,10),]
samples_info<-sample_info_RNAseq[c(1,4,5,8,9,10),]
sample_info<-sample_info[c(1,4,5,8,9,10),]
svacomRe8<-svacom(norm_samples, sample_info, treatment,replicate)
saveRDS(svacomRe8,"/home/ben/Desktop/clean_AllRna_8")

# 24 hours
norm_samples<-norm_samples[,c(2,3,6,7,11,12)]
y$samples<-y$samples[c(2,3,6,7,11,12),]
samples_info<-sample_info_RNAseq[c(2,3,6,7,11,12),]
sample_info<-sample_info[c(2,3,6,7,11,12),]
svacomRe24<-svacom(norm_samples, sample_info, treatment,replicate)
saveRDS(svacomRe24,"/home/ben/Desktop/clean_AllRna_24")

### keep clean matrices with DE genes
de8<-read.table("/home/ben/Downloads/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt",sep="\t")
de8<-rownames(de8)
clean_de_8<-svacomRe8[rownames(svacomRe8) %in% de8,]
saveRDS(clean_de_8,"/home/ben/Desktop/clean_DeRna_8")

de24<-read.table("/home/ben/Downloads/rnaseq_rankprod_results_24.txt",sep="\t")
de24<-rownames(de24)
genes <- getBM(attributes = attributes, filters=filters,values = de24, mart = homo.anno)

clean_de_24<-svacomRe24[rownames(svacomRe24) %in% genes$hgnc_symbol,]
saveRDS(clean_de_24,"/home/ben/Desktop/clean_DeRna_24")

### keep de with filtered and normalized matrices
filt_norm_de_8<-norm_samples[rownames(norm_samples) %in% de8,]
saveRDS(filt_norm_de_8,"/home/ben/Desktop/rna_filt_norm_de_8")

filt_norm_de_24<-norm_samples[rownames(norm_samples) %in% genes$hgnc_symbol,]
saveRDS(filt_norm_de_24,"/home/ben/Desktop/rna_filt_norm_de_24")
### OTHER NORMALIZATION METHODS

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards

y$counts <- cpm(y, log=F) 
#Normalize 
y <- calcNormFactors(y)




