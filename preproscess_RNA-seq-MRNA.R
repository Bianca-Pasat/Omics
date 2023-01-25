library("edgeR")
library("limma")
library("sva")
library("biomaRt")



raw_counts <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/countsaitor.txt", sep="", quote = "") 
raw_counts <- read.delim("/home/bianca/Downloads/countsaitor.txt", sep="", quote = "") 
rownames(raw_counts)<-raw_counts[,1]
exp <- raw_counts[,c(7:18)] # first 6 columns is annotation
sample_info_RNAseq <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/sample_info_RNAseq.txt", sep="", row.names=1) 
sample_info_RNAseq <- read.delim("/home/bianca/Downloads/sample_info_RNAseq.txt", sep="", row.names = 1) 
colnames(sample_info_RNAseq)<-c("sample", "replicate", "times","treatment","treat")
colnames(sample_info_RNAseq)<-c("sample", "replicate", "times","treatment")
sample_info_RNAseq$condition<-c(rep(0,6), rep(1,6))
colnames(exp)<-sample_info_RNAseq$sample
samples_info<-sample_info_RNAseq
sample_info<-sample_info_RNAseq

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

y <- DGEList(counts=exp, samples=sample_info_RNAseq , genes=rownames(exp))


# REMOVE LOW EXPRESSED GENES
keep.exprs <- aveLogCPM(y) > 0
y <- y[keep.exprs, keep.lib.sizes=F]

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

# 8 hours
norm_samples<-norm_samples[,c(1,4,5,8,9,10)]
y$samples<-y$samples[c(1,4,5,8,9,10),]
samples_info<-sample_info_RNAseq[c(1,4,5,8,9,10),]
sample_info<-samples_info


# 24 hours
norm_samples<-norm_samples[,c(2,3,6,7,11,12)]
y$samples<-y$samples[c(2,3,6,7,11,12),]
samples_info<-sample_info_RNAseq[c(2,3,6,7,11,12),]
sample_info<-samples_info


### OTHER NORMALIZATION METHODS

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards

y$counts <- cpm(y, log=F) 
#Normalize 
y <- calcNormFactors(y)




