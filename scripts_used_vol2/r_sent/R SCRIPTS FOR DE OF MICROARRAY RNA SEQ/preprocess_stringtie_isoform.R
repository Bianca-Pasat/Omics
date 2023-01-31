library(edgeR)
library("preprocessCore")
# isoformSwitchAnalyzer on all
samples<-as.data.frame(stringtie$counts)
rownames(samples)<-samples[,1]
samples<-samples[,-1]
samples<-data.frame(samples[,c(4,5,6,10,11,12,1,2,3,7,8,9)])

sample_info<-as.data.frame(sampleInfo)
sample_info$replicate<-c(paste("R",1:3,sep=""))
sample_info$treatment<-c(rep("DMSO",6),rep("MKC",6))
sample_info$time<-rep(c(rep(8,3),rep(24,3)),2)
colnames(sample_info)<-c("sampleId", "sample","replicate","treatment","time")
sample_info<-sample_info[c(1,2,3,7,8,9,4,5,6,10,11,12),]

### OR ADD 1 CONSTANT
samples<-samples+1


# REMOVE LOW EXPRESSED GENES
y <- DGEList(counts=samples, samples=sample_info , genes=rownames(samples))
keep.exprs <- aveLogCPM(y) > 0
y <- y[keep.exprs, keep.lib.sizes=F]
### Normalization
#TMM
y <- calcNormFactors(y)

samples<-y$counts


# or quantile
function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  
  return(norm)
}

samples<-function.norm(y$counts)
rm(y)
