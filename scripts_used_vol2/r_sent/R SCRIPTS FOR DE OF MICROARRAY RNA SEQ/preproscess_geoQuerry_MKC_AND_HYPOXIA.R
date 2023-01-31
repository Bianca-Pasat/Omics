library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library(GEOquery)


### inputs microarray
gse <- getGEO("GSE99766", GSEMatrix = TRUE)
columns<-pData(gse[[1]])
#fData(gse[[1]])
exp<-exprs(gse[[1]])


### keep mkc dmso y and normal at 4 hours
exp<-exp[,c(1:6,10:15)]
colnames(exp)<-c(paste("DMSO_4",1:3,sep="_"), paste("MKC_4",1:3,sep="_"),paste("DMSO_4",1:3,sep="_Y_"), paste("MKC_4",1:3,sep="_Y_"))

# make samples info
samples_info<-data.frame(treatment=rep(c(rep(1,3), rep(0,3)),2))
samples_info$replicate<-c(rep(c("1","2","3"), 4))
samples_info$sample <-colnames(exp)
samples_info$times <- rep(4,12)
samples_info$condition<- c(rep("Normal",6),rep("Hypoxia",6))
samples_info$lim1<-c(rep(1,6),rep(0,6))
samples_info$lim2<-rep(c(rep(1,3),rep(0,3)),2)
# ANNOTATION
httr::set_config(httr::config(ssl_verifypeer = FALSE))

values=rownames(exp)
attributes = c("affy_hugene_2_0_st_v1", "ensembl_gene_id")
attributes = c("affy_hugene_2_0_st_v1", "hgnc_symbol")
filters="affy_hugene_2_0_st_v1"


homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "useast")
#listAttributes(homo.anno)

genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

exp<-as.data.frame(exp)
exp$names<-rownames(exp)

te<-merge(exp, genes, by.x = "names", by.y=1)
te<-te[!duplicated(te$hgnc_symbol),]
rownames(te)<-te$hgnc_symbol
exp<-te[,-c(1,14)]

hif<-read.table("/home/bianca/Downloads/hif_signature.txt")
expHif<-exp[rownames(exp) %in% hif$V1,]
# create DGEobject 
y <- DGEList(counts=exp, samples=samples_info)
y <- DGEList(counts=expHif, samples=samples_info)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y, log=TRUE) 

y <- calcNormFactors(y)

### OUR DIFF with hif signature ATTENTION WHICH RP_TOP!!
confirmed_hif_signature<-RP_TOP[rownames(RP_TOP) %in% hif$V1,]

## rowk but not significant genes
sample_info$treatment<-relevel(sample_info$treatment,ref="0")
sample_info$lim<-c(rep(c(rep(0,3),rep(1,3)),2))
sample_info$lim2<-c(rep(0,3),rep(-1,3),rep(1,3),rep(-1,3))
sample_info$lim3<-c(rep(2,3),rep(0,3),rep(2,3),rep(1,3))
design<- model.matrix(~0+treatment:condition, as.data.frame(sample_info))
design<- model.matrix(~0+lim, as.data.frame(sample_info))
design
y<-estimateDisp(y)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit,contrast = c(-1,1))
qlf <- glmQLFTest(fit,coef=2)
interconnected<-topTags(qlf,n=Inf, p.value = 0.05)
