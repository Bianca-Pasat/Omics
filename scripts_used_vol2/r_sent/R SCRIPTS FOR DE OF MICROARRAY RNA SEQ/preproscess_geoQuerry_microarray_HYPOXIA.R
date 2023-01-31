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

# make samples info BECAUSE OF RANKPROD 0=CONDITION
samples_info<-data.frame(treatment=rep(c(rep(1,3), rep(0,3)),2))
samples_info$replicate<-c(rep(c("1","2","3"), 4))
samples_info$sample <-colnames(exp)
samples_info$times <- rep(4,12)
samples_info$condition<- c(rep("Normal",6),rep("Hypoxia",6))


# ANNOTATION
httr::set_config(httr::config(ssl_verifypeer = FALSE))

values=rownames(exp)
attributes = c("affy_hugene_2_0_st_v1", "ensembl_gene_id")
attributes = c("affy_hugene_2_0_st_v1", "hgnc_symbol")
filters="affy_hugene_2_0_st_v1"


homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")
#listAttributes(homo.anno)

genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

mat <- match(rownames(exp), genes$gene.id)
genes <- genes[mat,]

# create DGEobject 
y <- DGEList(counts=exp, samples=samples_info, genes=genes)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y, log=TRUE) 

#Normalize
y <- calcNormFactors(y)


#### dmso 4 h with hypoxia 10, 11, 12 columns
expIn<-exp[,1:6]
expInrev<-exp[,c(4,5,6,1,2,3)]
#colnames(expIn)<-c(rep(0,3), rep(1,3))


samples_vector_file<-data.frame(c(rep(0,3), rep(1,3)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 2))
#samples_vector_file$replicate<-c(rep(c("a","b","c"), 2))
samples_vector_file$sample <-c(paste("DMSO", 1:3, sep=""), paste("MKC",1:3, sep=""))
samples_vector_file$treatment <-c(rep("DMSO", 3), rep("MKC",3))
colnames(expIn)<-samples_vector_file$sample


######## FOR ALL OF THE MKC treatments
expMKC<-exp[,c(1:9)]
samples_vector_file<-data.frame(c(rep(0,3), rep(1,6)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 3))
#samples_vector_file$replicate<-c(rep(c("a","b","c"), 2))
samples_vector_file$times<-c(rep(c(8,24),c(6,3)))
samples_vector_file$sample <-c(paste("DMSO", 1:3, sep=""), paste("MKC",1:6, sep=""))
samples_vector_file$treatment <-c(rep("DMSO", 3), rep("MKC",6))
colnames(expMKC)<-samples_vector_file$sample


### OR READ FRP, FILE



#annotation of genes
httr::set_config(httr::config(ssl_verifypeer = FALSE))

values=rownames(expIn)
attributes = c("affy_hugene_2_0_st_v1", "ensembl_gene_id")
attributes = c("affy_hugene_2_0_st_v1", "hgnc_symbol")
filters="affy_hugene_2_0_st_v1"

# affy_hc_g110
# 107                                    affy_hg_focus
# 108                                    affy_hg_u133a
# 109                                  affy_hg_u133a_2
# 110                                    affy_hg_u133b
# 111                              affy_hg_u133_plus_2
# 112                                     affy_hg_u95a
# 113                                   affy_hg_u95av2
# 114                                     affy_hg_u95b
# 115                                     affy_hg_u95c
# 116                                     affy_hg_u95d
# 117                                     affy_hg_u95e
# 118                                     affy_hta_2_0
# 119                              affy_huex_1_0_st_v2
# 120                                    affy_hugenefl
# 121                            affy_hugene_1_0_st_v1
# 122                            affy_hugene_2_0_st_v1
# 123                            affy_hugene_2_1_st_v1
# 124                                   affy_primeview
# 125                                    affy_u133_x3p


homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")
#listAttributes(homo.anno)

genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

mat <- match(rownames(expIn), genes$gene.id)
genes <- genes[mat,]


# create DGEobject 
y <- DGEList(counts=expIn, samples=samples_vector_file, genes=genes)
#y <- DGEList(counts=expIn, samples=samples_vector_file, genes=rownames(expIn))

## FOR ALL THE MKC
y <- DGEList(counts=expMKC, samples=samples_vector_file, genes=genes)


# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y, log=TRUE) 

#Normalize
y <- calcNormFactors(y)

hyp4<-read.table("/home/bianca/Desktop/microarray/Hypoxia-DMSO/microarray_rankprod_DEregulatedBio_4.txt",sep="")
mkc4<-read.table("/home/bianca/Downloads/ANALYSIS_APRIL/DE/MICROARRAY/microarray_rankprod_DE_Bio_4.txt",sep="")