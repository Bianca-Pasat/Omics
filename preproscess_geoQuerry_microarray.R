library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library(GEOquery)


### inputs microarray
gse <- getGEO("GSE99766", GSEMatrix = TRUE)
#columns<-pData(gse[[1]])
#fData(gse[[1]])
exp<-exprs(gse[[1]])


# gse$GSE99766_series_matrix.txt.gz@phenoData@data$`treatment:ch1`
# [1] "DMSO treated, 4 hours"             "DMSO treated, 4 hours"             "DMSO treated, 4 hours"             "MKC-8866 treated, 4 hours"        
# [5] "MKC-8866 treated, 4 hours"         "MKC-8866 treated, 4 hours"         "MKC-8866 treated, 24 hours"        "MKC-8866 treated, 24 hours"       
# [9] "MKC-8866 treated, 24 hours"        "DMSO treated, 1% O2, 4 hours"      "DMSO treated, 1% O2, 4 hours"      "DMSO treated, 1% O2, 4 hours"     
# [13] "MKC-8866 treated, 1% O2, 4 hours"  "MKC-8866 treated, 1% O2, 4 hours"  "MKC-8866 treated, 1% O2, 4 hours"  "MKC-8866 treated, 1% O2, 24 hours"
# [17] "MKC-8866 treated, 1% O2, 24 hours" "MKC-8866 treated, 1% O2, 24 hours"


#### dmso  D_N_4 vs D_Y_4

expY<-as.data.frame(exp[,c(1:3,10:12)])

samples_vector_file<-data.frame(c(rep(0,3), rep(1,3)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 2))
#samples_vector_file$replicate<-c(rep(c("a","b","c"), 2))
samples_vector_file$sample <-c(paste("DMSO_N_4_", 1:3, sep=""), paste("DMSO_Y_4_",1:3, sep=""))
samples_vector_file$treatment <-c(rep("DMSO_N_4", 3), rep("DMSO_Y_4", 3))
colnames(expY)<-samples_vector_file$sample

#### dmso  D_Y_4 vs M_Y_4

expMD<-as.data.frame(exp[,10:16])


#### samples vector file for DMSO vs MKC at 4 hours
samples_vector_file<-data.frame(c(rep(0,3), rep(1,3)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 2))
#samples_vector_file$replicate<-c(rep(c("a","b","c"), 2))
samples_vector_file$sample <-c(paste("DMSO_Y_4_", 1:3, sep=""), paste("MKC_Y_4_",1:3, sep=""))
samples_vector_file$treatment <-c(rep("DMSO_Y_4", 3), rep("MKC_Y_4",3))
colnames(expIn)<-samples_vector_file$sample




#### samples vector file for DMSO vs MKC HYPOXIA at 24 hours
ex<-as.data.frame(exp[,c(10:13,17:19)])
samples_vector_file<-data.frame(c(rep(0,3), rep(1,3), rep(1,3)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 3))
samples_vector_file$sample <-c(paste("DMSO_Y_4_", 1:3, sep=""),paste("MKC_Y_24_",1:3, sep=""))
samples_vector_file$treatment <-c(rep("DMSO_Y_4", 3),  rep("MKC_Y_24",3))
colnames(ex)<-samples_vector_file$sample



# #annotation of genes
exp<-expY

httr::set_config(httr::config(ssl_verifypeer = FALSE))

homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl",
                        mirror = "asia")
# listAttributes(homo.anno)
values=rownames(expY)
attributes = c("affy_hugene_2_0_st_v1", "ensembl_gene_id")
filters="affy_hugene_2_0_st_v1"

genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)


# genetable<-data.frame(gene.id=rownames(exp))
# dup <- genes$gene.id[!duplicated(genes$affy_hugene_2_0_st_v1)]
# mat <- match(genetable$gene.id, genes$affy_hugene_2_0_st_v1)
# genes1 <- genes[mat,]

############## This computer can't establ;ish secure connection so this is why you read the files
# geneseh<-read.table("/home/ben/Downloads/genes_ens_hugo_micro.txt", sep=",")
# genes1<-read.table("/home/ben/Downloads/genes_affy_ens_micro.txt", sep=",")
# 
# 
### 4 hours HYPOXIA
expY$names<-rownames(expY)
te<-merge(expY, genes, by.x = "names", by.y=1)
te<-te[!duplicated(te$ensembl_gene_id),]
rownames(te)<-te$ensembl_gene_id

# ### 4 hours DMSO VS MKC HYPOXIA 
# expIn$names<-rownames(expIn)
# te<-merge(expIn, genes1, by.x = "names", by.y=1)
# te<-te[!duplicated(te$V2),]
# rownames(te)<-te$V2
# 
# 
# ### 24 hours DMSO VS MKC HYPOXIA 
# ex$names<-rownames(ex)
# te<-merge(ex, genes1, by.x = "names", by.y=1)
# te<-te[!duplicated(te$V2),]
# rownames(te)<-te$V2

# ### THIS WORKED BUT WE DO NOT HAVE MATRICES
# ordered_genes1<-genes1[order(genes1),]
# uog1<-as.character(unique(ordered_genes1$V1))
# ordered_expIn<-expIn[order(rownames(expIn)),]
# uoeIn<-as.character(unique(rownames(ordered_expIn)))
# now<-uoeIn[uoeIn %in% uog1]
# 
# counts<-expIn[now,]
# 
# ######### Now we have matrices. The reason we repeat ourselves is because something weird happens with %in% and reduces the rows. CHECK IT
# og1<-genes1[order(genes1$V1),]
# uog1<-og1[!duplicated(og1$V1),]
# uog1<-uog1[-40813,]
# 
# # oex<-expIn[order(rownames(expIn)),]
# # moex<-expIn[rownames(oex) %in% uog1$V1,]
# # rownames(moex)==as.character(uog1$V1)


# ##### First check 
# genes=uog1
# y <- DGEList(counts=counts, samples=samples_vector_file, genes=genes)
# 
# ######## Then keep only what you want
# 
# rownames(counts)=uog1$V2
# 
# counts1<-counts[!duplicated(rownames(counts)),]
# genes=uog1$V2


counts=expY
counts=te[,c(2:7)]
genes<-te$V2


# create DGEobject 
y <- DGEList(counts=counts, samples=samples_vector_file, genes=genes)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y, log=TRUE) 

#Normalize
y <- calcNormFactors(y)
