## writer bi
library(GEOquery)
library(biomaRt)

## get the data
gse <- getGEO("GSE47959", GSEMatrix = TRUE)
#columns<-pData(gse[[1]])
#fData(gse[[1]])
exp<-exprs(gse[[1]])

## annotation
mouse <- useEnsembl("ensembl",dataset = "mmusculus_gene_ensembl")
#listAttributes(mouse) 
values=rownames(exp)
attributes = c("affy_mouse430_2", "ensembl_gene_id")
filters="affy_mouse430_2"
genes <- getBM(attributes = attributes, filters=filters,values = values, mart = mouse)

genetable<-data.frame(gene.id=rownames(exp))
mat <- match(genetable$gene.id, genes$affy_mouse430_2)
genes1 <- genes[mat,]
rownames(exp)<-genes1$ensembl_gene_id

### make right column names
colnames(exp)<-c("CONDI","cond2")
# choose columns (not necessary)
#new<-exp[,c(1,2,5,6)]

## save table to load to galaxy
write.table(exp,"/home/bianca/Downloads/part1/table1.txt",sep="\t",quote=FALSE, row.names = FALSE)
