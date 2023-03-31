library(miRNAmeConverter)

now<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_bianca/correlation_8_rna_WITH_BIG_FC")
result1=miRNA_MatureToPrecursor(now$miRNA)



homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "uswest")
attributes <- c("mirbase_id", "hgnc_symbol")
filters<-"mirbase_id"
values=result1$Precursor
mirnas<-getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
mirnas<-merge(mirnas,result1, by.x=1, by.y=2)
now3<-merge(now, mirnas, by.x=1, by.y=3)
now3<-unique.array(now3)

cor8bio<-rbind(as.data.frame(cbind(as.character(now3$hgnc_symbol), as.numeric(now3$logratio.miRNA.), as.numeric(now3$P.adjust.miRNA.))), as.data.frame(cbind(as.character(now3$Gene), as.numeric(now3$logratio.gene.),as.numeric(now3$P.adjust.gene.))))
write.table(unique.data.frame(cor8bio), "/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_bianca/cor8bioinfo.txt")


# test<-for (i in 1:length(now$miRNA)){
#   now$newmi[i]<-ifelse(mirnas2$OriginalName[i] %in% now$miRNA, mirnas$hgnc_symbol[i],"NA")
# }
