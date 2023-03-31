library(biomaRt)
library(miRBaseConverter)
ensembl<-useMart("ensembl")
ensembl<-useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters="mirbase_id"
attributes=c("mirbase_id","ensembl_gene_id", "hgnc_symbol")
mirnas<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/rankprodDE_valuesfilterednorm_mirna.txt")
mirnas<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/TRUEfiltered_norm_all.txt")

mirnas$names<-rownames(mirnas)
mirnas$names<-gsub("-5p","",rownames(mirnas))
mirnas$names<-gsub("-3p","",mirnas$names)

result1=miRNA_NameToAccession(rownames(mirnas),version="v22")
result2 = getMiRNASequence(result1$Accession)

mirna_ids<-getBM(filters=filters, attributes=attributes, values=values, mart = ensembl)
mirnas$names<-gsub("R","r", mirnas$names)
mirnas_new<-merge(mirna_ids, mirnas, by.x=1, by.y=13)
dub<-duplicated(mirnas_new$hgnc_symbol)
mirnas_new[dub, "hgnc_symbol"]<-paste0(mirnas_new[dub,"hgnc_symbol"],"_",seq_along(dub)[dub])
rownames(mirnas_new)<-mirnas_new$hgnc_symbol
write.table(mirnas_new[,c(4:15)],"/home/ben/Desktop/analysis_march23/de_mirnas_hgnc.txt", sep="\t", quote=F, row.names = T)
ridd_mirnas<-c("MIMAT0019076","MIMAT0020957", "MIMAT0000080","MIMAT0005885")
# hsa-miR-1295a MIMAT0005885
# hsa-miR-24-3p MIMAT0000080
# hsa-miR-548ah-3p MIMAT0020957
# hsa-miR-548am-3p MIMAT0019076

wooo<-c("hsa-miR-548am","hsa-miR-548ah")
rownames(mirnas)<-sub("R","r",rownames(mirnas))
mirnas<-mirnas[rownames(mirnas) %in% mirna_ids$mirbase_id,]
