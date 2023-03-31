##### LOAD MATRICES THAT DON'T CONTAIN 0: REMOVED LOWLY EXPRESSED VALUES
library(omicade4)
mrna<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/rankprodDE_valuesfilterednorm_mrna.txt", sep="\t")
prot<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/rankprodDE_valuesfilterednorm_prot.txt", sep="\t")
mirna<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/rankprodDE_valuesfilterednorm_mirna.txt", sep="\t")
lipids<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/lipids_short.csv", sep="\t")

### all mirnas
mirna<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/filtered_norm_all.txt")
si_mirna<-readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/sample_info")
library(edgeR)
library(tidyverse)
library(biomaRt)
library(miRBaseConverter)
ensembl<-useMart("ensembl")
ensembl<-useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters="mirbase_id"
attributes=c("mirbase_id","ensembl_gene_id", "hgnc_symbol")
#mirnas<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_values_for_integration/rankprodDE_valuesfilterednorm_mirna.txt")
mirnas<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/TRUEfiltered_norm_all.txt")

mirnas$names<-rownames(mirnas)
mirnas$names<-gsub("-5p","",rownames(mirnas))
mirnas$names<-gsub("-3p","",mirnas$names)

result1=miRNA_NameToAccession(rownames(mirnas),version="v22")
result2 = getMiRNASequence(result1$Accession)
values=mirnas$names
mirna_ids<-getBM(filters=filters, attributes=attributes, values=values, mart = ensembl)
mirnas$names<-gsub("R","r", mirnas$names)
mirnas_new<-merge(mirna_ids, mirnas, by.x=1, by.y=13)
dub<-duplicated(mirnas_new$hgnc_symbol)
mirnas_new[dub, "hgnc_symbol"]<-paste0(mirnas_new[dub,"hgnc_symbol"],"_",seq_along(dub)[dub])
rownames(mirnas_new)<-mirnas_new$hgnc_symbol
write.table(mirnas_new[,c(4:15)],"/home/ben/Desktop/analysis_march23/filtered_mirnas_hgnc.txt", sep="\t", quote=F, row.names = T)

mirnas<-mirnas_new[,c(4:15)]
preproc<-function(counts,si,design,rpnorm){
  ### make object, filter and normalize
  counts <- mutate_all(counts, function(x) as.numeric(as.character(x)))
  y <- DGEList(counts=counts, samples=si)
  # REMOVE LOW EXPRESSED GENES
  keep.exprs <- filterByExpr(y,design=design)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  ####### normalization
  if(rpnorm){
    y$counts<-function.norm(y$counts)
  }
  else if(!rpnorm){
    y <- calcNormFactors(y, method = "upperquartile")
  }
  return(y)
}
design<-model.matrix(~0+joint, si_mirna)
mirnas<-preproc(mirna, si_mirna, design,FALSE)

mirnas <- TCGAanalyze_Filtering(mirna, method = "quantile")

l<-lipids
l$D8R2<-rowMeans(l[,c('D8R2', 'D8R3')])
l$D8R3<-rowMeans(l[,c('D8R4', 'D8R5')])
l$M8R2<-rowMeans(l[,c('M8R2', 'M8R3')])
l$M8R3<-rowMeans(l[,c('M8R4', 'M8R5')])
l$D24R2<-rowMeans(l[,c('D24R2', 'D24R3')])
l$D24R3<-rowMeans(l[,c('D24R4', 'D24R5')])
l$M24R2<-rowMeans(l[,c('M24R2', 'M24R3')])
l$M24R3<-rowMeans(l[,c('M24R4', 'M24R5')])
l<-l[,-c(4,5,9,10,14,15,19,20)]

l$names<-rownames(l)
lipid_name<-read.table("/home/ben/Desktop/lipids_names_metasyx_aitorfile.csv", sep=",")



mprot<-prot[,colnames(mrna)]
mmi<-mirna[,colnames(mrna)]
mmi<-mirnas[,colnames(mrna)]
lipids<-l[,colnames(mrna)]

mda<-list(mRNA=na.omit(mrna), miRNA=na.omit(mmi), proteins=na.omit(mprot), lipids=na.omit(lipids))
mcoin<-mcia(mda)
si<-as.data.frame(colnames(mrna),row.names = colnames(mrna))
si$condition<-gsub("R.*","", si$`colnames(mrna)`, perl=T)
plot(mcoin, axes=1:2, phenovec=si$condition, sample.lab=FALSE, df.color=1:4)
### olds
# integrated_featuresM24 <- selectVar(mcoin, a1.lim=c(0.5, Inf), a2.lim=c(-Inf, 0.5))
# integrated_featuresM24I <- selectVar(mcoin, a1.lim=c(0.5, 1.5), a2.lim=c(-1.5, 0.3))
# integrated_featuresM8 <- selectVar(mcoin, a1.lim=c(0.5, Inf), a2.lim=c(-Inf, 0.5))
# integrated_featuresM8I <- selectVar(mcoin, a1.lim=c(-2,-0.5 ), a2.lim=c(-0.2, 0.5))

plot(mcoin, axes=1:2, phenovec=si$condition, sample.lab=FALSE, df.color=1:4)
integrated_m8_mirnas<-selectVar(mcoin, a1.lim=c(-Inf, 0), a2.lim=c(-0.2, 1))
geneStat <- plotVar(mcoin, var=integrated_m8_mirnas$var, var.lab=TRUE)
plot(mcoin, axes=1:2, phenovec=si$condition, sample.lab=FALSE, df.color=1:4)
integrated_m24_mirnas<-selectVar(mcoin, a1.lim=c(0.2, Inf), a2.lim=c(-Inf, 0.2))
geneStat <- plotVar(mcoin, var=integrated_m24_mirnas$var, var.lab=TRUE)

# ints<-rbind(integrated_featuresM8, integrated_featuresM24)
# intsI<-rbind(integrated_featuresM8I, integrated_featuresM24I)
ints<-rbind(integrated_m8_mirnas, integrated_m24_mirnas)

inlipids<-ints[ints$lipids!=FALSE,]
inlipidsI<-ints[intsI$lipids!=FALSE,]
intlipidnames<-merge(inlipids,lipid_name, by.x="var", by.y="V1")
intlipidnamesI<-merge(inlipidsI,lipid_name, by.x="var", by.y="V1")

intsnolip<-ints[(ints$lipids==FALSE),]
intsnolipI<-intsI[(intsI$lipids==FALSE),]


integrated_variants<-c(intlipidnames$V2, intsnolip$var)
integrated_variantsI<-c(intlipidnamesI$V2, intsnolipI$var)
write.table(integrated_variants,"/home/ben/Desktop/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_allmirnas", sep="\t", quote=F, row.names = F)
write.table(integrated_variantsI,"/home/ben/Desktop/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_more_selected", sep="\t", quote=F, row.names = F)


##3 SANITY CHECK
geneStat <- plotVar(mcoin, var=integrated_featuresM24$var, var.lab=TRUE)
geneStat <- plotVar(mcoin, var=integrated_featuresM8$var, var.lab=TRUE)

write.table(integrated_featuresM8$var,"/home/ben/Desktop/MCIA-INTEGRATED/all_de8_integrated_mcia_features", sep="\t", quote=F, row.names = F)
write.table(integrated_featuresM24$var,"/home/ben/Desktop/MCIA-INTEGRATED/all_de24_integrated_mcia_features", sep="\t", quote=F, row.names = F)
write.table(ints$var,"/home/ben/Desktop/MCIA-INTEGRATED/all_de_integrated_mcia_features", sep="\t", quote=F, row.names = F)

save.image("/home/ben/Desktop/MCIA-INTEGRATED/mcia_on_rankprod_filtNorm_values.RData")


geneStat <- plotVar(mcoin, var=rownames(mrna), var.lab=TRUE)
geneStat1 <- plotVar(mcoin, var=rownames(mprot), var.lab=TRUE)
geneStat2 <- plotVar(mcoin, var=rownames(mmi), var.lab=TRUE)
geneStat3 <- plotVar(mcoin, var=rownames(lipids), var.lab=TRUE)
# x_b=rownames(prot)
# prot1=as.data.frame(sapply(prot, function(x) as.numeric(x)))
# rownames(prot)=x_b







