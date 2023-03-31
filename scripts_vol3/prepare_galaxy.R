# for galaxy

mrna8<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_AllRna_8")
mrna8<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_AllRna_24")

hey<-grep("^MIR", rownames(mrna8))
mirna8<-mrna8[hey,]
hey2<-grep("^LINC*", rownames(mrna8))
lncrna8<-mrna8[hey2,]

mrna88<-mrna8[!(rownames(mrna8) %in% rownames(mirna8)),]
mrna888<-mrna88[!(rownames(mrna88) %in% rownames(lncrna8)),]


## biomart
mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "uswest")
attributes = c("ensembl_gene_id", "hgnc_symbol")
filters="hgnc_symbol"


# #listAttributes(homo.anno)
values=rownames(mirna8)
mirna8_ens_tab <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
mirna88<-mirna8[mirna8_ens_tab$hgnc_symbol,]
rownames(mirna88)==mirna8_ens_tab$hgnc_symbol
rownames(mirna88)=mirna8_ens_tab$ensembl_gene_id
write.table(mirna88, "/home/ben/Desktop/galaxy_mirna8.txt", sep="\t", quote=F)
write.table(mirna88, "/home/ben/Desktop/galaxy_mirna24.txt", sep="\t", quote=F)


values=rownames(lncrna8)
lncrna8_ens_tab <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
lncrna88<-lncrna8[lncrna8_ens_tab$hgnc_symbol,]
rownames(lncrna88)==lncrna8_ens_tab$hgnc_symbol
rownames(lncrna88)=lncrna8_ens_tab$ensembl_gene_id
write.table(lncrna88, "/home/ben/Desktop/galaxy_lncrna8.txt", sep="\t", quote=F)
write.table(lncrna88, "/home/ben/Desktop/galaxy_lncrna24.txt", sep="\t", quote=F)

values=rownames(mrna888)
mrna8_ens_tab <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
mrna188<-mrna888[mrna8_ens_tab$hgnc_symbol,]
rownames(mrna188)==mrna8_ens_tab$hgnc_symbol
rownames(mrna188)=mrna8_ens_tab$ensembl_gene_id
write.table(mrna188, "/home/ben/Desktop/galaxy_mrna8_noNc.txt", sep="\t", quote=F)
write.table(mrna188, "/home/ben/Desktop/galaxy_mrna24_noNc.txt", sep="\t", quote=F)

