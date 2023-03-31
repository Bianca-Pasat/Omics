library(mixOmics)
library(biomaRt)


# sep 22
rna<-readRDS("/home/bianca/Desktop/rna_filt_norm.RDS")
r1=rna[,c(1,2,3,7,8,9)]
rn<-as.matrix(log2(r1+1))
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
rownames(samples_info)=colnames(rn)
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
colnames(samples_info)=c("condition")
rownames(samples_info)=colnames(r1)
der_8<-read.table("/home/bianca/Desktop/mRNA/RNAseq-avelogcpm_upper/8hours/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt",sep="\t")
#der_8<-read.table("/home/bianca/Desktop/mRNA/RNAseq-avelogcpm_upper/8hours/rnaseq_rankprod_DEregulatedBio_8.txt",sep="\t")
#der_8$ensembl_gene_id<-rownames(der_8)
# stupid annotation of esnembl ids
# httr::set_config(httr::config(ssl_verifypeer = FALSE))
# homo.anno = useEnsembl(biomart = "ensembl", 
#                        dataset = "hsapiens_gene_ensembl", 
#                        mirror = "useast")
# values<-rownames(der_8)
# attributes <- c("ensembl_gene_id", "hgnc_symbol")
# filters="ensembl_gene_id"
# der_8_name <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
# der_8_1<-merge(der_8,der_8_name, by="ensembl_gene_id")
# write.table(der,"/home/bianca/Desktop/mRNA/RNAseq-avelogcpm_upper/8hours/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt",sep="\t")
#der<-rn[rownames(rn) %in% der_8_1$hgnc_symbol,]

der<-rn[rownames(rn) %in% rownames(der_8),]


## sep 22
prot<-readRDS("/home/bianca/Desktop/prot.RDS")
x_b=rownames(prot)
prot1=as.data.frame(sapply(prot, function(x) as.numeric(x)))
rownames(prot)=x_b
pr1=prot[,c(1,2,3,7,8,9)]
pr1$name=rownames(pr1)
pr1<-as.matrix(log2(pr1+1))

de_8<-read.table("/home/bianca/Desktop/proteomics_all/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
# de_8$name=rownames(de_8)
# pr1$name=rownames(pr1)
# library(dplyr)
# dep<-left_join(as.data.frame(pr1),de_8, by="name")
# dep<-na.omit(dep)
dep=pr1[rownames(pr1) %in% rownames(de_8),]
dep=dep[,-7]
#############################################################################################################################
dep<-readRDS("/home/bianca/Desktop/de_p_8_normal_values")
der<-readRDS("/home/bianca/Desktop/de_r_8_normal_values")

dep=demi8
der=mrna8

samples_info<-data.frame(condition=c(rep("DMSO",3),rep("MKC",3)),row.names = colnames(dep))



X<-t(as.matrix(der))
Y<-t(as.matrix(dep))
splsRes<-spls(X,Y)
splsRes<-spls(X,Y, keepX = c(100,10), keepY = c(100,10))
plotVar(splsRes)
plotIndiv(splsRes)
samples_info<-pheno.mrna
plotIndiv(splsRes, group=samples_info$condition, 
          rep.space = "XY-variate", legend=TRUE, legend.title = "Condition",ind.names = rownames(samples_info), 
          title="miRNA and RNAseq together")

plotIndiv(splsRes, group=samples_info$condition, 
          rep.space = "XY-variate", legend=TRUE, legend.title = "Condition",ind.names = rownames(samples_info), 
          title="Lipids and proteins together")

#plotIndiv(splsRes, group=samples_info$condition, pch=colnames(pr1),
#          rep.space = "XY-variate", legend=TRUE, legend.title = "Condition",ind.names =F,
#          title="Proteomics and RNAseq together" )

plotVar(splsRes, cex=c(3,2), legend=T)
coordinates<-plotVar(splsRes, plot=F)
coordinatesnegative<-coordinates[coordinates$x < 0 ,]


X11()
cim(splsRes, comp=1)
#cim(splsRes, comp=1, save="jpeg", name.save="sPLScim_8")
network(splsRes, comp=1, cutoff = 0.6)
#network(splsRes, comp=1, cutoff = 0.6, save='jpeg', name.save = "sPLS_network_8")

# save it for cytoscape
myNetwork<-network(splsRes, comp=1, cutoff = 0.6)$gR

X11()
plotArrow(splsRes, group=samples_info$condition, legend=T, X.label="PLS comp 1", Y.label="PLS comp 2", legend.title = "Condition")


# variable selection outputs
mySelectedVariables=selectVar(splsRes, comp=1)
mySelectedVariables$X$name
mySelectedVariables$Y$name

X11()
plotLoadings(splsRes, comp=1, size.name = rel(0.5))


## tuning parameters
perf.spls<-perf(splsRes, validation="Mfold", folds=5, progressBar = T, nrepeat=10)
plot(perf.spls, criterion='Q2.total')


#### keep features with values Q2 less than 0.0975

# keep some features based on correlation
tune.spls.cor<-tune.spls(X,Y,ncomp=2, 
                         test.keepX = c(7,12),
                         test.keepY = c(12,17),
                         validation="Mfold", 
                         folds=5, nrepeat=2,measure='cor')

plot(tune.spls.cor, measure='cor')
plot(tune.spls.cor, measure='RSS')
