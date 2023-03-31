library(anamiR)
library(biomaRt)
library("mirTarRnaSeq")

# mrna8<-read.table("/home/ben/Downloads/rnaseq_rankprod_DEregulatedBio_8.txt")
# genes <- getBM(attributes = attributes, filters=filters,values = rownames(mrna8), mart = homo.anno)
# mrna8$ens=rownames(mrna8)
# mrna8<-merge(genes, mrna8, by.x=1,by.y=3)
# mrna8<-mrna8[,-1]
# colnames(mrna8)<-c("HGNC","logratio","Pvalue")
# 
# library(plyr)
# exp2<-ddply(mrna8, .(HGNC), numcolwise(mean))
# rownames(exp2)=exp2$HGNC
# mrna8<-exp2[,-1]
# write.table(mrna8,"/home/ben/Desktop/rna_bianca/rnaseq_rankprod_DEregulatedBio_8.txt",sep="\t",quote=F, row.names = T)
demrna8<-read.table("/home/ben/Desktop/rna_bianca/rnaseq_rankprod_DEregulatedBio_8.txt",sep="\t")
demrna8=demrna8[-1,]
demrna8$Gene=rownames(demrna8)
dem8=demrna8[,c(3,1,2)]


# colnames(demrna8)<-demrna8[1,]
# demrna8<-demrna8[-1,]
# ro<-rownames(demrna8)
# demrna8=as.data.frame(lapply(demrna8, function(x) as.numeric(x)))
# rownames(demrna8)<-ro

demirna8<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
colnames(demirna8)<-c("log-ratio","Pvalue")
demirna8$Gene=rownames(demirna8)
demirna8=demirna8[,c(3,1,2)]

mirna8<-readRDS("/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_8")
mrna8<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")

demi8<-mirna8[rownames(mirna8) %in% demirna8$Gene, ]
dedemi8<-demirna8[demirna8$Gene %in% rownames(mirna8) , ]

decor<-cor(t(dem8$logratio), t(demi8$logratio))

mirna8<-mirna8 + abs(min(mirna8))
pheno.mrna<-data.frame(condition=c(rep("DMSO",3),rep("MKC",3)), row.names = colnames(mrna8) )
pheno.mirna<-data.frame(condition=c(rep("DMSO",3),rep("MKC",3)),row.names = colnames(mirna8))

mrna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=mrna8),
  colData = pheno.mrna)
mrna_se2<-mrna_se

mirna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=as.matrix(demi8)),
  colData = pheno.mirna
)

# mrna_d <- differExp_discrete(se = mrna_se,
#                              class = "DMSO", method = "limma",
#                              t_test.var = FALSE, log2 =FALSE,
#                              p_value.cutoff = 1,  logratio = 0.0
# )

##
mrna_d <- differExp_continuous(se = mrna_se,
                                class = "condition", log2 = FALSE,
                                p_value.cutoff = 1)

mrna_d[,7]<-demrna8$logratio
mrna_d[,9]<-demrna8$Pvalue

mirna_d <- differExp_continuous(se = mirna_se,
                               class = "condition", log2 = FALSE,
                               p_value.cutoff = 1)
mirna_d<-mirna_d[rownames(dedemi8),]
dede8<-dedemi8[rownames(mirna_d),]
mirna_d[,7]<-dede8$`log-ratio`
mirna_d[,9]<-dede8$Pvalue

colnames(mrna_d)=c("D8R1","D8R3","D8R2","M8R2","M8R1","M8R3","log-ratio","P-Value","P-adjust","mean_case","mean_control")
colnames(mirna_d)=c("D8R1","D8R2","D8R3","M8R1","M8R2","M8R3","log-ratio","P-Value","P-adjust","mean_case","mean_control")

mirna_d <- differExp_continuous(se = mirna_se,
                              class = "condition", log2 = FALSE,
                              p_value.cutoff = 1)

cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_d,
                    method = "pearson", cut.off = -0.5)

cor8 <- negative_cor(mrna_data = demrna8, mirna_data = demirna8,
                     method = "pearson", cut.off = -0.5)
corr=list()
cor_cut=1
heydy<-function(data_1, data_2, cor_cut) {
  n <- 1
  corr <- list()
  for (i in seq_len(nrow(data_1))) {
    for (j in seq_len(nrow(data_2))) {
      mrna <- as.numeric(data_1[i, 1:(ncol(mrna_data) - 5)])
      mirna <- as.numeric(data_2[j, 1:(ncol(mirna_data) - 5)])
      tmp <- stats::cor(mrna, mirna, method = method)
      if (tmp < cor_cut) {
        corr[[n]] <- row.names(data_2)[j]
        corr[[n]][2] <- row.names(data_1)[i]
        corr[[n]][3] <- tmp
        corr[[n]][4] <- data_2[["logratio"]][j]
        corr[[n]][5] <- data_2[["Pvalue"]][j]
        corr[[n]][6] <- data_1[["logratio"]][i]
        corr[[n]][7] <- data_1[["Pvalue"]][i]
        n <- n + 1
      }
    }
  }
  return(corr)
}



cor8 <- negative_cor(mrna_data = t(demrna8[,c(1,2)]), mirna_data = t(demirna8[,c(1,2)]),
                    method = "pearson", cut.off = -0.5)

demrna8$Gene<-rownames(demrna8)
m8<-as.data.frame(demrna8)
m8=m8[-1,]
lm8<-list(m8)
mrna <- one2OneRnaMiRNA(lm8, gene_colname = "Gene", fc_colname = "logratio", pval_colname = "Pvalue")
demirna8$miRNA<-rownames(demirna8)
mi8<-as.data.frame(demirna8)
mi8=mi8[-1,]
lmi8<-list(mi8)
lmi8<-list(as.data.frame(demirna8))
mirna <- one2OneRnaMiRNA(lmi8, gene_colname = "miRNA", fc_colname = "logratio", pval_colname = "Pvalue")

ldem8<-list(dem8)
ldemi8<-list(demi8)
mrna <- one2OneRnaMiRNA(ldem8, gene_colname = "Gene", fc_colname = "logratio", pval_colname = "Pvalue")$foldchanges
mirna <- one2OneRnaMiRNA(ldemi8, gene_colname = "Gene", fc_colname = "logratio", pval_colname = "Pvalue")$fodchanges


mrna$pvalues<-demrna8$Pvalue
mirna$pvalues<-demirna8$Pvalue

corr_0 <- corMirnaRna(demrna8$logratio,demirna8$logratio,method="pearson")
corr_0 <- corMirnaRna(mrna$foldchanges,mirna$foldchanges,method="pearson")
corr_0 <- corMirnaRna(mrna$FC1,mirna$FC1,method="pearson")
corr_0 <- corMirnaRna(dem8$logratio,demi8$logratio,method="pearson")
corr_0 <- corMirnaRna(mrna,mirna,method="pearson")

outs <- sampCorRnaMirna(mrna$foldchanges, mirna$foldchanges,method="pearson",
                        Shrounds = 100, Srounds = 1000)



# Estimate the miRNA mRNA FC differences for your dataset
inter0 <- twoTimePoint(mrna, mirna)

#Make a background distribution for your miRNA mRNA FC differences
outs <- twoTimePointSamp(mrna, mirna,Shrounds = 10 )

#Import concordant miRanda file
miRanda <- getInputSpecies("Human", threshold = 140)

# 3.6 Plot density plots
# Density plot for background and corrs in our data. Note grey is the background distribution
# and red is the actual data.
# #Draw density plot
mirRnaDensityCor(corr_0, outs)



################# STACKOVERFLOW LALALALA 
# You should start by merging the two data frames:
#   

combined <- merge(Cali_Income, Cali_Asthma_Rates,
                         by = `County Name`,
                         suffixes = c(".Income", ".Asthma_Rate"))
# head(Cali_combined)
# (note: fix the by to be your header for "County Name" if there is an underscore or something that I didn't see)

#Then, you can do your correlation on the pairwise complete observations (for example)

with(Cali_combined, 
  cor(Rank.Income, Rank.Asthma_Rate,
    use = "pairwise.complete.obs",
    method = "spearman")
)


res <- cor.test(demrna8$logratio, demirna8$FC..class1.class2., 
                method = "kendal")

