library(anamiR)
library(biomaRt)
library("mirTarRnaSeq")


#### read all files in a list
dataFiles <- lapply(Sys.glob("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/*.txt"), read.table)
# transfom into dataframe
df<-data.table::rbindlist(dataFiles)

#### write multiple files into a folder
for(i in 1:length(data_names)) {                              # Head of for-loop
  write.csv2(get(data_names[i]),                              # Write CSV files to folder
             paste0("/home/ben/Desktop/",
                    data_names[i],
                    ".csv"),
             row.names = FALSE)
}

############ using loops
data_files <- list.files("/home/ben/Desktop/")  # Identify file names

for(i in 1:length(data_files)) {                              # Head of for-loop
    assign(paste0("data", i),                                   # Read and store data frames
           read.csv2(paste0("/home/ben/Desktop/",
                            data_files[i])))
}



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
dem8<-dem8[abs(dem8$logratio) > 0.5,]

## integrated rnas-proteins
Ximportant<-readRDS("/home/ben/Desktop/mrna_proteins_spls_on_de8_rnasonly")
Yimportant<-readRDS("/home/ben/Desktop/mrna_proteins_spls_on_de8_proteinsonly")

dem8=demrna8[Ximportant,]
mrna8<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")
mrna8=mrna8[Ximportant,]
deprot8<-read.table("/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8_DEregulatedBio.txt", sep="\t")
colnames(deprot8)<-c("logratio","Pvalue")
deprot8$Gene=rownames(deprot8)
prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot8")
deprot8=deprot8[Yimportant,]
prot8=prot8[Yimportant,]

rnamat<-rbind(mrna8, prot8)
rnademat<-rbind(dem8,deprot8)
rnademat<-rnademat[,c(3,1,2)]

pheno.mrna<-data.frame(condition=c(rep("DMSO",3),rep("MKC",3)), row.names = colnames(mrna8) )

mrna_se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(counts=as.matrix(rnamat)),
  colData = pheno.mrna)

mrna_d <- differExp_continuous(se = mrna_se,
                               class = "condition", log2 = FALSE,
                               p_value.cutoff = 1)

####  keep important features from integration
mrna_d<-mrna_d[rownames(rnademat),]
mrna_d[,7]<-rnademat$logratio
mrna_d[,9]<-rnademat$Pvalue

# colnames(demrna8)<-demrna8[1,]
# demrna8<-demrna8[-1,]
# ro<-rownames(demrna8)
# demrna8=as.data.frame(lapply(demrna8, function(x) as.numeric(x)))
# rownames(demrna8)<-ro

demirna8<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
colnames(demirna8)<-c("logratio","Pvalue")
demirna8$Gene=rownames(demirna8)
demirna8=demirna8[,c(3,1,2)]
demirna8=demirna8[abs(demirna8$logratio) > 1,]
demirna8=demirna8[abs(demirna8$logratio) > 0.5,]

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

#### keep only FC > 0.5
mrna_d<-mrna_d[rownames(dem8),]
mrna_d[,7]<-dem8$logratio
mrna_d[,9]<-dem8$Pvalue

mrna_d[,7]<-demrna8$logratio
mrna_d[,9]<-demrna8$Pvalue

mirna_d <- differExp_continuous(se = mirna_se,
                               class = "condition", log2 = FALSE,
                               p_value.cutoff = 1)
## keep only FC >1 
mirna_d[,7]<-demirna8$logratio
mirna_d[,9]<-demirna8$Pvalue


keep<-mirna_d[,ncol(mirna_d)-5] > 1
mirna_d1<-mirna_d[keep,]

## regular
mirna_d<-mirna_d[rownames(dedemi8),]

dede8<-dedemi8[rownames(mirna_d),]
mirna_d[,7]<-dede8$logratio
mirna_d[,9]<-dede8$Pvalue

colnames(mrna_d)=c("D8R1","D8R3","D8R2","M8R2","M8R1","M8R3","log-ratio","P-Value","P-adjust","mean_case","mean_control")
colnames(mirna_d)=c("D8R1","D8R2","D8R3","M8R1","M8R2","M8R3","log-ratio","P-Value","P-adjust","mean_case","mean_control")

mirna_d <- differExp_continuous(se = mirna_se,
                              class = "condition", log2 = FALSE,
                              p_value.cutoff = 1)


corr_0 <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_d,
                    method = "pearson", cut.off = -0.5)

corr_0 <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_d,
                       method = "pearson", cut.off = -0.9)

corr_0 <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_d,
                       method = "pearson", cut.off = -0.9)

mirnaFC1_rnaFeautres<-rbind(unique(corr_0[,c(1,4,5)]),unique(corr_0[,c(2,8,9)]))
write.table(mirnaFC1_rnaFeautres,"/home/ben/Desktop/cor_mirna_rnaFeautres.txt",sep="\t",quote=F, row.names = F)

corr_0<-readRDS("/home/ben/Desktop/cor_mirnaFC05_integrated_features_spls")
corr_0<-as.data.frame(corr_0[,c(1,2,3)])
corr_1<-corr_0[,c(2,1,3)]
colnames(corr_1)<-colnames(corr_0)
corr_0<-rbind(corr_0, corr_1)
corr_0$`Correlation( -0.9 )`<-as.numeric(corr_0$`Correlation( -0.9 )`)
cor_0<-corr_0 %>% pivot_wider(names_from = Gene, values_from = `Correlation( -0.9 )`)
cor_0[is.na(cor_0)]<-0



corr_0$`Correlation( -0.9 )`<-as.numeric(corr_0$`Correlation( -0.9 )`)
corr_0a<-corr_0[corr_0$Gene %in% m1,]
corr_0ab<-corr_0[corr_0$Gene %in% m2,]
png("/home/ben/Desktop/correlated_mirnasfc05_features_spls.png", width = 1200)
ggplot(as.data.frame(corr_0a), aes(Gene, miRNA, fill= `Correlation( -0.9 )`)) + 
    geom_tile() + ggtitle("Correlated miRNAs and selected features from sPLS") 
dev.off()




library(GGally)

custom_neg_cor<-function(mrna_data,mirna_data,method = c("pearson", "kendall","spearman"), cut.off = -0.5){
    method <- match.arg(method)
    print(method)
    data_1<-mrna_data
    data_2<-mirna_data
    cal_cor <- function(data_1, data_2, cor_cut) {
      n <- 1
      corr <- list()
      for (i in seq_len(nrow(data_1))) {
        for (j in seq_len(nrow(data_2))) {
          mrna <- as.numeric(data_1[i, 1:(ncol(mrna_data) - 1)])
          mirna <- as.numeric(data_2[j, 1:(ncol(mirna_data) - 1)])
          tmp <- stats::cor(mrna, mirna, method = method)
          if (tmp < cor_cut) {
            corr[[n]] <- row.names(data_2)[j]
            corr[[n]][2] <- row.names(data_1)[i]
            corr[[n]][3] <- tmp
            corr[[n]][4] <- data_2[["logratio"]][j]
            corr[[n]][5] <- data_2[["Pvalue"]][j]
            corr[[n]][8] <- data_1[["logratio"]][i]
            corr[[n]][9] <- data_1[["Pvalue"]][i]
            n <- n + 1
          }
        }
      }
      return(corr)
    }
    corr <- cal_cor(mrna_data, mirna_data, cut.off)
    corr <- do.call(rbind, corr)
    if (is.null(corr)) {
      cut.off <- cut.off + 0.2
      corr <- cal_cor(mrna_data, mirna_data, cut.off)
      corr <- do.call(rbind, corr)
    }
    last_column <- paste("Correlation(", cut.off, ")")
    colnames(corr) <- c("miRNA", "Gene", last_column, "logratio(miRNA)", 
                        "P-adjust(miRNA)","logratio(gene)", "P-adjust(gene)")
    return(corr)
}


# test2<-function (mrna_data, mirna_data, method = c("pearson", "kendall", 
#                                             "spearman"), cut.off = -0.5) 
# {
#   method <- match.arg(method)
#   print(method)
#   common_column <- intersect(colnames(mrna_data), colnames(mirna_data))
#   mrna_data <- as.data.frame(mrna_data[, common_column])
#   mirna_data <- as.data.frame(mirna_data[, common_column])
#   cal_cor <- function(data_1, data_2, cor_cut) {
#     n <- 1
#     corr <- list()
#     for (i in seq_len(nrow(data_1))) {
#       for (j in seq_len(nrow(data_2))) {
#         mrna <- as.numeric(data_1[i, 1:(ncol(mrna_data) - 
#                                           5)])
#         mirna <- as.numeric(data_2[j, 1:(ncol(mirna_data) - 
#                                            5)])
#         tmp <- stats::cor(mrna, mirna, method = method)
#         if (tmp < cor_cut) {
#           corr[[n]] <- row.names(data_2)[j]
#           corr[[n]][2] <- row.names(data_1)[i]
#           corr[[n]][3] <- tmp
#           corr[[n]][4] <- data_2[["log-ratio"]][j]
#           corr[[n]][5] <- data_2[["P-adjust"]][j]
#           corr[[n]][6] <- data_2[["mean_case"]][j]
#           corr[[n]][7] <- data_2[["mean_control"]][j]
#           corr[[n]][8] <- data_1[["log-ratio"]][i]
#           corr[[n]][9] <- data_1[["P-adjust"]][i]
#           corr[[n]][10] <- data_1[["mean_case"]][i]
#           corr[[n]][11] <- data_1[["mean_control"]][i]
#           n <- n + 1
#         }
#       }
#     }
#     return(corr)
#   }
#   corr <- cal_cor(mrna_data, mirna_data, cut.off)
#   corr <- do.call(rbind, corr)
#   if (is.null(corr)) {
#     cut.off <- cut.off + 0.2
#     corr <- cal_cor(mrna_data, mirna_data, cut.off)
#     corr <- do.call(rbind, corr)
#   }
#   last_column <- paste("Correlation(", cut.off, ")")
#   colnames(corr) <- c("miRNA", "Gene", last_column, "logratio(miRNA)", 
#                       "P-adjust(miRNA)", "mean_case(miRNA)", "mean_control(miRNA)", 
#                       "logratio(gene)", "P-adjust(gene)", "mean_case(gene)", 
#                       "mean_control(gene)")
#   return(corr)
# }

corr_0<-custom_neg_cor(rnademat,dedemi8, "pearson",-0.9)

heat_vis(corr_0, mrna_d, mirna_d)
mirna <- corr_0[, 1]
gene <- corr_0[, 2]
mirna_da <- mirna_21[, 1:(ncol(mirna_21) - 5)]
mrna_da <- mrna_d[, 1:(ncol(mrna_d) - 5)]
mirna_exp <- mirna_da[mirna,]
mrna_exp <- mrna_da[gene, ]

mrna_exp<-as.matrix(as.numeric(rnamat))
hmcols <- rev(gplots::redgreen(100))

mirna_exp<-as.matrix(demi8)
mrna_exp<-as.matrix(rnamat)
mrna_exp=mrna_exp[(unique(corr_0[,2])),]
mirna_exp=mirna_exp[(unique(corr_0[,1])),]

library(ggpubr)
theme_set(theme_pubr())

gplots::heatmap.2(mirna_exp, trace = "none", col = hmcols, 
                  Rowv = FALSE, scale = "row", Colv = "Rowv", margins = c(9,9), dendrogram = "none")
gplots::heatmap.2(mrna_exp, trace = "none", col = hmcols, 
                  Rowv = FALSE, scale = "row", Colv = "Rowv", margins = c(9,9), dendrogram = "none")

fig<-ggarrange(m,mi, ncol=2, nrow=1)
sup <- database_support(cor_data = corr_0,
                        org = "hsa", Sum.cutoff = 3)


# saveRDS(corr_0,"/home/ben/Desktop/mirna_bianca/correlation_8_rna")
saveRDS(corr_0,"/home/ben/Desktop/mirna_bianca/correlation_8_rna_WITH_BIG_FC") #### with FC for MRNA > 0.5 AND for MIRNA > 1
corr_01<-corr_0[,c(1,2,3)]
colnames(corr_01)=c("miRNA","Gene","Correlation")
hey<-as.data.frame(corr_01) %>% pivot_wider(names_from = Gene, values_from = Correlation)

cor8 <- negative_cor(mrna_data = demrna8, mirna_data = demirna8,
                     method = "pearson", cut.off = -0.5)
d[is.na(d)] <- 0


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
