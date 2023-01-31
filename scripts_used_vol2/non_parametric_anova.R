# cal_cor <- function(data_1) {
#   n <- 1
#   corr<-list()
#   for (i in seq_len(nrow(data_1))) {
#     for (j in seq_len(nrow(data_2))) {
#       mrna <- as.numeric(data_1[i, 1:(ncol(data_1) - 1)])
#       mirna <- as.numeric(data_2[j, 1:(ncol(data_2) - 1)])
#       tmp <- stats::cor(mrna, mirna, method = method)
#       if (tmp < cor_cut) {
#         corr[[n]] <- row.names(data_2)[j]
#         corr[[n]][2] <- row.names(data_1)[i]
#         corr[[n]][3] <- tmp
#         corr[[n]][4] <- data_2[["logFC"]][j]
#         corr[[n]][5] <- data_2[["Pval"]][j]
#         corr[[n]][6] <- data_1[["logFC"]][i]
#         corr[[n]][7] <- data_1[["Pval"]][i]
#         n <- n + 1
#       }
#     }
#   }
#   return(corr)
# }

# read samples info and samples
x_info<-read.table("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/all_sample_info.txt", sep="\t")
colnames(x_info)<-x_info[1,]
x_info<-x_info[-1,]

#proteins
x_info<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/proteins_sample_info")
x_info$condition<-x_info$treatment
x_info$treatment<-gsub("_.*","",x_info$treatment)
colnames(x_info)[3]<-"hours"

# lipids
x_info<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/lipidsb/lipids_sample_info")
colnames(x_info)<-c("sample","hours","treatment","replicate")
x_info$hours<-gsub(" .*","",x_info$hours)

# miRNA
x_info<-readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_bianca/sample_info")


#x_info<-rbind(x_info,xx_info)

# xx<-as.data.frame(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_AllRna_8"))
# xx<-xx[,x_info$sample]
# x24<-as.data.frame(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_AllRna_24"))
# x24<-x24[,xx_info$sample]
# 
# x_matrix<-cbind(xx,x24)

x_matrix<-read.table("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/all_separately_cleaned.txt", sep="\t")
x_matrix<-read.table("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/rna_normalized_final.txt", sep="\t")

# proteins
x_matrix<-as.data.frame(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/removed_low_counts_proteins_QUANTILE_NORM_log2"))

# lipids
x_matrix<-as.data.frame(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/lipidsb/mstus_lipids_all"))

# mirna
x_matrix<-as.data.frame(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_bianca/filtNormAll"))

## make sure they are the same order
x_matrix<-x_matrix[,x_info$sample]


## make a list with parametes you want to check
pars<-list(hours=c(8,24,48,72),
           treatment="DMSO")
pars<-list(hours=c(8,24),
           treatment="DMSO")
## choose parameter of interest from list. then you want to do this for all parameters and store them and then put them into functions
x2<-data.frame(row.names = rownames(x_matrix))
for (i in 1:length(pars$hours)){
  col1=which(x_info$hours==pars$hours[i])
  col2=which(x_info$treatment=="DMSO")
  col3=which(x_info$treatment!="DMSO")
  x1<-x_matrix %>% mutate(meanBase=rowMeans(x_matrix[,col2]))
  x1<-sapply(x1[,col3], function(x) {x - x1$meanBase})
  x2<-cbind(x2,x1)
  return(x2)
}


# col1=which(x_info$hours==8)
# 
# col2=which(x_info$treatment=="DMSO")
# col3=which(x_info$treatment!="DMSO")

# cols=eval(parse(text=opt$cols)) 
# cols=c(1,2,3)
# cols2=c(4,5,6)

# # make function
# x1<-xx %>% select(xx[,col1]) %>% mutate(meanBase=rowMeans(xx[,col2]))
# x1<-sapply(x1[,col3], function(x) {x - x1$meanBase})
# 
# # make function
# x1<-xx %>% mutate(meanBase=rowMeans(xx[,col2]))
# x1<-sapply(x1[,col3], function(x) {x - x1$meanBase})
# 
# X3<-get_de_genes(x2)
vec=as.factor(colnames(x1))
# rna
vec=as.factor(c(rep("MKC8",3), rep("MKC24",3)))
# proteins
vec=as.factor(c( rep("MKC24",3),rep("MKC48",3), rep("MKC72",3),rep("MKC8",3)))
# lipids
vec=as.factor(c( rep("MKC6",5),rep("MKC12",5), rep("MKC24",5),rep("MKC48",5),rep("MKC72",5)))
# mirna
vec=as.factor(c(rep("MKC8",3), rep("MKC24",3)))
x2<-na.omit(x2)


# below is the get_de_genes
# tmp <- apply(dataset, 1, kruskal.test, g = factor(labels))
tmp <- apply(x2, 1, kruskal.test, g =factor(vec))
ps <- unlist(lapply(tmp, "[[", "p.value"))

ps1<-ps
ps <- p.adjust(ps,method = "BH")
ps <- p.adjust(ps,method = "BY")
ps <- p.adjust(ps)
return(ps)

# rna
de_kw<-as.data.frame(ps[(which(ps < 0.5))])

# proteins
de_kw<-as.data.frame(ps[(which(ps < 0.69))])

# lipids
de_kw<-as.data.frame(ps[(which(ps < 0.05))])

# mirna
de_kw<-as.data.frame(ps[(which(ps < 0.8))])

dim(de_kw)
ps2 <- as.data.frame(unlist(lapply(tmp, "[[", "statistic")))
first.word <- function(my.string){
  unlist(strsplit(my.string, "[.]"))[1]
}
ps2$name<-sapply(rownames(ps2), first.word)
ps2<-ps2[(ps2$name %in% rownames(de_kw)),]

new<-cbind(de_kw,ps2)
new<-new[,c(1,2)]
colnames(new)<-c("pval","Fstat")

# rna
saveRDS(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/de_kw_onFilNorm_0.5cutof")
write.table(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/de_kw_onFilNorm_0.5cutof.txt",sep="\t",quote=F, row.names=T)

# proteins
saveRDS(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/de_kw_onFilNorm_0.69cutof")
write.table(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/de_kw_onFilNorm_0.69cutof.txt",sep="\t",quote=F, row.names=T)

# lipids
saveRDS(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/lipidsb/de_kw_onFilNorm_0.05cutof")
write.table(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/lipidsb//de_kw_onFilNorm_0.05cutof.txt",sep="\t",quote=F, row.names=T)

# miRNA
saveRDS(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_biancade_kw_onFilNorm_0.8cutof")
write.table(new,"/home/ben/Desktop/OCT_22_DIF_OMICS/mirna_bianca/de_kw_onFilNorm_0.8cutof.txt",sep="\t",quote=F, row.names=T)

############################
KW<-data.frame()
## wrong 
for (i in 1:length(rownames(x2))){
 wooo<-kruskal.test(x2[i,])
 KW[i,1]<-wooo$statistic
 KW[i,2]<-wooo$p.value
 return(KW)
}

## try anova!

## what you do for mkc do for dmso as well and then do rankprod and keep de genes 
