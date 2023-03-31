library(optparse)
library(stringr)


option_list = list(
  make_option(
    c("-m", "--demRNA"),
    action = "store",
    default = "/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/test2/*.txt",
    type = 'character',
    help = 'Path to mRNA with logFC AS SECOND and Pval as third column. For multiple files, gather all files in one path and take them all with *.txt'
  ),
  make_option(
    c("-i", "--demiRNA"),
    action = "store",
    default = "/home/ben/Desktop/mirna_bianca/test/*.txt",
    type = 'character',
    help = 'Path to miRNA with logFC and Pval as columns'
  ),
  make_option(
    c("-a", "--cleanmRNA"),
    action = "store",
    default = "/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/test2/test3/*",
    type = 'character',
    help = 'Path to RDS object. Same name as de matrix if you have multiple files'
  ),
  make_option(
    c("-b", "--cleanmiRNA"),
    action = "store",
    default = "/home/ben/Desktop/mirna_bianca/test/test2/*",
    type = 'character',
    help = 'Path to RDS object. Same name as de matrix if you have multiple files'
  ),
  make_option(
    c("-t", "--method"),
    action = "store",
    default = "pearson",
    type = 'character',
    help = 'Method to use for correlation. Available: pearson, spearman, kendall'
  ),
  make_option(
    c("-c", "--cutOff"),
    action = "store",
    default = 0.99,
    type = 'double',
    help = 'Cut off for correlation. Default to -0.9'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

############# keep name as variable
ohmyname<-function(mrna){
  whoa<-str_extract(mrna, regex("[^/]+$*[\\.]"))
  whoa2<-sub("\\.","",whoa)
  return(whoa2)
}


print("Importing")
mRNAfiles <- lapply(Sys.glob(opt$demRNA), read.table)
miRNAfiles <- lapply(Sys.glob(opt$demiRNA), read.table)
mRNAnames<-lapply(Sys.glob(opt$demRNA), ohmyname)
miRNAnames<-lapply(Sys.glob(opt$demiRNA), ohmyname)
mRNAc <- lapply(Sys.glob(opt$cleanmRNA), readRDS)
miRNAc <- lapply(Sys.glob(opt$cleanmiRNA), readRDS)
method=opt$method
cut.off=opt$cutOff




### functions #########################################################################################################

### FIRST two functions were based on negative_cor from anamiR package
cal_cor <- function(data_1, data_2, cor_cut) {
  n <- 1
  corr<-list()
  for (i in seq_len(nrow(data_1))) {
    for (j in seq_len(nrow(data_2))) {
      mrna <- as.numeric(data_1[i, 1:(ncol(data_1) - 1)])
      mirna <- as.numeric(data_2[j, 1:(ncol(data_2) - 1)])
      tmp <- stats::cor(mrna, mirna, method = method)
      if (tmp < cor_cut) {
        corr[[n]] <- row.names(data_2)[j]
        corr[[n]][2] <- row.names(data_1)[i]
        corr[[n]][3] <- tmp
        corr[[n]][4] <- data_2[["logFC"]][j]
        corr[[n]][5] <- data_2[["Pval"]][j]
        corr[[n]][6] <- data_1[["logFC"]][i]
        corr[[n]][7] <- data_1[["Pval"]][i]
        n <- n + 1
      }
    }
  }
  return(corr)
}

  
neg_cor_b<-function (mrna_data, mirna_data, method = c("pearson", "kendall","spearman"), cut.off = -0.1) {
  method <- match.arg(method)
  corr <- cal_cor(mrna_data, mirna_data, cut.off)
  corr <- do.call(rbind, corr)
  if (is.null(corr)) {
    cut.off <- cut.off + 0.2
    corr <- cal_cor(mrna_data, mirna_data, cut.off)
    corr <- do.call(rbind, corr)
  }
  last_column <- paste("Correlation(",cut.off,")")
  print(head(corr))
  colnames(corr) <- c("miRNA", "Gene", last_column, "logFC(miRNA)", 
                      "Pval(miRNA)","logFC(gene)", "Pval(gene)")
  return(corr)
}


# run analysis
whole_analysis<-function(mRNAfiles,miRNAfiles,mRNAc,miRNAc,method,cut.off){
  for(i in 1:length(mRNAfiles)){
    for (j in 1:length(miRNAfiles)){
      for(x in 1:length(mRNAc)){
        for (y in 1:length(miRNAc)){
          print("Preparing matrices")
          mRNA<-mRNAfiles[[i]]
          miRNA<-miRNAfiles[[j]]
          t1<-mRNAc[[x]]
          t2<-miRNAc[[y]]
          tt1<-t1[mRNAfiles[[i]]$V1,]
          tt2<-t2[miRNAfiles[[j]]$V1,]
          mrnaf<-cbind(tt1,mRNAfiles[[i]][,c(2,3)])
          colnames(mrnaf)=c(colnames(tt1), "logFC","Pval")
          mirnaf<-cbind(tt2,miRNAfiles[[j]][,c(2,3)])
          colnames(mirnaf)=c(colnames(tt2), "logFC","Pval")
          res[[n]]<-list(neg_cor_b(mrnaf,mirnaf,method,cut.off))
          n=n+1
          do.call(rbind,res)
          return(res)
        }
      }
    }
  }
}


#### write multiple files into a folder
writocsv<-function(res){
  for(i in 1:length(res)) {                              
    write.csv2(as.data.frame(res[[i]]),                              
               paste0("/home/ben/Desktop/",
                      i,
                      ".csv"),
               row.names = FALSE)
  }
}

writotxt<-function(res){
  for(i in 1:length(res)) {                              
    write.table(as.data.frame(res[[i]]),                              
               paste0("/home/ben/Desktop/",
                      i,
                      ".csv"),
               row.names = FALSE, quote=F, sep="\t")
  }
}
####################################################################################################################################

res=list()
print("Correlating..")
n=1
res<-whole_analysis(mRNAfiles, miRNAfiles, mRNAc, miRNAc, method, cut.off)
print("Printing results")
#writocsv(res)
writotxt(res)

