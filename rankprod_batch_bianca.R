library("sva")
library(tidyverse)
library(ggplot2)
library(preprocessCore)
library(RankProd)
library(biomaRt)
library("edgeR")
library("limma")
library("RColorBrewer")
library("stringr")
library("reshape2")
library(optparse)
library(stringr)

########################################### make options python wise ###############################################################
option_list = list(
  make_option(
    "--runRank",
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Run differential expression analysis with rankprod?'
  ),
  make_option(
    "--samples",
    action = "store",
    default = "/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/test2/detest/*.txt",
    type = 'character',
    help = 'Path to samples files, for multiple files add *.txt'
  ),
  make_option(
    "--samplesInfo",
    action = "store",
    default = "/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/test2/detest/samples_info/*.txt",
    type = 'character',
    help = 'Path to samplesInfo files, for multiple files add *.txt'
  ),
  make_option(
    "--cols",
    action = "store",
    default = "c(1,2,3,4,5,6)",
    type = 'character',
    help = 'log base, 2 or 10'
  ),
  make_option(
    "--rows",
    action = "store",
    default = "c(1,2,3,4,5,6)",
    type = 'character',
    help = 'log base, 2 or 10'
  ),
  make_option(
    "--logged",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Are values logged?'
  ),
  make_option(
    "--logbase",
    action = "store",
    default = 2,
    type = 'double',
    help = 'log base, 2 or 10'
  ),
  make_option(
    "--norm",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Are data normalized ?'
  ),
  make_option(
    "--wantgenes",
    action = "store",
    default = FALSE,
    type ='logical',
    help = 'Want specific number of genes instead of pval or pfp?'
  ),
  make_option(
    "--numgenes",
    action = "store",
    default = 2000,
    type = 'double',
    help = 'Number of genes'
  ),
  make_option(
    "--method",
    action = "store",
    default = "pval",
    type = 'character',
    help = 'Choose method between pvalue or pfp'
  ),
  make_option(
    "--cutoff",
    action = "store",
    default = 0.05,
    type = 'double',
    help = 'Cutoff of pval or pfp'
  ),
  make_option(
    "--outdir",
    action = "store",
    default = "/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/test2/detest/",
    type = 'character',
    help = 'Path to output directory'
  ),
  
  make_option(
    "--runBatch",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Run batch effect correction?'
  ),
  make_option(
    "--sva",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Sva for batch effect correction'
  ),
  make_option(
    "--svacom",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Sva and Combat or batch effect correction'
  ),
  make_option(
    "--comsva",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Combat and sva for batch effect correction'
  ),
  make_option(
    "--com",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Combat for batch effect correction'
  ),
  make_option(
    "--condition",
    action = "store",
    default = 'treatment',
    type ='character',
    help = 'Specify the condition of interest. Same as one of your column names of your samplesInfo!. If you have both drug and time make a common column for them. Design matrix is: model.matrix(~0 + condition)'
  ),
  make_option(
    "--batcheffect",
    action = "store",
    default = 'replicate',
    type = 'character',
    help = 'Specify the batch effect. Same as one of your column names of your samplesInfo!'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
#######################################################################################################################################

############################################### RANKPROD parameters ######################################################################
### # RP parameters:
cols=eval(parse(text=opt$cols)) # which cols from samples (==same as rows of sample info)
rows=eval(parse(text=opt$rows)) # which rows from samples info (==same as cols of sample info)
is_logged<- opt$logged

logbase <- opt$logbase
is_normalized <- opt$norm
data_size=opt$wantgenes
cl=c(0,1,1,1,0,0)
name="name"



if (data_size){
  num.gene=opt$numgenes 
}else{
  
  method <- opt$method
  #message("method can be 'pval', 'pfp'")
  cutoff <- opt$cutoff
}

################################################ back end functions ############################################################


# keep name as variable
ohmyname<-function(mrna){
  whoa<-str_extract(mrna, regex("[^/]+$*[\\.]"))
  whoa2<-sub("\\.","",whoa)
  return(whoa2)
}

#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  return(norm)
}


#RP function:  ###  BECAUSE RP DOES 0 VS 1 0 = condition!!!
function.rp <-function(norm_samples){
  ### TODO: here you can add parameter for which rankprod should run
  RP=RankProd::RankProducts(norm_samples, cl=cl, logged =is_logged, na.rm = TRUE, 
                            gene.names = rownames(norm_samples),
                            plot = FALSE, rand = 123, calculateProduct = TRUE, MinNumOfValidPairs = NA,
                            RandomPairs = NA, huge = FALSE)
  if (data_size){
    RP_top=RankProd::topGene(RP,gene.names=rownames(samples), num.gene = num.gene, logged = is_logged, logbase = logbase) 
    # logbase=2 is the default
  }else{
    RP_top=RankProd::topGene(RP,gene.names=rownames(norm_samples), method = method, cutoff = cutoff, logged = is_logged, logbase = logbase) 
    if (is.null(rownames(RP_top$Table1)) || is.null(rownames(RP_top$Table2)))
      RP_top1=matrix(rep("NAN",5),nrow=1,ncol=5)
    RP_top2=matrix(rep("NAN",5),nrow=1,ncol=5)
    #print(RP_top2)
    print("Rank Products did not identify differentially expressed features, either up or downregulated")
    if (!is.null(rownames(RP_top$Table1)))
      RP_top1=as.data.frame(RP_top$Table1) 
    RP_top1[,3]<- (log2(RP_top1[,3]))  ###  BECAUSE RP DOES 0 VS 1 0 = condition!!!
    colnames(RP_top2)<-colnames(RP_top1)
    if (!is.null(rownames(RP_top$Table2))){
      RP_top2=as.data.frame(RP_top$Table2)
      RP_top2[,3]<- (log2(RP_top2[,3])) ### BECAUSE RP DOES 0 VS 1 
      colnames(RP_top1)<-colnames(RP_top2)
    }
    RP_TOP<-rbind(RP_top1, RP_top2)
    return(list(RP_TOP, RP))
  }
}

# batch effect corrections

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

#sampleInfo = DATA FRAME samples treatment replicate
comsva<- function(samples,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  adjustedc <- ComBat(samples, batch=sampleInfo$replicate)
  svaob<-adjustedc - min(adjustedc)
  svobj <- svaseq(svaob, mod1, mod0) 
  comsvare <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(comsvare)
}

svacom<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(as.matrix(svaob), mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv[c(1,2)])) 
  svacomre <- ComBat(cleaned_count, batch=sampleInfo$replicate)
  return(svacomre)
}

svAll<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(svaob, mod1, mod0) 
  svAllre <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(svAllre)
}

com<- function(samples,sampleInfo,replicate){
  comre <- ComBat(samples, batch=sampleInfo$replicate)
  return(comre)
}


################################################################################################################################

print("Importing")

allsamples <- lapply(Sys.glob(opt$samples), read.table)
allsamples_info <- lapply(Sys.glob(opt$samplesInfo), read.table)
allnames<-lapply(Sys.glob(opt$samples), ohmyname)

################################################### rankprod ###################################################################
### check if this works and TODO: add all parameters here in order to make it fully dynamical
runRankProdBianca<-function(allsamples,allsamples_info,allnames){
  for(i in 1:length(allsamples)){
    for (j in 1:length(allsamples_info)){
      
      #parameters
      name=allnames[i]
      samples=as.data.frame(allsamples[i])
      samples=data.frame(lapply(samples,as.numeric), row.names = rownames(samples))
      samples=samples[,cols]
      cl=as.data.frame(allsamples_info[j])
      cl=cl[rows,1]
      out_dir=paste(opt$outdir,  "all_",allnames[i],".txt",sep="")
      upsave=paste(opt$outdir,  "upsaveBio_",name,".txt",sep="")
      downsave=paste(opt$outdir,  "downsaveBio_",name,".txt",sep="")
      allsave=paste(opt$outdir,  "allBio_",name,".txt",sep="")
      
      ### analysis
      if (!is_normalized){
        # normalize data
        print("Normalizing Data")
        norm_samples = function.norm(samples)
        # RankProducts run
        print("RankProducts run")
        RP_TOP=function.rp(norm_samples)[1]
        print="\n\ntopGene function run,\nT1: downregulated genes in condition vs control\nT2: upregulated genes in condition vs control\n\n"
        cat(print)
        # Saving outputs 
        print("Saving outputs")
        RP_TOP<-as.data.frame(RP_TOP)
        write.table(RP_TOP, out_dir,quote=FALSE, sep="\t", dec=".")
        upregulated<-RP_TOP[RP_TOP$FC..class1.class2. > 0 & RP_TOP$P.value < 0.05,]
        downregulated<-RP_TOP[RP_TOP$FC..class1.class2. < 0 & RP_TOP$P.value < 0.05,]
        write.table(upregulated[,c(3,5)], upsave, sep="\t", quote=F, row.names = rownames(upregulated))
        write.table(RP_TOP[,c(3,5)], allsave, sep="\t", quote=F, row.names = rownames(RP_TOP))
        write.table(downregulated[,c(3,5)], downsave, sep="\t", quote=F, row.names = rownames(downregulated))
        print("Worflow Finished")
      }else{
        print("RankProducts run with normalized data")
        RP_TOP=function.rp(samples)[1]
        print="\n\ntopGene function run,\nT1: downregulated genes in condition vs control\nT2: upregulated genes in condition vs control\n\n"
        cat(print)
        # Saving outputs 
        print("Saving outputs")
        RP_TOP<-as.data.frame(RP_TOP)
        write.table(RP_TOP, out_dir, quote=FALSE, sep="\t", dec=".")
        upregulated<-RP_TOP[RP_TOP$FC..class1.class2. > 0 & RP_TOP$P.value < 0.05,]
        downregulated<-RP_TOP[RP_TOP$FC..class1.class2. < 0 & RP_TOP$P.value < 0.05,]
        write.table(upregulated[,c(3,5)], upsave, sep="\t", quote=F, row.names = rownames(upregulated))
        write.table(RP_TOP[,c(3,5)], allsave, sep="\t", quote=F, row.names = rownames(RP_TOP))
        write.table(downregulated[,c(3,5)], downsave, sep="\t", quote=F, row.names = rownames(downregulated))
        print("Worflow Finished")
      }
    }
  }
}


# run rankprod 
if (opt$runRank){
  runRankProdBianca(allsamples,allsamples_info, allnames)
}

###################################################################################################################################

############################################## batch effect #######################################################################
runBatchCorrectorBianca<-function(allsamples,allsamples_info,allnames){
  for(i in 1:length(allsamples)){
    for (j in 1:length(allsamples_info)){
      
      #parameters
      #name=allnames[i]
      samples=as.data.frame(allsamples[i])
      samples=data.frame(lapply(samples,as.numeric), row.names = rownames(samples))
      sampleInfo=as.data.frame(allsamples_info[j])
      
      if(opt$sva){
        svallre<-svAll(samples, sampleInfo, opt$condition, opt$batcheffect)
        # save file 
        out_dir=paste(opt$outdir,  "cleanedWithSva_",allnames[i],".txt",sep="")
        write.table(svallre, out_dir,sep="\t", quote=FALSE)
      }
      if(opt$svacom){
        svacomre<-svacom(samples, sampleInfo, opt$condition, opt$batcheffect)
        out_dir=paste(opt$outdir,  "cleanedWithSvaCombat_",allnames[i],".txt",sep="")
        write.table(svacomre, out_dir,sep="\t", quote=FALSE)
      }
      if(opt$comsva){
        comsvare<-comsva(samples, sampleInfo, opt$condition, opt$batcheffect)
        out_dir=paste(opt$outdir,  "cleanedWithCombatSva_",allnames[i],".txt",sep="")
        write.table(comsvare, out_dir,sep="\t", quote=FALSE)
      }
      if(opt$com){
        comre<-com(samples, sampleInfo,opt$batcheffect)
        out_dir=paste(opt$outdir,  "cleanedWithCombat_",allnames[i],".txt",sep="")
        write.table(comre, out_dir,sep="\t", quote=FALSE)
      }
    }
  }
}

# run batch effect
if(opt$runBatch){
  runBatchCorrectorBianca(allsamples,allsamples_info,allnames)
}

###################################################################################################################################