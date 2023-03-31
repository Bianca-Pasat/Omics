library("sva")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("limma")
library("pheatmap")
library("GGally")

######################## LIPIDS ###################################################################################
# Load data and metadata 
P00994_Inspired_SampleInfo <- read.csv("~/aitor/METASYX_2.0/P00994_Inspired_SampleInfo.csv")
Raw_data_MNAR_annot <- read.csv("~/aitor/METASYX_2.0/Raw_data_MNAR_annot.csv")
Raw_data_MNAR <- Raw_data_MNAR_annot[,13:72]
row.names(Raw_data_MNAR) <- Raw_data_MNAR_annot$Peak_ID

# CELL COUNTS
probando <- P00994_Inspired_SampleInfo[order(P00994_Inspired_SampleInfo$replicate),]
p <- ggplot(probando, aes(probando$Sample_ID, probando$Cell_Number))
p + geom_bar(stat="identity", fill = probando$replicate) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# REMOVE REPEAT1
ind <- P00994_Inspired_SampleInfo$replicate != "1"
Raw_data_MNAR <- Raw_data_MNAR[,ind]
P00994_Inspired_SampleInfo <- P00994_Inspired_SampleInfo[ind,]
sampleInfo<-P00994_Inspired_SampleInfo
colnames(sampleInfo)[2]=c("sample")
colnames(sampleInfo)[12]=c("condition")
sampleInfo$treatment=paste(sampleInfo$condition, sampleInfo$Time_point, sep="_")
rm(P00994_Inspired_SampleInfo)


# MSTUS NORMALIZATION (normalized to sum of peak areas)
mstus_lipids <- sweep(Raw_data_MNAR ,2,colSums(Raw_data_MNAR)/100000000,`/`)
row.names(mstus_lipids) <- row.names(Raw_data_MNAR)
rm(Raw_data_MNAR)

# batch exploration
## on all
PreparePlotpca(mstus_lipids, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(mstus_lipids, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

comsvaRe<-comsva(mstus_lipids, sampleInfo, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(comsvaRe, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

svacomre<-svacom(as.matrix(mstus_lipids),sampleInfo,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(svacomre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

comre<-com(mstus_lipids,sampleInfo,replicate)
PreparePlotpca(comre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(comre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

savre<-svAll(as.matrix(mstus_lipids),sampleInfo,treatment,replicate)
PreparePlotpca(savre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(savre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

## time points

# 6 hours
sampleInfo6<-sampleInfo[which(sampleInfo$Time_point =="6 hour"),]
samples6<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo6$sample]

PreparePlotpca(samples6, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

comsvaRe<-comsva(samples6, sampleInfo6, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo6,sampleInfo6$condition,as.factor(sampleInfo6$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples6),sampleInfo6,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

comre<-com(samples6,sampleInfo6,replicate)
PreparePlotpca(comre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples6),sampleInfo6,treatment,replicate)
PreparePlotpca(savre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

clean6<-comsvaRe
saveRDS(clean6, "/home/ben/Desktop/lipidsb/clean6")

## 12
sampleInfo12<-sampleInfo[which(sampleInfo$Time_point =="12 hour"),]
samples12<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo12$sample]

PreparePlotpca(samples12, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

comsvaRe<-comsva(samples12, sampleInfo12, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo12,sampleInfo12$condition,as.factor(sampleInfo12$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples12),sampleInfo12,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

comre<-com(samples12,sampleInfo12,replicate)
PreparePlotpca(comre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples12),sampleInfo12,treatment,replicate)
PreparePlotpca(savre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

clean12a<-comsvaRe
clean12b<-svacomre
saveRDS(clean12b, "/home/ben/Desktop/lipidsb/clean12_svacom")
saveRDS(clean12a, "/home/ben/Desktop/lipidsb/clean12_comsva")

### 24

sampleInfo24<-sampleInfo[which(sampleInfo$Time_point =="24 hour"),]
samples24<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo24$sample]
PreparePlotpca(samples24, sampleInfo24,sampleInfo24$treatment,as.factor(sampleInfo24$replicate),values=cbp_12)

comsvaRe<-comsva(samples24, sampleInfo24, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo24,sampleInfo24$condition,as.factor(sampleInfo24$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples24),sampleInfo24,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo24,sampleInfo24$treatment,as.factor(sampleInfo24$replicate),values=cbp_12)

comre<-com(samples24,sampleInfo24,replicate)
PreparePlotpca(comre, sampleInfo24,sampleInfo24$treatment,as.factor(sampleInfo24$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples24),sampleInfo24,treatment,replicate)
PreparePlotpca(savre, sampleInfo24,sampleInfo24$treatment,as.factor(sampleInfo24$replicate),values=cbp_12)

clean24<-svacomre # *** consider removing replicate 6 !
saveRDS(clean24, "/home/ben/Desktop/lipidsb/clean24")

## 48
sampleInfo48<-sampleInfo[which(sampleInfo$Time_point =="48 hour"),]
samples48<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo48$sample]
PreparePlotpca(samples48, sampleInfo48,sampleInfo48$treatment,as.factor(sampleInfo48$replicate),values=cbp_12)

comsvaRe<-comsva(samples48, sampleInfo48, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo48,sampleInfo48$condition,as.factor(sampleInfo48$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples48),sampleInfo48,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo48,sampleInfo48$treatment,as.factor(sampleInfo48$replicate),values=cbp_12)

comre<-com(samples48,sampleInfo48,replicate)
PreparePlotpca(comre, sampleInfo48,sampleInfo48$treatment,as.factor(sampleInfo48$replicate),values=cbp_12)


savre<-svAll(as.matrix(samples48),sampleInfo48,treatment,replicate)
PreparePlotpca(savre, sampleInfo48,sampleInfo48$treatment,as.factor(sampleInfo48$replicate),values=cbp_12)


# remove replicate 3
s3=samples48[,c(1,3,4,5,6,8,9,10)]
s3i=sampleInfo48[c(1,3,4,5,6,8,9,10),]
comre<-com(s3,s3i,replicate)
svacomre<-svacom(as.matrix(s3),s3i,treatment,replicate)
PreparePlotpca(comre, s3i,s3i$treatment,as.factor(s3i$replicate),values=cbp_12)
PreparePlotpca(svacomre, s3i,s3i$treatment,as.factor(s3i$replicate),values=cbp_12)

clean48<-svacomre
saveRDS(clean48, "/home/ben/Desktop/lipidsb/clean48")

## 72
sampleInfo72<-sampleInfo[which(sampleInfo$Time_point =="72 hour"),]
samples72<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo72$sample]
PreparePlotpca(samples72, sampleInfo72,sampleInfo72$treatment,as.factor(sampleInfo72$replicate),values=cbp_12)

comsvaRe<-comsva(samples72, sampleInfo72, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo72,sampleInfo72$condition,as.factor(sampleInfo72$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples72),sampleInfo72,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo72,sampleInfo72$treatment,as.factor(sampleInfo72$replicate),values=cbp_12)

comre<-com(samples72,sampleInfo72,replicate)
PreparePlotpca(comre, sampleInfo72,sampleInfo72$treatment,as.factor(sampleInfo72$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples72),sampleInfo72,treatment,replicate)
PreparePlotpca(savre, sampleInfo72,sampleInfo72$treatment,as.factor(sampleInfo72$replicate),values=cbp_12)

clean72<-comsvaRe
saveRDS(clean72, "/home/ben/Desktop/lipidsb/clean72")

### whatif we removed bad replicates?
new<-cbind(samples6[,c(1,2,3,4,6,8,9)],samples12[,c(2,3,4,5,7,8,9,10)],samples24[,c(1,2,3,4,6,7,8,9)],samples48[,c(1,3,4,5,6,8,9,10)],samples72[,c(2,4,5,6,7,8,9,10)])
newI<-rbind(sampleInfo6[c(1,2,3,4,6,8,9),],sampleInfo12[c(2,3,4,5,7,8,9,10),],sampleInfo24[c(1,2,3,4,6,7,8,9),],sampleInfo48[c(1,3,4,5,6,8,9,10),],sampleInfo72[c(2,4,5,6,7,8,9,10),])
PreparePlotpca(new, newI,newI$treatment,as.factor(newI$replicate),values=cbp_12)

comsvaRe<-comsva(new,newI, treatment, replicate)
PreparePlotpca(comsvaRe, newI,newI$condition,as.factor(newI$replicate),values=cbp_12)
PreparePlotpca(comsvaRe, newI,newI$treatment,as.factor(newI$Time_point),values=cbp_12)

svacomre<-svacom(as.matrix(new),newI,treatment,replicate)
PreparePlotpca(svacomre, newI,newI$treatment,as.factor(newI$replicate),values=cbp_12)

clean_all<-comsvaRe
saveRDS(clean_all, "/home/ben/Desktop/lipidsb/comsvare_all_removed_replicates")

### tries to put svacom into play instead of comsva for 6 and 72 hours. RESULT = NOT BETTER
#6
# sa6<-samples6[,c(1,2,3,4,6,8,9)]
# sa6i<-sampleInfo6[c(1,2,3,4,6,8,9),]
# svacomre<-svacom(as.matrix(sa6),sa6i,treatment,replicate)
# PreparePlotpca(svacomre, sa6i,sa6i$treatment,as.factor(sa6i$replicate),values=cbp_12)
# comsvaRe<-comsva(sa6, sa6i, treatment, replicate)
# PreparePlotpca(comsvaRe, sampleInfo6,sampleInfo6$condition,as.factor(sampleInfo6$replicate),values=cbp_12)
# 
# # 72
# sa72<-samples72[,c(2,4,5,6,7,8,9,10)]
# sa72i<-sampleInfo72[c(2,4,5,6,7,8,9,10),]
# svacomre<-svacom(as.matrix(sa72),sa72i,treatment,replicate)
# PreparePlotpca(svacomre, sa72i,sa72i$treatment,as.factor(sa72i$replicate),values=cbp_12)
# comsvaRe<-comsva(sa72, sa72i, treatment, replicate)
# PreparePlotpca(comsvaRe, sampleInfo72,sampleInfo72$condition,as.factor(sampleInfo72$replicate),values=cbp_12)
# 
# tries for a better 12 hour removed asmples 2,6 from D AND M
s12<-samples12[,c(2,3,4,7,8,9)]
s12i<-sampleInfo12[c(2,3,4,7,8,9),]
PreparePlotpca(samples12, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

comsvaRe<-comsva(s12, s12i, treatment, replicate)
PreparePlotpca(comsvaRe, s12i,s12i$condition,as.factor(s12i$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(s12),s12i,treatment,replicate)
PreparePlotpca(svacomre, s12i,s12i$treatment,as.factor(s12i$replicate),values=cbp_12)

clean12_ultra<-svacomre
saveRDS(clean12_ultra, "/home/ben/Desktop/lipidsb/clean12_svacom_removed2_6_rep")



################## DE WITH RANKPROD ON TIME POINTS #################################
### # RP parameters:

is_logged <- as.logical("TRUE")
logbase <- 2
is_normalized <- as.logical("TRUE")
data_size=as.logical("FALSE")


# 6 HOURS
samples=clean6
cl=c(rep(1,5),rep(0,5))

out_dir="/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_6.txt"
upsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_6_upregulatedBio.txt"
downsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_6_downregulatedBio.txt"
allsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_6_DEregulatedBio.txt"

# 12
samples=clean12_ultra
cl=c(rep(1,3),rep(0,3))
out_dir="/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_12.txt"
upsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_12_upregulatedBio.txt"
downsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_12_downregulatedBio.txt"
allsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_12_DEregulatedBio.txt"

# 24
samples=clean24
cl=c(rep(1,5),rep(0,5))
out_dir="/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_24.txt"
upsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_24_upregulatedBio.txt"
downsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_24_downregulatedBio.txt"
allsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_24_DEregulatedBio.txt"

# 48
samples=clean48
cl=c(rep(1,4),rep(0,4))
out_dir="/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_48.txt"
upsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_48_upregulatedBio.txt"
downsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_48_downregulatedBio.txt"
allsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_48_DEregulatedBio.txt"

# 72
samples=clean72
cl=c(rep(1,5),rep(0,5))
out_dir="/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_72.txt"
upsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_72_upregulatedBio.txt"
downsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_72_downregulatedBio.txt"
allsave= "/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_72_DEregulatedBio.txt"

### run analysis
if (!is_normalized){
  # normalize data
  print("Normalizing Data")
  function.norm <- function(counts){
    norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
    rownames(norm)=rownames(counts)
    colnames(norm)=colnames(counts)
    return(norm)
  }
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
  print("Worflow Finished")
  print("Saving results")
  upregulated<-RP_TOP[RP_TOP$FC..class1.class2. > 0 & RP_TOP$P.value < 0.05,]
  downregulated<-RP_TOP[RP_TOP$FC..class1.class2. < 0 & RP_TOP$P.value < 0.05,]
  write.table(upregulated[,c(3,5)], upsave, sep="\t", quote=F, row.names = rownames(upregulated))
  write.table(RP_TOP[,c(3,5)], allsave, sep="\t", quote=F, row.names = rownames(RP_TOP))
  write.table(downregulated[,c(3,5)], downsave, sep="\t", quote=F, row.names = rownames(downregulated))
}


#### LOAD DE AND SUBSTRACT FROM EACH TIME POINT!
de6<-read.table("/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_6.txt",sep="\t")
de6<-rownames(de6)
de12<-read.table("/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_12.txt",sep="\t")
de12<-rownames(de12)
de24<-read.table("/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_24.txt",sep="\t")
de24<-rownames(de24)
de48<-read.table("/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_48.txt",sep="\t")
de48<-rownames(de48)
de72<-read.table("/home/ben/Desktop/lipidsb/rp_lip_mkc_vs_dmso_72.txt",sep="\t")
de72<-rownames(de72)

deAll<-c(de6,de12,de24,de48,de72)

lipids6<-clean6[rownames(clean6) %in% de6,]
lipids12<-clean12_ultra[rownames(clean12_ultra) %in% de12,]
lipids24<-clean24[rownames(clean24) %in% de24,]
lipids48<-clean48[rownames(clean48) %in% de48,]
lipids72<-clean72[rownames(clean72) %in% de72,]

saveRDS(lipids6,"/home/ben/Desktop/lipidsb/lipids6")
saveRDS(lipids12,"/home/ben/Desktop/lipidsb/lipids12")
saveRDS(lipids24,"/home/ben/Desktop/lipidsb/lipids24")
saveRDS(lipids48,"/home/ben/Desktop/lipidsb/lipids48")
saveRDS(lipids72,"/home/ben/Desktop/lipidsb/lipids72")

############# keep de lipids with original mstus matrix
mstus6<-mstus_lipids[rownames(mstus_lipids) %in% de6,]
mstus6=mstus6[,c(1:12)]
mstus12<-mstus_lipids[rownames(mstus_lipids) %in% de12,]
mstus12=mstus12[,c(13:24)]
mstus24<-mstus_lipids[rownames(mstus_lipids) %in% de24,]
mstus24=mstus24[,c(25:36)]
mstus48<-mstus_lipids[rownames(mstus_lipids) %in% de48,]
mstus48=mstus48[,c(37:48)]
mstus72<-mstus_lipids[rownames(mstus_lipids) %in% de72,]
mstus72=mstus72[,c(49:60)]

saveRDS(mstus6,"/home/ben/Desktop/lipidsb/mstus6")
saveRDS(mstus12,"/home/ben/Desktop/lipidsb/mstus12")
saveRDS(mstus24,"/home/ben/Desktop/lipidsb/mstus24")
saveRDS(mstus48,"/home/ben/Desktop/lipidsb/mstus48")
saveRDS(mstus72,"/home/ben/Desktop/lipidsb/mstus72")

### ANNOTATE LIPIDS : WHAT CLSAS THEY ARE IN AND WITH LIPIDR!!

## MOVE TO INTEGRATION

######################################## PROBABLY NO
### LOG2 TRANSFORMATION

log2_mstus_lipids<-log2(mstus_lipids)

# FILTRATION OF LOW COUNTS
# row has at least two columns with a count of over 3
# keep=apply(Y, 1, function(x) { length(x[x>3])>=2 } )
# Y= Y[keep,]

keep=apply(as.matrix(log2_mstus_lipids), 1, function(x) { length(x[x>5])>=25 } )
Y= log2_mstus_lipids[keep,]
#############################################################
# FILTRATION OF LOW COUNTS
# row has at least two columns with a count of over 3
# keep=apply(Y, 1, function(x) { length(x[x>3])>=2 } )
# Y= Y[keep,]

### LOG2 TRANSFORMATION

## QUANTILE NORMALIZATION (or whichever)

### SVA CORRECTION
## first you investigate and then add sv into design. Then you can correct with removeBatchEffect or ComBat
# groups <- as.factor(rep(c("DMSO","MKC"), each=3))
# test_model <- model.matrix(~groups)
# null_model <- test_model[,1]
# svar <- svaseq(arab_filtered, test_model, null_model, n.sv=1)
# design <- cbind(test_model, svar$sv) THIS DEISGN CAN EITHER BE USED DIRECTLY TO LIMMA OR WITH REMOVEBATCHEFFECT

### DIFFERENTIAL EXPRESSION