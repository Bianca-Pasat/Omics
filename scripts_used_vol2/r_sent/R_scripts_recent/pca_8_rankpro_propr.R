library("sva")
library(tidyverse)
library(ggplot2)


### proteins full and info
X<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final.txt",sep="\t")
X_info<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final_info_correct_order.txt",sep="\t")
colnames(X_info)=X_info[1,]
X_info=X_info[-1,]
X_info$sample=colnames(X)
x_b=X
X=as.data.frame(lapply(X, function(x) as.numeric(x)))
rownames(X)=rownames(x_b)

# ## exclude low counts
# keep=apply(X,1,function(x) sum(x>=0) >= 3)
# # 
# # keep=apply(X,1,function(x) sum(x>=15) >= 3)
# 
# X1=X[keep,]


## another exculsion method
X1=X[!(rowSums(X==0)==2),]

## rankprod differentially expressed
de_prot_8<-read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
trans<-read.table("/home/bianca/Desktop/transcripts_8.txt",sep="\t")


a=rownames(de_prot_8)
a=c(rownames(de_prot_8), trans$V1)

de_8_ex=X1[a,]
de_8_ex=de_8_ex[,c(10:12,22:24)]
# 
# all_8=X[,c(10:12,22:24)]
# X_info=X_info[c(10:12,22:24),]
# 
# keep=apply(de_8_ex,1,function(x) sum(x>=15) >=2)


de_8_ex1=de_8_ex[!(rowSums(de_8_ex==0)),]
de_8_ex1=de_8_ex[!(rowSums(de_8_ex==0)>=2),]

b=rownames(de_8_ex1)
wo<-de_prot_8[rownames(de_prot_8) %in% b,]
write.table(wo,"/home/bianca/Desktop/rankprod_propr_8_no0s.txt",sep="\t",quote=F, row.names = F)

trans_cts=as.data.frame(na.omit(de_8_ex1))
trans_cts=na.omit(trans_cts)
trans_cts=na.omit(de_8_ex)



#8 hours
de_8_ex=de_8_ex + 1
X3=as.data.frame(lapply(de_8_ex, function(x) as.numeric(x)))
rownames(X3)=rownames(de_8_ex)

trans_cts=as.data.frame(X3)

trans_cts=de_prot_8
X_info=X_info[c(10:12,22:24),]
sample_info=X_info
sample_info$sample=sample_info$short_name


### try
X1=X1[,c(10:12,22:24)]
X_info=X_info[c(10:12,22:24),]

trans_cts=as.data.frame(all_8)

trans_cts=as.data.frame(X1)
sample_info=X_info
sample_info$sample=sample_info$short_name



#### heatmap of common 

de8=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
de24=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_24_DEregulatedBio.txt",sep="\t")
de48=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_48_DEregulatedBio.txt",sep="\t")
de72=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_72_DEregulatedBio.txt",sep="\t")

de8$name=rownames(de8)
de24$name=rownames(de24)
de48$name=rownames(de48)
de72$name=rownames(de72)


de_all1=merge(de8,de24, by=3)
de_all2=merge(de_all1,de48, by="name")
de_all=merge(de_all2,de72, by="name")

de_all=de_all[,c(1,2,4,6,8)]
rownames(de_all)=de_all$gene_name
de_all=de_all[,-1]
colnames(de_all)=c("FC_08", "FC_24", "FC_48", "FC_72")
colnames(de_all)=c(8,24,48,72)


### technical replicates
library(limma)
pro=read.table("/home/bianca/Desktop/proteins_inputs/protein_matrix_gene_names_zero.txt", sep=",")
colnames(pro)=pro[1,]
pro=pro[-1,]
cl<-colnames(pro)
pro<-as.data.frame(lapply(pro, function(x) as.numeric(x)))
colnames(pro)=cl
# pro<-pro[,c(1:48)]


plotScatter <- function(df, bait, title = ''){
  
  # inital checks of input
  stopifnot(any(grepl('rep', colnames(df))))
  colRep <- as.vector(grepl('rep', colnames(df)) & unlist(lapply(df, is.numeric)))
  reps = regmatches(colnames(df), regexpr('rep[0-9]',colnames(df)))
  nRep <- sum(as.numeric(colRep))
  
  require(ggplot2)
  require(ggrepel)
  
  for (i in 1:(nRep-1)) {
    for (j in (i+1):nRep) {
      
      # set columns of replicates
      col1 <- reps[i]#paste("rep",i,sep="")
      col2 <- reps[j]#paste("rep",j,sep="")
      r <- cor(df[,col1],df[,col2]) # Pearson correlation
      temp_df <- data.frame(gene=df$gene,temp_rep1=df[,col1],temp_rep2=df[,col2],significant=df$significant)
      
      # start scatterplot
      p <- ggplot(temp_df, aes(x=temp_rep1, y=temp_rep2)) +
        
        # plot all proteins (green = significant, blue = not significant)
        geom_point(alpha=0.5, size=1.5, color=ifelse(df$significant, "springgreen3", "royalblue2")) +
        
        # label bait (red = signficant, orange = not significant)
        geom_point(subset(temp_df, gene==bait & significant), mapping=aes(x=temp_rep1, y=temp_rep2),size=2, color="red") + 
        geom_point(subset(temp_df, gene==bait & !significant), mapping=aes(x=temp_rep1, y=temp_rep2),size=2, color="orange") +
        geom_point(subset(temp_df, gene==bait), mapping=aes(x=temp_rep1, y=temp_rep2), size=2, color="black", shape=1) +	
        geom_text_repel(subset(temp_df, gene==bait), mapping=aes(label=gene),
                        arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.15, "lines"),
                        point.padding=unit(0.2, "lines"), color="black", size=3) +
        
        # identity line, title (with correlation), theme
        geom_abline(intercept=0, slope=1, linetype="longdash", size=0.2) +
        labs(title = title, subtitle = paste("correlation:",format(r,digits=3))) + xlab(col1) + ylab(col2) +
        theme_bw() + theme(axis.line=element_line(color="grey")) + ggstamp()
      
      print(p)
      
    }
  }
  return(list(r=r))
}

de=data.frame(unique(colnames(pro)))

duplicateCorrelation(pro)

