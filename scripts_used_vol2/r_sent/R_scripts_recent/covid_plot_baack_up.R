## COVID-19 genome mutation annotator ###
##########################################
## By Federico M. Giorgi


#################################
###### BASH section
### Download SARS-CoV-2 sequences in FASTA format
# input=input.fasta
### Run nucmer to obtain variant file
# ref=NC_045512.2.fa # The reference SARS-CoV-2 Wuhan Genome
# dos2unix $input
# nucmer --forward -p nucmer $ref $input
# show-coords -r -c -l nucmer.delta > nucmer.coords
# show-snps nucmer.delta -T -l > nucmer.snps



#################################
###### R section
library(seqinr)
library(Biostrings)

# Load variant list
nucmer<-read.delim("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/after_june_numcer/nucmer.snps",as.is=TRUE,skip=4,header=FALSE)

##until june 1
nucmer<-read.delim("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/until_june1/nucmer.snps",as.is=TRUE,skip=4,header=FALSE)
nucmer1<-nucmer[1:(1/2*(length(rownames(nucmer)))),]
nucmer1<-nucmer[(1/2*(length(rownames(nucmer)))):(3/4*(length(rownames(nucmer)))),]
nucmer1<-nucmer[(3/4*(length(rownames(nucmer)))):(length(rownames(nucmer))),]
nucmer=nucmer1

##until june 2
nucmer<-read.delim("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/until_june2/formatted/until_june2.snps",as.is=TRUE,skip=4,header=FALSE)

colnames(nucmer)<-c("rpos","rvar","qvar","qpos","","","","","rlength","qlength","","","rname","qname")
rownames(nucmer)<-paste0("var",1:nrow(nucmer))
nrow(nucmer) # 69989

# Fix IUPAC codes
table(nucmer$qvar)
nucmer<-nucmer[!nucmer$qvar%in%c("B","D","H","K","M","N","R","S","V","W","Y"),]
nrow(nucmer) # 69980


### Aminoacid variant list ----
# Load reference sequence

refseq<-read.fasta("/home/bianca/Downloads/ml_covid/first_inputs/wuhan.fasta",forceDNAtolower=FALSE)[[1]]

# Load GFF3
gff3<-read.delim("/home/bianca/Downloads/ml_covid/covid_wuhan.gff3",as.is=TRUE,skip=2,header=FALSE)
annot<-setNames(gff3[,10],gff3[,9])

### Merge neighboring events ----
samples<-unique(nucmer$qname)
length(samples) # 9884
pb<-txtProgressBar(0,length(samples),style=3)
for (pbi in 1:length(samples)){ # This will update the nucmer object
  sample<-samples[pbi]
  allvars<-nucmer[nucmer$qname==sample,]
  snps<-allvars[(allvars[,"rvar"]!=".")&(allvars[,"qvar"]!="."),]
  inss<-allvars[(allvars[,"rvar"]=="."),]
  dels<-allvars[(allvars[,"qvar"]=="."),]
  # Merge insertions
  prevqpos<-0
  prevrowname<-NULL
  remove<-c()
  i<-1
  corrector<-0
  while(i<=nrow(inss)){
    rpos<-inss[i,"rpos"]
    rvar<-inss[i,"rvar"]
    qvar<-inss[i,"qvar"]
    qpos<-inss[i,"qpos"]
    if((qpos!=1)&(qpos==(prevqpos+1+corrector))){
      inss<-inss[-i,]
      inss[prevrowname,"qvar"]<-paste0(inss[prevrowname,"qvar"],qvar)
      corrector<-corrector+1
      i<-i-1
    } else {
      corrector<-0
      prevrowname<-rownames(inss)[i]
      prevqpos<-qpos
    }
    i<-i+1
  }
  # Merge deletions
  prevqpos<-0
  prevrowname<-NULL
  remove<-c()
  i<-1
  while(i<=nrow(dels)){
    rpos<-dels[i,"rpos"]
    rvar<-dels[i,"rvar"]
    qvar<-dels[i,"qvar"]
    qpos<-dels[i,"qpos"]
    
    if((qpos!=1)&(qpos==(prevqpos))){
      dels<-dels[-i,]
      dels[prevrowname,"rvar"]<-paste0(dels[prevrowname,"rvar"],rvar)
      i<-i-1
    } else {
      prevrowname<-rownames(dels)[i]
      prevqpos<-qpos
    }
    i<-i+1
  }
  # Merge SNPs
  prevqpos<-0
  prevrowname<-NULL
  remove<-c()
  i<-1
  corrector<-0
  while(i<=nrow(snps)){
    rpos<-snps[i,"rpos"]
    rvar<-snps[i,"rvar"]
    qvar<-snps[i,"qvar"]
    qpos<-snps[i,"qpos"]
    
    if((qpos!=1)&(qpos==(prevqpos+1+corrector))){
      snps<-snps[-i,]
      snps[prevrowname,"rvar"]<-paste0(snps[prevrowname,"rvar"],rvar)
      snps[prevrowname,"qvar"]<-paste0(snps[prevrowname,"qvar"],qvar)
      corrector<-corrector+1
      i<-i-1
    } else {
      corrector<-0
      prevrowname<-rownames(snps)[i]
      prevqpos<-qpos
    }
    i<-i+1
  }
  
  # Remerge back
  allvars2<-rbind(snps,inss,dels)
  remove<-setdiff(rownames(allvars),rownames(allvars2))
  nucmer<-nucmer[setdiff(rownames(nucmer),remove),]
  nucmer[rownames(allvars2),]<-allvars2
  setTxtProgressBar(pb,pbi)
}


### Provide effect of each SNP and indel ----
header<-c("sample","refpos","refvar","qvar","qpos","qlength","protein","variant","varclass","annotation")
results<-matrix(NA,ncol=length(header),nrow=0)
colnames(results)<-header

refirst<-results

samples<-unique(nucmer$qname)
pb<-txtProgressBar(0,length(samples),style=3)
for (pbi in 1:length(samples)){ # This will update the nucmer object
  sample<-samples[pbi]
  allvars<-nucmer[nucmer$qname==sample,]
  # Check changes in query protein sequence according to variants
  for(i in 1:nrow(allvars)){ # Assuming they are sorted numerically
    nucline<-allvars[i,]
    rpos<-nucline[1,"rpos"]
    rvar<-nucline[1,"rvar"]
    qvar<-nucline[1,"qvar"]
    qpos<-nucline[1,"qpos"]
    qlength<-nucline[1,"qlength"]
    
    # Match over GFF3 annotation
    a<-rpos-gff3[,4]
    b<-rpos-gff3[,5]
    signs<-sign(a)*sign(b)
    w<-which(signs==-1)
    
    # Outside genes scenarios
    if(length(w)==0){
      if(rpos<gff3[1,4]){
        protein<-"5'UTR";output<-c(rpos,"extragenic")
      } else if(rpos>gff3[1,5]){
        protein<-"3'UTR";output<-c(rpos,"extragenic")
      } else {
        protein<-"intergenic";output<-c(rpos,"extragenic")
      }
      
    } else{ # Inside genes scenario
      start<-gff3[w,4]
      end<-gff3[w,5]
      protein<-gff3[w,9]
      refdnaseq<-DNAString(paste0(refseq[start:end],collapse=""))
      refpepseq<-Biostrings::translate(refdnaseq)
      refpepseq<-strsplit(as.character(refpepseq),"")[[1]]
      if(qvar=="."){ # Deletion scenario
        if((nchar(rvar)%%3)!=0){ # Deletion frameshift scenario
          mutpos<-round((rpos-start+1)/3)
          output<-c(paste0(refpepseq[mutpos],mutpos),"deletion_frameshift")
        } else { # In-frame deletion
          varseq<-refseq
          varseq<-varseq[-(rpos:(rpos+nchar(rvar)-1))]
          varseq<-varseq[start:(end-nchar(rvar))]
          vardnaseq<-DNAString(paste0(varseq,collapse=""))
          varpepseq<-Biostrings::translate(vardnaseq)
          varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
          
          for(j in 1:length(refpepseq)){
            refj<-refpepseq[j]
            varj<-varpepseq[j]
            if(refj!=varj){
              if(varj=="*"){
                output<-c(paste0(refj,j),"deletion_stop")
              } else {
                output<-c(paste0(refj,j),"deletion")
              }
              break()
            }
          }
        }
      } else if(rvar=="."){ # Insertion scenario
        if((nchar(qvar)%%3)!=0){ # Insertion frameshift scenario
          mutpos<-round((rpos-start+1)/3)
          output<-c(paste0(refpepseq[mutpos],mutpos),"insertion_frameshift")
        } else { # In-frame insertion
          varseq<-c(refseq[1:rpos],strsplit(qvar,"")[[1]],refseq[(rpos+1):length(refseq)])
          varseq<-varseq[start:(end+nchar(qvar))]
          vardnaseq<-DNAString(paste0(varseq,collapse=""))
          varpepseq<-Biostrings::translate(vardnaseq)
          varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
          
          for(j in 1:length(refpepseq)){
            refj<-refpepseq[j]
            varj<-varpepseq[j]
            if(refj!=varj){
              nr_aa_inserted<-nchar(qvar)/3
              multivarj<-varpepseq[j:(j+nr_aa_inserted-1)]
              if(any(multivarj=="*")){
                multivarj<-paste0(multivarj,collapse="")
                output<-c(paste0(multivarj,j),"insertion_stop")
              } else{
                multivarj<-paste0(multivarj,collapse="")
                output<-c(paste0(multivarj,j),"insertion")
              }
              break()
            }
          }
        }
      } else { # SNP scenario
        if(nchar(qvar)==1){ # Single nucleotide scenario
          varseq<-refseq
          varseq[rpos]<-qvar
          varseq<-varseq[start:end]
          vardnaseq<-DNAString(paste0(varseq,collapse=""))
          varpepseq<-Biostrings::translate(vardnaseq)
          varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
          mutpos<-which(varpepseq!=refpepseq)
          if(length(mutpos)==0){ # Silent SNP scenario
            mutpos<-round((rpos-start+1)/3)
            refaa<-refpepseq[mutpos]
            varaa<-varpepseq[mutpos]
            output<-c(paste0(refaa,mutpos,varaa),"SNP_silent")
          } else { # Changed aa scenario
            refaa<-refpepseq[mutpos]
            varaa<-varpepseq[mutpos]
            if(varaa=="*"){
              output<-c(paste0(refaa,mutpos,varaa),"SNP_stop")
            } else {
              output<-c(paste0(refaa,mutpos,varaa),"SNP")
            }
          }
        } else { # Multiple neighboring nucleotides
          varlength<-nchar(qvar)
          varseq<-refseq
          varseq[rpos:(rpos+varlength-1)]<-strsplit(qvar,"")[[1]]
          varseq<-varseq[start:end]
          vardnaseq<-DNAString(paste0(varseq,collapse=""))
          varpepseq<-Biostrings::translate(vardnaseq)
          varpepseq<-strsplit(as.character(varpepseq),"")[[1]]
          mutpos<-which(varpepseq!=refpepseq)
          if(length(mutpos)==0){ # Silent SNP scenario
            mutpos<-round((rpos-start+1)/3)
            refaa<-refpepseq[mutpos]
            varaa<-varpepseq[mutpos]
            output<-c(paste0(refaa,mutpos,varaa),"SNP_silent")
          } else { # Changed aa scenario
            refaa<-paste0(refpepseq[mutpos],collapse="")
            varaa<-paste0(varpepseq[mutpos],collapse="")
            if(any(varaa=="*")){
              output<-c(paste0(refaa,mutpos[1],varaa),"SNP_stop")
            } else {
              output<-c(paste0(refaa,mutpos[1],varaa),"SNP")
            }
          }
        }
      }
    }
    results<-rbind(results,c(sample,rpos,rvar,qvar,qpos,qlength,protein,output,annot[protein]))
  }
  setTxtProgressBar(pb,pbi)
}

results<-as.data.frame(results,stringsAsFactors=FALSE)

results_after_june<-as.data.frame(results,stringsAsFactors=FALSE) ##** this is AFTER june. it a mistake
results_until_june1a<-as.data.frame(results,stringsAsFactors=FALSE)
results_until_june1b<-as.data.frame(results,stringsAsFactors=FALSE)
results_until_june1c<-as.data.frame(results,stringsAsFactors=FALSE)
results_until_june2<-as.data.frame(results,stringsAsFactors=FALSE)

# results=results_after_june
# results=results_until_june1
# results=results_until_june2

results_all=rbind(results_after_june,results_until_june1a,results_until_june1b,results_until_june1c,results_until_june2)




write.csv(results,file="/home/bianca/Downloads/ml_covid/inputs_for_nextclade/after_june_numcer/covid_annot_after_june.csv",row.names = FALSE)

#until june 1
write.csv(results,file="/home/bianca/Downloads/ml_covid/inputs_for_nextclade/covid_annot_until_june1.csv",row.names = FALSE)

#until june 2
write.csv(results,file="/home/bianca/Downloads/ml_covid/inputs_for_nextclade/covid_annot_until_june2.csv",row.names = FALSE)

# all results
write.csv(results_all,file="/home/bianca/Downloads/ml_covid/inputs_for_nextclade/covid_annot_all_together.csv",row.names = FALSE)
save(results,file="covid_annot.rda")
results<-read.csv("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/covid_annot_all_together.csv")
## read metadata and join with results

metadata<-read.table("/home/bianca/Downloads/ml_covid/first_inputs/all_covid.tsv",sep="\t")
colnames(metadata)=metadata[1,]
metadata=metadata[-1,]
metadata$date2=as.Date(metadata$date)
colnames(metadata)=c("sample","country","division","region","date","subdivision","date2")

library(tidyverse)
resmeta=full_join(results,metadata,by=("sample"))

# exclude Alpha mutations
alpha_mutations<-read.table("/home/bianca/Downloads/ml_covid/alpha_mutations.tsv",sep="\t")
results1=resmeta[!(resmeta$variant %in% alpha_mutations$V1),]
rm(resmeta)
rm(results)
rm(alpha_mutations)
rm(metadata)
write.table(resmeta,"/home/bianca/Downloads/ml_covid/vcf_covid_and_metadata_all.txt", sep="\t",quote=F,row.names = T)

# Attica
attica_b=results1[results1$date2 >= "2020-12-01" & results1$date2 <= "2021-06-01" & results1$division=='Attica',]
attica_a=results1[results1$date2 >= "2021-06-01" & results1$date2 <= "2021-10-30" & results1$division=='Attica',]
rm(results1)



###reduce results 
attica_b_nona=na.omit(attica_b)
attica_a_nona=na.omit(attica_a)
rm(attica_a)
rm(attica_b)

# "brown3"= resistant, "deepskyblue4"=susceptible

results_before=attica_b_nona
vec1<-as.factor(attica_a_nona$variant %in% attica_b_nona$variant)

myColors1<-apply(attica_b_nona, 1,function(x) ifelse(attica_b_nona$variant %in% attica_a_nona$variant, "brown3","deepskyblue4"))
myColors1=myColors1[,8]

results_after=attica_a_nona
vec2<-as.factor(attica_b_nona$variant %in% attica_a_nona$variant)

myColors<-apply(attica_a_nona, 2,function(x) ifelse(attica_b_nona$variant %in% attica_a_nona$variant, "brown3","deepskyblue4"))
myColors=myColors[,8]
myColors<-lapply(vec1,function(x) ifelse(vec1=="TRUE", "brown3","deepskyblue4"))

data=results_before
data=results
data=attica_b_nona
data=attica_a_nona
myColors <- ifelse(levels(vec1)=="TRUE" , "brown3" ,
                   ifelse(levels(vec2)=="FALSE", "deepskyblue4",
                          "grey90" ) )
png("/home/bianca/Downloads/ml_covid/attica_b_a.png",w=3000,h=2000,res=300)
par(mfrow=c(2,2))

# Nucleotide events
occ<-sort(table(apply(results_before[,c("refvar","refpos","qvar")],1,paste0,collapse="")),dec=TRUE)[1:10]
par(las=2,mar=c(8,5,5,1))
barplot(occ,ylab="nr of samples",main="Most frequent events (nucleotide) 
        for Attica before vaccination",col=myColors1)

# 
# Protein events

occ<-sort(table(apply(results_before[,c("protein","variant")],1,paste0,collapse=":")),dec=TRUE)[1:10]
par(las=2,mar=c(8,5,5,1))
barplot(occ,ylab="nr of samples",main="Most frequent events (protein) 
        for Attica before vaccination",col=myColors1)

# nuc
#png("/home/bianca/Downloads/ml_covid/test.png",w=3000,h=2000,res=300)
occ<-sort(table(apply(results_after[,c("refvar","refpos","qvar")],1,paste0,collapse="")),dec=TRUE)[1:10]
#occ<-occ[-7]
par(las=2,mar=c(8,5,5,1))
barplot(occ,ylab="nr of samples",main="Most frequent events (nucleotide) 
        for Attica after vaccination",col=myColors)
#heat.colors(length(occ)
# Protein events
occ<-sort(table(apply(results_after[,c("protein","variant")],1,paste0,collapse=":")),dec=TRUE)[1:10]
par(las=2,mar=c(8,5,5,1))
barplot(occ,ylab="nr of samples",main="Most frequent events (protein) for Attica 
        after vaccination",col=myColors)

dev.off()


# Prepare a vector of colors with specific color rsistant sus

# ### Describe Results ----
# png("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/after_june_numcer/covid_annot.png",w=3000,h=2000,res=300)
# par(mfrow=c(2,3))
# 
# # Most mutated samples
# occ<-sort(table(results$sample),dec=TRUE)[1:20]
# par(las=2,mar=c(15,5,5,1))
# barplot(occ,ylab="nr of mutations",main="Most mutated samples",col=heat.colors(length(occ)))
# 
# # Mutations per sample
# occ<-table(table(results$sample))
# par(las=2,mar=c(5,5,5,1))
# barplot(occ,xlab="nr of mutations",main="Overall mutations per sample",col="cornflowerblue")
# 
# # Variant classes
# occ<-sort(table(results$varclass),dec=TRUE)
# par(las=2,mar=c(8,5,5,1))
# barplot(occ,ylab="nr of events",main="Most frequent events per class",col=heat.colors(length(occ)))
# 
# # Variant class (A/T, etc)
# occ<-sort(table(apply(results[,c("refvar","qvar")],1,paste0,collapse=">")),dec=TRUE)[1:20]
# par(las=2,mar=c(5,5,5,1))
# barplot(occ,ylab="nr of samples",main="Most frequent per type",col=heat.colors(length(occ)))



##histograms

na156s<-read.table("/home/bianca/Downloads/ml_covid/metadata_clades/na156s.tsv",sep="\t", row.names = 1)
na156s1<-read.table("/home/bianca/Downloads/ml_covid/metadata_clades/ng204p.tsv",sep="\t", row.names = 1)
na156s2<-read.table("/home/bianca/Downloads/ml_covid/metadata_clades/sd138h.tsv",sep="\t", row.names = 1)
na156s3<-read.table("/home/bianca/Downloads/ml_covid/metadata_clades/sg1219v.tsv",sep="\t", row.names = 1)
na156s4<-read.table("/home/bianca/Downloads/ml_covid/metadata_clades/sy508h.tsv",sep="\t", row.names = 1)

fc_covid<-read.table("/home/bianca/Downloads/ml_covid/fc_covid_division_time.tsv",sep="\t")
colnames(fc_covid)<-fc_covid[1,]
fc_covid=fc_covid[-1,]

fc_covid1<-read.table("/home/bianca/Downloads/ml_covid/longer_fc_covid.tsv",sep="\t")
colnames(fc_covid1)<-fc_covid1[1,]
fc_covid1=fc_covid1[-1,]


fc_covid2<-read.table("/home/bianca/Downloads/ml_covid/only_fc_covid.csv",sep="\t")
colnames(fc_covid2)<-fc_covid2[1,]
fc_covid2=fc_covid2[-1,]


colnames(na156s)<-na156s[1,]
na156s<-na156s[-1,]
na156s$date2<-as.Date(na156s$date, format="%Y-%m-%d")

colnames(na156s1)<-na156s1[1,]
na156s1<-na156s1[-1,]
na156s1$date2<-as.Date(na156s1$date, format="%Y-%m-%d")

colnames(na156s2)<-na156s2[1,]
na156s2<-na156s2[-1,]
na156s2$date2<-as.Date(na156s2$date, format="%Y-%m-%d")


colnames(na156s3)<-na156s3[1,]
na156s3<-na156s3[-1,]
na156s3$date2<-as.Date(na156s3$date, format="%Y-%m-%d")

colnames(na156s4)<-na156s4[1,]
na156s4<-na156s4[-1,]
na156s4$date2<-as.Date(na156s4$date, format="%Y-%m-%d")



greek<-c("Western Macedonia","Thessaly","Attica","Epirus","Central Macedonia","Western Greece","South Aegean","Ionian Islands","Crete","Eastern Macedonia and Thrace","Peloponnese","Central Greece","North Aegean")

now<-na156s[na156s$division %in% greek,]
now1<-na156s1[na156s1$division %in% greek,]
now2<-na156s2[na156s2$division %in% greek,]
now3<-na156s3[na156s3$division %in% greek,]
now4<-na156s4[na156s4$division %in% greek,]

breaks<-c("2021-01-1","2021-03-1","2021-05-1","2021-07-1","2021-09-1")
cbp_12 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#882255","#661100","#999933","#44AA99")

cbp_13 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#882255","#661100","#999933","#44AA99","#999999")

### SAVE EACH ONE IN A DIFFERENT NOW

one<-ggplot(now,aes(x=date2,fill=division)) +
  geom_histogram() + 
  scale_fill_manual(values=cbp_12) + 
  scale_x_date(breaks = date_breaks("1 month"),labels = date_format("%m")) + 
  labs(y="Frequency",x="Months",title="Histogram of N:A156S sublineage") + 
  theme_bw()

two<-ggplot(now1,aes(x=date2,fill=division)) +
  geom_histogram() + 
  scale_fill_manual(values=cbp_12) + 
  scale_x_date(breaks = date_breaks("1 month"),labels = date_format("%m")) + 
  labs(y="Frequency",x="Months",title="Histogram of N:G204P sublineage") + 
  theme_bw()

three<-ggplot(now2,aes(x=date2,fill=division)) +
  geom_histogram() + 
  scale_fill_manual(values=cbp_13) + 
  scale_x_date(breaks = date_breaks("1 month"),labels = date_format("%m")) + 
  labs(y="Frequency",x="Months",title="Histogram of S:D138H sublineage") + 
  theme_bw()

four<-ggplot(now3,aes(x=date2,fill=division)) +
  geom_histogram() + 
  scale_fill_manual(values=cbp_13) + 
  scale_x_date(breaks = date_breaks("1 month"),labels = date_format("%m")) + 
  labs(y="Frequency",x="Months",title="Histogram of S:G1219V sublineage") + 
  theme_bw()


five<-ggplot(now4,aes(x=date2,fill=division)) +
  geom_histogram() + 
  scale_fill_manual(values=cbp_13) + 
  scale_x_date(breaks = date_breaks("1 month"),labels = date_format("%m")) + 
  labs(y="Frequency",x="Months",title="Histogram of S:Y508H sublineage") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggarrange(one,two,three, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 2)

png("/home/bianca/Downloads/ml_covid/metadata_clades/histograms2.png",w=3000,h=2000,res=1024)
par(mfrow=c(2,3))
ggarrange(four,five + rremove("x.text"),
          labels = c("D","E"),
          ncol = 2, nrow = 2)
dev.off()


now$sublineage<-"N:A156S"
now1$sublineage<-"N:G204P"
now2$sublineage<-"S:D138H"
now3$sublineage<-"S:G1219V"
now4$sublineage<-"S:Y508H"

now_all<-rbind(now[,c(3,6,7)],now1[,c(6,12,13)],now2[,c(6,11,12)],now3[,c(6,11,12)],now4[,c(6,11,12)])

library(dplyr)
now_all$values<-1/ave(seq_along(now_all$sublineage), now_all$sublineage, FUN = length)
wow<-now_all %>% group_by(division,date2) %>% mutate( percentage = values/sum(values))




wow %>%
  ggplot(aes(date2,as.character(division),fill = division)) + 
  geom_boxplot() + scale_fill_manual(values=cbp_13) + 
  facet_grid(cols = vars(sublineage)) + 
  labs (y = "Division", x = "Months") +
  ggtitle("Distribution of important sublineages") +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))



fc_covid %>%
  ggplot(aes(date2,as.character(division),fill = division)) + 
  geom_boxplot() + scale_fill_manual(values=cbp_13) + 
  facet_grid(cols = vars(sublineage)) + 
  labs (y = "Division", x = "Months") +
  ggtitle("Distribution of important sublineages") +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fc_covid1 %>%
  ggplot(aes(as.numeric(cases),as.character(Division),col = waves)) + 
  geom_point() + scale_fill_manual(values=cbp_13) + 
  #facet_grid(cols = vars(waves)) + 
  labs (x = "Cases", y = "Division") +
  ggtitle("Number of cases for each period and division") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


fc_covid2 %>%
  ggplot(aes(log2(as.numeric(FC)),as.character(Division),col = waves)) + 
  geom_point() + scale_fill_manual(values=cbp_13) + 
  #facet_grid(cols = vars(waves)) + 
  labs (x = "Cases", y = "Division") +
  ggtitle("Number of cases for each period and division") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

