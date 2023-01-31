#### Proteomics analysis with rankprod and normalize with 


proteins<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final.txt",sep="\t")
proteins_info<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final_info_correct_order.txt",sep="\t")
colnames(proteins)=proteins[1,]
proteins=proteins[-1,]

proteins8=proteins[,c(10:12,22:24)]
proteins24=proteins[,c(10:12,13:15)]
proteins48=proteins[,c(10:12,16:18)]
proteins72=proteins[,c(10:12,19:21)]


norm_proteins=function.norm(proteins8)
norm_proteins=function.norm(proteins24)
norm_proteins=function.norm(proteins48)
norm_proteins=function.norm(proteins72)



#### vs DMSO 8 ALWAYS
cl=c(rep(1,3),rep(0,3))
## dmso 8 vs mkc 8
samples=function.norm(proteins8)
out_dir="/home/bianca/Desktop/rp_prot_mkc_8_vs_dmso_8.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_8_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_8_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_8_vs_dmso_8_DEregulatedBio.txt"

## dmso 8 vs mkc 24
samples=function.norm(proteins24)


out_dir="/home/bianca/Desktop/rp_prot_mkc_24_vs_dmso_8.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_24_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_24_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_24_vs_dmso_8_DEregulatedBio.txt"
## dmso 8 vs mkc 48
samples=function.norm(proteins48)

out_dir="/home/bianca/Desktop/rp_prot_mkc_48_vs_dmso_8.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_48_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_48_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_48_vs_dmso_8_DEregulatedBio.txt"
## dmso 8 vs mkc 72
samples=function.norm(proteins72)

out_dir="/home/bianca/Desktop/rp_prot_mkc_72_vs_dmso_8.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_72_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_72_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_72_vs_dmso_8_DEregulatedBio.txt"


dmso72=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/proteomics_time_series_72_dmso_bio.txt",sep="\t")
dmso72=rownames(dmso72)

### EBSeqHMM
library(EBSeq)
