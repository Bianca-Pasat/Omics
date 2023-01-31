de_prot_8<-read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
X2=t(X1)
de_8_ex=X2[rownames(X2) %in% rownames(de_prot_8),]

#8 hours
trans_cts=as.data.frame(de_8_ex)
trans_cts=trans_cts +1
sample_info=X_info
sample_info$sample=sample_info$short_name


trans_cts=as.data.frame(t(X1))
sample_info=X_info
sample_info$sample=sample_info$short_name
