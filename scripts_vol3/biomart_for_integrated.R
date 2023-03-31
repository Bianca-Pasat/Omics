library(biomaRt)

httr::set_config(httr::config(ssl_verifypeer = FALSE))
genes<-read.table("/home/ben/Desktop/MCIA-INTEGRATED/first2000_go_integrated.tsv", sep="\t")[,2]
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "useast")
seq<-getSequence(id=genes, type="hgnc_symbol", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/MCIA-INTEGRATED/first2000_go_integrated.fasta")



attributes <- c("entrezgene", "hgnc_symbol")
filters="hgnc_symbol"
values_ridd<-genes
genes_ridd <- getBM(attributes = attributes, filters=filters,values = values_ridd, mart = homo.anno)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)




transcriptCoordsByGene.GRangesList <-
  transcriptsBy (TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene") [e2f3]
# a GrangesList of length one, describing three transcripts

promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList,
                                 Hsapiens, upstream=1000, downstream=20)

ftable<-read.table("/home/ben/Downloads/sample.tab", sep="\t")
consensus_seq <- as.character("TGACGTCA")

for(i in 1:nrow(ftable)) {
  ftable$V3[i] <- str_detect(ftable$V2[i], consensus_seq)
}

integrated_ridds<-read.table("/home/ben/Desktop/MCIA-INTEGRATED/first2000_go_integrated_out2.txt")
