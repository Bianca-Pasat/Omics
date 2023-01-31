library(biomaRt)
# list<- read.table("/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/RIDD_EXON_SHARED_RESULTS/ridd/results/RIDD_candidate_genes_crosstalk_between_two_methods")
# head(list)
list<-c("APP","MACF1","JUP","MYH10","FLNA","MEGF8","ABL1","ATM","CTNNB1","PKD1","SPTBN4")
listconfr<-c("ANXA8","UNC13A","CD44","TNFRSF14","ADAMTS14","CDK11A","ATAD3C","GPR39","ZNF23","SLC16A1","USP45","CDH1","NDUFA3")
listEns<-read.table("/home/ben/Desktop/exon_coordinates/hg_enst.csv", sep="\t")
list<-"SCD"
list<-c("SCD","VEGFA","LRP5","G0S2")

mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))

homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")

enst_to_gene<-getBM(attributes = c("ensembl_transcript_id","hgnc_symbol"), 
                    filters="ensembl_transcript_id", values=listEns, mart=mart)


enst_to_gene_exons_8<-getBM(attributes = c("ensembl_transcript_id","hgnc_symbol"), 
                    filters="ensembl_transcript_id", values=new_group, mart=homo.anno)

goids = getBM(attributes = c('hgnc_symbol','go_id','name_1006', 'definition_1006'), 
              filters = 'hgnc_symbol', 
              values = listconfr, 
              mart = mart)

write.table(goids,"/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/RIDD_EXON_SHARED_RESULTS/ridd/results/selected_conf_goid",sep="\t", quote=F)

write.table(goids,"/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/RIDD_EXON_SHARED_RESULTS/ridd/results/crosstalk_go_terms_exact_seq",sep="\t", quote=F)

idridd<-read.table("/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/PRIMARY_INPUTS/ridd/all_ridd_targets.txt")
id<-unique(idridd)
seq<-getSequence(id=id, type="hgnc_symbol", seqType="gene_exon",mart=homo.anno)
seq<-getSequence(id=listconfr, type="hgnc_symbol", seqType="cdna",mart=homo.anno)
seq<-getSequence(id=list, type="hgnc_symbol", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/ridd/primary_inputs/ridd_input.fasta")
exportFASTA(seq, file="/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/ridd_selected.fasta")
exportFASTA(seq, file="/home/bianca/Desktop/common_ridd.fasta")

list<-c("CDH1", "CDK11A", "USP45","ZNF23","SCD","VEGFA","LRP5","G0S2")
seq<-getSequence(id=new_group, type="ensembl_transcript_id", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/bianca/Desktop/alice_selected.fasta")

idexon<-read.table("/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/PRIMARY_INPUTS/exon_usage/exon_all_ids.csv")
id<-unique(idexon)
seq<-getSequence(id=id, type="hgnc_symbol", seqType="gene_exon",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/exon_usage/primary_input/exons_all_input.fasta")



options(tibble.pri_max=Inf)

Go_tbl<-biomartr::getGO(organism="Homo sapiens",
                        genes=list,
                        filters="hgnc_symbol")
Go_tbl
idridd2<-read.table("/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/ridd/check_with_old/exact_seq", sep="\t")
idridd3<-read.table("/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/ridd/check_with_old_2/exact_seq", sep="\t")
