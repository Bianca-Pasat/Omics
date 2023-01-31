
library(data.table)
library(tidyverse)
library(jsonlite)
library(RJSONIO)
library(rjson)

drugs_genes<-fread("/home/bianca/Downloads/interactions.tsv", data.table=F, na.strings = "")
drugs_genes<-fread("/home/bianca/Downloads/interactions.tsv", data.table=)
which(is.na(drugs_genes$drug_claim_name))
which(is.na(drugs_genes$drug_name))
which(is.na(drugs_genes$gene_name))
which(is.na(drugs_genes$gene_claim_name))
# nala2<-drugs_genes %>%  mutate(gene_name = replace_na('none')) 
# nala3<-nala2%>% pivot_wider(names_from = gene_name, values_from = drug_claim_name)
# nala2<-drugs_genes%>% pivot_wider(names_from=drug_claim_name, values_from = gene_name)
# 
# nala<-as.data.frame(drugs_genes[,c(1,8,10)])
nala<-as.data.frame(drugs_genes[,c(1,8)])
nala<-as.data.frame(drugs_genes)
# nala2<-nala %>% pivot_wider(names_from = drug_name, values_from = interaction_group_score)
# nala2<-nala %>% pivot_wider(names_from=drug_name, values_from = gene_name, values_fn = list)
nala2<-nala %>% pivot_wider(names_from=gene_name, values_from = drug_name, values_fn = list)
nala2<-nala %>% pivot_wider(names_from=gene_name, values_from = he, values_fn = list)
nala2<-nala %>% pivot_wider(values_from=gene_name, names_from = drug_name)
write.table(nala2,"/home/bianca/Desktop/nala2.tsv", sep=",",quote=F,row.names = F)
nalajson<-toJSON(nala2,indent = 1)
nalajson<-toJSON(nala2,indent = 0, "R") #  this was last used
write(nalajson,file="nalajson2.json")
write(nalajson,file="DnG_whole.json")

### without tidyr blah blah
dat <- drugs_genes
colnames(drugs_genes)
by(dat, list(dat$gene_claim_name), function(x) {
  
  outer_template <- '"gene_claim_name":"%s", 
      [ 
%s 
      ]'
  
  inner_template <- '        [ "drug_claim_name":"%s",
          "entrez_id":"%s", 
          "interaction_claim_source":"%s",
        ]'
  
  condition <- paste0(apply(x, 1, function(y) {
    sprintf(inner_template, y["drug_claim_name"],
            tolower(y["entrez_id"]), y["interaction_claim_source"])
  }), collapse=",\n")
  
  sprintf(outer_template, x$gene_claim_name[1], condition)
  
}) -> genes

genes_json <- paste0(genes, collapse=",\n")
cat(genes_json)
write(genes_json,file="DnG_more_columns.json")

by(dat, list(dat$gene_claim_name), function(x) {
  
  outer_template <- '{ "gene_claim_name":"%s","attributes":[%s]}'
  inner_template <- '{ "drug_claim_name":"%s","entrez_id":"%s","interaction_claim_source":"%s"}'

  condition <- paste0(apply(x, 1, function(y) {
    sprintf(inner_template, y["drug_claim_name"],
            tolower(y["entrez_id"]), y["interaction_claim_source"])
  }), collapse=",\n")
  
  sprintf(outer_template, x$gene_claim_name[1], condition)
  
}) -> genes

## alll
dat<-drugs_genes[,c(2,1,3:11)]
by(dat, list(dat$gene_claim_name), function(x) {
  
  outer_template <- '{ "gene_claim_name":"%s","attributes":[%s]}'
  inner_template <- '{ "gene_name":"%s","entrez_id":"%s","interaction_claim_source":"%s","interaction_types":"%s","drug_claim_name":"%s","drug_claim_primary_name":"%s","drug_name":"%s","drug_concept_id":"%s","interaction_group_score":"%s","PMIDs":"%s"}'
  
  condition <- paste0(apply(x, 1, function(y) {
    sprintf(inner_template,y["gene_name"],y["entrez_id"],y["interaction_claim_source"],y["interaction_types"],y["drug_claim_name"],y["drug_claim_primary_name"],y["drug_name"],y["drug_concept_id"],y["interaction_group_score"],y["PMIDs"])
  }), collapse=",\n")
  
  sprintf(outer_template, x$gene_claim_name[1], condition)
  
}) -> genes

## all sligthly different THIS ONE! CHECK HOW YOU COULD PUT OUTER BRACKETS
by(dat, list(dat$gene_claim_name), function(x) {
  
  outer_template <- '{ "gene_claim_name":"%s","gene_name":"%s","entrez_id":"%s","attributes":[%s]}'
  inner_template <- '{ "interaction_claim_source":"%s","interaction_types":"%s","drug_claim_name":"%s","drug_claim_primary_name":"%s","drug_name":"%s","drug_concept_id":"%s","interaction_group_score":"%s","PMIDs":"%s"}'
  
  condition <- paste0(apply(x, 1, function(y) {
    sprintf(inner_template,y["interaction_claim_source"],y["interaction_types"],y["drug_claim_name"],y["drug_claim_primary_name"],y["drug_name"],y["drug_concept_id"],y["interaction_group_score"],y["PMIDs"])
  }), collapse=",\n")
  
  sprintf(outer_template, x$gene_claim_name[1], x$gene_name[1],x$entrez_id[1],condition)
  
}) -> genes

genes_json <- paste0(genes, collapse=",\n")
cat(genes_json)
write(genes_json,file="DnG_grouped_on_genes.json")

#### in relation to drugs uncomplete!
by(dat, list(dat$drug_claim_name), function(x) {
  
  outer_template <- '{ "drug_claim_name":"%s","drug_name":"%s","entrez_id":"%s","attributes":[%s]}'
  inner_template <- '{ "interaction_claim_source":"%s","interaction_types":"%s","drug_claim_name":"%s","drug_claim_primary_name":"%s","drug_name":"%s","drug_concept_id":"%s","interaction_group_score":"%s","PMIDs":"%s"}'
  
  condition <- paste0(apply(x, 1, function(y) {
    sprintf(inner_template,y["interaction_claim_source"],y["interaction_types"],y["drug_claim_name"],y["drug_claim_primary_name"],y["drug_name"],y["drug_concept_id"],y["interaction_group_score"],y["PMIDs"])
  }), collapse=",\n")
  
  sprintf(outer_template, x$gene_claim_name[1], x$gene_name[1],x$entrez_id[1],condition)
  
}) -> genes

genes_json <- paste0(genes, collapse=",\n")
cat(genes_json)
write(genes_json,file="DnG_all_columns_lighter.json")
# json_data_raw<-fromJSON("nalajson.json")
# json_data_raw<-nalajson
# 
# json_file <- lapply(json_data_raw, function(x) {
#   x[sapply(x, is.null)] <- NA
#   unlist(x)
# })
# 
# output <- do.call("rbind", json_file)
# write.csv(output, file="genedrugjson.csv",row.names = FALSE,quote=F)
# file.show("genedrugjson.csv")

for (i in 1:dim(nala2)[2])
{
  lapply(nala2[[i]], function(x) write.table( data.frame(x), '/home/bianca/Downloads/test.csv'  , append= T, sep=',' ))
}
hah<-data.frame()

for (i in 1:dim(nala2)[2]){
  hahahah<-do.call(rbind, nala2[[i]])
  hah<-bind_rows(hah,hahahah)
}


listnala<-{nala} %>%
  dplyr::group_by(drug_name) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
length(nala2[[1]])
