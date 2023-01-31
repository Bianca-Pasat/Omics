drugs_genes<-fread("/home/bianca/Downloads/interactions.tsv", data.table=F, na.strings = "")


nala2<-drugs_genes %>%  mutate(gene_name = replace_na('none')) 
nala3<-nala2%>% pivot_wider(names_from = gene_name, values_from = drug_claim_name)
nala2<-drugs_genes%>% pivot_wider(names_from=drug_claim_name, values_from = gene_name)

nala<-as.data.frame(drugs_genes[,c(1,8,10)])
nala<-as.data.frame(drugs_genes[,c(1,8)])
nala2<-nala %>% pivot_wider(names_from = drug_name, values_from = interaction_group_score)
nala2<-nala %>% pivot_wider(names_from=drug_name, values_from = gene_name, values_fn = list)
nala2<-nala %>% pivot_wider(names_from=gene_name, values_from = drug_name, values_fn = list)
write.table(nala2,"/home/bianca/Desktop/nala2.tsv", sep=",",quote=F,row.names = F)
nalajson<-toJSON(nala2,indent = 1)
nalajson<-toJSON(nala2,indent = 0)
write(nalajson,file="nalajson.json")

library(jsonlite)
library(RJSONIO)
library(rjson)
# json_data_raw<-fromJSON("mydata.txt")
json_data_raw<-nalajson

json_file <- lapply(json_data_raw, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

output <- do.call("rbind", json_file)
write.csv(output, file="genedrugjson.csv",row.names = FALSE,quote=F)
file.show("json.csv")

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
