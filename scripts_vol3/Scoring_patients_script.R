#IRE1sign38 <- c("ASS1","C3","CCL20","COL4A6","CXCL2","CXCL5","CXCL8","IFI44L","IL1B","IL6","KCNN2","MMP1","MMP12","MMP3","PLA2G4A","PPP4R4","SERPINB2","TFPI2","ZNF804A","ANGPT1","CFH","CFI","CLEC3B","COL3A1","COL8A1","DACH1","DCN","FHL1","GAS1","LUM","OXTR","PLAC8","RGS4","TAGLN","TGFB2","THBS1","TIMP3","TMEM255A")

XBP1sign <- c("ASS1","C3","CCL20","COL4A6","CXCL2","CXCL5","CXCL8","IFI44L","IL1B","IL6","KCNN2","MMP1","MMP12","MMP3","PLA2G4A","PPP4R4","SERPINB2","TFPI2","ZNF804A")
RIDDsign <- c("ANGPT1","CFH","CFI","CLEC3B","COL3A1","COL8A1","DACH1","DCN","FHL1","GAS1","LUM","OXTR","PLAC8","RGS4","TAGLN","TGFB2","THBS1","TIMP3","TMEM255A")

index.rown.df <- data.frame("signature_id"=rownames(aux),"label" = rep("unlabel"))
df.xbp1 <- match(XBP1sign,index.rown.df$signature_id)
df.xbp1 <- na.omit(df.xbp1)

index.rown.df$label[df.xbp1] <- "XBP1"
df.ridd <- match(RIDDsign,index.rown.df$signature_id)
df.ridd <- na.omit(df.ridd)
index.rown.df$label[df.ridd] <- "RIDD"

#fill in the cells of aux matrix
for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(rownames(aux)[i] %in% XBP1sign){
        if (logCPM.genesign[i,j] <= q[i,1]){ 
          aux[i,j] <- 1
        }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
          aux[i,j] <- 2
        }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
          aux[i,j] <- 3
        }else if (logCPM.genesign[i,j] >= q[i,3]){
          aux[i,j] <- 4
        }
      }else if(rownames(aux)[i] %in% RIDDsign){
        if (logCPM.genesign[i,j] <= q[i,1]){ 
          aux[i,j] <- 4
        }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
          aux[i,j] <- 3
        }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
          aux[i,j] <- 2
        }else if (logCPM.genesign[i,j] >= q[i,3]){
          aux[i,j] <- 1
        }
    }else {
      next
    }
  }
}

score.xbp1 <- c(rep(0,ncol(aux)))
score.ridd <- c(rep(0,ncol(aux)))

#calculate XBP1 and RIDD score
for(k in 1:ncol(aux)){
  
  score <- aggregate(x = aux[,k],                 
                     by = list(index.rown.df$label),       
                     FUN = sum)                       
  
  score.xbp1[k] <- score[score$Group.1=="XBP1",2] / length(which(index.rown.df$label == "XBP1"))
  score.ridd[k] <- score[score$Group.1=="RIDD",2] / length(which(index.rown.df$label == "RIDD"))
}
