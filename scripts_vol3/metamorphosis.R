library(dplyr)

dt<-read.table("/home/ben/Desktop/exon_coordinates/test_cor.csv", sep="\t", header=T)
dt1<-do.call(rbind, lapply(dt, function(x) as.data.frame(t(x))))
ref<-read.table("/home/ben/Desktop/exon_coordinates//test_gene_ref.csv", sep="\t", header=T,fill=TRUE)
ref<-as.data.frame(ref)
ref1<-ref %>% mutate_all(~replace(., is.na(.), 1000000))


found_genes<-ref1[ref1$gene_name %in% dt$gene_name,]


df1<-data.frame()

dts<-dt[1:3,1:3]
fgs<-found_genes[1:3,]

# return(df1<-cbind(dt$gene_name[x],colnames(found_genes)[j]))
# 
# #dt[x,4] <- colnames(found_genes)[j]
# #return(dt)


### MY AMAZING RECURSIVE ALGORITHM ***** ALGORITHM WORKS BUT LOOPING DOES NOT BECAUSE THE PBJECTS ARE LISTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
my_fun3<- function(found_genes,dt){
  for (i in 1:nrow(found_genes)){
    for (j in 2:ncol(found_genes)){
      for (x in 1:nrow(dt)){
        for (y in 2:ncol(dt)){
          g<-ncol(found_genes)
          #print(dt[i,j])
          found_genes1<-found_genes
          if (dt[i,2] >= found_genes[i,g]){
            print("bigger")
            dt[i,4]<-colnames(found_genes)[g]
          }
          else if(ncol(found_genes) == 2) {
              my_fun3(found_genes[,-g],dt)}
          else {
            print("else")
            dt[i,4]<-"exon_1"
          }
      }
    }
  }
  }
  #found_genes<-found_genes1
  return(dt)
}



df1<-data.frame()
my_fun3<- function(found_genes,dt){
  for (i in 1:nrow(found_genes)){
          g<-ncol(found_genes)
          #found_genes1<-found_genes
          if (dt[i,2] >= found_genes1[i,g]){
            #dt[i,4]<-colnames(found_genes1)[g]
            print("bigger")
            df1[i,1]<-colnames(found_genes)[g]
          }
          else if(ncol(found_genes1) > 2) {
            my_fun3(found_genes1[,-g],dt)
          }
          else if(ncol(found_genes1)==2){
            print("else")
            #dt[i,4]<-"exon1"
            df1[i,1]<-"exon1"
          }
          #found_genes1<-found_genes
  }
  return(df1)
  #return(dt)
}


##THIS ONE
my_fun3<- function(found_genes,dt){
  for (i in 1:nrow(found_genes)){
    g<-ncol(found_genes)
    if( ncol(found_genes) > 2){
    if (dt[i,2] >= found_genes[i,g]){
      print(colnames(found_genes)[g])
      df1[i,1]<-colnames(found_genes)[g]
      #dt[i,4]<-colnames(found_genes1)[g]
    }
    else  {
      my_fun3(found_genes[,-g],dt)
    }
      }
    else if(ncol(found_genes)==2){
      print("exon1")
      df1[i,1]<-"exon1"
      #dt[i,4]<-"exon1"
    }
  }
  return(df1)
  #return(dt)
}



# this loop displays one item per column goes through all the columns and then prints
a<-for (i in 1:nrow(found_genes)){
      for (x in 2:ncol(found_genes)){
      print(found_genes[i,x])
    }
  }

############# MAKE FIRST A STRING AND THEN ADD IT TO THE MATRIX ****************************
my_fun3<- function(found_genes,dt){
  for (i in 1:nrow(found_genes)){
    for (j in 2:ncol(found_genes)){
      for (x in 1:nrow(dt)){
        for (y in 2:ncol(dt)){
          g<-ncol(found_genes)
          if (dt[x,2] >= found_genes[x,g]){
            line<-cbind(line, colnames(found_genes)[g] )
            #dt[x,4]<-colnames(found_genes)[g]
          }
          else if(ncol(found_genes) > 2) {
            return (my_fun3(found_genes[,-g],dt))
          }
          else {
            line<-cbind(line, colnames(found_genes)[g] )
            #dt[x,4]<-"exon_1"
          }
          return(line)
        }
      }
    }
  }
}


####RECURSIVE BUT ON ONE MATRIX

test1<-merge(dt,found_genes)
dt1<-data.frame()
my_fun4<- function(matrix){
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
          g<-ncol(matrix)
          print(matrix[i,j])
          if (matrix[i,2] >= matrix[i,g]){
            line<-cbind(line, colnames(matrix)[g] )
            #matrix[i,ncol(matrix) + 1]<-colnames(matrix)[g]
          }
          else if(ncol(matrix) > 3) {
            return (my_fun4(matrix[,-g]))
          }
          else {
            #matrix[i,ncol(matrix) + 1]<-"exon_1"
            line<-cbind(line,"exon_1")
          }
          
        }
  }
  #return(matrix)
  return(line)
}


### THIS FUNCTION NEEDS TO BE TRANSFOEMED INTO RECURSIVE
my_fun2<- function(found_genes,dt){
  for (i in 1:nrow(found_genes)){
    for (j in 2:ncol(found_genes)){
      for (x in 1:nrow(dt)){
        for (y in 2:ncol(dt)){
          g=ncol(found_genes)
          dt$winner <- ifelse(dt[,y] > found_genes[,g], 'Yes',
                              ifelse(my_fun2(found_genes[,g-1]), 'No', 'exon_1'))
          return(dt)
        }
      }
    }
  }
}



#BACKUP
# my_fun2<- function(found_genes,dt){
#   for (i in 1:nrow(found_genes)){
#     for (j in 2:ncol(found_genes)){
#       for (x in 1:nrow(dt)){
#         for (y in 2:ncol(dt)){
#           if (dt[x,y] >= found_genes[i,j]){
#             if (dt[x,y] <= found_genes[i,(j+1)]){
#               return(df1<-cbind(dt$gene_name[x],colnames(found_genes)[j]))
#             }
#             else {
#               return(df1<-cbind(dt$gene_name[x],colnames(found_genes)[j+1]))
#             }
#           }
#           return(dt)
#         }
#       }
#     }
#   }
# }


# dt1<-as.numeric(dt[,c(2,3)])
# my_fun<- function(dt1,dt2){
#   for (i in 1:nrow(found_genes)){
#     for (j in 2:ncol(found_genes)){
#       for (x in 1:nrow(dt)){
#         for (y in 2:ncol(dt)){
#          # if (dt[x,y] > found_genes[i,j] && dt[x,(y+1)] < found_genes[i,(j+1)])
#           #{
#           #dt$V4[x]<- colnames(found_genes)[j+1]
#           #}
#         }
#       }
#     }
#   }
# }

# my_fun2<- function(found_genes,dt){
#   for (i in 1:nrow(found_genes)){
#     for (j in 2:ncol(found_genes)){
#       for (x in 1:nrow(dt)){
#         for (y in 2:ncol(dt)){
#           if (dt[x,y] >= found_genes[i,j]){
#             if (dt[x,y] <= found_genes[i,(j+1)]){
#               dt[x,4]<-colnames(found_genes)[j]
#             }
#             else {
#               dt[x,4]<-colnames(found_genes)[j+1]
#             }
#           }
#           return(dt)
#         }
#       }
#     }
#   }
# }