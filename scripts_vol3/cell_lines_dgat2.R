mD231<-read.table("/home/ben/Downloads/MDAMB231.csv", sep=",")
colnames(mD231)=mD231[1,]
mD231=mD231[-1,]
rownames(mD231)=mD231[,1]
mD231=mD231[,-1]
md231=t(mD231)
rm(mD231)
md1231<-apply(md231,2,function(x) as.numeric(x))
rownames(md1231)=rownames(md231)
md231=md1231

mD468<-read.table("/home/ben/Downloads/MDAMB468.csv", sep=",")
colnames(mD468)=mD468[1,]
mD468=mD468[-1,]
rownames(mD468)=mD468[,1]
mD468=mD468[,-1]
md468=t(mD468)
rm(mD468)
md1468<-apply(md468,2,function(x) as.numeric(x))
rownames(md1468)=rownames(md468)
md468=md1468
rm(md1231)
rm(md1468)

cell_all=cbind(md231, md468)
cell_all=cell_all[!(rowSums(cell_all==0)),]


md231=md231[,1:100]
md468=md468[,1:100]


cl=c(rep(1,length(colnames(md231))), rep(0, length(colnames(md468))))
samples=cbind(md231, md468)

cl=cl=c(rep(1,50), rep(0,50)) 
samples=cell_all[,c(1:50,1407:1457)]
