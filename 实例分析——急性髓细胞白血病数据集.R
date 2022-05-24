###白血病基因数据集###
setwd("E:/GEO基因数据库")
a=read.table("GSE34860_series_matrix.txt.gz",row.names = 1,header=T,sep="\t",quote="",fill=T,comment.char = "!")
X=t(a)
Y=as.matrix(data.frame(npm_status=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0)))
s.X=scale(X)#标准化变量X

#探针--基因对应表
library(plyr)
GPL<-read.delim("GPL96.annot",stringsAsFactors=FALSE,skip = 27 )
symbola<-data.frame(probe=GPL$ID,symbol=GPL$Gene.symbol,stringsAsFactors = FALSE)
symbola<-symbola[which(symbola$symbol!= ""),]#删去symbol为空的
symbolb<-strsplit(symbola$symbol,split = "///")# 分割字符串
names(symbolb)<-symbola$probe
symbolc<-ldply(symbolb,data.frame)#list转换为data.frame
colnames(symbolc) <- c("probe","symbol")#最终生成探针和基因的对应表

###SIS
SIS_X=SIS(s.X,Y)#最终选出的变量
SIS_X=as.matrix(gsub('[""]','',SIS_X))#最终选出的变量ID
###MMLE-SIS
MMLE_SIS_X=MMLE_SIS(s.X,Y)#最终选出的变量
MMLE_SIS_X=as.matrix(gsub('[""]','',MMLE_SIS_X))#最终选出的变量ID
###CSIS
gene_k<-as.matrix(c("NPM1","DNMT3A","FLT3","IDH1","IDH2","CEBPA","ASXL1","RUNX1","GATA2"))#已知致病基因
colnames(gene_k)='symbol'
work_gene_k<-merge(gene_k,symbolc,by='symbol')
Xc_name<-work_gene_k$probe#已知致病基因的对应探针
s.X_name<-as.matrix(gsub('[""]','',colnames(s.X)))#s.X的列变量探针名
loca<-NA
for (i in Xc_name){
  loca[i]<-grep(i, s.X_name)#每个已知探针所在s.X的位置
}
loca<-as.matrix(na.omit(loca))#已知探针所在s.X的位置
loca<-as.vector(loca)#存为向量
s.Xc=s.X[,loca]#已知探针
s.Xp<-s.X[ , !colnames(s.X) %in% colnames(as.data.frame(s.Xc))]#未知探针
CSIS_X=CSIS(s.Xc,s.Xp,Y)#最终选出的探针
CSIS_X=as.matrix(gsub('[""]','',CSIS_X))#最终选出的探针ID

#找到筛选出的探针对应的基因名(symbol)
#SIS
prb_SIS<-SIS_X
colnames(prb_SIS)='probe'
work_SIS<-merge(prb_SIS,symbolc,by='probe')
#MMLE-SIS
prb_MMLE_SIS<-MMLE_SIS_X
colnames(prb_MMLE_SIS)='probe'
work_MMLE_SIS<-merge(prb_MMLE_SIS,symbolc,by='probe')
#CSIS
prb_CSIS<-CSIS_X
colnames(prb_CSIS)='probe'
work_CSIS<-merge(prb_CSIS,symbolc,by='probe')