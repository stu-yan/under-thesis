####函数SIS
#点二列相关系数
rpb<-function(X,Y){
  y1<-which(Y==1,arr.ind=TRUE)[,1]
  y0<-which(Y==0,arr.ind=TRUE)[,1]
  X1<-0
  for (i in y1){
    X1<-X1+X[i]
  }
  m_X1<-X1/length(y1)
  z1<-length(y1)/length(Y)
  X0<-0
  for (i in y0){
    X0<-X0+X[i]
  }
  m_X0<-X0/length(y0)
  z0<-length(y0)/length(Y)
  St<-sd(X)
  r_pb<-((m_X1-m_X0)/St)*sqrt(z1*z0)
  return(r_pb)
}

#通过SIS函数计算点二列相关系数绝对值并排序选出重要变量
SIS<-function(X,Y){
  library(Matrix)
  library(glmnet)
  n=nrow(X)#样本量
  p=ncol(X)#变量
  c=as.data.frame(matrix(data = NA,nrow = p,ncol = 1))
  for(i in 1:p)
  {
    c[i,]=abs(rpb(X[,i],Y))#计算点二列相关系数的绝对值
    rownames(c)[i]<-colnames(X)[i]
  }
  abs_c<-as.matrix(c)
  r=n-1#计算选择变量数
  c_order=sort(abs_c[,1],decreasing=T)#从大到小排序
  v=c_order[1:r]#选取排在前r个变量，即重要变量
  v_names<-names(v)#输出重要变量名称
  Xr=X[,v_names]#选出X中的重要变量
  X_name=colnames(Xr)#最终选出的变量
  return(X_name)
}

####函数MMLE_SIS
#通过MMLE-SIS函数计算最大似然估计绝对值并排序选出重要变量 
MMLE_SIS<-function(X,Y){
  library(Matrix)
  library(glmnet)
  n=nrow(X)#样本量
  p=ncol(X)#变量数
  beta=as.data.frame(matrix(data = NA,nrow = p,ncol = 1))
  #计算beta最大似然估计值
  f1<-function(xx){
    a<-which(names(X)==colnames(xx))
    rownames(beta)[a]<-colnames(xx)
    Y_logistic=glm(Y~xx,data=data.frame(cbind(xx,Y)),family = binomial(link="logit"))
    beta[a]<-abs(Y_logistic$coeff[2])
  }
  beta0<-apply(X,2,f1)#向量化计算
  abs_beta<-as.matrix(beta0)
  beta_order=sort(abs_beta[,1],decreasing=T)#从大到小排序
  r=n-1#计算选择变量数，这里选择的是SIS选择的变量个数
  v=beta_order[1:r]#选取排在前r个变量，即重要变量
  v_names<-names(v)#输出重要变量名称
  Xr=X[,v_names]#选出X中的重要变量
  X_name=colnames(Xr)
  return(X_name)
}

####函数CSIS
#通过CSIS函数计算最大似然估计绝对值并排序选出重要变量 
CSIS<-function(Xc,Xp,Y){
  library(Matrix)
  library(glmnet)
  library(dplyr)
  n=nrow(Xp)#样本量
  p=ncol(Xp)#未知变量数
  c=ncol(Xc)#已知变量数
  beta=as.data.frame(matrix(data = NA,nrow = p,ncol = 1))
  #计算已知Xc为重要变量的情况下的其他beta估计值
  f2<-function(xx){
    a<-which(names(Xp)==colnames(xx))
    rownames(beta)[a]<-colnames(xx)
    Y_logistic=glm(Y~Xc+xx,data=data.frame(cbind(cbind(Xc,xx),Y)),family = binomial(link="logit"))
    beta[a]<-abs(Y_logistic$coeff[c+2])
  }
  beta0<-apply(Xp,2,f2)#向量化计算
  abs_beta<-as.matrix(beta0)
  beta_order=sort(abs_beta[,1],decreasing=T)#从大到小排序
  r=n-1#计算选择变量数
  v=beta_order[1:r]#选取排在前r个变量，即重要变量
  v_names<-names(v)#输出重要变量名称
  Xr=X[,v_names]#选出X中的重要变量
  X_name=colnames(Xr)
  return(X_name)
}

####模拟数据集函数
generate_X<-function(q,b,n,r){
  library(MASS)
  library(Matrix)
  mu<-rep(0,q*b)#多元高斯分布的均值向量（均值都为0）# 100维
  Sigma<-matrix(r,ncol=q,nrow = q)#初始化协方差阵为全部元素r的矩阵
  diag(Sigma)<-1#协方差阵的对角线更正为1
  Sigmas<-as.matrix(bdiag(lapply(1:b,function(o) Sigma)))#生成大协方差矩阵，由Sigma组成的对角矩阵块
  X<-mvrnorm(n=n,mu=mu,Sigma=Sigmas)#产生服从N(0,Sigmas)的随机数
  return(X)
}
generate_Y<-function(X){
  library(MASS)
  #X:设计矩阵
  #Y=-x1+2x21-3x41+4x61-5x81+e
  #根据上面公式提取出beta系数
  beta<-rep(0,ncol(X))#初始化beta系数为与X相同列数的0向量
  idx<-seq(1,100,20)#生成1,21,41,61,81
  beta[idx]<-c(-1,2,-3,4,-5)#更正每组第一个变量的系数
  Y<-X%*%beta+rnorm(nrow(X))#X.beta+E即为得到Y，E为标准正态随机数(rnorm()生成)
  return(Y)
}

###模拟###
TPR_SIS<-0
TPR_MMLE_SIS<-0
TPR_CSIS<-0
FPR_SIS<-0
FPR_MMLE_SIS<-0
FPR_CSIS<-0
M<-c("V1","V21","V41","V61","V81")

memory.limit(100000)
##100次模拟##
for (j in 1:100){
  ###模拟数据集
  b<-100# b:同分布的随机变量组数100
  q<-40# q:每组随机变量个数10/20/40
  p<-100
  n<-400# n:样本量100,200,400
  r<-0.6# r:同组随机变量间的任意两个随机变量的相关系数0.4/0.6
  X<-generate_X(q,b,n,r)
  y<-generate_Y(X)
  p=1/(1+exp(-y))
  Y=rep(0,nrow(p))
  for (i in 1:nrow(p)){
    Y[i]=rbinom(1,1,p[i])
  }
  Y=as.matrix(Y)
  colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  s.X=scale(X)#标准化变量X
  
  ###SIS
  SIS_X=SIS(s.X,Y)#最终选出的变量
  t=0
  for (i in SIS_X){
    if(i%in%M){
      t=t+1
    }
  }
  TPR_S=t/length(M)#真正率
  TPR_SIS=TPR_SIS+TPR_S
  FPR_S=(length(SIS_X)-t)/(ncol(X)-length(M))#假阳性率
  FPR_SIS=FPR_SIS+FPR_S
  
  ###MMLE-SIS
  MMLE_SIS_X=MMLE_SIS(s.X,Y)#最终选出的变量
  t=0
  for (i in MMLE_SIS_X){
    if(i%in%M){
      t=t+1
    }
  }
  TPR_MMLE_S=t/length(M)#真正率
  TPR_MMLE_SIS=TPR_MMLE_SIS+TPR_MMLE_S
  FPR_MMLE_S=(length(MMLE_SIS_X)-t)/(ncol(X)-length(M))#假阳性率
  FPR_MMLE_SIS=FPR_MMLE_SIS+FPR_MMLE_S
  
  ###CSIS
  s.Xc=s.X[,c('V1','V2','V3')]#已知变量
  s.Xp<-s.X[ , !colnames(s.X) %in% colnames(as.data.frame(s.Xc))]#未知变量
  CSIS_X=CSIS(s.Xc,s.Xp,Y)#最终选出的变量
  t=0
  for (i in CSIS_X){
    if(i%in%M){
      t=t+1
    }
  }
  TPR_CS=t/(length(M)-1)#真正率
  TPR_CSIS=TPR_CSIS+TPR_CS
  FPR_CS=(length(CSIS_X)-t)/(ncol(X)-length(M))#假阳性率
  FPR_CSIS=FPR_CSIS+FPR_CS
}


#计算100次模拟后的平均真正率TPR
T_S=TPR_SIS/100
T_M=TPR_MMLE_SIS/100
T_C=TPR_CSIS/100
#计算100次模拟后的平均假阳性率FPR
F_S=FPR_SIS/100
F_M=FPR_MMLE_SIS/100
F_C=FPR_CSIS/100






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



