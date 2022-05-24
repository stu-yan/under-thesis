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