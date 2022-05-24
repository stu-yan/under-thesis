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