#######################################################################
######### The code below is for our proposed method OWC################
#########in our paper "An Optimally Weighted Combination Method########
#########to Detect Novel Disease Associated Genes Using Publicly####### 
#########Available GWAS Summary Data". If you have any questions,###### 
#########please contact Xuexia Wang at xuexia.wang@unt.edu#############
#######################################################################

####Zs denotes Z-statistics, R denotes LD matrix among SNPs, 
####n.perm denotes the number of permutation, 
####rho denotes the range of rho (ncol=4), W denotes the weightes of WSS. 
library(MASS)
OWC=function(Zs,R,n.perm,rho,W){
  k <- length(Zs)
  Z=matrix(Zs,ncol=1)
  
  eS <- eigen(R, symmetric = TRUE)
  ev <- eS$values
  CovSsqrt <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), k)
  
  a1=matrix(rep(1,k),ncol=1)%*%matrix(rep(1,k),nrow=1)
  a2=matrix(W,ncol=1)%*%matrix(W,nrow=1)
  a3=ginv(R)
  a4=diag(1,k) 
  
  Pva=matrix(0,nrow=nrow(rho),ncol=n.perm+1)
  s <- sample(1:10^5, 1)
  for(i in 1:nrow(rho)){
    set.seed(s)
    A=rho[i,1]*a1+rho[i,2]*a2+rho[i,3]*a3+rho[i,4]*a4
    ev <- eigen(R%*%A, symmetric = TRUE)$values
    a=sum(ev^3)/sum(ev^2)
    b=sum(ev)-(sum(ev^2))^2/sum(ev^3)
    d=(sum(ev^2))^3/(sum(ev^3))^2
    
    T0=t(Z)%*%A%*%Z
    Pva[i,1]=1-pchisq((T0-b)/a,df=d)
    for (bb in 1:n.perm){
      U00 <- rnorm(k, 0, 1)
      U0 <- CovSsqrt %*% U00
      T0=t(U0)%*%A%*%U0
      Pva[i,bb+1]=1-pchisq((T0-b)/a,df=d)
    }
  }
  Ts=apply(Pva,2,min)
  Pvs <- (sum(Ts[-1] <= Ts[1]) + 1)/(n.perm + 1)
  return(Pvs)
}  
