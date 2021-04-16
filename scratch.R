# source("robust_stat.R")
library(robustbase)
library(huge)
library(Rcpp)
library(RcppArmadillo)
# library(Matrix)
sourceCpp("test.cpp")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("golubEsets")
# 
# BiocManager::install("multtest")

library(golubEsets)
# library("multtest")

data("Golub_Merge")

data<-t(Golub_Merge@assayData$exprs)
rm(Golub_Merge)
gc()


easy.psd<-function(sigma, method="perturb"){
  if (method=="perturb"){
    p<-ncol(sigma)
    eig<-eigen(sigma, symmetric = T, only.values = T) #only.values returns only eigenvalues
    const<-abs(min(eig$values,0))
    sigma.psd<-sigma+diag(p)*const
  }
  if(method=="npd"){
    eig<-eigen(sigma,symmetric = T)
    print(8)
    d<-pmax(eig$values,0)
    print(9)
    sigma.psd=eig$vectors%*%diag(d)%*%t(eig$vectors)
  }
  return (sigma.psd)
}

theta.sparse<-function(sigma.psd,n)
{
  # print(1)
  huge.out<-huge(sigma.psd,method="glasso",verbose=T)
  # print(2)
  my.bic= rep(0,10)
  for (i in 1:10) {
    print(i)
    Theta.hat<-huge.out$icov[[i]]
    ID.hat<-eigenMapMatMult2(Theta.hat, sigma.psd, n_cores = 12)
    my.bic[i] <- -log(det(Theta.hat)) + sum(diag(ID.hat)) +sum(Theta.hat[upper.tri(Theta.hat,diag=TRUE)]!=0)*log(n)/n
    # my.bic[i] = (-huge.out$loglik[[i]]) +sum(A[upper.tri(A,diag=TRUE)]!=0)*log(n)/n
  }
  # print(3)
  opt.i=which.min(my.bic)
  # print(4)
  return(huge.out$icov[[opt.i]])
}

# Function that computes the robust sample covariance matrix with quadrant rank
# correlation and made positive semidefinite.
quadrant.transformed<-function(x,method="perturb"){
  print(1)
  x.m=apply(x,2,median)
  print(2)
  x=sweep(x,2,x.m)
  print(3)
  x.s=sign(x)
  print(4)
  x.q=apply(x,2,Qn)# here Qn is from the package robustbase
  print(5)
  cor.quadrant=sin(pi*cor(x.s)/2)
  print(6)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  print(7)
  return(easy.psd(sigma.quadrant,method))
}



# quadrant 3-step
threestep.quadrant <- quadrant.transformed(data, method="npd")
write.table( threestep.quadrant,"golub_quadrant_S.csv" , sep=",",  row.names=FALSE, col.names=FALSE)
nrows=nrow(data)
rm(data)
gc()
threestep.quadrant<-as.matrix(scan("golub_quadrant_S2.csv",sep=","))
threestep.quadrant<-threestep.quadrant[,-1]
ths.quad.theta <- theta.sparse(threestep.quadrant, n=nrows)
write.table( ths.quad.theta,"golub_quadrant_Theta.csv" , sep=",",  row.names=FALSE, col.names=FALSE)
image(ths.quad.theta[1:250,1:250])

# spearman 3-step
spearman.transformed<-function(x,method="perturb")
{
  # print(1)
  x.r=apply(x,2,rank)
  # print(2)
  x.q=apply(x,2,Qn)
  # print(3)
  cor.sp=2*sin(pi*cor(x.r)/6)
  # print(4)
  sigma.sp=diag(x.q)%*%cor.sp%*%diag(x.q)
  # print(5)
  return(easy.psd(sigma.sp,method))
}

threestep.spearman <- spearman.transformed(data, method="npd")
write.table(threestep.spearman, "golub_spearman_S.csv", sep=",",  row.names=FALSE, col.names=FALSE)
ths.sp.theta <- theta.sparse(threestep.spearman, n=nrows)
write.table(ths.sp.theta, "golub_spearman_theta.csv", sep=",",  row.names=FALSE, col.names=FALSE)

# Gauss
Grank<-function(x)
{
  print(1)
  n=nrow(x)
  print(2)
  x.q=apply(x,2,Qn)
  print(3)
  x.r=apply(x,2,rank)
  print(4)
  cor.Grank=cor(qnorm(x.r/(n+1)))
  print(5)
  sigma.quadrant=diag(x.q)%*%cor.Grank%*%diag(x.q)
  print(6)
  return(sigma.quadrant)
}

twostep.gaussian <- Grank(data)
write.table(twostep.gaussian, "golub_gauss_S.csv", sep=",",  row.names=FALSE, col.names=FALSE)
tws.gauss.theta <- theta.sparse(twostep.gaussian, n=nrows)
write.table(tws.gauss.theta, "golub_gauss_theta.csv", sep=",",  row.names=FALSE, col.names=FALSE)
