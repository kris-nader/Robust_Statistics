## Computation using the Robust and sparse estimation of the inverse
## covariance matrix using rank correlation measures as reference

#if you have errors, uncomment and install these
#install.packages("robustbase")
#install.packages("huge")

library("robustbase")
library("huge")
library(Rcpp)
library(RcppArmadillo)
# library(Matrix)
sourceCpp("test.cpp")

# Function to compute the perturbed (method='perturb') sample covariance matrix and
# the nearest positive definite matrix by eliminating the negative eigenspaces (method='npd')
easy.psd<-function(sigma, method="perturb"){
  if (method=="perturb"){
    p<-ncol(sigma)
    eig<-eigen(sigma, symmetric = T, only.values = T) #only.values returns only eigenvalues
    const<-abs(min(eig$values,0))
    sigma.psd<-sigma+diag(p)*const
  }
  if(method=="npd"){
    eig<-eigen(sigma,symmetric = T)
    # print(8)
    d<-pmax(eig$values,0)
    # print(9)
    sigma.psd=eig$vectors%*%diag(d)%*%t(eig$vectors)
  }
  return (sigma.psd)
}

# Function that computes the robust sample covariance matrix with quadrant rank
# correlation and made positive semidefinite.
quadrant.transformed<-function(x,method="perturb"){
  # print(1)
  x.m=apply(x,2,median)
  # print(2)
  x=sweep(x,2,x.m)
  # print(3)
  x.s=sign(x)
  # print(4)
  x.q=apply(x,2,Qn)# here Qn is from the package robustbase
  # print(5)
  cor.quadrant=sin(pi*cor(x.s)/2)
  # print(6)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  # print(7)
  return(easy.psd(sigma.quadrant,method))
}

# Function that computes the robust sample covariance matrix with quadrant rank
# correlation.
quadrant.untransformed<-function(x){
  x.m=apply(x,2,median)
  x=sweep(x,2,x.m)
  x.s=sign(x)
  x.q=apply(x,2,Qn)
  cor.quadrant=cor(x.s)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  return(sigma.quadrant)
}

# Function that computes the robust sample covariance matrix with spearman rank
# correlation and made positive semidefinite.
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

# Function that computes the robust sample covariance matrix with spearman rank
# correlation.
spearman.untransformed<-function(x)
{
  x.r=apply(x,2,rank)
  x.q=apply(x,2,Qn)
  cor.sp=cor(x.r)
  sigma.sp=diag(x.q)%*%cor.sp%*%diag(x.q)
  return(sigma.sp)
}

# Function that computes the robust sample covariance matrix with Gaussian rank
# correlation.
Grank<-function(x)
{
  # print(1)
  n=nrow(x)
  # print(2)
  x.q=apply(x,2,Qn)
  # print(3)
  x.r=apply(x,2,rank)
  # print(4)
  cor.Grank=cor(qnorm(x.r/(n+1)))
  # print(5)
  sigma.quadrant=diag(x.q)%*%cor.Grank%*%diag(x.q)
  # print(6)
  return(sigma.quadrant)
}

# Function that computes the precision matrix given a (robust) sample covariance
# matrix.
theta.sparse<-function(sigma.psd,n)
{
  # print(1)
  huge.out<-huge(sigma.psd,method="glasso",verbose=F)
  # print(2)
  my.bic= rep(0,10)
  for (i in 1:10) {
    # print(i)
    Theta.hat<-huge.out$icov[[i]]
    ID.hat<-eigenMapMatMult2(Theta.hat, sigma.psd, n_cores = 12)
    my.bic[i] <- -log(det(Theta.hat)) + sum(diag(ID.hat)) +sum(Theta.hat[upper.tri(Theta.hat,diag=TRUE)]!=0)*log(n)/n
    # my.bic[i] = (-huge.out$loglik[[i]]) +sum(A[upper.tri(A,diag=TRUE)]!=0)*log(n)/n
  }
  # print(3)
  print(my.bic)
  opt.i=which.min(my.bic)
  print(opt.i)
  # print(4)
  return(huge.out$icov[[opt.i]])
}

# theta.sparse<-function(sigma.psd,n)
# {
#   # print(1)
#   huge.out<-huge(sigma.psd,method="glasso",verbose=T)
#   # print(2)
#   my.bic= rep(0,10)
#   for (i in 1:10) {
#     print(i)
#     Theta.hat<-huge.out$icov[[i]]
#     my.bic[i] <- -log(det(Theta.hat)) + sum(diag(Theta.hat%*%sigma.psd)) +sum(Theta.hat[upper.tri(Theta.hat,diag=TRUE)]!=0)*log(n)/n
#     # my.bic[i] = (-huge.out$loglik[[i]]) +sum(A[upper.tri(A,diag=TRUE)]!=0)*log(n)/n
#   }
#   # print(3)
#   opt.i=which.min(my.bic)
#   # print(4)
#   return(huge.out$icov[[opt.i]])
# }

kullback.leibler<-function(Theta.inv, Theta.est) {
  #  Theta.inv <- solve(Theta) ##I decided to get rid of this line as its probably computationally expensive
  ID<-Theta.inv%*%Theta.est
  return(sum(diag(ID))-log(det(ID))-ncol(Theta.inv))
}
