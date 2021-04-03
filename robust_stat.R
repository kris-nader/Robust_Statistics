## Computation using the Robust and sparse estimation of the inverse
## covariance matrix using rank correlation measures as reference

install.packages("robustbase")
install.packages("huge")
library("robustbase")
library("huge")

easy.psd<-function(sigma, method="perturb"){
  if (method=="perturb"){
    p==ncol(sigma)
    eig=eigen(sigma, symmetric = T, only.values = T)
    const=abs(min(eig$values,0))
    sigma.psd=sigma+diag(p)*const
  }
  if(method=="npd"){
    eig=eigen(sigma,symmetric = T)
    d=pmax(eig$values,0)
    sigma.psd=eig$vectors%*%
      diag(d)%*%
      t(eig$vectors)
  }
  return (sigma.psd)
}

quadrant.transformed<-function(x,method="perturb"){
  x.m=apply(x,2,median)
  x=sweep(x,2,x,m)
  x.s=sign(x)
  x.q=apply(x,2,Qn)# here Qn is from the package robustbase
  cor.quadrant=sin(pi*cor(x.s)/2)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  return(easy.psd(sigma.quadrant,method))
}

spearman.transformed<-function(x,method="perturb")
{
  x.r=apply(x,2,rank)
  x.q=apply(x,2,Qn)
  cor.sp=2*sin(pi*cor(x.r)/6)
  sigma.sp=diag(x.q)%*%cor.sp%*%diag(x.q)
  return(easy.psd(sigma.sp,method))
}

Grank<-function(x)
{
  n=nrow(x)
  x.q=apply(x,2,Qn)
  x.r=apply(x,2,rank)
  cor.Grank=cor(qnorm(x.r/(n+1)))
  sigma.quadrant=diag(x.q)%*%cor.Grank%*%diag(x.q)
  return(sigma.quadrant)
}

theta.sparse<-function(sigma.psd,n)
{
  huge.out<-huge(sigma.psd,method="glasso",verbose=F)
  my.bic=-huge.out$loglik+huge.out$df*log(n)/n
  opt.i=which.min(my.bic)
  return(huge.out$icov[[opt.i]])
}


