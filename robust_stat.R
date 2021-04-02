## Computation using the Robust and sparse estimation of the inverse
## covariance matrix using rank correlation measures as reference

install.packages("robustbase")
library("robustbase")

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
  x.q=apply(x,2,Qn)
  cor.quadrant=sin(pi*cor(x.s)/2)
  sigma.quadrant=diag(x.q)%*%cor.quadrant%*%diag(x.q)
  return(easy.psd(sigma.quadrant,method))
}




