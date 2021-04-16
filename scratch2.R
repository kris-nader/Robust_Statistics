n<-50
p<-15

Theta1 <- matrix(0, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    Theta1[i,j] = 0.6^(abs(i-j))
  }
}

rm(i,j)
Sigma<-chol2inv(chol(Theta1))
set.seed(1)

U <- chol(Theta1) # By default R's chol returns upper Cholesky factor
Z <- matrix(rnorm(p*n),p,n) #p*n standarnormal samples in a pxn matrix
# The following is one sample (we need M)
X <- t(backsolve(U,Z)) # more efficient and stable than actually inverting  

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

S.hat <- quadrant.transformed(X,"npd")

huge.out<-huge(S.hat,method="glasso",verbose=F)
# print(2)
my.bic= rep(0,10)
for (i in 1:10) {
  Theta.hat<-huge.out$icov[[i]]
  my.bic[i] = -log(det(Theta.hat)) + sum(diag(Theta.hat%*%S.hat)) +sum(Theta.hat[upper.tri(Theta.hat,diag=TRUE)]!=0)*log(n)/n
}
my.bic
# print(3)
opt.i=which.min(my.bic)
# print(4)
A<-huge.out$icov[[opt.i]]
B<-huge.out$icov[[2]]
