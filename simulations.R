source("robust_stat.R")
library(MASS)

set.seed(1)
#We will try to estimate a banded precision matrix.
n <- 50
M <- 100
#Theta0 will be a p x p = 5 x 5 matrix with Theta_{i,j} = 0.6^abs(i-j)
p <- 5
Theta0 <- matrix(0, nrow = p, ncol = p) 
for (i in 1:p) {
  for (j in 1:p) {
    Theta0[i,j] = 0.6^(abs(i-j)) # Makes it banded
  }
} 
# We need to sample N(0, S)
# Our covariance matrix S is Theta0 inverted.
Theta0.inv <- solve(Theta0) # but this is not very robust
# We sample by using X = mu + (U_Theta)^-1 Z
# Here mu stands for the mean and will be zero, Z stands for the standard
# multivariate normal distribution and U_Theta stands for the
# Upper triangular matrix in the Cholesky decomposition of Theta:
# Theta = L_Theta%*%U_Theta (L_Theta = t(U_Theta))
sample.results<-matrix(0,M,7)
colnames(sample.results) <- c("Glasso", "Sample Covariance", "Twostep Quadrant", "Threestep Quadrant", "Twostep Spearman", "Threestep Spearman", "Gaussian")
bool.contaminated = 1
for (i in 1:M) {
  Z <- matrix(rnorm(p*n),p,n) #p*n standarnormal samples in a pxn matrix
  if (bool.contaminated == 1) {
    indices = sample.int(250,25)
    for (i in indices) {
      Z[i] = rnorm(1, mean = 10, sd = sqrt(0.2))
    }
  }
  U <- chol(Theta0) # By default R's chol returns upper Cholesky factor
  # The following is one sample (we need M)
  X <- t(backsolve(U,Z)) # more efficient and stable than actually inverting  
  
  #glasso based on sample covariance matrix
  S <- cov(X)
  sample.theta <- theta.sparse(S, n)
  sample.theta.KL <- kullback.leibler(Theta0.inv,sample.theta)
  sample.results[i,1] = sample.theta.KL
  #precision matrix obtained by inverting S
  inv.theta <- solve(S)
  inv.theta.KL <- kullback.leibler(Theta0,inv.theta)
  sample.results[i,2] = inv.theta.KL
  #twostep quadrant
  twostep.quadrant <- quadrant.untransformed(X)
  tws.quad.theta <- theta.sparse(twostep.quadrant,n)
  tws.quad.theta.KL <- kullback.leibler(Theta0.inv,tws.quad.theta)
  sample.results[i,3] = tws.quad.theta.KL
  #threestep quadrant npd
  threestep.quadrant <- quadrant.transformed(X, method="npd")
  ths.quad.theta <- theta.sparse(threestep.quadrant, n)
  ths.quad.theta.KL <- kullback.leibler(Theta0.inv,ths.quad.theta)
  sample.results[i,4] = ths.quad.theta.KL
  #twostep spearman
  twostep.spearman <- spearman.untransformed(X)
  tws.sp.theta <- theta.sparse(twostep.spearman, n)
  tws.sp.theta.KL <- kullback.leibler(Theta0.inv,tws.sp.theta)
  sample.results[i,5] = tws.sp.theta.KL
  #threestep spearman
  threestep.spearman <- spearman.transformed(X, method="npd")
  ths.sp.theta <- theta.sparse(threestep.spearman, n)
  ths.sp.theta.KL <- kullback.leibler(Theta0.inv,ths.sp.theta)
  sample.results[i,6] = ths.sp.theta.KL
  #Gaussian
  twostep.gaussian <- Grank(X)
  tws.gauss.theta <- theta.sparse(twostep.gaussian, n)
  tws.gauss.theta.KL <- kullback.leibler(Theta0.inv,tws.gauss.theta)
  sample.results[i,7] = tws.gauss.theta.KL
}


p <- 100
#Theta1 will be an n x p = 50 x 100 matrix
Theta1 <- matrix(0, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    Theta1[i,j] = 0.6^(abs(i-j))
  }
}
rm(i,j)
# Our covariance matrix is Theta1 inverted.
Theta1.inv <- solve(Theta1) #inv doesn't actually exist when p>n, but whatever.
# We sample by using X = mu + (U_Theta)^-1 Z
# Here mu stands for the mean and will be zero, Z stands for the standard
# multivariate normal distribution and U_Theta stands for the
# Upper triangular matrix in the Cholesky decomposition of Theta:
# Theta = L_Theta%*%U_Theta (L_Theta = t(U_Theta))
sample.results100<-matrix(0,M,6) #6 instead of 7 as we don't need to invert sample covariance.
colnames(sample.results100) <- c("Glasso", "Twostep Quadrant", "Threestep Quadrant", "Twostep Spearman", "Threestep Spearman", "Gaussian")
A <- Sys.time()
bool.contaminated = 1
for (i in 1:M) {
  Z <- matrix(rnorm(p*n),p,n) #p*n standarnormal samples in a pxn matrix
  if (bool.contaminated == 1) {
    indices <- sample.int(500,50)
    for (j in indices) {
      Z[j] <- rnorm(1, mean = 10, sd = sqrt(0.2))
    }
  }
  U <- chol(Theta1) # By default R's chol returns upper Cholesky factor
  # The following is one sample (we need M)
  X <- t(backsolve(U,Z)) # more efficient and stable than actually inverting  
  
  
  #glasso based on sample covariance matrix
  S <- cov(X)
  sample.theta <- theta.sparse(S, n)
  sample.theta.KL <- kullback.leibler(Theta1.inv,sample.theta)
  sample.results100[i,1] = sample.theta.KL
  #twostep quadrant
  twostep.quadrant <- quadrant.untransformed(X)
  tws.quad.theta <- theta.sparse(twostep.quadrant,n)
  tws.quad.theta.KL <- kullback.leibler(Theta1.inv,tws.quad.theta)
  sample.results100[i,2] = tws.quad.theta.KL
  #threestep quadrant npd
  threestep.quadrant <- quadrant.transformed(X, method="npd")
  ths.quad.theta <- theta.sparse(threestep.quadrant, n)
  ths.quad.theta.KL <- kullback.leibler(Theta1.inv,ths.quad.theta)
  sample.results100[i,3] = ths.quad.theta.KL
  #twostep spearman
  twostep.spearman <- spearman.untransformed(X)
  tws.sp.theta <- theta.sparse(twostep.spearman, n)
  tws.sp.theta.KL <- kullback.leibler(Theta1.inv,tws.sp.theta)
  sample.results100[i,4] = tws.sp.theta.KL
  #threestep spearman
  threestep.spearman <- spearman.transformed(X, method="npd")
  ths.sp.theta <- theta.sparse(threestep.spearman, n)
  ths.sp.theta.KL <- kullback.leibler(Theta1.inv,ths.sp.theta)
  sample.results100[i,5] = ths.sp.theta.KL
  #Gaussian
  twostep.gaussian <- Grank(X)
  tws.gauss.theta <- theta.sparse(twostep.gaussian, n)
  tws.gauss.theta.KL <- kullback.leibler(Theta1.inv,tws.gauss.theta)
  sample.results100[i,6] = tws.gauss.theta.KL
  if (i%%5==0) {
    print(i)
  }
}
B <- Sys.time()
time = B-A
time
avg.results100 <- colMeans(sample.results100)
avg.results100


