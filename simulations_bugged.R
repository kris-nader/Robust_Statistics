source("robust_stat.R")
library(MASS)

set.seed(1)

simulate_precision<-function(n,p,M, bool.contaminated = 1) {
  #Theta1 will be an n x p = 50 x 100 matrix
  Theta1 <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Theta1[i,j] = 0.6^(abs(i-j))
    }
  }
  rm(i,j)
  # Our covariance matrix is Theta1 inverted.
  Theta1.inv <- chol2inv(chol(Theta1)) #inv doesn't actually exist when p>n, but whatever.
  # We sample by using X = mu + (U_Theta)^-1 Z
  # Here mu stands for the mean and will be zero, Z stands for the standard
  # multivariate normal distribution and U_Theta stands for the
  # Upper triangular matrix in the Cholesky decomposition of Theta:
  # Theta = L_Theta%*%U_Theta (L_Theta = t(U_Theta))
  if (p>n) {
    sample.results<-matrix(0,M,6) #6 instead of 7 as we don't need to invert sample covariance.    
    colnames(sample.results) <- c("Glasso", "Twostep Quadrant", "Threestep Quadrant", "Twostep Spearman", "Threestep Spearman", "Gaussian")
  } else {
    sample.results<-matrix(0,M,7)
    colnames(sample.results) <- c("Glasso", "Twostep Quadrant", "Threestep Quadrant", "Twostep Spearman", "Threestep Spearman", "Gaussian", "Sample Covariance")
  }
  
  for (i in 1:M) {
    print("one")
    Z <- matrix(rnorm(p*n),p,n) #p*n standarnormal samples in a pxn matrix
    if (bool.contaminated == 1) {
      indices <- sample.int(p*n,p*n*0.1)
      for (j in indices) {
        Z[j] <- rnorm(1, mean = 10, sd = sqrt(0.2))
      }
    }
    U <- chol(Theta1) # By default R's chol returns upper Cholesky factor
    # The following is one sample (we need M)
    X <- t(backsolve(U,Z)) # more efficient and stable than actually inverting  
    
    print("two")
    #glasso based on sample covariance matrix
    S <- cov(X)
    sample.theta <- theta.sparse(S, n)
    sample.theta.KL <- kullback.leibler(Theta1.inv,sample.theta)
    sample.results[i,1] = sample.theta.KL
    print("three")
    #twostep quadrant
    print("three.one")
    twostep.quadrant <- quadrant.untransformed(X)
    print("three.two")
    tws.quad.theta <- theta.sparse(twostep.quadrant,n)
    print("three.three")
    tws.quad.theta.KL <- kullback.leibler(Theta1.inv,tws.quad.theta)
    print("three.four")
    sample.results[i,2] = tws.quad.theta.KL
    print("four")
    #threestep quadrant npd
    threestep.quadrant <- quadrant.transformed(X, method="npd")
    ths.quad.theta <- theta.sparse(threestep.quadrant, n)
    ths.quad.theta.KL <- kullback.leibler(Theta1.inv,ths.quad.theta)
    sample.results[i,3] = ths.quad.theta.KL
    print("five")
    #twostep spearman
    twostep.spearman <- spearman.untransformed(X)
    tws.sp.theta <- theta.sparse(twostep.spearman, n)
    tws.sp.theta.KL <- kullback.leibler(Theta1.inv,tws.sp.theta)
    sample.results[i,4] = tws.sp.theta.KL
    print("six")
    #threestep spearman
    threestep.spearman <- spearman.transformed(X, method="npd")
    ths.sp.theta <- theta.sparse(threestep.spearman, n)
    ths.sp.theta.KL <- kullback.leibler(Theta1.inv,ths.sp.theta)
    sample.results[i,5] = ths.sp.theta.KL
    print("seven")
    #Gaussian
    twostep.gaussian <- Grank(X)
    tws.gauss.theta <- theta.sparse(twostep.gaussian, n)
    tws.gauss.theta.KL <- kullback.leibler(Theta1.inv,tws.gauss.theta)
    sample.results[i,6] = tws.gauss.theta.KL
    print("eight")
    if (p<=n) {
      #precision matrix obtained by inverting S
      inv.theta <- chol2inv(chol(S))
      inv.theta.KL <- kullback.leibler(Theta1,inv.theta)
      sample.results[i,7] = inv.theta.KL
    }
    print("nine")
    if (i%%5==0) {
      print(i)
    }
  }
  avg.results <- colMeans(sample.results)  
  return(avg.results)
}

avg.results<-simulate_precision(50,100,10,1)





