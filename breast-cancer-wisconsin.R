source("robust_stat.R")
library(MASS)
data = read.csv("breast-cancer-wisconsin.data",header=F, na="?")

# See if it's safe with NAs
NAs = is.na(data)
count = 0
for (i in 1:699) {
  if (is.element(TRUE,NAs[i,])) {
    count = count + 1
  }
} # count = 16 observations with an NA out of 699 observations.
# Let's just omit these observations
data <- na.omit(data)
data <- data[,-1]
# We can now continue
data.cov <- cov(data)
data.prec <- solve(data.cov) #not singular

# quadrant 3-step
threestep.quadrant <- quadrant.transformed(data, method="npd")
ths.quad.theta <- theta.sparse(threestep.quadrant, n=nrow(data))

# spearman 3-step
threestep.spearman <- spearman.transformed(data, method="npd")
ths.sp.theta <- theta.sparse(threestep.spearman, n=nrow(data))

# Gauss
twostep.gaussian <- Grank(data)
tws.gauss.theta <- theta.sparse(twostep.gaussian, n=nrow(data))

write.table(ths.quad.theta, file="breastW_quadrant.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(ths.sp.theta, file="breastW_spearman.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(tws.gauss.theta, file="breastW_gauss.csv", sep=",", row.names=FALSE, col.names=FALSE)
