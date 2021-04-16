source("robust_stat.R")
library(MASS)
data = read.csv("covtype.data",header=F)

data.cov <- cov(data)
data.prec <- solve(data.cov) #error: singular


# quadrant 3-step
threestep.quadrant <- quadrant.transformed(data, method="npd")
ths.quad.theta <- theta.sparse(threestep.quadrant, n=nrow(data))
A <- zapsmall(ths.quad.theta)
# spearman 3-step
threestep.spearman <- spearman.transformed(data, method="npd")
ths.sp.theta <- theta.sparse(threestep.spearman, n=nrow(data))

# Gauss
twostep.gaussian <- Grank(data)
tws.gauss.theta <- theta.sparse(twostep.gaussian, n=nrow(data))

write.table(ths.quad.theta, file="covtype_quadrant.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(ths.sp.theta, file="covtype_spearman.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(tws.gauss.theta, file="covtype_gauss.csv", sep=",", row.names=FALSE, col.names=FALSE)
