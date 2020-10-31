path <- "I:/JULIA_PC_GSK/julia-7076ab06f1"
setwd(path)

ewcdf <- spatstat::ewcdf


dd <- read.csv("juliavertex.csv", header=FALSE)
VERTEX <- as.matrix(dd)
dd <- read.csv("juliaweight.csv", header=FALSE)
WEIGHT <- as.matrix(dd)[1,]


v <- 2
vertJ <- VERTEX[v,]
weightJ <- WEIGHT
	final.ecdfJ <- ewcdf(vertJ,weightJ)
plot(final.ecdfJ, xlim=c(-8,5), main=paste("Vertex",v))

