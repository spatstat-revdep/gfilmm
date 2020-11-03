N <- 10000
WEIGHT <- readBin("inst/julia/weight.bin", "double", n = N)
VERTEX <- matrix(readBin("inst/julia/vertex.bin", "double", n = 3*N), 3, N)

gfi <- list(
  VERTEX = setNames(as.data.frame(t(VERTEX)), c("mu", "sigma_group", "sigma_error")),
  WEIGHT = WEIGHT
)

library(gfilmm)

gfiConfInt(~ mu, gfi, conf = 0.5)
gfiConfInt(~ sigma_error, gfi, conf = 0.5)

dat <- data.frame(
  ylwr = c(2.0, 2.0, 3.0, 4.0, 4.0, 6.0),
  yupr = c(3.0, 3.0, 4.0, 5.0, 5.0, 7.0),
  group = gl(3,2)
)

gfiR <- gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 15000)
gfiSummary(gfiR, conf = 0.5)

