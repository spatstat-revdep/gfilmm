library(gfilmm)

nsims <- 2L
confint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_CV           <- matrix(NA_real_, nrow = nsims, ncol = 2L)

mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
n            <- 8L # sample size per group

for(i in 1L:nsims){
  means <- rnorm(2L, mu, sigmaBetween)
  y1    <- rnorm(n, means[1L], sigmaWithin) 
  y2    <- rnorm(n, means[2L], sigmaWithin) 
  y     <- round(c(y1, y2), digits = 1)
  dat   <- data.frame(
    ylwr = y - 0.05,
    yupr = y + 0.05,
    group <- gl(2L, n)
  )
  gfi <- gfilmm(~ cbind(ylwr, yupr), fixed = ~ 1, random = ~ group, 
                data = dat, N = 200)
  confint_grandMean[i,]    <- gfiConfInt(~ `(Intercept)`, gfi)
  confint_sigmaBetween[i,] <- gfiConfInt(~ sigma_group, gfi)
  confint_sigmaWithin[i,]  <- gfiConfInt(~ sigma_error, gfi)
  confint_sigmaTotal[i,]   <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2), gfi)
  confint_CV[i,]           <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2)/`(Intercept)`, gfi)
}
