library(gfilmm)

nsims <- 1000L
confint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_CV           <- matrix(NA_real_, nrow = nsims, ncol = 2L)

mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
sigmaTotal   <- sqrt(sigmaBetween^2 + sigmaWithin^2)
CV           <- sigmaTotal / mu
n            <- 8L # sample size per group

set.seed(666L)
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
                data = dat, N = 10000)
  confint_grandMean[i,]    <- gfiConfInt(~ `(Intercept)`, gfi)
  confint_sigmaBetween[i,] <- gfiConfInt(~ sigma_group, gfi)
  confint_sigmaWithin[i,]  <- gfiConfInt(~ sigma_error, gfi)
  confint_sigmaTotal[i,]   <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2), gfi)
  confint_CV[i,]           <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2)/`(Intercept)`, gfi)
}


results <- list(
  grandMean = confint_grandMean,
  sigmaWithin = confint_sigmaWithin,
  sigmaBetween = confint_sigmaBetween,
  sigmaTotal = confint_sigmaTotal,
  CV = confint_CV
)

saveRDS(results, "~/Work/R/gfilmm/inst/essais/simulations01.rds")

################################################################################
results <- readRDS("~/Work/R/gfilmm/inst/essais/simulations01.rds")

cvrg_twoSided <- c(
  sum(results$grandMean[,1] <= mu & mu <= results$grandMean[,2]),
  sum(results$sigmaWithin[,1] <= sigmaWithin & sigmaWithin <= results$sigmaWithin[,2]),
  sum(results$sigmaBetween[,1] <= sigmaBetween & sigmaBetween <= results$sigmaBetween[,2]),
  sum(results$sigmaTotal[,1] <= sigmaTotal & sigmaTotal <= results$sigmaTotal[,2]),
  sum(results$CV[,1] <= CV & CV <= results$CV[,2]) 
) / nsims

cvrg_leftSided <- c(
  sum(mu <= results$grandMean[,2]),
  sum(sigmaWithin <= results$sigmaWithin[,2]),
  sum(sigmaBetween <= results$sigmaBetween[,2]),
  sum(sigmaTotal <= results$sigmaTotal[,2]),
  sum(CV <= results$CV[,2]) 
) / nsims

cvrg_rightSided <- c(
  sum(results$grandMean[,1] <= mu),
  sum(results$sigmaWithin[,1] <= sigmaWithin),
  sum(results$sigmaBetween[,1] <= sigmaBetween),
  sum(results$sigmaTotal[,1] <= sigmaTotal),
  sum(results$CV[,1] <= CV) 
) / nsims

