library(gfilmm)
library(AOV1R)
library(rstanarm)
options(mc.cores = parallel::detectCores())

mu           <- 10000 # grand mean
sigmaWithin  <- 1
ratio        <- 50
sigmaBetween <- sigmaWithin * ratio
I            <- 10L # number of groups
J            <- 5L # sample size per group

nsims <- 15L
fidIntervals <- freqIntervals <- stanIntervals <- list(
  within = `colnames<-`(matrix(NA_real_, nrow = nsims, ncol = 3), 
                      c("estimate", "lwr", "upr")),
  between = `colnames<-`(matrix(NA_real_, nrow = nsims, ncol = 3), 
                      c("estimate", "lwr", "upr"))
)

set.seed(666L)
for(i in 1L:nsims){
  groupmeans <- rnorm(I, mu, sigmaBetween)
  y          <- c(
    vapply(groupmeans, function(gmean) rnorm(J, gmean, sigmaWithin), numeric(J))
  )
  y_rounded  <- round(y, digits = 4L)
  dat        <- data.frame(
    y = y,
    ylwr = y_rounded - 5e-5,
    yupr = y_rounded + 5e-5,
    group = gl(I, J)
  )
  x <- confint(aov1r(y ~ group, data = dat))
  freqIntervals$within[i, ] <- unlist(x["within", ])  
  freqIntervals$between[i, ] <- unlist(x["between", ])  

  N <- 10000L
  gfi <- gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = N)
  x <- gfiSummary(gfi)
  fidIntervals$within[i, ] <- x["sigma_error", 2L:4L]  
  fidIntervals$between[i, ] <- x["sigma_group", 2L:4L]  
  
  stan <- stan_lmer(
    y ~ (1|group), data = dat, 
    prior_aux = cauchy(0, 5),
    prior_covariance = decov(1, 1, 1, 1), 
    iter = 5000L
  )
  pstrr <- as.data.frame(
    stan, 
    pars = c(
      "(Intercept)", 
      "sigma", 
      "Sigma[group:(Intercept),(Intercept)]"
    )
  )
  names(pstrr)[2L:3L] <- c("sigma_error",  "sigma_group")
  pstrr$sigma_total <- with(pstrr, sqrt(sigma_group + sigma_error^2))
  pstrr[["sigma_group"]] <- sqrt(pstrr[["sigma_group"]])
  x <- t(vapply(pstrr, quantile, numeric(3L), probs = c(50, 2.5, 97.5)/100))
  stanIntervals$within[i, ] <- x["sigma_error", ]  
  stanIntervals$between[i, ] <- x["sigma_group", ]  
}
#table(attr(gfi, "ESS")/(N/2) < 1)

Results <- list(
  freqIntervals = freqIntervals,
  fidIntervals = fidIntervals,
  stanIntervals = stanIntervals
)

saveRDS(Results, "~/Work/R/gfilmm/inst/stan/AOV1R_highRatio_simulations.rds")

stop()

# ####
Results <- readRDS("~/Work/R/gfilmm/inst/stan/AOV1R_highRatio_simulations.rds")

nsims <- 15L
freqWithin <- as.data.frame(Results$freqIntervals$within)
fidWithin <- as.data.frame(Results$fidIntervals$within)
stanWithin <- as.data.frame(Results$stanIntervals$within)
freqWithin$simulation <- fidWithin$simulation <- 
  stanWithin$simulation <- factor(1L:nsims)
freqWithin$inference <- "frequentist"
fidWithin$inference <- "fiducial"
freqAndFid <- rbind(freqWithin, fidWithin)

library(ggplot2)
ggplot(
  freqAndFid, 
  aes(
    x = simulation, y = estimate, ymin = lwr, ymax = upr, 
    group = simulation, color = inference
  )
) + geom_pointrange(position = position_dodge2(width = 0.5)) + 
  scale_discrete_manual("colour", values = c("red", "blue")) + 
  ylab("interval") + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  ggtitle("Confidence intervals about the within standard deviation", 
          subtitle = "Fiducial and frequentist")

ggplot(
  stanWithin, 
  aes(
    x = simulation, y = estimate, ymin = lwr, ymax = upr
  )
) + geom_pointrange() + 
  ylab("interval") + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  ggtitle("Confidence intervals about the within standard deviation", 
          subtitle = "STAN")
