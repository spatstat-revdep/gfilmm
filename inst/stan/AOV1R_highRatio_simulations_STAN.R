library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# compile Stan model ####
stanmodel <- stan_model(
  file = "~/Work/R/gfilmm/inst/stan/AOV1R.stan", 
  model_name = "AOV1R")

# parameters ####
mu           <- 10000 # grand mean
sigmaWithin  <- 1
ratio        <- 50
sigmaBetween <- sigmaWithin * ratio
I            <- 10L # number of groups
J            <- 5L # sample size per group

# Stan initial values ####
estimates <- function(dat, perturb=FALSE){
  if(perturb) dat$y <- dat$y + rnorm(length(dat$y), 0, .1)
  x <- confint(AOV1R::aov1r(y ~ group, data = dat))
  return(list(mu = x["grandMean","estimate"], sigma_error = x["within","estimate"], sigma_group = x["between","estimate"]))
}
inits <- function(chain_id){
  values <- estimates(dat, perturb = chain_id > 1)
  return(values)
}

# simulations ####
nsims <- 15L
Intervals <- list(
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
  dat        <- data.frame(
    y = y,
    group = gl(I, J)
  )

  # Stan data ####
  standata <- list(
    y=dat$y, N=nrow(dat), I = I,
    groupID = as.integer(dat$group)
  )
  
  # run Stan ####
  stansamples <- sampling(stanmodel, data = standata, init=inits, 
                          iter = 10000, warmup = 5000, chains = 4)

  # outputs ###
  codasamples <- do.call(
    mcmc.list, 
    plyr::alply(
      rstan::extract(stansamples, permuted=FALSE, 
                     pars = c("mu", "sigma_error", "sigma_group", "sigma_total")), 
      2, mcmc))
  x <- summary(codasamples)$quantiles[, c(3L, 1L, 5L)]
  
  Intervals$within[i, ] <- x["sigma_error", ]
  Intervals$between[i, ] <- x["sigma_group", ]  

}

saveRDS(Intervals, "~/Work/R/gfilmm/inst/stan/AOV1R_highRatio_simulations_STAN.rds")

stop()

# ####
STAN <- readRDS("~/Work/R/gfilmm/inst/stan/AOV1R_highRatio_simulations_STAN.rds")
Results <- readRDS("~/Work/R/gfilmm/inst/stan/AOV1R_highRatio_simulations.rds")

nsims <- 15L
freqWithin <- as.data.frame(Results$freqIntervals$within)
STANWithin <- as.data.frame(STAN$within)
freqWithin$simulation <- STANWithin$simulation <- factor(1L:nsims)
freqWithin$inference <- "frequentist"
STANWithin$inference <- "Bayesian"
freqAndSTAN <- rbind(freqWithin, STANWithin)

library(ggplot2)
ggplot(
  freqAndSTAN, 
  aes(
    x = simulation, y = estimate, ymin = lwr, ymax = upr, 
    group = simulation, color = inference
  )
) + geom_pointrange(position = position_dodge2(width = 0.5)) + 
  scale_discrete_manual("colour", values = c("red", "blue")) + 
  ylab("interval") + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  ggtitle("Confidence intervals about the within standard deviation", 
          subtitle = "Frequentist and Bayesian with Cauchy priors")

# between ####
freqBetween <- as.data.frame(Results$freqIntervals$between)
STANBetween <- as.data.frame(STAN$between)
freqBetween$simulation <- STAN$simulation <- factor(1L:nsims)
freqBetween$inference <- "frequentist"
STANBetween$inference <- "Bayesian"
freqAndSTAN <- rbind(freqBetween, STANBetween)

library(ggplot2)
ggplot(
  freqAndSTAN, 
  aes(
    x = simulation, y = estimate, ymin = lwr, ymax = upr, 
    group = simulation, color = inference
  )
) + geom_pointrange(position = position_dodge2(width = 0.5)) + 
  scale_discrete_manual("colour", values = c("red", "blue")) + 
  ylab("interval") + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  ggtitle("Confidence intervals about the between standard deviation", 
          subtitle = "Frequentist and Bayesian with Cauchy priors")

