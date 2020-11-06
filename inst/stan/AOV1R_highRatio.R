library(gfilmm)
library(AOV1R)
library(rstanarm)
options(mc.cores = parallel::detectCores())

mu           <- 10000 # grand mean
sigmaWithin  <- 1e-3
ratio        <- 50
( sigmaBetween <- sigmaWithin * ratio )
I            <- 10L # number of groups
J            <- 5L # sample size per group

set.seed(31415926L)
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

confint(aov1r(y ~ group, data = dat))

N <- 5000L
double = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = N)#, thresh = 5000)
gfiSummary(double)
table(attr(double, "ESS")/(N/2) < 1)

long = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 6000L, long = TRUE)
gfiSummary(long)


stan <- stan_lmer(
  y ~ (1|group), data = dat, 
  prior_aux = cauchy(0, 5),
  prior_covariance = decov(1, 1, 1, 1)
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
t(vapply(pstrr, quantile, numeric(2), probs = c(2.5,97.5)/100))

#####
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### compile Stan model
stanmodel <- stan_model(
  file = "~/Work/R/gfilmm/inst/stan/AOV1R.stan", 
  model_name = "AOV1R")

### Stan data
standata <- list(
  y=dat$y, N=nrow(dat), I = I,
  groupID = as.integer(dat$group)
)

### Stan initial values
estimates <- function(dat, perturb=FALSE){
  if(perturb) dat$y <- dat$y + rnorm(length(dat$y), 0, 1)
  x <- confint(aov1r(y ~ group, data = dat))
  return(list(mu = x["grandMean","estimate"], sigma_error = x["within","estimate"], sigma_group = x["between","estimate"]))
}
inits <- function(chain_id){
  values <- estimates(dat, perturb = chain_id > 1)
  return(values)
}

### run Stan
stansamples <- sampling(stanmodel, data = standata, init=inits, 
                        iter = 5000, warmup = 1000, chains = 4)
#                        control=list(adapt_delta=0.999, max_treedepth=15))

### outputs
library(coda)
codasamples <- do.call(
  mcmc.list, 
  plyr::alply(
    rstan::extract(stansamples, permuted=FALSE, 
                   pars = c("mu", "sigma_error", "sigma_group", "sigma_total")), 
    2, mcmc))
summary(codasamples)

