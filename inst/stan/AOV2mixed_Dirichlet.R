library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

SimAV2mixed <- function(I, J, Kij, mu=0, alphai, sigmaO=1,
                        sigmaPO=1, sigmaE=1, factor.names=c("Part","Operator"),
                        resp.name="y", keep.intermediate=FALSE){
  if(length(Kij) == 1L){
    Kij <- rep(Kij, I*J)
  }
  Operator <- factor(
    rep(paste0("B", sprintf(paste0("%0", floor(log10(J)) + 1, "d"), 1:J)), each = I),
  )
  Part <- factor(
    rep(paste0("A", sprintf(paste0("%0", floor(log10(I)) + 1, "d"), 1:I)), times = J),
  )
  Oj <- rep(rnorm(J, 0, sigmaO), each=I)
  Pi <- rep(alphai, times=J)
  POij <- rnorm(I*J, 0, sigmaPO)
  simdata0 <- data.frame(Part, Operator, Pi, Oj, POij)
  simdata <- 
    as.data.frame(
      sapply(simdata0, function(v) rep(v, times=Kij), simplify=FALSE))
  Eijk <- rnorm(sum(Kij), 0, sigmaE)
  simdata <- cbind(simdata, Eijk)
  simdata[[resp.name]] <- mu + with(simdata, Oj+Pi+POij+Eijk)
  names(simdata)[1:2] <- factor.names
  if(!keep.intermediate) simdata <- simdata[,c(factor.names,resp.name)]
  simdata
}

set.seed(666)  
I = 2; J = 6; Kij = rpois(I*J, 1) + 3
alphai <- c(10, 20)
sigmaO <- 1
sigmaPO <- 0.5
sigmaE <- 2
dat <- SimAV2mixed(I, J, Kij, mu = 0, alphai = alphai, 
                   sigmaO = sigmaO, sigmaPO = sigmaPO, sigmaE = sigmaE)

### compile Stan model
stanmodel <- stan_model(
  file = "~/Work/R/gfilmm/inst/stan/AOV2mixed_Dirichlet.stan", 
  model_name = "AOV2mixed_Dirichlet")

### Stan data
standata <- list(
  y=dat$y, N=nrow(dat), I = I, J = J,
  PartID = as.integer(dat$Part),
  OperatorID = as.integer(dat$Operator),
  InteractionID = as.integer(dat$Part:dat$Operator),
  shape = 1e6, scale = 1/1e-7, concentration = 1
)

### Stan initial values
inits <- function(chain_id){
  u <- runif(1,0.3,0.7)
  list(
    PartA = alphai + chain_id,
    sigmaE = chain_id/2 + 1,
    tau = rgamma(1, shape = 10, scale = 1),
    theta = c(u, 1-u)
  )
}

### run Stan
stansamples <- sampling(stanmodel, data = standata, init=inits, 
                        iter = 10000, warmup = 5000, chains = 4)
#                        control=list(adapt_delta=0.999, max_treedepth=15))

### outputs
library(coda)
codasamples <- do.call(
  mcmc.list, 
  plyr::alply(
    rstan::extract(stansamples, permuted=FALSE, 
                   pars = c("PartA", "sigmaE", "sigmaO", "sigmaOP", "sigmaTotal")), 
    2, mcmc))
x <- summary(codasamples)
x$quantiles

# RSTANARM ####
library(rstanarm)
rstanarm <- stan_lmer(
  y ~  Part + (1|Operator) + (1|Operator:Part), data = dat,
  prior = normal(0, 100),
  prior_aux = cauchy(0, 5),
  prior_covariance = decov(concentration = 1, shape = 1e6, scale = 1e-7/sqrt(2)),
  iter = 5000
)
posterior_interval(rstanarm, prob = 95/100)
