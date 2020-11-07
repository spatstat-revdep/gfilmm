library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

SimAV2mixed <- function(I, J, Kij, mu=0, alphai, sigmaO=1,
                        sigmaPO=1, sigmaE=1, factor.names=c("Part","Operator"),
                        resp.name="y", keep.intermediate=FALSE){
  Operator <- rep(1:J, each=I)
  Oj <- rep(rnorm(J, 0, sigmaO), each=I)
  Part <- rep(1:I, times=J)
  Pi <- rep(alphai, times=J)
  POij <- rnorm(I*J, 0, sigmaPO)
  simdata0 <- data.frame(Part, Operator, Pi, Oj, POij)
  simdata0$Operator <- factor(simdata0$Operator)
  levels(simdata0$Operator) <- 
    sprintf(paste0("%0", floor(log10(J))+1, "d"), 1:J)
  simdata0$Part <- factor(simdata0$Part)
  levels(simdata0$Part) <- sprintf(paste0("%0", floor(log10(I))+1, "d"), 1:I)
  simdata <- 
    as.data.frame(
      sapply(simdata0, function(v) rep(v, times=Kij), simplify=FALSE))
  Eijk <- rnorm(sum(Kij), 0, sigmaE)
  simdata <- cbind(simdata, Eijk)
  simdata[[resp.name]] <- mu + with(simdata, Oj+Pi+POij+Eijk)
  levels(simdata[,1]) <- paste0("A", levels(simdata[,1]))
  levels(simdata[,2]) <- paste0("B", levels(simdata[,2]))
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
  shape = 10, scale = 1, concentration = 5
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
                        iter = 5000, warmup = 2500, chains = 4)
#                        control=list(adapt_delta=0.999, max_treedepth=15))

### outputs
library(coda)
codasamples <- do.call(
  mcmc.list, 
  plyr::alply(
    rstan::extract(stansamples, permuted=FALSE, 
                   pars = c("PartA", "sigmaE", "sigmaO", "sigmaOP", "sigmaTotal")), 
    2, mcmc))
summary(codasamples)


# RSTANARM ####
library(rstanarm)
rstanarm <- stan_lmer(
  y ~  0 + Part + (1|Operator) + (1|Operator:Part), data = dat,
  prior_aux = cauchy(0, 5),
  prior_covariance = decov(concentration = 5, shape = 10, scale = 1),
  iter = 5000
)
