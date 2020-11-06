library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(666)  
I = 2; J = 6; Kij = rpois(I*J, 1) + 3#rep(3, I*J) #c(2,3,4,5,6,7)#rep(5,I*J)
alphai <- c(10, 20)
sigmaO <- 1
sigmaPO <- 0.5
sigmaE <- 2
dat <- SimAV2mixed(I, J, Kij, mu = 0, alphai = alphai, 
                   sigmaO = sigmaO, sigmaPO = sigmaPO, sigmaE = sigmaE)

### compile Stan model
stanmodel <- stan_model(
  file = "~/Work/R/gfilmm/inst/stan/AOV2mixed.stan", 
  model_name = "AOV2mixed")

### Stan data
standata <- list(
  y=dat$y, N=nrow(dat), I = I, J = J,
  PartIndicator = as.integer(dat$Part),
  OperatorIndicator = as.integer(dat$Operator),
  InteractionIndicator = as.integer(dat$Part:dat$Operator)
)

### Stan initial values
estimates <- function(dat, perturb=FALSE){
  if(perturb) dat$y <- dat$y + rnorm(length(dat$y), 0, 1)
  PartA <- aggregate(y ~ Part, data=dat, FUN=mean)$y
  sd(subset(aggregate(y ~ Part:Operator, data = dat, FUN=mean), Part == "A1")$y)
  sigma2O_plus_sigma2OP <- 
    mean(aggregate(y ~ Part, data = aggregate(y ~ Part:Operator, data = dat, FUN=mean), FUN = var)$y)
  sigmaO <- mean(aggregate(y ~ Part:Operator, data = dat, FUN=sd)$y * (I-1)/I)
  sigmaOP <- sqrt(sigma2O_plus_sigma2OP - sigmaO^2)
  sigmaE <- mean(aggregate(y ~ Part:Operator, data = dat, FUN=sd)$y) # non,c'est le même truc que sigmaO !
  return(list(PartA = PartA, sigmaE = sigmaE, sigmaO = sigmaO, sigmaOP = sigmaOP))
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
                   pars = c("PartA", "sigmaE", "sigmaO", "sigmaOP", "sigmaTotal")), 
    2, mcmc))
summary(codasamples)
