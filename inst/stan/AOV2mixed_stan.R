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
