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
I = 2; J = 6; Kij = rep(5,I*J) #c(2,3,4,5,6,7)#rep(5,I*J)
alphai <- c(3, -3)
sigmaO <- 2
sigmaPO <- 1
sigmaE <- sqrt(3)
dat <- SimAV2mixed(I, J, Kij, mu = 0, alphai = alphai, 
                   sigmaO = sigmaO, sigmaPO = sigmaPO, sigmaE = sigmaE)

library(lme4)
lmer(y ~ 0 + Part + (1|Operator) + (1|Operator:Part), data = dat)

library(rstanarm)
options(mc.cores = parallel::detectCores())
stanfit <- stan_lmer(
  y ~ 0 + Part + (1|Operator) + (1|Operator:Part), data = dat,
  prior_aux = cauchy(0, 5)
)
pstrr <- as.data.frame(
  stanfit, 
  pars = c(
    "PartA1", 
    "PartA2", 
    "sigma", 
    "Sigma[Operator:(Intercept),(Intercept)]",
    "Sigma[Operator:Part:(Intercept),(Intercept)]"
  )
)
names(pstrr) <- 
  c("PartA1", "PartA2", "sigma_error", "sigma_Operator", "sigma_Operator:Part")
pstrr[["sigma_Operator"]] <- sqrt(pstrr[["sigma_Operator"]])
pstrr[["sigma_Operator:Part"]] <- sqrt(pstrr[["sigma_Operator:Part"]])


library(gfilmm)
h <- 0.01
gfi <- gfilmm(~ cbind(y-h,y+h), ~ Part, ~ Operator, data = dat, 
              N = 100, thresh=0)
