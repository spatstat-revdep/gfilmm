library(gfilmm)
library(rstanarm)
options(mc.cores = parallel::detectCores())

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
sigmaTotal <- sqrt(sigmaO^2 + sigmaPO^2 + sigmaE^2)

nsims <- 30
fidResults <- stanResults <- vector("list", nsims)

for(i in 1:nsims){
  dat <- SimAV2mixed(I, J, Kij, mu = 0, alphai = alphai, 
                     sigmaO = sigmaO, sigmaPO = sigmaPO, sigmaE = sigmaE)
  
  stanfit <- stan_lmer(
    y ~  0 + Part + (1|Operator) + (1|Operator:Part), data = dat,
    prior_aux = cauchy(0, 5),
    prior_covariance = decov(1, 1, 1, 100),
    iter = 2500
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
  pstrr$sigma_total <- 
    with(pstrr, sqrt(sigma_Operator + `sigma_Operator:Part` + sigma_error^2))
  pstrr[["sigma_Operator"]] <- sqrt(pstrr[["sigma_Operator"]])
  pstrr[["sigma_Operator:Part"]] <- sqrt(pstrr[["sigma_Operator:Part"]])
  
  stan <- t(vapply(pstrr, quantile, numeric(3), probs = c(50,2.5,97.5)/100))
  parms <- c(alphai, sigmaE, sigmaO, sigmaPO, sigmaTotal)
  stanResults[[i]] <- cbind(
    stan, 
    catch_both = stan[,2] < parms & parms < stan[,3],
    catch_left = stan[,3] > parms,
    catch_right = stan[,2] < parms
  )
  
  h <- 0.01
  gfi <- gfilmm(~ cbind(y-h,y+h), ~ 0 + Part, ~ Operator + Operator:Part, 
                data = dat, N = 5000)
  fid <- gfiSummary(gfi)
  fid <- rbind(
    fid, 
    sigma_total = c(
      NA, 
      spatstat::quantile.ewcdf(gfiCDF(~ sqrt(sigma_Operator^2 + `sigma_Operator:Part`^2 + sigma_error^2), gfi), probs = 0.5), 
      gfiConfInt(~ sqrt(sigma_Operator^2 + `sigma_Operator:Part`^2 + sigma_error^2), gfi), 
      0
    )
  )
  parms <- c(alphai, sigmaO, sigmaPO, sigmaE, sigmaTotal)
  fidResults[[i]] <- cbind(
    fid,
    catch_both = fid[,3] < parms & parms < fid[,4],
    catch_left = fid[,4] > parms,
    catch_right = fid[,3] < parms
  )
  cat("\n**********", i, "**********\n")
}

Results <- list(stanResults = stanResults, fidResults = fidResults)
saveRDS(Results, "~/Work/R/gfilmm/inst/stan/Results.rds")

stop()

################################################################################
Results <- readRDS("~/Work/R/gfilmm/inst/stan/Results.rds")
fidResults <- Results$fidResults
stanResults <- Results$stanResults

fidCoverage <- t(vapply(rownames(fidResults[[1]]), function(i){
  vapply(
    c("catch_both", "catch_left", "catch_right"),
    function(y){
      mean(vapply(fidResults, function(x) x[i,][y], numeric(1)))
    },
    numeric(1)
  )
}, numeric(3))) * 100


stanCoverage <- t(vapply(rownames(stanResults[[1]]), function(i){
  vapply(
    c("catch_both", "catch_left", "catch_right"),
    function(y){
      mean(vapply(stanResults, function(x) x[i,][y], numeric(1)))
    },
    numeric(1)
  )
}, numeric(3))) * 100

colnames(fidCoverage) <- colnames(stanCoverage) <- 
  c("two-sided", "left-sided", "right-sided")


fidMedians <- as.data.frame(vapply(rownames(fidResults[[1]]), function(i){
  vapply(fidResults, function(x) x[i,]["median"], numeric(1))
}, numeric(100))) 

stanMedians <- as.data.frame(vapply(rownames(stanResults[[1]]), function(i){
  vapply(stanResults, function(x) x[i,]["50%"], numeric(1))
}, numeric(100))) 

library(kde1d)
fmedian <- kde1d(fidMedians$sigma_Operator, xmin = 0, mult = 4)
plot(fmedian)
smedian <- kde1d(stanMedians$sigma_Operator, xmin = 0, mult = 4)
plot(smedian)
summary(fidMedians$sigma_Operator)
summary(stanMedians$sigma_Operator)

