library(rgr)

SimAOV2R <- function(I, J, Kij, mu = 0, sigmaP = 1, sigmaO = 1, sigmaPO = 1,
                     sigmaE = 1, factor.names = c("Part", "Operator"),
                     resp.name = "y", keep.intermediate = FALSE) {
  if (length(Kij) == 1L) {
    Kij <- rep(Kij, I * J)
  }
  Operator <- factor(
    rep(sprintf(paste0("B%0", floor(log10(J)) + 1, "d"), 1:J), each = I),
  )
  Part <- factor(
    rep(sprintf(paste0("A%0", floor(log10(I)) + 1, "d"), 1:I), times = J),
  )
  Oj <- rep(rnorm(J, 0, sigmaO), each = I)
  Pi <- rep(rnorm(I, 0, sigmaP), times = J)
  POij <- rnorm(I * J, 0, sigmaPO)
  simdata0 <- data.frame(Part, Operator, Pi, Oj, POij)
  simdata <- simdata0[rep(1:nrow(simdata0), times = Kij), ]
  Eijk <- rnorm(sum(Kij), 0, sigmaE)
  simdata <- cbind(simdata, Eijk)
  simdata[[resp.name]] <- mu + with(simdata, Oj + Pi + POij + Eijk)
  names(simdata)[1:2] <- factor.names
  if (!keep.intermediate) simdata <- simdata[, c(factor.names, resp.name)]
  simdata
}


confintAOV2R <- function(dat, alpha = 0.05) {
  dat <- gx.sort.df(~ Operator + Part, dat)
  I <- length(levels(dat$Part))
  J <- length(levels(dat$Operator))
  Kij <- aggregate(y ~ Part:Operator, FUN = length, data = dat)$y
  Khmean <- I * J / (sum(1 / Kij))
  ag <- aggregate(y ~ Part:Operator, FUN = mean, data = dat)
  Ybarij <- ag$y
  Ybari <- aggregate(y ~ Part, FUN = mean, data = ag)$y
  Ybarj <- aggregate(y ~ Operator, FUN = mean, data = ag)$y
  Ybar <- mean(Ybarij)
  S2.P <- J * Khmean * crossprod(Ybari - Ybar) / (I - 1) # I=p J=o
  S2.O <- I * Khmean * crossprod(Ybarj - Ybar) / (J - 1)
  Ybari <- rep(Ybari, times = J)
  Ybarj <- rep(Ybarj, each = I)
  S2.PO <- Khmean * crossprod(Ybarij - Ybari - Ybarj + Ybar) / (I - 1) / (J - 1) #
  S2.E <- crossprod(dat$y - rep(Ybarij, times = Kij)) / (sum(Kij) - I * J) # je crois qu'on peut aussi faire anova() pour S2.PO et S2.E
  USS <- c(S2.P, S2.O, S2.PO, S2.E)
  gammaP <- (S2.P - S2.PO) / J / Khmean
  gammaM <- (S2.O + (I - 1) * S2.PO + I * (Khmean - 1) * S2.E) / I / Khmean
  gammaT <- (I * S2.P + J * S2.O + (I * J - I - J) * S2.PO + I * J * (Khmean - 1) * S2.E) / I / J / Khmean
  G <- 1 - sapply(
    c(I - 1, J - 1, (I - 1) * (J - 1), sum(Kij) - I * J),
    function(d) {
      d / qchisq(1 - alpha / 2, d)
    }
  ) # n/qchisq(alpha,n) = qf(1-alpha,Inf,n)
  H <- sapply(
    c(I - 1, J - 1, (I - 1) * (J - 1), sum(Kij) - I * J),
    function(d) {
      d / qchisq(alpha / 2, d)
    }
  ) - 1
  F <- rep(NA, 4)
  F[c(1, 3)] <- qf(1 - alpha / 2, df1 = I - 1, df2 = c((I - 1) * (J - 1), J - 1))
  F[c(2, 4)] <- qf(alpha / 2, df1 = I - 1, df2 = c((I - 1) * (J - 1), J - 1))
  G13 <- ((F[1] - 1)^2 - G[1]^2 * F[1]^2 - H[3]^2) / F[1]
  H13 <- ((1 - F[2])^2 - H[1]^2 * F[2]^2 - G[3]^2) / F[2]
  VLP <- G[1]^2 * S2.P^2 + H[3]^2 * S2.PO^2 + G13 * S2.P * S2.PO
  VUP <- H[1]^2 * S2.P^2 + G[3]^2 * S2.PO^2 + H13 * S2.P * S2.PO
  wUSS2 <- (c(1, I - 1, I * (Khmean - 1)) * USS[2:4])^2
  VLM <- crossprod(G[2:4]^2, wUSS2)
  VUM <- crossprod(H[2:4]^2, wUSS2)
  wUSS2 <- (c(I, J, I * J - I - J, I * J * (Khmean - 1)) * USS)^2
  VLT <- crossprod(G^2, wUSS2)
  VUT <- crossprod(H^2, wUSS2)
  out <- data.frame(
    Parameter = c("gammaP", "gammaM", "gammaT"),
    Estimate = c(gammaP, gammaM, gammaT),
    low.bound = c(gammaP - sqrt(VLP) / J / Khmean, gammaM - sqrt(VLM) / I / Khmean, gammaT - sqrt(VLT) / I / J / Khmean),
    up.bound = c(gammaP + sqrt(VUP) / J / Khmean, gammaM + sqrt(VUM) / I / Khmean, gammaT + sqrt(VUT) / I / J / Khmean)
  )
  attributes(out) <- c(attributes(out), alpha = alpha)
  return(out)
}


set.seed(3141)
dat <- SimAOV2R(4, 3, 2)
confintAOV2R(dat)

library(rstanarm)
options(mc.cores = parallel::detectCores())
rsa <- stan_lmer(y ~ (1|Part) + (1|Operator) + (1|Part:Operator), data = dat, 
                 iter = 10000,
                 prior_aux = cauchy(0, 5),
                 prior_covariance = decov(shape = 1/15, scale = 15))
tail(posterior_interval(rsa, prob = 0.95))

library(gfilmm)
library(doParallel)
cl <- makePSOCKcluster(3L)
registerDoParallel(cl)
dat <- gx.sort.df(~ Part + Operator, dat)
gfs <- foreach(i = c(3L,4L,5L), .combine=list, .multicombine = TRUE, .export = "gfilmm") %dopar% 
  gfilmm(~ cbind(y-0.01, y+0.01), ~ 1, ~ Part+Operator, data = dat, N = 10000*i, long = FALSE)
stopCluster(cl)
lapply(gfs, gfiSummary)



