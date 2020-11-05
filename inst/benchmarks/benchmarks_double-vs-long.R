library(gfilmm)
library(microbenchmark)

mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
I            <- 6L # number of groups
J            <- 5L # sample size per group

set.seed(31415926L)
groupmeans <- rnorm(I, mu, sigmaBetween)
y          <- c(
  vapply(groupmeans, function(gmean) rnorm(J, gmean, sigmaWithin), numeric(J))
)
y_rounded  <- round(y, digits = 1L)
dat        <- data.frame(
  ylwr = y_rounded - 0.05,
  yupr = y_rounded + 0.05,
  group = gl(J, I)
)

microbenchmark(
  double = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 20000L),
  long = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 20000L, long = TRUE),
  times = 3
)

################################################################################
library(gfilmm)

mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
I            <- 20L # number of groups
J            <- 15L # sample size per group

set.seed(31415926L)
groupmeans <- rnorm(I, mu, sigmaBetween)
y          <- c(
  vapply(groupmeans, function(gmean) rnorm(J, gmean, sigmaWithin), numeric(J))
)
y_rounded  <- round(y, digits = 1L)
dat        <- data.frame(
  ylwr = y_rounded - 0.05,
  yupr = y_rounded + 0.05,
  group = gl(J, I)
)


double = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 200L)
gfiSummary(double)
long = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 200L, long = TRUE)
gfiSummary(long)

library(AOV1R)
aovfit <- aov1r(y ~ dat$group)
confint(aovfit)
