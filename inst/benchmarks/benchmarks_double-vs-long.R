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
  group = gl(I, J)
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
sigmaWithin  <- 1e-3
I            <- 10L # number of groups
J            <- 5L # sample size per group

set.seed(31415926L)
groupmeans <- rnorm(I, mu, sigmaBetween)
y          <- c(
  vapply(groupmeans, function(gmean) rnorm(J, gmean, sigmaWithin), numeric(J))
)
y_rounded  <- round(y, digits = 4L)
dat        <- data.frame(
  ylwr = y_rounded - 5e-5,
  yupr = y_rounded + 5e-5,
  group = gl(I, J)
)


double = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 30000L)
gfiSummary(double)
long = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 6000L, long = TRUE)
gfiSummary(long)

library(AOV1R)
aovfit <- aov1r(y ~ dat$group)
confint(aovfit)
