library(gfilmm)
library(ggplot2)

dat <- AOV1R::simAOV1R(6, 5, 0, 3, 4)

gfi <- gfilmm(~ cbind(y - 0.05, y + 0.05), ~ 1, ~ group, dat, 5000)

indcs <- sample.int(5000, replace = TRUE, prob = gfi$WEIGHT)
sims <- data.frame(
  between = gfi$VERTEX[["sigma_group"]][indcs]^2,
  within = gfi$VERTEX[["sigma_error"]][indcs]^2
)

ggplot(sims, aes(between, within)) + 
  geom_point(alpha = 0.2) +
  geom_abline(slope = -1, intercept = 25)

################################################################################
library(GGally)

ggpairs(
  gfi$VERTEX[indcs,],
  #mapping = ggplot2::aes(),
  upper = NULL,#list(continuous = ggally_density),
  lower = list(continuous = wrap("points", alpha = 0.2)),
  diag = NULL,
  title = "Diamonds"
)

#############################################################
set.seed(666L)
n <- 30L
x <- 1L:n
y <- rnorm(n, mean = x, sd = 2)
y_rounded <- round(y, digits = 1L)
dat <- data.frame(
  ylwr = y_rounded - 0.05,
  yupr = y_rounded + 0.05,
  x = x
)

library(gfilmm)
fidSims <- gfilmm(
  y = ~ cbind(ylwr, yupr), # interval data
  fixed = ~ x,             # fixed effects
  random = NULL,           # random effects
  data = dat,              # data
  N = 20000L               # number of simulations
)

indcs <- sample.int(20000L, replace = TRUE, prob = fidSims$WEIGHT)

ggpairs(
  fidSims$VERTEX[indcs,],
  #mapping = ggplot2::aes(),
  upper = list(continuous = ggally_density),
  lower = list(continuous = wrap("points", alpha = 0.1))
  #diag = NULL,
)

library(car)
lmfit <- lm(y ~ x)

plot(x ~ `(Intercept)`, data = fidSims$VERTEX[indcs,])
confidenceEllipse(lmfit, 1:2, add = TRUE)
