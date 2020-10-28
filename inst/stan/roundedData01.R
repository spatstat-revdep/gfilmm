y0 <- c(rep(10,10),rep(9,7),rep(11,8),8,8,8,12,12)
y <- sample(y0, 30)
hist(y, breaks = 7:12+0.5)

# 6 op, each 3 measurements
set.seed(3141592L)
mu <- 10
opmeans <- rnorm(6L, mean = mu, sd = 1)
y <- c(vapply(opmeans, function(m) rnorm(3L, mean = m, sd = 1), numeric(3L)))
y_rounded <- round(y, digits = 0L)
dat <- data.frame(
  ylwr = y_rounded - 0.05,
  yupr = y_rounded + 0.05,
  operator = gl(6L, 3L)
)

library(gfilmm)
gfi <- gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ operator, dat, N = 10000L)
gfiSummary(gfi)
gfiConfInt(~ sigma_operator^2 + sigma_error^2, gfi)




