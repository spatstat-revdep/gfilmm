library(gfilmm)
library(ggplot2)

dat <- AOV1R::simAOV1R(6, 5, 0, 3, 4)

gfi <- gfilmm(~ cbind(y - 0.05, y + 0.05), ~ 1, ~ group, dat, 5000)

indcs <- sample.int(5000, replace = TRUE, prob = gfi$WEIGHT)
sims <- data.frame(
  between = gfi$VERTEX[["sigma_group"]][indcs],
  within = gfi$VERTEX[["sigma_error"]][indcs]
)

ggplot(sims, aes(log(between), log(within))) + 
  geom_point(alpha = 0.2) 
#  geom_abline(slope = -1, intercept = 25)

