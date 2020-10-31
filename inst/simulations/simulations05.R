library(gfilmm)
library(AOV1R)

pconfint <- function(fit){
  piv <- rGPQ(fit, 50000L)
  piv$GPQ_sigmatot <- with(piv, sqrt(GPQ_sigma2b + GPQ_sigma2w))
  cis <- lapply(piv, quantile, probs = c(2.5, 97.5)/100)
  cis$GPQ_sigma2b <- sign(cis$GPQ_sigma2b) * sqrt(abs(cis$GPQ_sigma2b))
  cis$GPQ_sigma2w <- sqrt(cis$GPQ_sigma2w)
  predi <- with(
    piv,
    quantile(rnorm(50000L, GPQ_mu, GPQ_sigmatot), probs = c(0.025, 0.975))
  )
  list(conf = cis, pred = predi)
}

nsims <- 20L
confint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_CV           <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
pconfint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
pconfint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
pconfint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
pconfint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
predint  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fpredint <- matrix(NA_real_, nrow = nsims, ncol = 2L)
gpredint <- matrix(NA_real_, nrow = nsims, ncol = 2L)


mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
sigmaTotal   <- sqrt(sigmaBetween^2 + sigmaWithin^2)
CV           <- sigmaTotal / mu
I            <- 3L # number of groups
J            <- 2L # sample size per group

set.seed(666L)
for(i in 1L:nsims){
  means <- rnorm(I, mu, sigmaBetween)
  y <- c(vapply(means, function(m) rnorm(J, m, sigmaWithin), numeric(J)))
  y_rounded <- round(y, digits = 1L)
  group <- gl(I, J)
  dat   <- data.frame(
    ylwr = y_rounded - 0.05,
    yupr = y_rounded + 0.05,
    group = group
  )
  gfi <- gfilmm(~ cbind(ylwr, yupr), fixed = ~ 1, random = ~ group, 
                data = dat, N = 3000L)
  confint_grandMean[i,]    <- gfiConfInt(~ `(Intercept)`, gfi)
  confint_sigmaBetween[i,] <- gfiConfInt(~ sigma_group, gfi)
  confint_sigmaWithin[i,]  <- gfiConfInt(~ sigma_error, gfi)
  confint_sigmaTotal[i,]   <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2), gfi)
  confint_CV[i,]           <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2)/`(Intercept)`, gfi)
  gfipred <- gfilmmPredictive(gfi, as.data.frame(matrix(nrow=1L,ncol=0L)))
  predint[i,] <- gfiConfInt(~ y1, gfipred)
  # frequentist
  aovfit <- aov1r(y ~ group)
  cis <- confint(aovfit)
  fconfint_grandMean[i,]    <- unlist(cis["grandMean", c("lwr", "upr")])
  fconfint_sigmaBetween[i,] <- unlist(cis["between", c("lwr", "upr")])
  fconfint_sigmaWithin[i,]  <- unlist(cis["within", c("lwr", "upr")])
  fconfint_sigmaTotal[i,]   <- unlist(cis["total", c("lwr", "upr")])
  fpredint[i,] <- predict(aovfit)
  cis <- pconfint(aovfit)
  pconfint_grandMean[i,]    <- cis[["cis"]]$GPQ_mu
  pconfint_sigmaBetween[i,] <- cis[["cis"]]$GPQ_sigma2b
  pconfint_sigmaWithin[i,]  <- cis[["cis"]]$GPQ_sigma2w
  pconfint_sigmaTotal[i,]   <- cis[["cis"]]$GPQ_sigmatot
  gpredint[i,] <- cis[["pred"]]
}


results <- list(
  grandMean = confint_grandMean,
  sigmaWithin = confint_sigmaWithin,
  sigmaBetween = confint_sigmaBetween,
  sigmaTotal = confint_sigmaTotal,
  CV = confint_CV,
  fgrandMean = fconfint_grandMean,
  fsigmaWithin = fconfint_sigmaWithin,
  fsigmaBetween = fconfint_sigmaBetween,
  fsigmaTotal = fconfint_sigmaTotal,
  pgrandMean = pconfint_grandMean,
  psigmaWithin = pconfint_sigmaWithin,
  psigmaBetween = pconfint_sigmaBetween,
  psigmaTotal = pconfint_sigmaTotal
)

saveRDS(results, "~/Work/R/gfilmm/inst/simulations/simulations04.rds")

################################################################################
results <- readRDS("~/Work/R/gfilmm/inst/simulations/simulations04.rds")

cvrg_twoSided <- c(
  sum(results$grandMean[,1] <= mu & mu <= results$grandMean[,2]),
  sum(results$sigmaWithin[,1] <= sigmaWithin & sigmaWithin <= results$sigmaWithin[,2]),
  sum(results$sigmaBetween[,1] <= sigmaBetween & sigmaBetween <= results$sigmaBetween[,2]),
  sum(results$sigmaTotal[,1] <= sigmaTotal & sigmaTotal <= results$sigmaTotal[,2]),
  sum(results$CV[,1] <= CV & CV <= results$CV[,2]) 
) / nsims

cvrg_leftSided <- c(
  sum(mu <= results$grandMean[,2]),
  sum(sigmaWithin <= results$sigmaWithin[,2]),
  sum(sigmaBetween <= results$sigmaBetween[,2]),
  sum(sigmaTotal <= results$sigmaTotal[,2]),
  sum(CV <= results$CV[,2]) 
) / nsims

cvrg_rightSided <- c(
  sum(results$grandMean[,1] <= mu),
  sum(results$sigmaWithin[,1] <= sigmaWithin),
  sum(results$sigmaBetween[,1] <= sigmaBetween),
  sum(results$sigmaTotal[,1] <= sigmaTotal),
  sum(results$CV[,1] <= CV) 
) / nsims

fcvrg_twoSided <- c(
  sum(results$fgrandMean[,1] <= mu & mu <= results$fgrandMean[,2]),
  sum(results$fsigmaWithin[,1] <= sigmaWithin & sigmaWithin <= results$fsigmaWithin[,2]),
  sum(results$fsigmaBetween[,1] <= sigmaBetween & sigmaBetween <= results$fsigmaBetween[,2]),
  sum(results$fsigmaTotal[,1] <= sigmaTotal & sigmaTotal <= results$fsigmaTotal[,2])
) / nsims

fcvrg_leftSided <- c(
  sum(mu <= results$fgrandMean[,2]),
  sum(sigmaWithin <= results$fsigmaWithin[,2]),
  sum(sigmaBetween <= results$fsigmaBetween[,2]),
  sum(sigmaTotal <= results$fsigmaTotal[,2])
) / nsims

fcvrg_rightSided <- c(
  sum(results$fgrandMean[,1] <= mu),
  sum(results$fsigmaWithin[,1] <= sigmaWithin),
  sum(results$fsigmaBetween[,1] <= sigmaBetween),
  sum(results$fsigmaTotal[,1] <= sigmaTotal)
) / nsims


pcvrg_twoSided <- c(
  sum(results$pgrandMean[,1] <= mu & mu <= results$pgrandMean[,2]),
  sum(results$psigmaWithin[,1] <= sigmaWithin & sigmaWithin <= results$psigmaWithin[,2]),
  sum(results$psigmaBetween[,1] <= sigmaBetween & sigmaBetween <= results$psigmaBetween[,2]),
  sum(results$psigmaTotal[,1] <= sigmaTotal & sigmaTotal <= results$psigmaTotal[,2])
) / nsims

pcvrg_leftSided <- c(
  sum(mu <= results$pgrandMean[,2]),
  sum(sigmaWithin <= results$psigmaWithin[,2]),
  sum(sigmaBetween <= results$psigmaBetween[,2]),
  sum(sigmaTotal <= results$psigmaTotal[,2])
) / nsims

pcvrg_rightSided <- c(
  sum(results$pgrandMean[,1] <= mu),
  sum(results$psigmaWithin[,1] <= sigmaWithin),
  sum(results$psigmaBetween[,1] <= sigmaBetween),
  sum(results$psigmaTotal[,1] <= sigmaTotal)
) / nsims

############################################################################
library(rAmCharts4)

dat <- data.frame(
  x = 1:nsims,
  x1 = results$grandMean[,1],
  x2 = results$grandMean[,2],
  y1 = results$fgrandMean[,1],
  y2 = results$fgrandMean[,2],
  z1 = results$pgrandMean[,1],
  z2 = results$pgrandMean[,2]
)

dat <- data.frame(
  x = 1:nsims,
  x1 = results$sigmaWithin[,1],
  x2 = results$sigmaWithin[,2],
  y1 = results$fsigmaWithin[,1],
  y2 = results$fsigmaWithin[,2],
  z1 = results$psigmaWithin[,1],
  z2 = results$psigmaWithin[,2]
)

dat <- data.frame(
  x = 1:nsims,
  x1 = results$sigmaBetween[,1],
  x2 = results$sigmaBetween[,2],
  y1 = results$fsigmaBetween[,1],
  y2 = results$fsigmaBetween[,2],
  z1 = results$psigmaBetween[,1],
  z2 = results$psigmaBetween[,2]
)

dat <- data.frame(
  x = 1:nsims,
  x1 = results$sigmaTotal[,1],
  x2 = results$sigmaTotal[,2],
  y1 = results$fsigmaTotal[,1],
  y2 = results$fsigmaTotal[,2],
  z1 = results$psigmaTotal[,1],
  z2 = results$psigmaTotal[,2]
)

amDumbbellChart(
  width = NULL,
  data = dat[1:20,],
  draggable = FALSE,
  category = "x",
  values = rbind(c("x1","x2"), c("y1","y2"), c("z1","z2")),
  seriesNames = c("Fiducial", "Frequentist", "Pivotal"),
  chartTitle = "Confidence intervals about the total standard deviation",
#  yLimits = c(-10, 200),
  segmentsStyle = list(
    "Fiducial" = amSegment(width = 2, color = "red"),
    "Frequentist" = amSegment(width = 2, color = "green"),
    "Pivotal" = amSegment(width = 2, color = "blue")
  ),
  bullets = list(
    x1 = amTriangle(strokeWidth = 0, color = "red"),
    x2 = amTriangle(rotation = 180, strokeWidth = 0, color = "red"),
    y1 = amTriangle(strokeWidth = 0, color = "green"),
    y2 = amTriangle(rotation = 180, strokeWidth = 0, color = "green"),
    z1 = amTriangle(strokeWidth = 0, color = "blue"),
    z2 = amTriangle(rotation = 180, strokeWidth = 0, color = "blue")
  ),
  tooltip = amTooltip("upper: {openValueY}\nlower: {valueY}", scale = 0.75),
  xAxis = list(
    title = amText(
      "simulation",
      fontSize = 17, fontWeight = "bold", fontFamily = "Helvetica"
    ),
    labels = amAxisLabels(fontSize = 12)
  ),
  yAxis = list(
    title = amText(
      "interval",
      fontSize = 17, fontWeight = "bold", fontFamily = "Helvetica"
    ),
    gridLines = amLine("silver", width = 1, opacity = 0.4)
  ),
  legend = amLegend(position = "right", itemsWidth = 15, itemsHeight = 15),
  backgroundColor = "lightyellow",
  theme = "dataviz",
  export = TRUE
)
