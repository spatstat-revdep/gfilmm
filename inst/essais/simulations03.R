library(gfilmm)
library(AOV1R)

nsims <- 100L
confint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)
confint_CV           <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_grandMean    <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaWithin  <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaBetween <- matrix(NA_real_, nrow = nsims, ncol = 2L)
fconfint_sigmaTotal   <- matrix(NA_real_, nrow = nsims, ncol = 2L)

mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
sigmaTotal   <- sqrt(sigmaBetween^2 + sigmaWithin^2)
CV           <- sigmaTotal / mu
n            <- 8L # sample size per group

set.seed(666L)
for(i in 1L:nsims){
  means <- rnorm(2L, mu, sigmaBetween)
  y1    <- rnorm(n, means[1L], sigmaWithin) 
  y2    <- rnorm(n, means[2L], sigmaWithin) 
  y     <- round(c(y1, y2), digits = 1)
  group <- gl(2L, n)
  dat   <- data.frame(
    ylwr = y - 0.05,
    yupr = y + 0.05,
    group = group
  )
  gfi <- gfilmm(~ cbind(ylwr, yupr), fixed = ~ 1, random = ~ group, 
                data = dat, N = 5000L)
  confint_grandMean[i,]    <- gfiConfInt(~ `(Intercept)`, gfi)
  confint_sigmaBetween[i,] <- gfiConfInt(~ sigma_group, gfi)
  confint_sigmaWithin[i,]  <- gfiConfInt(~ sigma_error, gfi)
  confint_sigmaTotal[i,]   <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2), gfi)
  confint_CV[i,]           <- 
    gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2)/`(Intercept)`, gfi)
  aovfit <- aov1r(c(y1,y2) ~ group)
  cis <- confint(aovfit)
  fconfint_grandMean[i,]    <- unlist(cis["grandmean", c("lwr", "upr")])
  fconfint_sigmaBetween[i,] <- unlist(cis["between", c("lwr", "upr")])
  fconfint_sigmaWithin[i,]  <- unlist(cis["within", c("lwr", "upr")])
  fconfint_sigmaTotal[i,]   <- unlist(cis["total", c("lwr", "upr")])
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
  fsigmaTotal = fconfint_sigmaTotal
)

saveRDS(results, "~/Work/R/gfilmm/inst/essais/simulations03.rds")

################################################################################
results <- readRDS("~/Work/R/gfilmm/inst/essais/simulations03.rds")

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

############################################################################
library(rAmCharts4)

dat <- data.frame(
  x = 1:nsims,
  y1 = results$grandMean[,1],
  y2 = results$grandMean[,2],
  z1 = results$fgrandMean[,1],
  z2 = results$fgrandMean[,2]
)

dat <- data.frame(
  x = 1:nsims,
  y1 = results$sigmaBetween[,1],
  y2 = results$sigmaBetween[,2],
  z1 = results$fsigmaBetween[,1],
  z2 = results$fsigmaBetween[,2]
)

amDumbbellChart(
  width = NULL,
  data = dat[1:20,],
  draggable = FALSE,
  category = "x",
  values = rbind(c("y1","y2"), c("z1","z2")),
  seriesNames = c("Fiducial", "Frequentist"),
  chartTitle = "Confidence intervals about the between standard deviation",
  yLimits = c(-10, 200),
  segmentsStyle = list(
    "Fiducial" = amSegment(width = 2),
    "Frequentist" = amSegment(width = 2)
  ),
  bullets = list(
    y1 = amTriangle(strokeWidth = 0),
    y2 = amTriangle(rotation = 180, strokeWidth = 0),
    z1 = amTriangle(strokeWidth = 0),
    z2 = amTriangle(rotation = 180, strokeWidth = 0)
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
