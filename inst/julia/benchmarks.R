library(XRJulia)
library(gfilmm)
library(microbenchmark)

juliaSource(normalizePath("inst/julia/forbenchmark.jl"))
x <- juliaEval("f(10)")
juliaGet(x)

dat <- data.frame(
  ylwr = c(2.0, 2.0, 3.0, 4.0, 4.0, 6.0),
  yupr = c(3.0, 3.0, 4.0, 5.0, 5.0, 7.0),
  group = gl(3,2)
)

microbenchmark(
  R = gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 30000),
  Julia = juliaEval("f(30000)"),
  times = 2
)



