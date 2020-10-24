#' Generalized fiducial inference
#' @description Samples the fiducial distributions.
#'
#' @param y a formula of the form \code{~ cbind(lower,upper)} for the interval 
#'   data 
#' @param fixed formula for the fixed effects
#' @param random formula for the random effects
#' @param N number of simulations
#' @param thresh threshold, default \code{N/2}; for experts only
#'
#' @return A list with two components: \code{VERTEX} and \code{WEIGHT}.
#' 
#' @importFrom lazyeval f_eval_rhs as.lazy lazy_eval
#' @export
#' 
#' @references J. Cisewski and J.Hannig. 
#' \emph{Generalized fiducial inference for normal linear mixed models}. 
#' The Annals of Statistics 2012, Vol. 40, No. 4, 2102–2127.
#'
#' @examples
#' h <- 0.01
#' gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, N=500)
#' library(spatstat) # to use ewcdf 
#' f <- ewcdf(gfi$VERTEX[1,], gfi$WEIGHT)
#' plot(f, xlim = c(40,65))
#' quantile(f, c(0.025,0.975))
gfilmm <- function(y, fixed, random, data, N, thresh=N/2){
  data <- droplevels(data)
  Y <- f_eval_rhs(y, data = data)
  yl <- Y[,1L]; yu <- Y[,2L]
  FE <- model.matrix(fixed, data = data)
  tf <- terms.formula(random)
  factors <- rownames(attr(tf, "factors"))
  tvars <- attr(tf, "variables")
  tlabs <- attr(tf, "term.labels")
  rdat <- lapply(eval(tvars, envir = data), as.factor)
  names(rdat) <- factors
  # laz <- as.lazy(tlabs[1L], env = rdat) 
  # eval(laz$expr, envir = laz$env)
  # lazy_eval(as.lazy(tlabs[3L]), data = rdat)
  RE <- lapply(tlabs, function(tlab){
    droplevels(lazy_eval(as.lazy(tlab), data = rdat))
  })
  
  n <- nrow(Y)
  re <- length(RE)+1L 

  E <- c(vapply(RE, nlevels, integer(1L)), n)
  
  RE2 <- c(RE, list(error = factor(1L:n))) #Adds the error effect 
  RE <- NULL 
  for(i in seq_along(E)){ # Builds an indicator RE matrix for the effects
    re_levels <- levels(RE2[[i]])
    for(j in 1L:E[i]){
      temp1 <- which(RE2[[i]] == re_levels[j]) 
      temp2 <- integer(n) 
      temp2[temp1] <- 1L 
      RE <- cbind(RE, temp2)
    }
  } 
  RE2 <- vapply(RE2, as.integer, integer(n)) - 1L
  gfilmm_(yl, yu, FE, RE, RE2, E, N, thresh)
}

######################
# dat <- data.frame(
#   A = c("a", "b", "c"),
#   B = c("x", "y", "z"),
#   NotUsed = c(1, 2, 3)
# )
# 
# frml <- ~ A + B + A:B
# 
# # [[1]]
# # [1] a b c
# # Levels: a b c
# # 
# # [[2]]
# # [1] x y z
# # Levels: x y z
# # 
# # [[3]]
# # [1] a:x b:y c:z
# # Levels: a:x b:y c:z
# 
# library(lazyeval) # to use 'as.lazy' and 'lazy_eval'
# tf <- terms.formula(frml)
# factors <- rownames(attr(tf, "factors"))
# tvars <- attr(tf, "variables")
# tlabs <- attr(tf, "term.labels")
# used <- lapply(eval(tvars, envir = dat), as.factor)
# names(used) <- factors
# lapply(tlabs, function(tlab){
#   droplevels(lazy_eval(as.lazy(tlab), data = used))
# })
