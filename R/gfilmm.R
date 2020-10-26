#' Generalized fiducial inference
#' @description Samples the fiducial distributions.
#'
#' @param y a right-sided formula of the form \code{~ cbind(lower,upper)} for 
#'   the interval data 
#' @param fixed a right-sided formula for the fixed effects
#' @param random a right-sided formula for the random effects, or \code{NULL} 
#'   for no random effect
#' @param data the data, a dataframe
#' @param N desired number of simulations
#' @param thresh threshold, default \code{N/2}; for experts only
#'
#' @return A list with three components: \code{VERTEX}, \code{WEIGHT} and 
#'   \code{ESS}.
#' 
#' @importFrom lazyeval f_eval_rhs as.lazy lazy_eval
#' @importFrom stats terms.formula model.matrix
#' @export
#' 
#' @references J. Cisewski and J.Hannig. 
#'   \emph{Generalized fiducial inference for normal linear mixed models}. 
#'   The Annals of Statistics 2012, Vol. 40, No. 4, 2102–2127.
#'
#' @examples h <- 0.01
#' gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N=5000)
#' # fiducial cumulative distribution function of the intercept:
#' library(spatstat) # to use ewcdf 
#' f <- ewcdf(gfi$VERTEX["(Intercept)",], gfi$WEIGHT)
#' plot(f, xlim = c(40,65))
#' # fiducial confidence interval of the intercept:
#' quantile(f, c(0.025,0.975))
#' # fiducial density function of the intercept:
#' library(kde1d)
#' kfit <- kde1d(gfi$VERTEX["(Intercept)",], weights = gfi$WEIGHT)
#' curve(dkde1d(x, kfit), from = 45, to = 65)
gfilmm <- function(y, fixed, random, data, N, thresh=N/2){
  data <- droplevels(data)
  Y <- f_eval_rhs(y, data = data)
  if(!is.matrix(Y) || ncol(Y) != 2L){
    stop(
      "`y` must be a right-sided formula of the form `~ cbind(lower, upper)`."
    )
  }
  if(!is.numeric(Y)){
    stop(
      "Invalid `y` argument: the response must be given by two numerical ",
      "columns of the data."
    )
  }
  n <- nrow(Y)
  yl <- Y[,1L]; yu <- Y[,2L]
  if(any(yl >= yu)){
    stop(
      "Invalid interval data: found some values in the first column higher ",
      "than the corresponding values in the second column."
    )
  }
  FE <- model.matrix(fixed, data = data)
  if(!is.null(random)){
    tf <- terms.formula(random)
    factors <- rownames(attr(tf, "factors"))
    tvars <- attr(tf, "variables")
    tlabs <- attr(tf, "term.labels")
    tvars <- eval(tvars, envir = data)
    numerical <- vapply(tvars, is.numeric, logical(1L))
    if(any(numerical)){
      warning(
        "Numeric random effects are not supported; converting to factors."
      )
    }
    rdat <- lapply(tvars, as.factor)
    names(rdat) <- factors
    # laz <- as.lazy(tlabs[1L], env = rdat) 
    # eval(laz$expr, envir = laz$env)
    # lazy_eval(as.lazy(tlabs[3L]), data = rdat)
    RE <- lapply(tlabs, function(tlab){
      droplevels(lazy_eval(as.lazy(tlab), data = rdat))
    })
  }else{
    tlabs <- NULL
    RE <- NULL
  }
  RE2 <- c(RE, list(error = factor(1L:n))) #Adds the error effect 
  E <- c(vapply(RE, nlevels, integer(1L)), n)
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
  gfi <- gfilmm_(yl, yu, FE, RE, RE2, E, N, thresh)
  rownames(gfi$VERTEX) <-
    c(colnames(FE), paste0("sigma_", c(tlabs, "error")))
  gfi
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
