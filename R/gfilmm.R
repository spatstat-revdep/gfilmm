#' Generalized fiducial inference
#' @description Samples the fiducial distributions.
#'
#' @param y a formula of the form \code{~ cbind(lower,upper)} for the interval 
#'   data 
#' @param fixed formula for the fixed effects
#' @param random formula for the random effects
#' @param data XXX
#' @param N number of simulations
#' @param thresh threshold, default \code{N/2}; for experts only
#'
#' @return A list with two components: \code{VERTEX} and \code{WEIGHT}.
#' 
#' @importFrom lazyeval f_eval_rhs as.lazy lazy_eval
#' @importFrom stats terms.formula model.matrix
#' @export
#' 
#' @references J. Cisewski and J.Hannig. 
#' \emph{Generalized fiducial inference for normal linear mixed models}. 
#' The Annals of Statistics 2012, Vol. 40, No. 4, 2102–2127.
#'
#' @examples
#' h <- 0.01
#' # gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N=500)
#' # library(spatstat) # to use ewcdf 
#' # f <- ewcdf(gfi$VERTEX[1,], gfi$WEIGHT)
#' # plot(f, xlim = c(40,65))
#' # quantile(f, c(0.025,0.975))
#' 
#' ##########################################################################
#' SimAV2mixed <- function(I, J, Kij, mu=0, alphai, sigmaO=1,
#'                         sigmaPO=1, sigmaE=1, factor.names=c("Part","Operator"),
#'                         resp.name="y", keep.intermediate=FALSE){
#'   Operator <- rep(1:J, each=I)
#'   Oj <- rep(rnorm(J, 0, sigmaO), each=I)
#'   Part <- rep(1:I, times=J)
#'   Pi <- rep(alphai, times=J)
#'   POij <- rnorm(I*J, 0, sigmaPO)
#'   simdata0 <- data.frame(Part, Operator, Pi, Oj, POij)
#'   simdata0$Operator <- factor(simdata0$Operator)
#'   levels(simdata0$Operator) <- 
#'     sprintf(paste0("%0", floor(log10(J))+1, "d"), 1:J)
#'   simdata0$Part <- factor(simdata0$Part)
#'   levels(simdata0$Part) <- sprintf(paste0("%0", floor(log10(I))+1, "d"), 1:I)
#'   simdata <- 
#'     as.data.frame(
#'       sapply(simdata0, function(v) rep(v, times=Kij), simplify=FALSE))
#'   Eijk <- rnorm(sum(Kij), 0, sigmaE)
#'   simdata <- cbind(simdata, Eijk)
#'   simdata[[resp.name]] <- mu + with(simdata, Oj+Pi+POij+Eijk)
#'   levels(simdata[,1]) <- paste0("A", levels(simdata[,1]))
#'   levels(simdata[,2]) <- paste0("B", levels(simdata[,2]))
#'   names(simdata)[1:2] <- factor.names
#'   if(!keep.intermediate) simdata <- simdata[,c(factor.names,resp.name)]
#'   simdata
#' }
#' 
#' 
#' set.seed(666)  
#' I = 2; J = 3; Kij = rep(5,I*J)
#' alphai <- c(3, -3)
#' sigmaO <- 2
#' sigmaPO <- 1
#' sigmaE <- sqrt(3)
#' dat <- SimAV2mixed(I, J, Kij, mu = 0, alphai = alphai, 
#'                    sigmaO = sigmaO, sigmaPO = sigmaPO, sigmaE = sigmaE)
#' 
#' 
#' library(gfilmm)
#' h <- 0.01
#' gfi <- gfilmm(~ cbind(y-h,y+h), ~ Part, ~ Operator, data = dat, 
#'               N = 100, thresh=0)
gfilmm <- function(y, fixed, random, data, N, thresh=N/2){
  data <- droplevels(data)
  Y <- f_eval_rhs(y, data = data)
  n <- nrow(Y)
  yl <- Y[,1L]; yu <- Y[,2L]
  FE <- model.matrix(fixed, data = data)
  if(!is.null(random)){
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
