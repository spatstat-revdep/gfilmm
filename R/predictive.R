#' Generalized fiducial predictive distributions
#' @description Simulations of the generalized fiducial predictive 
#'   distributions.
#'
#' @param gfi a \code{\link{gfilmm}} object
#' @param newdata dataframe in which to look for variables with which to predict
#'
#' @return A list with two fields: \code{FPD}, a dataframe containing the 
#'   simulations, and \code{WEIGHT}, their weight.
#' 
#' @importFrom stats model.matrix rnorm
#' @importFrom utils head tail
#' @export
#'
#' @examples gfi <- gfilmm(~ cbind(yield-0.1, yield+0.1), ~ N, ~ P, npk, 20)
#' gfiPredictive(gfi, data.frame(N = "0", P = "0"))
gfiPredictive <- function(gfi, newdata){
  if(anyDuplicated(newdata)){
    stop(
      "There are some duplicated rows in `newdata`."
    )
  }
  cvrts <- attr(gfi, "covariates")
  factors <- names(cvrts[["categorical"]])
  continuous <- names(cvrts[["continuous"]])
  nms   <- c(continuous, factors)
  if(any(!is.element(nms, names(newdata)))){
    stop(
      "`newdata` does not contain all the necessary variables."
    )
  }
  newdata <- droplevels(newdata)
  for(fact in factors){
    newdata[[fact]] <- as.factor(newdata[[fact]]) -> f
    levs <- cvrts[["categorical"]][[fact]]
    if(any(!is.element(levels(f), levs))){
      stop(
        "Found a factor in `newdata` with an unrecognized level."
      )
    }
    levels(newdata[[fact]]) <- levs
  }
  for(x in continuous){
    if(!is.numeric(newdata[[x]])){
      stop(
        sprintf("Variable `%s` should be numeric", x)
      )
    }
  }
  X <- model.matrix(attr(gfi, "fixed"), data = newdata) 
  RE2 <- getRE2(newdata, attr(gfi, "random"), check = FALSE)
  Z   <- getZ(RE2)
  N <- length(gfi[["WEIGHT"]])
  vertices <- t(as.matrix(gfi[["VERTEX"]]))
  neffects <- attr(gfi, "effects")
  fparams  <- head(vertices, neffects[["fixed"]])
  vars     <- tail(vertices, neffects[["random"]])
  E   <- vapply(RE2, nlevels, integer(1L))
  gauss <- matrix(rnorm(N*nrow(newdata)), nrow = N, ncol = nrow(newdata))
  out <- matrix(NA_real_, nrow = N, ncol = nrow(newdata))
  colnames(out) <- paste0("y", seq_len(nrow(newdata)))
  for(i in 1L:N){
    cholSigma <- chol(Z %*% (rep(vars[, i], times = E) * t(Z)))
    Mu <- X %*% fparams[, i]
#    out[i,] <- rmvn(1L, Mu, Sigma)
    # A <- matrix(NA_real_, nrow = 1L, ncol = neffects[["fixed"]])
    # .Call("rmvnCpp", n_ = 1L, mu_ = Mu, sigma_ = Sigma, ncores_ = 1L,
    #       isChol_ = FALSE, A_ = A, PACKAGE = "mvnfast")
    # out[i,] <- A
    out[i,] <- t(Mu) + t(gauss[i,]) %*% cholSigma # or gauss[i,,drop=FALSE] instead of t() => TODO: benchmarks
  }
  out <- list(FPD = as.data.frame(out), WEIGHT = gfi[["WEIGHT"]])
  class(out) <- "gfilmm.pred"
  out
}