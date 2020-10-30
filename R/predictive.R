#' Generalized fiducial predictive distributions
#' @description Simulations of the generalized fiducial predictive 
#'   distributions.
#'
#' @param gfi a \code{\link{gfilmm}} object
#' @param newdata dataframe in which to look for variables with which to predict
#'
#' @return The simulations in a dataframe.
#' 
#' @importFrom stats model.matrix rnorm
#' @importFrom utils head tail
#' @export
#'
#' @examples gfi <- gfilmm(~ cbind(yield-0.1, yield+0.1), ~ N, ~ P, npk, 20)
#' gfiPredictive(gfi, data.frame(N = "0", P = "0"))
gfiPredictive <- function(gfi, newdata){
  cvrts <- attr(gfi, "covariates")
  nms   <- c(names(cvrts[["continuous"]]), names(cvrts[["categorical"]]))
  if(any(!is.element(nms, names(newdata)))){
    stop(
      "`newdata` does not contain all the necessary variables."
    )
  }
  # TODO: check levels
  fixed  <- attr(gfi, "fixed")
  X <- model.matrix(fixed, data = newdata) # pb si nlevels=1 => ajouter les niveaux
  random <- attr(gfi, "random")
  Z <- getZ(getRE2(newdata, random, check = FALSE))
  N <- length(gfi[["WEIGHT"]])
  vertices <- t(as.matrix(gfi[["VERTEX"]]))
  neffects <- attr(gfi, "effects")
  fe       <- neffects[["fixed"]]
  fparams  <- head(vertices, fe)
  vars     <- tail(vertices, neffects[["random"]])
  E <- attr(gfi, "E")
  E[length(E)] <- nrow(newdata)
  gauss <- matrix(rnorm(N*fe), nrow = N, ncol = fe)
  out <- matrix(NA_real_, nrow = N, ncol = nrow(newdata))
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
  as.data.frame(out)
}