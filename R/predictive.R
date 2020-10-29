#' Generalized fiducial predictive distributions
#' @description Simulations of the generalized fiducial predictive 
#'   distributions.
#'
#' @param gfi a \code{\link{gfilmm}} object
#' @param newdata dataframe in which to look for variables with which to predict
#'
#' @return The simulations in a dataframe.
#' 
#' @importFrom mvnfast rmvn
#' @importFrom stats model.matrix
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
  X <- model.matrix(fixed, data = newdata)
  random <- attr(gfi, "random")
  Z <- getZ(getRE2(newdata, random, check = FALSE))
  N <- length(gfi[["WEIGHT"]])
  vertices <- t(as.matrix(gfi[["VERTEX"]]))
  neffects <- attr(gfi, "effects")
  fparams  <- head(vertices, neffects[["fixed"]])
  vars     <- tail(vertices, neffects[["random"]])
  E <- attr(gfi, "E")
  E[length(E)] <- nrow(newdata)
  out <- matrix(NA_real_, nrow = N, ncol = nrow(newdata))
  for(i in 1L:N){
    Sigma <- Z %*% diag(rep(vars[, i], times = E)) %*% t(Z)
    Mu <- X %*% fparams[, i]
    out[i,] <- rmvn(1L, Mu, Sigma)
    # A <- matrix(NA_real_, nrow = 1L, ncol = neffects[["fixed"]])
    # .Call("rmvnCpp", n_ = 1L, mu_ = Mu, sigma_ = Sigma, ncores_ = 1L, 
    #       isChol_ = FALSE, A_ = A, PACKAGE = "mvnfast")
    # out[i,] <- A
  }
  as.data.frame(out)
}