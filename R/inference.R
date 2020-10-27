inference <- function(gfi, v, alpha=0.05){ 
  out <- numeric(5L)
  names(out) <- c("mean", "median", "lwr", "upr", "Pr(=0)")
  vertices <- gfi$VERTEX[[v]]
  weights <- gfi$WEIGHT
  out[1L] <- sum(vertices*weights) # mean
  h <- cbind(vertices,weights)
  hsort <- h[order(h[,1L]),] # gx.sort(h,1L)
  hsum <- cumsum(hsort[,2L])
  ci_u <- min(which(hsum >= 1-alpha/2)) 
  ci_l <- min(which(hsum >= alpha/2))   
  ci_m <- min(which(hsum >= 0.5))
  out[3L] <- hsort[ci_l,1L] # lower bound
  out[4L] <- hsort[ci_u,1L] # upper bound
  out[2L] <- hsort[ci_m,1L] # estimate (median)
  fe <- attr(gfi, "effects")[["fixed"]]
  if(v <= fe){
    out[5L] <- NA_real_
  }else{
    zeros <- vertices == 0
    if(any(zeros)){
      out[5L] <- sum(weights[zeros])
    }else{
      out[5L] <- 0
    }
  }
  out
}

#' Summary of fiducial distributions
#'
#' @param gfi output of \code{\link{gfilmm}}
#' @param conf confidence level
#'
#' @return A matrix with summary statistics: means, medians, confidence 
#'   intervals, and probabilities that the standard deviations equal 0.
#' @export
#'
#' @examples
#' data(KM41)
#' h <- 0.005
#' gfi <- gfilmm(~ cbind(y-h, y+h), ~ 1, ~ Batch, data = KM41, N = 5000)
#' gfiSummary(gfi)
gfiSummary <- function(gfi, conf = 0.95){
  seq_ <- 1L:ncol(gfi$VERTEX)
  names(seq_) <- names(gfi$VERTEX)
  out <- t(vapply(seq_, function(v) inference(gfi, v, 1-conf), numeric(5L)))
  attr(out, "confidence level") <- conf
  out
}


#' Fiducial confidence interval
#' @description Fiducial confidence interval of a parameter of interest.
#'
#' @param parameter a right-sided formula defining the parameter of interest, 
#'   like \code{~ sigma_error/`(Intercept)`}
#' @param gfi the output of \code{\link{gfilmm}}
#' @param conf confidence level
#'
#' @return The fiducial confidence interval of the parameter.
#' 
#' @importFrom lazyeval f_eval_rhs
#' @importFrom spatstat ewcdf quantile.ewcdf
#' @export
#'
#' @examples h <- 0.01
#' gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N=5000)
#' gfiConfInt(~ sqrt(sigma_block^2 + sigma_error^2)/`(Intercept)`, gfi)
gfiConfInt <- function(parameter, gfi, conf = 0.95){#, side = "two-sided"){
  #side <- match.arg(side, c("two-sided", "left", "right"))
  fsims <- f_eval_rhs(parameter, data = gfi$VERTEX)
  fcdf <- ewcdf(fsims, weights = gfi$WEIGHT)
  alpha <- 1 - conf
  quantile.ewcdf(fcdf, c(alpha/2, 1-alpha/2))
}

#' Fiducial cumulative distribution function
#' @description Fiducial cumulative distribution function of a parameter of 
#'   interest.
#'
#' @param parameter a right-sided formula defining the parameter of interest, 
#'   like \code{~ sigma_error/`(Intercept)`}
#' @param gfi the output of \code{\link{gfilmm}}
#'
#' @return The fiducial cumulative distribution function of the parameter.
#' 
#' @importFrom lazyeval f_eval_rhs
#' @importFrom spatstat ewcdf 
#' @export
#'
#' @examples h <- 0.01
#' gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N=5000)
#' F <- gfiCDF(~ sqrt(sigma_block^2 + sigma_error^2)/`(Intercept)`, gfi)
#' plot(F, xlim = c(0, 0.3), main = "Coefficient of variation", 
#'      ylab = expression("Pr("<="x)"))
#' F(0.2)
gfiCDF <- function(parameter, gfi){
  fsims <- f_eval_rhs(parameter, data = gfi$VERTEX)
  ewcdf(fsims, weights = gfi$WEIGHT) 
}
