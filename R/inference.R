inference <- function(gfi, v, alpha=0.05){ 
  out <- numeric(5L)
  names(out) <- c("mean", "median", "lwr", "upr", "Pr(=0)")
  vertices <- gfi$VERTEX[v,]
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
  seq_ <- 1L:nrow(gfi$VERTEX)
  names(seq_) <- rownames(gfi$VERTEX)
  out <- t(vapply(seq_, function(v) inference(gfi, v, 1-conf), numeric(5L)))
  attr(out, "confidence level") <- conf
  out
}


#' Fiducial confidence interval
#' @description Fiducial confidence interval.
#'
#' @param parameter a right-sided formula defining the paramter of interest, 
#'   like \code{~ sigma_error/`(Intercept)`}
#' @param gfi the output of \code{\link{gfilmm}}
#' @param conf confidence level
#'
#' @return The fiducial confidence interval of the parameter.
#' 
#' @importFrom lazyeval f_eval_rhs
#' @export
#'
#' @examples h <- 0.01
#' gfi <- gfilmm(~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N=5000)
#' gfiConfInt(~ (sigma_block + sigma_error)/`(Intercept)`, gfi)
gfiConfInt <- function(parameter, gfi, conf = 0.95){#, side = "two-sided"){
  #side <- match.arg(side, c("two-sided", "left", "right"))
  vertices <- as.data.frame(t(gfi$VERTEX))
  fsims <- f_eval_rhs(parameter, data = vertices)
  fcdf <- ewcdf(fsims, weights = gfi$WEIGHT)
  alpha <- 1 - conf
  quantile(fcdf, c(alpha/2, 1-alpha/2))
}