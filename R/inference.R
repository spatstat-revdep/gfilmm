inference <- function(gfi, v, alpha=0.05){ 
  out <- numeric(4L)
  names(out) <- c("mean","median","low","up")
  vertex <- gfi$VERTEX[v,]
  weight <- gfi$WEIGHT
  out[1] <- sum(vertex*weight) # mean
  h <- cbind(vertex,weight)
  hsort <- h[order(h[,1L]),] # gx.sort(h,1L)
  hsum <- cumsum(hsort[,2L])
  ci_u <- min(which(hsum >= 1-alpha/2)) 
  ci_l <- min(which(hsum >= alpha/2))   
  ci_m <- min(which(hsum >= 0.5))
  out[3] <- hsort[ci_l,1L] # lower bound
  out[4] <- hsort[ci_u,1L] # upper bound
  out[2] <- hsort[ci_m,1L] # estimate (median)
  out
}

#' Summary of fiducial distributions
#'
#' @param gfi output of \code{\link{gfilmm}}
#' @param conf confidence level
#'
#' @return Summary statistics in a matrix.
#' @export
#'
#' @examples
#' data(KM41)
#' h <- 0.005
#' gfi <- gfilmm(~ cbind(y-h, y+h), ~ 1, ~ Batch, N = 500)
#' gfiSummary(gfi)
gfiSummary <- function(gfi, conf = 0.95){
  t(vapply(1L:nrow(gfi$VERTEX), 
           function(v) inference(gfi, v, 1-conf), numeric(4L)))
}
