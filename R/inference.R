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
#' @param gfi Output of \code{\link{gfimm}}
#' @param conf confidence level
#'
#' @return Summary statistics in a matrix.
#' @export
#'
#' @examples
#' data(KM41)
#' h <- 0.005
#' yl <- KM41$y-h
#' yu <- KM41$y+h
#' FE <- as.matrix(rep(1,nrow(KM41))) # intercept
#' RE <- data.frame(Batch = KM41$Batch)
#' gfi <- gfilmm(yl, yu, FE, RE, N=5000)
#' gfiSummary(gfi)
gfiSummary <- function(gfi, conf=0.95){
  t(vapply(1:nrow(gfi$VERTEX), 
           function(v) inference(gfi, v, 1-conf), numeric(4)))
}
