#' Generalized fiducial inference
#' @description Samples the fiducial distributions.
#'
#' @param yl,yu data 
#' @param FE fixed effects
#' @param RE random effects
#' @param N number of simulations
#' @param thresh threshold, default \code{N/2}
#'
#' @return A list with two components: \code{VERTEX} and \code{WEIGHT}.
#' @export
#'
#' @examples
#' h <- 0.01
#' yl <- npk$yield-h
#' yu <- npk$yield+h
#' FE <- as.matrix(rep(1,nrow(npk))) # intercept
#' RE <- data.frame(block = npk$block)
#' gfi <- gfilmm(yl, yu, FE, RE, N=5000)
#' library(spatstat) # to use ewcdf 
#' f <- ewcdf(gfi$VERTEX[1,], gfi$WEIGHT)
#' plot(f, xlim = c(40,65))
#' quantile(f, c(0.025,0.975))
gfilmm <- function(yl, yu, FE, RE, N, thresh=N/2){
  n <- length(yl)
  re <- ncol(RE)+1L 
  E <- integer(re) # E[i] = number of levels of i-th random effect
  E[re] <- n 
  for(i in 1L:(re-1L)){
    E[i] <- nlevels(RE[,i]) 
  }
  RE2 <- cbind(RE, factor(1L:n)) #Adds the error effect 
  RE <- NULL 
  for(i in 1L:re){ # Builds an indicator RE matrix for the effects
    re_levels <- levels(RE2[,i])
    for(j in 1L:E[i]){
      temp1 <- which(RE2[,i]==re_levels[j]) 
      temp2 <- integer(n) 
      temp2[temp1] <- 1L 
      RE <- cbind(RE,temp2)
    }
  } 
  RE2 <- sapply(RE2, function(x) as.integer(x)-1L)
  gfilmm_(yl, yu, FE, RE, RE2, E, N, thresh)
}
