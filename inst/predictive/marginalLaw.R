data <- data.frame(
  group = gl(2L, 4L),
  treatment = gl(2L, 2L, 8L)
)

data <- data.frame(
  group = c("1", "2"),
  treatment = c("1", "1")
)

#fixed <- ~ 1
random <- ~ group * treatment

library(lazyeval)
data <- droplevels(data)

tf <- terms.formula(random)
factors <- rownames(attr(tf, "factors"))
tvars <- attr(tf, "variables")
tlabs <- attr(tf, "term.labels")
tvars <- eval(tvars, envir = data)
rdat <- lapply(tvars, function(tvar) droplevels(as.factor(tvar)))
names(rdat) <- factors
RE <- lapply(tlabs, function(tlab){
  droplevels(lazy_eval(as.lazy(tlab), data = rdat))
})
# just to show:
setNames(as.data.frame(RE), tlabs)
#   group treatment group:treatment
# 1     1         1             1:1
# 2     1         2             1:2
# 3     1         3             1:3
# 4     2         1             2:1
# 5     2         2             2:2
# 6     2         3             2:3
n <- nrow(data)
RE2 <- c(RE, list(error = factor(1L:n))) # Adds the error effect 
E <- c(vapply(RE, nlevels, integer(1L)), n)
Z <- NULL 
for(i in seq_along(E)){ # Builds an indicator matrix for the effects
  re_levels <- levels(RE2[[i]])
  for(j in 1L:E[i]){
    temp1 <- which(RE2[[i]] == re_levels[j]) 
    temp2 <- integer(n) 
    temp2[temp1] <- 1L 
    Z <- cbind(Z, temp2)
  }
} 

variances <- c(1, 2, 3, 4)

Z %*% diag(rep(variances, times = E)) %*% t(Z)
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]   10    6    1    1    2    2    0    0
# [2,]    6   10    1    1    2    2    0    0
# [3,]    1    1   10    6    0    0    2    2
# [4,]    1    1    6   10    0    0    2    2
# [5,]    2    2    0    0   10    6    1    1
# [6,]    2    2    0    0    6   10    1    1
# [7,]    0    0    2    2    1    1   10    6
# [8,]    0    0    2    2    1    1    6   10
