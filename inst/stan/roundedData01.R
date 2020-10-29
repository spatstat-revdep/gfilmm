y0 <- c(rep(10,10),rep(9,7),rep(11,8),8,8,8,12,12)
y <- sample(y0, 30)
hist(y, breaks = 7:12+0.5)

# 6 op, each 3 measurements
set.seed(3141592L)
mu <- 10
# between 2 within 0.5
opmeans <- rnorm(4L, mean = mu, sd = 1)
y <- c(vapply(opmeans, function(m) rnorm(3L, mean = m, sd = 1), numeric(3L)))
y_rounded <- round(y, digits = 0L)
dat <- data.frame(
  ylwr = y_rounded - 0.5,
  yupr = y_rounded + 0.5,
  operator = gl(4L, 3L)
)

library(ggplot2)
lo <- 6:12 + 0.5
up <- 7:13 + 0.5
x <- (lo+up)/2 
cuts <- c(-Inf, 6:13 + 0.5, Inf)
classes <- cut(x, cuts, right=FALSE)
dd <- data.frame(class=classes, lo=lo, up=up)
dd <- transform(dd, Pr=pnorm(up, 10, 1)-pnorm(lo, 10, 1))
ggplot(data=dd, aes(x=class, y=Pr)) +
  geom_bar(stat="identity", colour="black", fill="pink") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Pr(class)") + xlab("class of y")

library(gfilmm)
gfi <- gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ operator, dat, N = 10000L)
gfiSummary(gfi)
gfiConfInt(~ sqrt(sigma_operator^2 + sigma_error^2), gfi)

library(AOV1R)
aovfit <- aov1r(y_rounded ~ dat$operator)
confint(aovfit)




