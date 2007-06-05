### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the risk theory facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)
if(dev.cur() <= 1) get(getOption("device"))()

op <- par(ask = interactive() &&
          (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")),
          col = "black")


###
### DISCRETIZATION OF CLAIM AMOUNT DISTRIBUTIONS
###

### Some numerical techniques to compute the aggregate claim amount
### distribution (see below) require a discrete arithmetic claim amount
### distribution. Function 'discretize' can be used to, well,
### discretize a continuous distribution. (The function can also be
### used to "discretize" [modify the support of] an already discrete
### distribution, but this requires additional care.)
###
### Currently, four discretization methods are supported: upper and
### lower discretization, rounding of the random variable (midpoint
### method), and local matching of the first moment. Usage is similar
### to 'curve' of package 'graphics'.

## Upper and lower discretization of a Gamma(2, 1) distribution with a
## step (or span, or lag) of 0.5. The value of 'to' is chosen so as to
## cover most of the distribution.
x <- seq(0, qgamma(1 - 1E-6, 2, 1), by = 0.5)
xu <- tail(x, 1)
fu <- discretize(pgamma(x, 2, 1), method = "upper",
                 from = 0, to = xu, step = 0.5)
fl <- discretize(pgamma(x, 2, 1), method = "lower",
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fu)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(x, diffinv(fl)), pch = 19, add = TRUE)
par(col = "black")

## Discretization with the rounding method, which has the true cdf
## pass through the midpoints of the intervals [x - step/2, x +
## step/2).
fr <- discretize(pgamma(x, 2, 1), method = "rounding",
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fr)), pch = 19, add = TRUE)
par(col = "black")

## Local matching of the first moment. This ensures that the total
## probability and the expected value on interval [from, to] match
## between the discretized and the true distributions. This requires a
## function to compute the limited expected value of the true
## distribution in any point.
fb <- discretize(pgamma(x, 2, 1), method = "unbiased",
                 lev = levgamma(x, 2, 1),
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(x, diffinv(fb)), pch = 19, add = TRUE)
par(col = "black")

all.equal(diff(pgamma(range(x), 2, 1)),
          sum(fb))                      # same total probability
all.equal(levgamma(xu, 2, 1) - xu * pgamma(xu, 2, 1, lower.tail = FALSE),
          drop(crossprod(x, fb)))       # same expected value

## Comparison of all four methods
fu <- discretize(plnorm(x), method = "upper", from = 0, to = 5)
fl <- discretize(plnorm(x), method = "lower", from = 0, to = 5)
fr <- discretize(plnorm(x), method = "rounding", from = 0, to = 5)
fb <- discretize(plnorm(x), method = "unbiased", from = 0, to = 5,
                 lev = levlnorm(x))
curve(plnorm(x), from = 0, to = 5, lwd = 2)
par(col = "blue")
plot(stepfun(0:4, diffinv(fu)), pch = 19, add = TRUE)
par(col = "red")
plot(stepfun(0:5, diffinv(fl)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(0:4, diffinv(fr)), pch = 19, add = TRUE)
par(col = "magenta")
plot(stepfun(0:5, diffinv(fb)), pch = 19, add = TRUE)
legend(4, 0.2, legend = c("upper", "lower", "rounding", "unbiased"),
       col = c("blue", "red", "green", "magenta"), lty = 1, pch = 19,
       text.col = "black")
par(col = "black")


###
### CALCULATION OF THE AGGREGATE CLAIM AMOUNT DISTRIBUTION
###

### Function 'aggregateDist' computes the aggregate claim amount
### distribution of a portfolio using of the following five methods:
### Panjer's recursive algorithm, convolutions, normal approximation,
### Normal Power II approximation, simulation. (More methods can be
### added in the future.)
###
### The function returns a function to compute the cdf of the
### distribution in any point, much like 'ecdf' of package
### 'stats'. Methods exist to summarize, plot, compute the mean and
### compute quantiles of the aggregate claim amount distribution.

## Recursive method. Requires a discrete claim amount distribution
## typically obtained from 'discretize'. Argument 'x.scale' is used to
## specify how much a value of 1 is really worth.
fx.b <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "unbiased", lev = levgamma(x, 2, 1))
Fs.b <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.b, lambda = 10, x.scale = 0.5)
summary(Fs.b)                           # summary method
knots(Fs.b)                             # support of Fs.b (knots)
Fs.b(knots(Fs.b))                       # evaluation at knots
plot(Fs.b, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs.b)                              # empirical mean
quantile(Fs.b)                          # quantiles

## Convolutions (exact calculation). Requires a vector of
## probabilities for the frequency model. This method can quickly
## become impractical for a large expected number of claims.
pn <- dpois(0:qpois(1-1E-6, 10), 10)
Fs <- aggregateDist("convolution", model.freq = pn, model.sev = fx.b,
                    x.scale = 0.5)
summary(Fs)                             # summary method
knots(Fs)                               # support of Fs (knots)
Fs(knots(Fs))                           # evaluation at knots
plot(Fs, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs)                                # empirical mean
quantile(Fs)                            # quantiles

## Normal approximation. Not hugely useful, but simple to implement...
Fs.n <- aggregateDist("normal", moments = c(20, 60))
summary(Fs.n)                           # summary method
plot(Fs.n, xlim = c(0, 60))             # graphic
mean(Fs.n)                              # true mean
quantile(Fs.n)                          # normal quantiles

## Normal Power II approximation. The approximation is valid for
## values above the expected value only.
Fs.np <- aggregateDist("npower", moments = c(20, 60, 0.516398))
summary(Fs.np)                          # summary method
plot(Fs.np, xlim = c(0, 60))            # truncated graphic

## Simulation method. Function 'simpf' is used to simulate the data
## (see the 'credibility' demo for examples). The result is a step
## function, but with very small steps if the number of simulations is
## large.
Fs.s <- aggregateDist("simulation",
                      model.freq = expression(y = rpois(10)),
                      model.sev = expression(y = rgamma(2, 1)),
                      nb.simul = 10000)
summary(Fs.s)                           # summary method
plot(Fs.s, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs.s)                              # empirical mean
quantile(Fs.s)                          # quantiles

## Graphic comparing the cdfs obtained by a few methods.
fx.u <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "upper")
Fs.u <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.u, lambda = 10, x.scale = 0.5)
fx.l <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "lower")
Fs.l <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.l, lambda = 10, x.scale = 0.5)
par(col = "black")
plot(Fs.b, do.points = FALSE, verticals = TRUE, xlim = c(0, 60), sub = "")
par(col = "blue")
plot(Fs.u, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "red")
plot(Fs.l, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "green")
plot(Fs.s, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "magenta")
plot(Fs.n, add = TRUE, sub = "")
legend(40, 0.2,
       legend = c("recursive + unbiased", "recursive + upper",
                  "recursive + lower", "simulation",
                  "normal approximation"),
       col = c("black", "blue", "red", "green", "magenta"),
       lty = 1, text.col = "black", cex = 1.2)

## Table of quantiles for the same methods as graphic above.
x <- knots(Fs.l)
m <- which.min(x[round(Fs.l(x), 6) > 0])
M <- which.max(x[round(Fs.u(x), 6) < 1])
x <- x[round(seq.int(from = m, to = M, length = 30))]
round(cbind(x = x,
            Lower = Fs.l(x),
            Unbiased = Fs.b(x),
            Upper = Fs.u(x),
            Simulation = Fs.s(x),
            Normal = Fs.n(x)), 6)


par(op)
