### ACT-66389
### Mathématiques actuarielles IARD
### Automne 2005
###
### Vincent Goulet, Hélène Cossette
### École d'actuariat, Université Laval
###
### Solutions de la série d'exercices 8


### Fonctions de calcul exact, Panjer et approximation Normal Power.
source("/home/vincent/travail/lib/exact.R")
source("/home/vincent/travail/lib/panjer.R")
source("/home/vincent/travail/lib/np.R")


### Vérification de la fonction 'exact': exemple 4.4 de Loss Models.
fx <- c(0, 15, 20, 25, 12.5, 7.5, 5, 5, 5, 2.5, 2.5)/100
pn <- c(5, 10, 15, 20, 25, 15, 6, 3, 1)/100
fs <- exact(fx, pn)
round(cbind(0:21, fs[1:22]), 5)


###
### EXERCICE 1
###

m <- qnbinom(1-1E-10, size=1.5, mu=100)
k <- qpois(1-1E-10, lambda=3)
fx <- dnbinom(1:m, size=1.5, mu=100)/(1 - dnbinom(0, size=1.5, mu=100))
pn <- dpois(0:k, lambda=3)
fs.conv <- exact(fx, pn)
cbind(0:98*20, round(cumsum(fs.conv), 6)[0:98 * 20 + 1])



###
### EXERCICE 2
###

fs.p <- panjer(fx, freq="poisson", par=list(lambda=3))
all.equal(fs.p, head(fs.conv, length(fs.p)))



###
### EXERCICE 3
###

### Fonction pour faire le calcul direct de F_S(x) dans le cas Poisson
### composée/gamma.
direct <- function(x, lambda, shape, rate)
{
    n <- 1:qpois(1-1E-10, lambda)
    sapply(x, function(x) dpois(0, lambda) + sum(pgamma(x, n * shape, rate) * dpois(n, lambda)))
}


### (a) E[N] = 10, X ~ Gamma(2, 1), h = 0.025

## Caractéristiques des distributions.
h <- 0.025
lambda <- 10; shape <- 2; rate <- 1
E.S <- lambda*shape/rate
Var.S <- lambda*shape*(shape + 1)/rate^2
g1.S <- (shape + 2)/sqrt(lambda * shape * (shape + 1))

## Discrétisation de la distribution de X.
Fx <- pgamma(seq(from=0, to=qgamma(1 - 1E-10, shape, rate), by=h), shape, rate)
fx.u <- diff(Fx)
fx.l <- c(0, fx.u)

## Bornes supérieure et inférieure.
fs.l.a <- panjer(fx.l, freq="poisson", par=list(lambda=lambda))
fs.u.a <- panjer(fx.u, freq="poisson", par=list(lambda=lambda))
Fs.l <- cumsum(fs.l.a)
Fs.u <- cumsum(c(exp(-lambda), fs.u.a))

## Simulation.
r <- (1.645/0.01)^2 * (1 + Var.S/E.S^2) # nombre de simulations
s <- sapply(rpois(r, lambda), function(n) sum(rgamma(n, shape, rate)))
Fn <- ecdf(s) # ogive

## Valeurs de x et F_S(x) à afficher.
m <- which.max(round(Fs.l, 6) > 0) - 1
M <- which.min(round(Fs.u, 6) < 1)
pas <- round((M - m)/30 * h)
x <- seq(from=floor(m * h), to=ceiling(M * h), by=pas)

## Calcul direct par convolutions.
Fs.ex <- direct(x, lambda, shape, rate)

## Approximations normale et Normal Power.
Fs.n <- pnorm((x - E.S)/sqrt(Var.S))
Fs.np <- np(x, mean=E.S, var=Var.S, skewness=g1.S)

## Tableau de résultats.
round(cbind(x=x,
            exact=Fs.ex,
            Fs.L=Fs.l[x/h + 1],
            Fs.U=Fs.u[x/h + 1],
            Normale=Fs.n,
            NP=Fs.np,
            Sim=Fn(x)), 6)

## (b) E[N] = 50,  X ~ Gamma(2, 1), h = 0.025

## Caractéristiques des distributions.
h <- 0.025
lambda <- 50; shape <- 2; rate <- 1
E.S <- lambda*shape/rate
Var.S <- lambda*shape*(shape + 1)/rate^2
g1.S <- (shape + 2)/sqrt(lambda * shape * (shape + 1))

## Discrétisation de la distribution de X.
Fx <- pgamma(seq(from=0, to=qgamma(1 - 1E-10, shape, rate), by=h), shape, rate)
fx.u <- diff(Fx)
fx.l <- c(0, fx.u)

## Bornes supérieure et inférieure.
fs.l.b <- panjer(fx.l, freq="poisson", par=list(lambda=lambda))
fs.u.b <- panjer(fx.u, freq="poisson", par=list(lambda=lambda))
Fs.l <- cumsum(fs.l.b)
Fs.u <- cumsum(c(exp(-lambda), fs.u.b))

## Simulation.
r <- (1.645/0.01)^2 * (1 + Var.S/E.S^2) # nombre de simulations
s <- sapply(rpois(r, lambda), function(n) sum(rgamma(n, shape, rate)))
Fn <- ecdf(s) # ogive

## Valeurs de x et F_S(x) à afficher.
m <- which.max(round(Fs.l, 6) > 0) - 1
M <- which.min(round(Fs.u, 6) < 1) + 1
pas <- round((M - m)/30 * h)
x <- seq(from=floor(m * h), to=ceiling(M * h), by=pas)

## Calcul direct par convolutions.
Fs.ex <- direct(x, lambda, shape, rate)

## Approximations normale et Normal Power.
Fs.n <- pnorm((x - E.S)/sqrt(Var.S))
Fs.np <- np(x, mean=E.S, var=Var.S, skewness=g1.S)

## Tableau de résultats.
round(cbind(x=x,
            exact=Fs.ex,
            Fs.L=Fs.l[x/h + 1],
            Fs.U=Fs.u[x/h + 1],
            Normale=Fs.n,
            NP=Fs.np,
            Sim=Fn(x)), 6)

### (c) E[N] = 50, X ~ LN(ln 2 - 0.1, 0.2), h = 0.025

## Caractéristiques des distributions.
h <- 0.025
lambda <- 50; s2 <- 0.2; mu <- log(2) - s2/2
E.S <- lambda * exp(mu + s2/2)
Var.S <- lambda * exp(2 * (mu + s2))
g1.S <- exp(1.5 * s2)/sqrt(lambda)

## Discrétisation de la distribution de X.
Fx <- plnorm(seq(from=0, to=qlnorm(1 - 1E-10, mu, sqrt(s2)), by=h),
             mu, sqrt(s2))
fx.u <- diff(Fx)
fx.l <- c(0, fx.u)

## Bornes supérieure et inférieure.
fs.l.c <- panjer(fx.l, freq="poisson", par=list(lambda=lambda))
fs.u.c <- panjer(fx.u, freq="poisson", par=list(lambda=lambda))
Fs.l <- cumsum(fs.l.c)
Fs.u <- cumsum(c(exp(-lambda), fs.u.c))

## Simulation.
r <- (1.645/0.01)^2 * (1 + Var.S/E.S^2) # nombre de simulations
s <- sapply(rpois(r, lambda), function(n) sum(rlnorm(n, mu, sqrt(s2))))
Fn <- ecdf(s) # ogive

## Valeurs de x et F_S(x) à afficher.
m <- which.max(round(Fs.l, 6) > 0) - 1
M <- which.min(round(Fs.u, 6) < 1) + 1
pas <- round((M - m)/30 * h)
x <- seq(from=floor(m * h), to=ceiling(M * h), by=pas)

## Approximations normale et Normal Power.
Fs.n <- pnorm((x - E.S)/sqrt(Var.S))
Fs.np <- np(x, mean=E.S, var=Var.S, skewness=g1.S)

## Tableau de résultats.
round(cbind(x=x,
            Fs.L=Fs.l[x/h + 1],
            Fs.U=Fs.u[x/h + 1],
            Normale=Fs.n,
            NP=Fs.np,
            Sim=Fn(x)), 6)
