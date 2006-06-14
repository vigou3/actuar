### Fonction pour faire le calcul direct de F_S(x) dans le cas Poisson
### composÅ√Å©e/gamma.
direct <- function(x, lambda, shape, rate)
{
    n <- 1:qpois(1-1E-10, lambda)
    sapply(x, function(x) dpois(0, lambda) + sum(pgamma(x, n * shape, rate) * dpois(n, lambda)))
}


### (a) E[N] = 10, X ~ Gamma(2, 1), h = 0.025

## CaractÅ√Å©ristiques des distributions.
h <- 0.025
lambda <- 10; shape <- 2; rate <- 1
E.S <- lambda*shape/rate
Var.S <- lambda*shape*(shape + 1)/rate^2
g1.S <- (shape + 2)/sqrt(lambda * shape * (shape + 1))

## DiscrÅ√Å©tisation de la distribution de X.
Fx <- pgamma(seq(from=0, to=qgamma(1 - 1E-10, shape, rate), by=h), shape, rate)
fx.u <- diff(Fx)
fx.l <- c(0, fx.u)

## Bornes supÅ√Å©rieure et infÅ√Å©rieure.
fs.l.a <- panjer(fx.l, freq="poisson", par=list(lambda=lambda))
fs.u.a <- panjer(fx.u, freq="poisson", par=list(lambda=lambda))
Fs.l <- cumsum(fs.l.a$fs)
Fs.u <- cumsum(c(exp(-lambda), fs.u.a$fs))

## Simulation.
r <- (1.645/0.01)^2 * (1 + Var.S/E.S^2) # nombre de simulations
s <- sapply(rpois(r, lambda), function(n) sum(rgamma(n, shape, rate)))
Fn <- ecdf(s) # ogive

## Valeurs de x et F_S(x) Å√Å† afficher.
m <- which.max(round(Fs.l, 6) > 0) - 1
M <- which.min(round(Fs.u, 6) < 1)
pas <- round((M - m)/30 * h)
x <- seq(from=floor(m * h), to=ceiling(M * h), by=pas)

## Calcul direct par convolutions.
Fs.ex <- direct(x, lambda, shape, rate)

## Approximations normale et Normal Power.
Fs.n <- pnorm((x - E.S)/sqrt(Var.S))
Fs.np <- np2(x, mean=E.S, var=Var.S, skewness=g1.S)

## Approximation WHgen
Fs.WH <- WHgen(x, mean = E.S, var = Var.S, skewness = g1.S)

## Tableau de rÅ√Å©sultats.
round(cbind(x=x,
            exact=Fs.ex,
            Fs.L=Fs.l[x/h + 1],
            Fs.U=Fs.u[x/h + 1],
            Normale=Fs.n,
            NP2=Fs.np$Fs,
            Sim=Fn(x)), 6)

##########
model.freq <- list(dist = "pois", par = c(lambda = 4))
model.sev <- list(dist = "norm", par = c(mean = 500, sd = 150))


