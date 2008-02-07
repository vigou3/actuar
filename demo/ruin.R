### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the ruin theory facilities provided by actuar
###
### AUTHOR: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)
#rm(list = ls(all=TRUE))
#library(actuar)

## demo of ruin

## Exponential and hyperexponential (mixture of exponential) distributions are used
## in the first serie of examples. By default, the safety loading is 20 percent.
## premium are always calculated according to the expected value principle.

## Case with an explicit formula: exponential claims and interarrival times.

psi <- ruin(claims = "e", par.claims = list(rate = 1),
            wait   = "e", par.wait   = list(rate = 1), 
			premium = 1.2)

psi(0:10)


## exponential claims and hyper-exponential interarrival times.
psi <- ruin(claims = "e", par.claims = list(rate = 2),
            wait   = "e", par.wait   = list(rate = c(1, 3/2, 1/2), w = c(1/3, 1/2, 1/6)), 
			premium = 1.2)

psi(0:10)

## hyper-exponential claims and interarrival times.
psi <- ruin(claims = "e", par.claims = list(rate = c(1, 3/2, 1/2), w = c(1/3, 1/2, 1/6)),
            wait   = "e", par.wait   = list(rate = c(1/2, 3/4, 1/4), w = c(1/3, 1/2, 1/6)), 
			premium = 0.6)

psi(0:10)

## hyper-exponential claims and exponential interarrival times.
psi <- ruin(claims = "e", par.claims = list(rate = c(1, 3/2, 1/2), w = c(1/3, 1/2, 1/6) ),
            wait   = "e", par.wait   = list(rate = 1/2), 
			premium = 0.6)

psi(0:10)

## 2nd serie : Erlang distribution and mixture of Erlang 
## distribution are now tested. 
## exponential claims and Erlang interarrival times
psi <- ruin(claims = "e", par.claims = list(rate = 2),
            wait   = "Erlang", par.wait   = list(shape = 2, rate = 1), 
			premium = 1.2)

psi(0:10)

## Erlang claims and exponential interarrival times
psi <- ruin(claims = "Erlang", par.claims = list(shape = 2, rate = 2),
            wait   = "e", par.wait   = list(rate = 1), 
			premium = 1.2)

psi(0:10)

## Erlang claims and interarrival times
psi <- ruin(claims = "Erlang", par.claims = list(shape = 2, rate = 2),
            wait   = "Erlang", par.wait   = list(shape = 2, rate = 1), 
			premium = 0.6)

psi(0:10)


## mixture of Erlang distribution for interarrival times and
## exponential claims
psi <- ruin(claims = "e", par.claims = list(rate = 3/2),
            wait   = "Erlang", par.wait   = list(shape = c(3,2,4), 
			rate = c(3,2,4), w=c(1/3, 1/3, 1/3)), 
			premium = 1.2*2/3)

psi(0:10)


## mixture of Erlang distribution for claims and Erlang waiting
## times

psi <- ruin(claims = "Erlang", par.claims = list(shape = c(2, 4), 
			rate = c(1, 3), w=c(1/3, 2/3)),
            wait   = "Erlang", par.wait   = list(shape = 2, rate = 1), 
			premium = 1.2)

psi(0:10)

## 3rd serie : phase-type distributed claims and interarrival times i.e.
## claims follow a generalized Erlang distribution and waiting times
## a mixture of two generalized Erlang distribution
prob.c <- c(1/3, 0, 2/3, 0)
subintens.c <- cbind(c(-1, 0, 0, 0), c(1, -3, 0, 0), c(0, 0, -2, 0), c(0, 0, 2, -3))
mean.c <- mphtype(1, prob.c, subintens.c)

prob.w <- c(1, 0, 0)
subintens.w <- cbind(c(-1, 0, 0), c(1, -2, 0), c(0, 2, -3))
mean.w <- mphtype(1, prob.w, subintens.w)

psi <- ruin(claims = "p", par.claims = list(prob = prob.c, 
			rate = subintens.c),
            wait   = "p", par.wait   = list(prob = prob.w,
			rate = subintens.w), 
			premium = 1.2*mean.c/mean.w)

psi(0:10)


## demo of adjCoef
## 1st no reinsurance : generalized Erlang claims and
## inverse gamma waiting times

claim <- function(x) 1/(1-x)*2/(2-x)*3/(3-x)
wait <- function(x) mgfinvgamma(x, 2, 6/11)

## claims and waiting times have the same mean.
## the adjustment coefficient is an increasing
## of the safety loading.
adjCoef(claim, wait, 1.1, 1)
adjCoef(claim, wait, 1.2, 1)
adjCoef(claim, wait, 1.3, 1)

## 2nd with proportional and excess of loss reinsurance
## examples when claims and waiting times are independent
## can be found in the help.

## we assume here, the dependence is structured by
## a Clayton copula with exponential marginals.
## code from the Copula package.

rclayton <- function(alpha, n) 
{
  val <- cbind(runif(n), runif(n))
  val[,2] <- (val[,1]^(-alpha) * (val[,2]^(-alpha/(alpha + 1)) - 1) + 1)^(-1/alpha)
  return(val)
}
## positive dependence
simu <- rclayton(2, 1000)
claim <- qexp(simu[,1])
wait <- qexp(simu[,2])

## proportional reinsurance with insurer's safety loading
## of 20% and reinsures's one of 30%. premium are calculated
## according to the expected value principle.
premProp <- function(a) mean(claim)/mean(wait)*(-.1 + a*1.3)

emph <- function(r, a) mean( exp( r*(a*claim - premProp(a)*wait) ) )

Rprop <- adjCoef(h = emph, upper = 1, reinsurance = "prop", from = 1/3, to = 1, n=101)

plot(Rprop)

## independence
simu <- rclayton(1, 1000)
claim <- qexp(simu[,1])
wait <- qexp(simu[,2])

RpropInd <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 1/3, to = 1, n=101)

plot(RpropInd, add=TRUE, col="green")

## excess-of-loss reinsurance with insurer's safety loading
## of 20% and reinsures's one of 30%. premium are calculated
## according to the expected value principle.

## positive dependence
simu <- rclayton(2, 1000)
claim <- qexp(simu[,1])
wait <- qexp(simu[,2])

premEOL <- function(L) mean(claim)/mean(wait) * (-.1) + 1.3 * mean(pmin(L, claim))/mean(wait)

h <- function(r, L) mean( exp( r*(pmin(L, claim) - premEOL(L)*wait) ) )

REOL <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 0, to = 10, n=101)

plot(REOL)

## independence
simu <- rclayton(1, 1000)
claim <- qexp(simu[,1])
wait <- qexp(simu[,2])

REOL2 <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 0, to = 10, n=101)

plot(REOL2, add=TRUE, col="green")
