\name{ruin}
\alias{ruin}
\title{Probability of Ruin}
\description{
  Calulation of infinite time probability of ruin in the models of
  \enc{Cramér}{Cramer}-Lundberg and Sparre Andersen, that is with
  exponential or phase-type (including mixtures of exponentials, Erlang
  and mixture of Erlang) claims interarrival time.
}
\usage{
ruin(claims = c("exponential", "Erlang", "phase-type"), par.claims,
     wait = c("exponential", "Erlang", "phase-type"), par.wait,
     premium.rate = 1, tol = sqrt(.Machine$double.eps),
     maxit = 200, echo = FALSE)
}
\arguments{
  \item{claims}{character; the type of claim severity distribution.}
  \item{wait}{character; the type of claim interarrival (wait) time
    distribution.}
  \item{par.claims, par.wait}{named list containing the parameters of
    the distribution (see details).}
  \item{premium.rate}{numeric vector of length 1; the premium rate.}
  \item{tol, maxit, echo}{respectively the tolerance level of the
    stopping criteria, the maximum number of iterations and whether or
    not to echo the procedure when the transition rates matrix is
    determined iteratively. Ignored if \code{wait = "exponential"}.}
}
\details{
  The names of the parameters in \code{par.claims} and \code{par.wait}
  must the same as in \code{\link[stats]{dexp}},
  \code{\link[stats]{dgamma}} or \code{\link{dphtype}}, as appropriate.
  A model will be a mixture of exponential or Erlang distributions (but
  not phase-type) when the parameters are vectors of length \eqn{> 1}
  and the parameter list contains a vector \code{weights} of the
  coefficients of the mixture.

  Parameters are recycled when needed. Their names can be abbreviated.

  Combinations of exponentials as defined in Dufresne and Gerber (1988)
  are \emph{not} supported.
  
  Ruin probabilities are evaluated using \code{\link{pphtype}} except
  when both distributions are exponential, in which case an explicit
  formula is used.
  
  When \code{wait != "exponential"} (Sparre Andersen model), the
  transition rate matrix \eqn{\boldsymbol{Q}}{Q} of the distribution of
  the probability of ruin is determined iteratively using a fixed
  point-like algorithm. The stopping criteria used is%
  \deqn{\max \left\{
    \sum_{j = 1}^n |\boldsymbol{Q}_{ij} - \boldsymbol{Q}_{ij}^\prime|
    \right\} < \code{tol},}{%
    max(rowSum(|Q - Q'|)) < tol,}%
  where \eqn{\boldsymbol{Q}}{Q} and \eqn{\boldsymbol{Q}^\prime}{Q'} are
  two successive values of the matrix.
}
\value{
  A function to compute the probability of ruin given initial surplus
  levels. The function has arguments:
  \item{u}{numeric vector of initial surplus levels;}
  \item{survival}{logical; if \code{FALSE} (default), probabilities are
    \eqn{\psi(u)}{psi(u)}, otherwise, \eqn{\phi(u) = 1 - \psi(u)}{phi(u)
      = 1 - psi(u)};}
  \item{lower.tail}{an alias for \code{!survival}.}
}
\references{
  Asmussen, S. and Rolski, T. (1991), Computational methods in risk theory:
  A matrix algorithmic approach, \emph{Insurance: Mathematics and
    Economics} \bold{10}, 259--274.

  Dufresne, F. and Gerber, H. U. (1988), Three methods to calculate
  the probability of ruin, \emph{Astin Bulletin} \bold{19}, 71--90.

  Gerber, H. U. (1979), \emph{An Introduction to Mathematical Risk
    Theory}, Huebner Foundation.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  and Christophe Dutang
}
\examples{
## Case with an explicit formula: exponential claims and exponential
## interarrival times.
psi <- ruin(claims = "e", par.claims = list(rate = 5),
            wait   = "e", par.wait   = list(rate = 3))
psi
psi(0:10)
curve(psi(x), xlim = c(0, 10))

## Example of Gerber (1979): mixture of two exponentials for claims,
## exponential interarrival times
psi <- ruin(claims = "e", par.claims = list(rate = c(3, 7), w = 0.5),
            wait   = "e", par.wait   = list(rate = 3), pre = 1)
u <- 0:10
psi(u)
(24 * exp(-u) + exp(-6 * u))/35	# same

## Figure 3 of Asmussen & Rolski (1992) - Phase-type claims
## and exponential interarrival times
##
## (Small differences within numerical error margin)
prob.B1 <- c(0.5614, 0.4386)
rates.B1 <- matrix(c(-8.64, 0.101, 1.997, -1.095), 2, 2)
lambda1 <- 1/(1.1 * mphtype(1, prob.B1, rates.B1))

prob.B2 <- c(0.9731, 0.0152, 0.0106, 0.001)
rates.B2 <- matrix(c(-28.648,  0.102,  0.113,  0.100,
                      28.532, -8.255,  0.107,  0.102,
                       0.089,  8.063, -5.807,  0.111,
                       0.027,  0.086,  5.296, -2.176), 4, 4)
lambda2 <- 1/(1.1 * mphtype(1, prob.B2, rates.B2))

f0 <- ruin(claims = "e", par.claims = list(rate = 1),
           wait   = "e", par.wait   = list(rate = 1/1.1))
f1 <- ruin(claims = "p", par.claims = list(prob = prob.B1, rates = rates.B1),
           wait   = "e", par.wait   = list(rate = lambda1))
f2 <- ruin(claims = "p", par.claims = list(prob = prob.B2, rates = rates.B2),
           wait   = "e", par.wait   = list(rate = lambda2))

curve(f0, xlim = c(0, 50), ylim = c(1E-4, 1), log = "y")
curve(f1, add = TRUE, col = "blue")
curve(f2, add = TRUE, col = "orange")

## Example of p. 266 of Asmussen & Rolski (1992) - Phase-type claims
## and mixture of two exponentials for interarrival times
##
## (Somewhat larger differences than above, but still reasonable to
## point at different numerical methods)
a <- (0.4/5 + 0.6) * lambda1
ruin(claims = "p", par.claims = list(prob = prob.B1, rates = rates.B1),
     wait   = "e", par.wait   = list(rate = c(5 * a, a), weights =
c(0.4, 0.6)), p = 1, maxit = 300, echo = TRUE)
}
\keyword{models}
