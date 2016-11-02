\name{ruin}
\alias{ruin}
\alias{plot.ruin}
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
     maxit = 200L, echo = FALSE)

\method{plot}{ruin}(x, from = NULL, to = NULL, add = FALSE,
     xlab = "u", ylab = expression(psi(u)),
     main = "Probability of Ruin", xlim = NULL, \dots)
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
  \item{x}{an object of class \code{"ruin"}.}
  \item{from, to}{the range over which the function will be plotted.}
  \item{add}{logical; if \code{TRUE} add to already existing plot.}
  \item{xlim}{numeric of length 2; if specified, it serves as default
    for \code{c(from, to)}.}
  \item{xlab, ylab}{label of the x and y axes, respectively.}
  \item{main}{main title.}
  \item{\dots}{further graphical parameters accepted by
    \code{\link[graphics]{curve}}.}
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
  A function of class \code{"ruin"} inheriting from the
  \code{"function"} class to compute the probability of ruin given
  initial surplus levels. The function has arguments:
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
plot(psi, from = 0, to = 10)

## Mixture of two exponentials for claims, exponential interarrival
## times (Gerber 1979)
psi <- ruin(claims = "e", par.claims = list(rate = c(3, 7), w = 0.5),
            wait   = "e", par.wait   = list(rate = 3), pre = 1)
u <- 0:10
psi(u)
(24 * exp(-u) + exp(-6 * u))/35	# same

## Phase-type claims, exponential interarrival times (Asmussen and
## Rolski 1991)
p <- c(0.5614, 0.4386)
r <- matrix(c(-8.64, 0.101, 1.997, -1.095), 2, 2)
lambda <- 1/(1.1 * mphtype(1, p, r))
psi <- ruin(claims = "p", par.claims = list(prob = p, rates = r),
            wait   = "e", par.wait   = list(rate = lambda))
psi
plot(psi, xlim = c(0, 50))

## Phase-type claims, mixture of two exponentials for interarrival times
## (Asmussen and Rolski 1991)
a <- (0.4/5 + 0.6) * lambda
ruin(claims = "p", par.claims = list(prob = p, rates = r),
     wait   = "e", par.wait   = list(rate = c(5 * a, a), weights =
                                     c(0.4, 0.6)),
     maxit = 225L)
}
\keyword{models}
