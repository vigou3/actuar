\name{adjCoef}
\alias{adjCoef}
\alias{plot.adjCoef}
\title{Adjustment Coefficient}
\description{
  Compute the adjustment coefficient in ruin theory, or return a
  function to compute the adjustment coefficient for various reinsurance
  retentions.
}
\usage{
adjCoef(mgf.claim, mgf.wait = mgfexp, premium.rate, upper.bound,
        h, reinsurance = c("none", "proportional", "excess-of-loss"),
        from, to, n = 101)

\method{plot}{adjCoef}(x, xlab = "x", ylab = "R(x)",
     main = "Adjustment Coefficient", sub = comment(x),
     type = "l", add = FALSE, \dots)
}
\arguments{
  \item{mgf.claim}{an expression written as a function of \code{x} or of
    \code{x} and \code{y}, or alternatively the name of a function,
    giving the moment generating function (mgf) of the claim severity
    distribution.}
  \item{mgf.wait}{an expression written as a function of \code{x}, or
    alternatively the name of a function, giving the mgf of the
    claims interarrival time distribution. Defaults to an exponential
    distribution with mean 1.}
  \item{premium.rate}{if \code{reinsurance = "none"}, a numeric value of
    the premium rate; otherwise, an expression written as a function of
    \code{y}, or alternatively the name of a function, giving the
    premium rate function.}
  \item{upper.bound}{numeric; an upper bound for the coefficient,
    usually the upper bound of the support of the claim severity mgf.}
  \item{h}{an expression written as a function of \code{x} or of
    \code{x} and \code{y}, or alternatively the name of a function,
    giving function \eqn{h} in the Lundberg equation (see below);
    ignored if \code{mgf.claim} is provided.}
  \item{reinsurance}{the type of reinsurance for the portfolio; can
     be abbreviated.}
  \item{from, to}{the range over which the adjustment coefficient will
    be calculated.}
  \item{n}{integer; the number of values at which to evaluate the
    adjustment coefficient.}
  \item{x}{an object of class \code{"adjCoef"}.}
  \item{xlab, ylab}{label of the x and y axes, respectively.}
  \item{main}{main title.}
  \item{sub}{subtitle, defaulting to the type of reinsurance.}
  \item{type}{1-character string giving the type of plot desired; see
    \code{\link[graphics]{plot}} for details.}
  \item{add}{logical; if \code{TRUE} add to already existing plot.}
  \item{\dots}{further graphical parameters accepted by
    \code{\link[graphics]{plot}} or \code{\link[graphics]{lines}}.}
}
\details{
  In the typical case \code{reinsurance = "none"}, the coefficient of
  determination is the smallest (strictly) positive root of the Lundberg
  equation%
  \deqn{h(x) = E[e^{x B - x c W}] = 1}{h(x) = E[exp(x B - x c W)] = 1}%
  on \eqn{[0, m)}, where \eqn{m =} \code{upper.bound}, \eqn{B} is the
  claim severity random variable, \eqn{W} is the claim interarrival
  (or wait) time random variable and \eqn{c =} \code{premium.rate}. The
  premium rate must satisfy the positive safety loading constraint
  \eqn{E[B - c W] < 0}.

  With \code{reinsurance = "proportional"}, the equation becomes
  \deqn{h(x, y) = E[e^{x y B - x c(y) W}] = 1,}{%
    h(x, y) = E[exp(x y B - x c(y) W)] = 1,}
  where \eqn{y} is the retention rate and \eqn{c(y)} is the premium rate
  function.

  With \code{reinsurance = "excess-of-loss"}, the equation becomes
  \deqn{h(x, y) = E[e^{x \min(B, y) - x c(y) W}] = 1,}{%
    h(x, y) = E[exp(x min(B, y) - x c(y) W)] = 1,}
  where \eqn{y} is the retention limit and \eqn{c(y)} is the premium rate
  function.

  One can use argument \code{h} as an alternative way to provide
  function \eqn{h(x)} or \eqn{h(x, y)}. This is necessary in cases where
  random variables \eqn{B} and \eqn{W} are not independent.

  The root of \eqn{h(x) = 1} is found by minimizing \eqn{(h(x) - 1)^2}.
}
\value{
  If \code{reinsurance = "none"}, a numeric vector of length one.
  Otherwise, a function of class \code{"adjCoef"} inheriting from the
  \code{"function"} class.
}
\references{
  Bowers, N. J. J., Gerber, H. U., Hickman, J., Jones, D. and Nesbitt, C.
  (1986), \emph{Actuarial Mathematics}, Society of Actuaries.

  Centeno, M. d. L. (2002), Measuring the effects of reinsurance by the
  adjustment coefficient in the Sparre-Anderson model, \emph{Insurance:
    Mathematics and Economics} \bold{30}, 37--49.

  Gerber, H. U. (1979), \emph{An Introduction to Mathematical Risk
    Theory}, Huebner Foundation.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2008),
  \emph{Loss Models, From Data to Decisions, Third Edition}, Wiley.
}
\author{
  Christophe Dutang, Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Basic example: no reinsurance, exponential claim severity and wait
## times, premium rate computed with expected value principle and
## safety loading of 20\%.
adjCoef(mgfexp, premium = 1.2, upper = 1)

## Same thing, giving function h.
h <- function(x) 1/((1 - x) * (1 + 1.2 * x))
adjCoef(h = h, upper = 1)

## Example 11.4 of Klugman et al. (2008)
mgfx <- function(x) 0.6 * exp(x) + 0.4 * exp(2 * x)
adjCoef(mgfx(x), mgfexp(x, 4), prem = 7, upper = 0.3182)

## Proportional reinsurance, same assumptions as above, reinsurer's
## safety loading of 30\%.
mgfx <- function(x, y) mgfexp(x * y)
p <- function(x) 1.3 * x - 0.1
h <- function(x, a) 1/((1 - a * x) * (1 + x * p(a)))
R1 <- adjCoef(mgfx, premium = p, upper = 1, reins = "proportional",
              from = 0, to = 1, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "p",
             from = 0, to = 1, n = 101)
R1(seq(0, 1, length = 10))	# evaluation for various retention rates
R2(seq(0, 1, length = 10))	# same
plot(R1)		        # graphical representation
plot(R2, col = "green", add = TRUE) # smoother function

## Excess-of-loss reinsurance
p <- function(x) 1.3 * levgamma(x, 2, 2) - 0.1
mgfx <- function(x, l)
    mgfgamma(x, 2, 2) * pgamma(l, 2, 2 - x) +
    exp(x * l) * pgamma(l, 2, 2, lower = FALSE)
h <- function(x, l) mgfx(x, l) * mgfexp(-x * p(l))
R1 <- adjCoef(mgfx, upper = 1, premium = p, reins = "excess-of-loss",
             from = 0, to = 10, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "e",
             from = 0, to = 10, n = 101)
plot(R1)
plot(R2, col = "green", add = TRUE)
}
\keyword{optimize}
\keyword{univar}
