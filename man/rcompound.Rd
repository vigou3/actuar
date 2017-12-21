\name{rcompound}
\alias{rcompound}
\alias{rcomppois}
\title{Simulation from Compound Models}
\description{
  \code{rcompound} generates random variates from a compound model.

  \code{rcomppois} is a simplified version for a common case.
}
\usage{
rcompound(n, model.freq, model.sev, SIMPLIFY = TRUE)

rcomppois(n, lambda, model.sev, SIMPLIFY = TRUE)}
\arguments{
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{model.freq, model.sev}{expressions specifying the frequency and
    severity simulation models with the number of variates omitted (see
    details).}
  \item{lambda}{Poisson parameter.}
  \item{SIMPLIFY}{boolean; if \code{FALSE} the frequency and severity
    variates are returned along with the aggregate variates.}
}
\details{
  \code{rcompound} generates variates from a random variable of the form
  \deqn{S = X_1 + ... X_N,}
  where \eqn{N} is the frequency random variable and \eqn{X_1, X_2,
    \dots} are the severity random variables. The latter are mutually
  independent, identically distributed and independent from \eqn{N}.

  \code{model.freq} and \code{model.sev} specify the simulation models
  for the frequency and the severity random variables, respectively. A
  model is a complete call to a random number generation function, with
  the number of variates omitted. This is similar to
  \code{\link{simul}}, but the calls need not be wrapped into
  \code{\link{expression}}. Either argument may also be the name of an
  object containing an expression, in which case the object will be
  evaluated in the evaluation frame to retrieve the expression.

  The argument of the random number generation functions for the number
  of variates to simulate \strong{must} be named \code{n}.

  \code{rcomppois} generates variates from the common Compound Poisson
  model, that is when random variable \eqn{N} is Poisson distributed
  with mean \code{lambda}.
}
\value{
  When \code{SIMPLIFY = TRUE}, a vector of aggregate amounts \eqn{S_1,
  \dots, S_n}. 

  When \code{SIMPLIFY = FALSE}, a list of three elements:
  \item{\code{aggregate}}{vector of aggregate amounts \eqn{S_1, \dots,
      S_n};}
  \item{\code{frequency}}{vector of frequencies \eqn{N_1, \dots,
      N_n};}
  \item{\code{severity}}{vector of severities \eqn{X_1, X_2, \dots}.}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{simul}} to simulate from compound hierarchical models.
}
\examples{
## Compound Poisson model with gamma severity.
rcompound(10, rpois(2), rgamma(2, 3))
rcomppois(10, 2, rgamma(2, 3))          # same

## Frequencies and individual claim amounts along with aggregate
## values.
rcomppois(10, 2, rgamma(2, 3), SIMPLIFY = FALSE)

## Wrapping the simulation models into expression() is allowed, but
## not needed.
rcompound(10, expression(rpois(2)), expression(rgamma(2, 3)))

\dontrun{## Speed comparison between rcompound() and simul().
## [Also note the simpler syntax for rcompound().]
system.time(rcompound(1e6, rpois(2), rgamma(2, 3)))
system.time(simul(1e6, expression(rpois(2)), expression(rgamma(2, 3))))}
## The severity can itself be a compound model. It makes sense
## in such a case to use a zero-truncated frequency distribution
## for the second level model.
rcomppois(10, 2,
          rcompound(rztnbinom(1.5, 0.7), rlnorm(1.2, 1)))
}
\keyword{datagen}
