\name{rcompound}
\alias{rcompound}
\title{Simulation from Compound Models}
\description{
  \code{rcompound} simulates variates from standard (non hierarchical)
  compound models.
}
\usage{
rcompound(n, model.freq, model.sev)
}
\arguments{
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{model.freq, model.sev}{expressions specifying the frequency and
    severity simulation models with the number of variates omitted (see
    details).}
}
\details{
  \code{rcompound} generate variates from a random variable of the form
  \deqn{S = X_1 + ... X_N,}
  where \eqn{N} is the frequency random variable and \eqn{X_1, X_2,
    \dots} are the severity random variables. The latter are mutually
  independent, identically distributed and independent from \eqn{N}.

  \code{model.freq} and \code{model.sev} specify the simulation models
  for the frequency and the severity random variables, respectively. A
  model is a complete call to a random number generation function, with
  the number of variates omitted. This is similar to
  \code{\link{simul}}, but a streamlined interface since the calls are
  not wrapped into \code{\link{expression}}.

  The argument of the random number generation functions for the number
  of variates to simulate \strong{must} be named \code{n}.
}
\value{
  A vector of aggregate amounts \eqn{S_1, \dots, S_n}.
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

\dontrun{## Speed comparison between rcompound() and simul().
## [Also note the simpler syntax for rcompound().]
system.time(rcompound(1e6, rpois(2), rgamma(2, 3)))
system.time(simul(1e6, expression(rpois(2)), expression(rgamma(2, 3))))}
## The severity can itself be a compound model. It makes sense
## in such a case to use a zero-truncated frequency distribution
## for the second level model.
rcompound(10, rpois(2),
          rcompound(rztnbinom(1.5, 0.7), rlnorm(1.2, 1)))
}
\keyword{datagen}
