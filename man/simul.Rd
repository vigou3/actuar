\name{simul}
\alias{simul}
\alias{simpf}
\alias{print.portfolio}
\title{Simulation of a portfolio of data}
\description{
  \code{simul} simulates data for insurance applications allowing
  hierarchical structures and separate models for the frequency and
  severity of claims distributions.
}
\usage{
simul(nodes, model.freq = NULL, model.sev = NULL, weights = NULL)

\method{print}{portfolio}(x, \dots)
}
\arguments{
  \item{nodes}{a named list giving the number of "nodes" at each level
    in the hierarchy of the portfolio. The nodes are listed from top
    (portfolio) to bottom (usually the years of experience).}
  \item{model.freq}{a named vector of expressions specifying the
    frequency of claims model (see details); if \code{NULL}, only claim
    amounts are simulated.}
  \item{model.sev}{a named vector of expressions specifying the severity
    of claims model (see details); if \code{NULL}, only claim numbers
    are simulated.}
  \item{weights}{a vector of weights.}
  \item{x}{a \code{portfolio} object.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  The order and the names of the elements in \code{nodes},
  \code{model.freq} and \code{model.sev} must match.  At least one of
  \code{model.freq} and \code{model.sev} must be non \code{NULL}.

  \code{nodes} specifies the hierarchical layout of the portfolio.  Each
  element of the list is a vector of the number of nodes at a given
  level.  Vectors are recycled as necessary.

  \code{model.freq} and \code{model.sev} specify the simulation models
  for claim numbers and claim amounts, respectively.  A model is
  expressed in a semi-symbolic fashion using an object of mode
  \code{\link[base]{expression}}.  Each element of the object
  must be named and should be a complete call to a random number
  generation function, with the number of variates omitted.  Hierarchical
  (or mixtures of) models are achieved by replacing one or more
  parameters of a distribution at a given level by any combination of
  the names of the levels above.  If no mixing is to take place at a
  level, the model for this level can be \code{NULL}.

  The argument of the random number generation functions for the number
  of variates to simulate \strong{must} be named \code{n}.

  Weights will be used wherever the name \code{"weights"} appears in a
  model.  It is the user's responsibility to ensure that the length of
  \code{weights} will match the number of nodes when weights are to be
  used.  Normally, there should be one weight per node at the lowest
  level of the model.

  Data is generated in lexicographic order, that is by row in the output
  matrix.
}
\value{
  An object of \code{\link[base]{class}} \code{"portfolio"}. A
  \code{print} method for this class displays the models used in the
  simulation as well as the frequency of claims for each year and entity
  in the portfolio.

  An object of class \code{"portfolio"} is a list containing the
  following components:
  \item{data}{a two dimension list where each element is a vector of
    claim amounts;}
  \item{weights}{the vector of weights given in argument reshaped as a
    matrix matching element \code{data}, or \code{NULL};}
  \item{classification}{a matrix of integers where each row is a unique
    set of subscripts identifying an entity in the portfolio
    (e.g. integers \eqn{i}, \eqn{j} and \eqn{k} for data
    \eqn{X_{ijkt}}{X[ijkt]});}
  \item{nodes}{the \code{nodes} argument, appropriately recycled;}
  \item{model.freq}{the frequency model as given in argument;}
  \item{model.sev}{the severity model as given in argument.}

  It is recommended to manipulate objects of class \code{"portfolio"} by
  means of the corresponding methods of functions \code{aggregate},
  \code{frequency} and \code{severity}.
}
\references{
  Goulet, V. (2006), \emph{Credibility for severity
    revisited}, North American Actuarial Journal \bold{10}, 49--62.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Sébastien Auclair and Louis-Philippe Pouliot
}
\seealso{
  \code{\link{simul.summaries}} for the functions to create the
  matrices of aggregate claim amounts, frequencies and individual
  claim amounts.
}
\examples{
## Simple two level (contracts and years) portfolio with frequency model
## Nit|Theta_i ~ Poisson(Theta_i), Theta_i ~ Gamma(2, 3) and severity
## model X ~ Lognormal(5, 1)
simul(nodes = list(contract = 10, year = 5),
      model.freq = expression(contract = rgamma(2, 3),
                              year = rpois(contract)),
      model.sev = expression(contract = NULL,
                             year = rlnorm(5, 1)))

## Model with weights and mixtures for both frequency and severity
## models
nodes <- list(entity = 8, year = c(5, 4, 4, 5, 3, 5, 4, 5))
mf <- expression(entity = rgamma(2, 3),
                 year = rpois(weights * entity))
ms <- expression(entity = rnorm(5, 1),
                 year = rlnorm(entity, 1))
wit <- sample(2:10, 35, replace = TRUE)
pf <- simul(nodes, mf, ms, wit)
pf 				# print method
weights(pf)			# extraction of weights
aggregate(pf)[, -1]/weights(pf)[, -1] # ratios

## Four level hierarchical model for frequency only
nodes <- list(sector = 3, unit = c(3, 4),
              employer = c(3, 4, 3, 4, 2, 3, 4), year = 5)
mf <- expression(sector = rexp(1),
                 unit = rexp(sector),
                 employer = rgamma(unit, 1),
                 year = rpois(employer))
pf <- simul(nodes, mf, NULL)
pf 				# print method
aggregate(pf) 			# aggregate claim amounts
frequency(pf)  			# frequencies
severity(pf)			# individual claim amounts
}
\keyword{datagen}
