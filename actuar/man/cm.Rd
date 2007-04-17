\name{cm}
\alias{cm}
\alias{print.cm}
\alias{predict.cm}
\alias{summary.cm}
\alias{print.summary.cm}
\title{Hierarchical Credibility Models}
\description{
  Fit a credibility model in the formulation of variance components as
  described in Dannenburg, Kaas and Goovaerts (1996). Models supported
  are part of a generalized hierarchical credibility theory as
  introduced in Dannenburg (1995).
}
\usage{
cm(formula, data, ratios, weights, subset, TOL = 1E-6, echo = FALSE)

\method{print}{cm}(x, \dots)

\method{predict}{cm}(object, levels = NULL, \dots)

\method{summary}{cm}(object, levels = NULL, \dots)

\method{print}{summary.cm}(x, \dots)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}
  \item{data}{a matrix or a data frame containing the portfolio
    structure, the ratios or claim amounts and their associated weights,
    if any.}
  \item{ratios}{expression indicating the columns of \code{data}
    containing the ratios or claim amounts.}
  \item{weights}{expression indicating the columns of \code{data}
    containing the weights associated with \code{ratios}.}
  \item{subset}{an optional logical expression indicating a subset of
    observations to be used in the modeling process. All observations
    are included by default.}
  \item{TOL}{maximum absolute error in the iterative calculation of the
    variance estimators.}
  \item{echo}{logical; whether to echo the iterative procedure or not}
  \item{x, object}{an object of class \code{"cm"}}
  \item{levels}{character vector indicating the levels to predict or to
    include in the summary; if \code{NULL} all levels are included.}
  \item{\dots}{additional attributes to attach to the result for the
    \code{predict} and \code{summary} methods; further arguments to
    \code{\link[base]{format}} for the \code{print.summary} method;
    unused for the \code{print} method.}
}
\details{
  The \code{formula} argument symbolically describes the structure of
  the portfolio in the form \eqn{~ terms}. Each term is an interaction
  between risk factors contributing to the total variance of the
  portfolio data. At the moment, \code{cm} can only handle hierarchical
  interactions, with any given number of levels.  The Bühlmann-Straub
  model is a special case with only one level.

  Terms are separated by \code{+} operators and interactions within each
  term by \code{:}. For a portfolio divided first into sectors, then
  units and finally contracts, \code{formula} would be \eqn{~ sector +
  sector:unit + sector:unit:contract}, where \code{sector}, \code{unit}
  and \code{contract} are column names in \code{data}. In general, the
  formula should be of the form \eqn{~ a + a:b + a:b:c + a:b:c:d +
    ...}.

  Arguments \code{ratios}, \code{weights} and \code{subset} are used
  like arguments \code{select}, \code{select} and \code{subset} of
  function \code{\link[base]{subset}}.

  Data does not have to be sorted by level. Nodes with no data (complete
  lines of \code{NA} except for the portfolio structure) are allowed.
  
  The credibility premium at one level is a convex combination between
  the linearly sufficient statistic of a node and the credibility
  premium of the level above. (For the first level, the complement of
  the credibility is given to the collective premium.) The linearly
  sufficient statistic one node is the credibility weighted average of
  the data of the node, except at the last level, where natural weights
  are used. The credibility factor of node \eqn{i} is equal to
  \deqn{\frac{w_i}{w_i + B/W},}{w[i]/(w[i] + B/W),}
  where \eqn{w_i}{w[i]} is the weight of the node used in the linearly
  sufficient statistic, \eqn{B} is the average between node variance and
  \eqn{W} is the average within node variance.

  Estimators of the variance structure parameters are pseudo-estimators
  of the form
  \deqn{\frac{1}{d} \sum_i w_i (X_i - \bar{X})^2,}{%
    sum(i; w[i] * (X[i] - Xb)^2)/d,}
  where \eqn{X_i}{X[i]} is the linearly sufficient statistic of one
  level, \eqn{\bar{X}}{Xb} is the linearly sufficient statistic of
  the level above and \eqn{d} is the effective number of nodes at one
  level minus the effective number of nodes of the level above.
}
\value{
  Function \code{cm} computes the structure parameters estimators of the
  model specified in \code{formula}. The value returned is an object of
  class \code{cm}.

  An object of class \code{"cm"} is a list with the following
  components:
  \item{means}{a list containing, for each level, the vector of linearly
    sufficient statistics.}
  \item{weights}{a list containing, for each level, the vector of total
    weights.}
  \item{variances}{a vector containing the variance structure parameters
    estimators.}
  \item{cred}{a list containining, for each level, the vector of
    credibility factors.}
  \item{levels}{the columns of \code{data} containing the portfolio
    structure.}
  \item{nodes}{a list containing, for each level, the vector of the
    number of nodes in the level.}
  \item{ordering}{a list containing, for each level, the affiliation of
    a node to the node of the level above.}
  \item{call}{the function call.}

  The method of \code{predict} for objects of class \code{"cm"} computes
  the credibility premiums for the nodes of every level included in
  argument \code{levels} (all by default). Result is always a list the
  same length as \code{levels} or the number of levels in
  \code{formula}.
}
\references{
  Goulet, V. (1998), \emph{Principles and Application of Credibility
    Theory}, Journal of Actuarial Practice, Volume 6, ISSN 1064-6647.
  
  Goovaerts, M. J. and Hoogstad, W. J. (1987), \emph{Credibility Theory},
  Surveys of Actuarial Studies, NO.4, Nationale-Nederlanden N.V.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  and Louis-Philippe Pouliot
}
\seealso{
  \code{\link{bstraub}} for faster calculations of simple Bühlmann or
  Bühlmann-Straub models; \code{\link[base]{subset}},
  \code{\link[base]{formula}}.
}
\examples{
data(hachemeister)

## Simple Bühlmann-Straub model
cm(~state, hachemeister,
   ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
bstraub(hachemeister[, 2:13], hachemeister[, 14:25]) # same, faster

## Two-level hierarchical model. Notice that data does not have
## to be sorted by level
X <- cbind(unit = c(1, 2, 1, 2, 2), hachemeister)
fit <- cm(~unit + unit:state, X, ratio.1:ratio.12, weight.1:weight.12)
fit				# print method
summary(fit)			# more details
predict(fit)			# credibility premiums
}
\keyword{models}
