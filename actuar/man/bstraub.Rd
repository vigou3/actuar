\name{bstraub}
\alias{bstraub}
\alias{print.bstraub}
\alias{predict.bstraub}
\alias{summary.bstraub}
\alias{print.summary.bstraub}
\title{Buhlmann-Straub Credibility Model}
\description{
  \code{bstraub} computes structure parameters estimators in the
  Bühlmann-Straub credibility model and \code{predict.bstraub} computes
  the credibility premiums.
}
\usage{
bstraub(ratios, weights,
        heterogeneity = c("iterative", "unbiased"),
        TOL = 1e-06, echo = FALSE)

\method{print}{bstraub}(x, \dots)

\method{predict}{bstraub}(object, \dots)

\method{summary}{bstraub}(object, \dots)

\method{print}{summary.bstraub}(x, \dots)
}
\arguments{
  \item{ratios}{a matrix of ratios (contracts in lines, years in columns).}
  \item{weights}{a matrix of weights corresponding to ratios.}
  \item{heterogeneity}{estimator of the between contract heterogeneity
    parameter used in premium calculation; \code{"iterative"} for the
    Bischel-Straub estimator; \code{"unbiased"} for the usual
    Bühlmann-Straub estimator (see below).}
  \item{TOL}{maximum relative error in the iterative procedure.}
  \item{echo}{logical; whether to echo iterative procedure or not.}
  \item{x, object}{an object of class \code{"bstraub"}.} 
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  The credibility premium of contract \eqn{i} is given by
  \deqn{z_i X_{iw} + (1 - z_i) X_{zw},}{z[i] X[iw] + (1 - z[i]) X[zw],}
  where
  \deqn{z_{i} = \frac{w_{i\cdot} \hat{a}}{w_{i\cdot} \hat{a} + \hat{s}^2},}{%
    z[i] = (w[i.] a)/(w[i.] a + s^2),}
  \eqn{X_{iw}}{X[iw]} is the weighted average of the ratios of contract
  \eqn{i}, \eqn{X_{zw}}{X[zw]} is the weighted average of the matrix of
  ratios using credibility factors and \eqn{w_{i\cdot}}{w[i.]} is
  the total weight of a contract. \eqn{\hat{s}^2}{s^2} is the estimator
  of the within contract heterogeneity and \eqn{\hat{a}}{a} is the
  estimator of the between contract heterogeneity.

  Missing data are represented by \code{NA} in both the matrix of ratios
  and the matrix of weights. The function can cope with complete lines
  of \code{NA} in case a contract has no experience.

  \code{bstraub} computes the structure parameters estimators and
  returns an object of class \code{"bstraub"}. The method of
  \code{\link[base]{summary}} for such objects displays further
  details and the method of \code{\link[stats]{predict}} computes the
  credibility premiums.
}
\section{Estimation of a}{
  The Bühlmann-Straub unbiaised estimator (\code{heterogeneity =
    "unbiased"}) of the between contracts heterogeneity parameter is
  \deqn{%
    \hat{a} = c \left( \sum_{i = 1}^I w_{i\cdot} (X_{iw} - X_{ww})^2 -
    (I - 1)\hat{s}^2 \right),}{%
    a = c sum(w[i.] * (X[iw] - X[ww])^2 - (I - 1) * s^2),}
  where \eqn{c = w_{\cdot\cdot}/(w_{\cdot\cdot}^2 - \sum_{i = 1}^I
    w_{i\cdot}^2)}{c = w[..]/(w[..]^2 - sum(w[i.]^2))} and \eqn{I} is the
  number of contracts.

  The Bishel-Straub pseudo-estimator (\code{heterogeneity =
    "iterative"}) is obtained recursively as the solution of
  \deqn{%
    \hat{a} = \frac{1}{I - 1} \sum_{i=1}^I z_i (X_{iw} - X_{zw})^2.}{%
    a = 1/(I - 1) sum(z[i] * (X[iw] - X[zw])^2).}
  The fixed point algorithm is used with a relative error of \code{TOL}
  as stopping criteria.
}
\value{
  For \code{bstraub}, an object of \code{\link[base]{class}}
  \code{"bstraub"}.

  An object of class \code{"bstraub"} is a list with the following
  components:
  \item{model}{the name of the model used (\code{"Buhlmann"} or
    \code{"Buhlmann-Straub"});}
  \item{individual}{vector of contract weighted averages;}
  \item{collective}{collective premium estimator;}
  \item{weights}{vector of contracts total weights, as used in
    credibility factors;}
  \item{s2}{estimator of the within contract heterogeneity parameter;}
  \item{unbiased}{unbiased estimator of the between contract
    heterogeneity parameter;}
  \item{iterative}{iterative estimator of the between contract
    heterogeneity parameter.}

  For \code{predict.bstraub}, a vector of credibility premiums.
}
\references{
  Goulet, V. (1998), \emph{Principles and Application of Credibility
    Theory}, Journal of Actuarial Practice \bold{6}, 5--62.

  Goovaerts, M. J. and Kaas, R. and van Heerwaarden, A. E. and
  Bauwelinckx, T. (1990), \emph{Effective actuarial methods},
  North-Holland.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Sébastien Auclair and Louis-Philippe Pouliot
}
\examples{
data(hachemeister)

## Credibility premiums calculated with the iterative estimator
fit <- bstraub(hachemeister[, 2:13], hachemeister[, 14:25])
fit 				# print method
summary(fit)			# more details
predict(fit)			# credibility premiums
}
\keyword{models}
