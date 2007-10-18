\name{bstraub}
\alias{bstraub}
\alias{predict.bstraub}
\alias{predict.bstraub.old}
\title{Buhlmann-Straub Credibility Model}
\description{
  \code{bstraub} computes structure parameters estimators in the
  Bühlmann-Straub credibility model and \code{predict.bstraub} computes
  the credibility premiums.
}
\usage{
bstraub(ratios, weights,
        method = c("unbiased", "iterative"),
        tol = sqrt(.Machine$double.eps), maxit = 100,
        echo = FALSE, old.format = TRUE)

\method{predict}{bstraub}(object, levels = NULL, newdata, \dots)
\method{predict}{bstraub.old}(object, \dots)
}
\arguments{
  \item{ratios}{a matrix of ratios (contracts in lines, years in columns).}
  \item{weights}{a matrix of weights corresponding to ratios.}
  \item{method}{estimator of the between contract heterogeneity
    parameter used in premium calculation; \code{"unbiased"} for the
    usual Bühlmann-Straub estimator, \code{"iterative"} for the
    Bischel-Straub estimator (see below).}
  \item{tol}{maximum relative error in the iterative procedure.}
  \item{maxit}{maximum number of iterations.}
  \item{echo}{logical; whether to echo iterative procedure or not.}
  \item{old.format}{logical; if \code{TRUE}, return results in the
    deprecated pre-0.9-4 format.}
  \item{object}{an object of class \code{"bstraub"}.}
  \item{levels, newdata, \dots}{unused arguments.}
}
\details{
  Direct usage of this function is deprecated. The function will not be
  exported in future versions of the package. Use
  \code{\link[actuar]{cm}} instead.
  
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
  returns an object of class \code{"bstraub"}. The methods of
  \code{\link[stats]{predict}} compute the credibility premiums.
}
\section{Estimation of a}{
  The Bühlmann-Straub unbiased estimator (\code{heterogeneity =
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
  The fixed point algorithm is used with a relative error of \code{tol}
  as stopping criteria.
}
\value{
  For \code{bstraub} with \code{old.format = TRUE} , an object of
  \code{\link[base]{class}} \code{"bstraub.old"}. This format is
  deprecated. An object of class \code{"bstraub.old"} is a list with the
  following components:
  \item{individual}{vector of contract weighted averages;}
  \item{collective}{collective premium estimator;}
  \item{weights}{vector of contracts total weights, as used in
    credibility factors;}
  \item{s2}{estimator of the within contract heterogeneity parameter;}
  \item{unbiased}{unbiased estimator of the between contract
    heterogeneity parameter;}
  \item{iterative}{iterative estimator of the between contract
    heterogeneity parameter;}
  \item{cred}{vector of credibility factors;}

  For \code{bstraub} with \code{old.format = FALSE}, an object of
  \code{\link[base]{class}} \code{"bstraub"}. An object of class
  \code{"bstraub"} is a list with the following components:
  \item{means}{a list containing the collective premium estimator and
    vector of contract weighted averages.}
  \item{weights}{a list containing the total portfolio weight and the
    vector of contracts total weights, as used in credibility factors;}
  \item{unbiased}{a vector containing the unbiased variance components
    estimators.}
  \item{iterative}{a vector containing the iterative variance components
    estimators, or \code{NULL}.}
  \item{cred}{vector of credibility factors.}
  \item{nodes}{a list containing the number of contracts in the
    portfolio.}
  
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
\seealso{
  \code{\link{cm}}
}
\examples{
data(hachemeister)

## Credibility premiums calculated with the iterative estimator
fit <- bstraub(hachemeister[, 2:13], hachemeister[, 14:25],
               old.format = FALSE)
fit 				# a list
predict(fit)			# credibility premiums
}
\keyword{models}
