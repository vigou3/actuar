\name{bstraub}
\alias{bstraub}
\title{Buhlmann-Straub Model}
\description{
  Calculate credibility premiums in the Bühlmann-Straub model.
}
\usage{
bstraub(ratios, weights, heterogeneity = c("iterative", "unbiased"), TOL = 1e-06, echo = FALSE)
}

\arguments{
  \item{ratios}{matrix of ratios (contracts in lines, years in columns)}
  \item{weights}{matrix of weights corresponding to ratios}
  \item{heterogeneity}{estimator of the between contract heterogeneity
    parameter used in premium calculation; "iterative" for the
    Bischel-Straub estimator; "unbiased" for the usual
    Bühlmann-Straub estimator (see below)}
  \item{TOL}{maximum relative error in the iterative procedure}
  \item{echo}{boolean, whether to print iterative procedure or not}
}

\details{
  The credibility premium of contract \eqn{i} is given by
  \deqn{z_i X_{iw} + (1 - z_i)X_{zw}}{z[i]X[iw] + (1-z[i])X[zw]}
  where
  \deqn{z_{i} = \frac{w_{i\cdot}\hat{a}}{w_{i\cdot}\hat{a} + \hat{s}^2}}{%
        z[i]=w[i.]a'/(w[i]a'+ s'^2)}
  \eqn{X_{iw}}{X[iw]} is the weighted average of the ratios of contract
  \eqn{i}, \eqn{X_{zw}}{X[zw]} is the weighted average using credibility
  factors of the matrix of ratios and \eqn{w_{i\cdot}}{w[i.]} is
  the total weight of a contract. \eqn{\hat{s}^2}{s'^2} is the estimator
  of the within contract heterogeneity and  \eqn{\hat{a}}{a'} the
  estimator of the between contract heterogeneity.

  Missing data are represent by NA's in both the matrix of ratios and
  the matrix of weights. The function can also cope with complete lines
  of NA's in case a contract has no experience.
}
\section{Estimation of a}{
  The Bühlmann-Straub unbiaised estimator of the between contracts
  heterogeneity parameter \eqn{a} is
  \deqn{
  \hat{a} = \frac{w_{\cdot \cdot}}{w_{\cdot \cdot}^2 - \sum_i w_{i\cdot}^2}
  \left( \sum_{i=1}^N w_{i\cdot} (X_{iw} - X_{ww})^2 - (N-1)\hat{s}^2
  \right)}{%
        a'= w[..]/(w[..]^2-sum(w[i.]^2)) sum(w[i.]*(X[iw]-X[ww])^2-(N-1)*s'^2)},
  where N is the number of contracts.

  The Bishel-Straub pseudo-estimator (``iterative'') is obtained recursively
  as the solution of
  \deqn{
    \hat{a} = \frac{1}{N-1} \sum_{i=1}^N z_i (X_{iw} - X_{zw})^2}{a'= 1/(N-1) sum(z[i]*(X[iw]-X[zw])^2)}.
  See the function vignette for details about the notation.
}

\value{
  A list with the following components:
  \item{premiums}{vector of credibility premiums}
  \item{individual}{vector of contract weightedaverages}
  \item{collective}{collective premium estimator}
  \item{weights}{vector of contracts total weights, as used in
    credibility factors}
  \item{s2}{estimator of the within contract heterogeneity parameter}
  \item{unbiased}{unbiased estimator of the between contract
    heterogeneity parameter}
  \item{iterative}{iterative estimator of the between contract
    heterogeneity parameter}
}

\references{   
  Goulet, V. (1998), \emph{Principles and Application of Credibility
    Theory}, Journal of Actuarial Practice, Volume 6, ISSN 1064-6647.

  Goovaerts, M. J. and Kaas, R. and van Heerwaarden, A. E. and
  Bauwelinckx, T. (1990), \emph{Effective actuarial methods},
  North-Holland.
}

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Sébastien Auclair
}

\examples{
data(hachemeister)
### Buhlmann-Straub calculated with the iterative estimator
bstraub(hachemeister$claims, hachemeister$weights)

### Buhlmann-Straub calculated with the unbiased estimator
bstraub(hachemeister$claims, hachemeister$weights, heterogeneity="unbiased")
}


\keyword{file}

