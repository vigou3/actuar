\name{ghcm}
\alias{ghcm}
\alias{print.ghcm}
\title{Generalized Hierarchical Credibility Models}
\description{Fit a credibility model in the formulation of variance components 
as described in Dannenburg, Kaas and Goovaerts (1996). Models supported
are part of a generalized hierarchical credibility theory as introduced
in Dannenburg (1995).}
\usage{
ghcm(formula, data, years, weights, subset, TOL = 1E-6, echo = FALSE)
\method{print}{ghcm}(x,\dots)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}
  \item{data}{a data frame containing the ratios, their associated
    weights and information on the hierarchical relationship between
    levels of contracts (See details below).}
  \item{years}{expression, indicating the columns of \code{data} that represent
    the years in which ratios were observed.}
  \item{weights}{expression, indicating the columns of \code{data} that represent
    the associated weights. Must be the same length as \code{years}}
  \item{subset}{an optional expression indicating the subset of the rows of
    'data' that should be used in the fit. This can be a logical
    vector, or a numeric vector indicating which observation
    numbers are to be included, or a  character  vector of the
    row names to be included.  All observations are included by
    default.}
  \item{TOL}{maximum relative error in the iterative procedure} 
  \item{echo}{boolean, whether to echo iterative procedure or not}
  \item{x}{a \code{"ghcm"} object}
  \item{\dots}{potential further arguments required by generic}
}
\details{
  The \code{formula} argument describes the contribution of a
  hierarchical interaction between risk factors to the total variance of
  the insured's claim ratio random variable. Typically, one could want
  to fit a model in which a portfolio is divided into sectors, which are
  further divided into classes, each consisting of homogeneous insureds.
  A formula describing this hierarchical interaction would
  then be \eqn{~ sector + sector:class + sector:class:insured}, given
  that the terms of the formula are column names in \code{data}. At the
  moment, \code{ghcm} can only handle a hierarchical model with any given number of
  levels, for which the Buhlmann-Straub model is a special case where
  there is only one level. Consequently, the formula has to be
  equivalent to \eqn{~ a + a:b + a:b:c + a:b:c:d + ...}. As there can be many
  equivalents to one formula, more formally, we expect the
  \code{order} attribute of \code{terms(formula)} to be the sequence of
  the \code{n} first integers, where \code{n} is the number of levels.

  \code{data} is a data frame where each row describes a particular
  contract. Columns then describe the contract's belonging at each
  level, they give the ratios for each year and their weights if any.
}
\value{
  Returns an object of class \code{ghcm}.
  The function \code{predict} is used to obtain the credibility premiums
  for each contract, or more straightforwardly, the claim ratios of year
  \eqn{t+1} given \eqn{t} years of experience.

  An object of class \code{"ghcm"} is a list containing the
  following components:
  \item{param}{a vector containing the structure parameters}
  \item{weights}{a list containing, for each level, the vector of weights}
  \item{means}{a list containing the linearly sufficient statistics of
    each level}
  \item{cred}{a list containing, for each level, the vector of
    credibility factors}
  \item{call}{the matched call}
  \item{data}{the specified data frame, or a modified data frame
    resulting of the use of argument \code{subset}, with an added column
    representing the affiliation to the one global portfolio. This
    addition is useful to simmetrize computations.}
  \item{levs}{the hierarchically ordered names of the levels}
  \item{aff}{a list expressing the affiliation structure of the
    portfolio}
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
\examples{
## Build a 'data' object
data(hachemeister)
contract <- c("C1", "C2" ,"C3", "C4", "C5")
sector <- LETTERS[c(1,2,1,2,2)]
(x <- data.frame(sector, contract, hachemeister$claims,
hachemeister$weights))

ghcm(~sector + sector:contract, x, X1:X12, X1.1:X12.1)
}
\keyword{models}
