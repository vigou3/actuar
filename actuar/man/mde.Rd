\name{mde}
\alias{Mde}
\alias{mde}
\title{Minimum Distance Estimation}
\description{
  Minimum distance estimators for individual and grouped data.
}
\usage{
  mde(x, fun, start, measure  = c("CvM", "chi-square", "LAS"), weights = NULL, ...)
}
\arguments{
  \item{x}{individual or grouped data.}
  \item{fun}{Either a character string or a function returning a
    cumulative distribution evaluated at its first argument.}
  \item{start}{A named list giving the parameters to be optimized with
    initial values }
  \item{measure}{Either \code{CvM} for Cramer von Mises' method,
    \code{chi-square} for modified chi-square's method, and \code{LAS}
    for layer average severity's method.}
  \item{weights}{weights. If arguments \code{weights} is missing,
    all weights are equal to 1 for \code{CvM} and \code{LAS}
    methods. For \code{chi-square}, weights are equal to \code{1 / nj}.}
  \item{\dots}{Additional parameters, either for \code{fun} or for \code{optim}. In particular, it can be used to specify bounds via \code{lower}
    or \code{upper} or both.  If arguments of \code{ecdf} are included they will be held fixed. }
}
\details{
  For \code{CvM}'s method direct optimization of the distance between
  theorical and empirical cumulative distribution function is performed
  using \code{optim}, for modified \code{chi-square}'s method
  direct optimization of the chi-square statistic is performed using
  \code{optim}, and for \code{LAS}'s direct optimization of distance
  between theorical and empirical limited distribution is performed
  using \code{optim}.
  
  For one-dimensional problems the Nelder-Mead method is used and for multi-dimensional problems the BFGS method, unless arguments named \code{lower} or \code{upper} are
  supplied when \code{L-BFGS-B} is used or \code{method} is supplied
  explicitly.
  
  There is \code{print} method for this class.
}
\value{
  An object od class \code{mde}, a list with tow components,
  
  estimate: the parameter estimates, and
  
  distance: the distance calculated.
  
  Invalid arguments will return warning.
}
\references{
  Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
x <- rpareto(10000, 8, 10)
d <- list(shape = 2, scale = 2)
mde(x, ppareto, d, "CvM")
t <- grouped(x = c(2, 4, 6, 10, 20), y = c(1, 4, 2, 1))
mde(t, pgamma, d, "chi-square")
}
\keyword{models}
