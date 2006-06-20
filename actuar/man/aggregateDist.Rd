\name{aggregateDist}
\alias{aggregateDist}
\alias{print.aggregateDist}
\title{Generate an Aggregate Distribution of Claim Amounts }
\description{
  Uses one of five methods to generate the distribution of the
  aggregate claims of a portfolio over one year given a frequency and a
  severity model or the true moments of the distribution.  
}
\usage{
aggregateDist(method = c("normal", "np2", "simulation", "recursive", "exact"),
      model.freq, model.sev, moments, h, p0, TOL = 1e-06, ...)
\method{print}{aggregateDist}(x,\dots)
}

\arguments{
  \item{method}{ method to be used}
  \item{model.freq}{ named list containing the frequency model (see details
    below). Ignored when using \code{normal} and \code{np2} methods.}
  \item{model.sev}{ named list containing the severity model (see details
    below). Ignored when using \code{normal} and \code{np2} methods.}
  \item{moments}{ named vector containing the moments required by the
    method chosen. Ignored unless \code{normal} or \code{np2} method is called.}
  \item{h}{ discretization step used when discretizing the severity
    distribution or when evaluating a continuous distribution over
    discrete values, depending on the chosen method}
  \item{p0}{ arbitrary amount of probability at zero given to the
    frequency distribution. It creates zero-modified or
    zero-truncated distributions. Ignored unless \code{recursive}
    method is called.}
  \item{TOL}{ passed on to a quantile function to limit the range of
    evaluation in the \code{normal} and \code{np2} methods, the inverse of
    the number of obsevations in the \code{simulation} method, and the level of
    precision in the \code{recursive} method.}
  \item{x}{ an \code{aggregateDist} object} 
  \item{...}{ arguments potentially required by additionnal methods}   	  
}
\details{
  \code{aggregateDist} essentiallly calls functions creating \code{aggregateDist}
  objects. \code{\link{normal}}, \code{\link{np2}}, and \code{\link{exact}} methods call
  functions of the same name.  \code{simulation} and \code{recursive}
  methods respectively call functions \code{simS} and \code{panjer}.

  \code{normal} and \code{np2} methods require only a discretization
  step ( \code{h} ) and a \code{moments} argument, which is a named
  vector that contains the moments required by the respective function,
  that is the true mean and the true variance in both cases, and an
  additionnal \code{skewness} element for \code{np2}. Other arguments
  will be ignored. With respect to \code{h}, \code{aggregateDist} is intelligent enough to evaluate the called
  function in an appropriate range of values. Argument \code{moments}
  will be read only if one of those two methods are called. 

  The \code{simulation} method directly calls simS(), with
  a number of observations equal to the inverse of \code{TOL} and
  \code{model.freq} and \code{model.sev} arguments being
  integrally passed on. Only those three arguments will be read. 

  The \code{recursive} method discretizes de distribution of the
  severity model before calling panjer().

  The \code{exact} method also discretizes the distribution of the
  severity model and calls the much slower function exact().  
}
\value{
  An list of class \code{aggregateDist}, containing the following
  objects:

  \item{fs}{ a vector of the discretized density of the total amount of
    claims distribution}
  \item{Fs}{ a vector of the discretized cumulative distribution of the
    total amount of claims}
  \item{call}{ the function call}
  \item{FUN}{ an approximated function returned by \code{\link{approxfun}}
    to evaluate the CDF at any point.}
  \item{h}{ the discretization step given in argument. NULL if
    \code{aggregateDist} object is not created by \code{aggregateDist()}.}
  
}
\seealso{
  \code{\link{panjer}} \code{\link{exact}}
  \code{\link{simS}} \code{\link{np2}} \code{\link{normal}}
  }

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot}

\examples{
### 'normal' method
moments <- c(mean = 5, var = 2.25)
h <- 0.001
(res <- aggregateDist("normal", moments = moments, h = h))

### 'np2' method
moments <- c(mean = 5, var = 2.25, skewness = 0.25)
(res <- aggregateDist("np2", moments = moments, h = h))

### 'sim' method
h <- 0.025
model.freq <- list(dist = "pois", par = c(lambda = 4))
model.sev <- list(dist = "norm", par = c(mean = 10, sd = 1.5))
(res <- aggregateDist("simulation", model.freq = model.freq, model.sev = model.sev, h =
h, TOL = 1e-05))

### 'recursive' method
(res <- aggregateDist("recursive", model.freq = model.freq, model.sev =
model.sev, h = h, TOL = 1e-08 ))

### 'exact' method
(res <- aggregateDist("exact", model.freq = model.freq, model.sev = model.sev, h = h, TOL = 1e-08 ))
}
\keyword{distribution}
\keyword{models}