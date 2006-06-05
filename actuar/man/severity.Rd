\name{severity}
\alias{severity}
\alias{severity.default}
\alias{severity.simpf}
\title{Display a Two-Dimension Version of a Matrix of Vectors}
\description{
  severity() displays all values of a matrix of vectors by
  \sQuote{unrolling} the object vertically or horizontally.   
}
\usage{
## default S3 method:
severity(x, byrow = TRUE, \dots)

\method{severity}{simpf}(x, byrow = TRUE, y.exclude = 0, \dots)
}
\arguments{
  \item{x}{matrix of vectors or a \code{\link{simpf}} object}
  \item{byrow}{boolean, whether to \sQuote{unroll} horizontally or vertically}
  \item{y.exclude}{integer, a number of years of observations at the tail of the
    portfolio to be displayed seperately}
  \item{\dots}{potential further arguments for future methods}

}
\details{
  \code{severity} gives the possibility to visualize all data contained
  in a matrix of elements of length greater than 1, that is,
  a matrix where each element can be a vector. The name \sQuote{severity}
  was chosen because of the actuarial context where the need of \sQuote{unrolling}
  a matrix appeared. The function was first designed to treat the matrix
  of data returned in a \code{link{simpf}} object, where each element is
  a vector of losses during a particular year for a particular
  contract. Because the \sQuote{unrolling} task is not specific to this
  context, a generic was created along with a method for class
  \code{"simpf"}.

  
}
\value{
  The object returned is a matrix where vector elements are concatenated
  by row if argument \code{byrow} is TRUE or by column in the opposite
  case, making every element of every vector visible. This operation can
  cause rows or columns to be of different lengths, so the matrix is
  filled with \code{NA}.
  
  If \code{y.exclude} is greater than 0, a second matrix is returned. It 
  displays the losses for the \dQuote{y.exclude} last years, with the same method as for
  the preceding years.  
  
  It is important to realize that the function brings a nice visualisation
  with detriment to precision. In an insurance context where we have a
  matrix of losses, for example, with vectors of claims for each
  contract and each year, the matrix returned by severity() will cause
  either years or contracts to be indistinguishable, depending on the
  value of \code{byrow}. 
}
\seealso{
  \code{\link{aggregate}} \code{\link{simpf}} \code{\link{matrix}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{

x <- list(c(1:3), c(1:8), c(1:4), c(1:3))
(mat <- matrix(x, 2, 2))             ## The matrix of vectors
severity(mat, byrow = FALSE)


## With a simpf object

modelfreq <- list(dist1 = "pois",
                  par1 = list(lambda = quote(Lambda * weights)),
                  dist2 = "gamma",
                  par2 = c(shape = 2, rate = 1))
modelsev<-list(dist1 = "lnorm",
               par1 = list(meanlog = quote(Theta), sdlog = 1),
               dist2 = "norm",
               par2 = c(mean = 5, sd = 1))
data(hachemeister)
weights <- hachemeister$weights/mean(hachemeister$weights)
x <- simpf(5, 12, modelfreq, modelsev, weights)           ## The simpf object

severity(x,                  ## A matrix of losses,
         byrow = TRUE,       ## 'unrolled' by row,
         y.exclude = 1)      ## with claims for the last year displayed separately.
                                          
}
\keyword{datagen}
\keyword{manip}

