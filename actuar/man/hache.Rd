\name{hache}
\alias{hache}
\title{Hachemeister Regression Model (1975)}
\description{
  The Hachemeister model extends the Bühlmann-Straub model by
  introducing a regression technique.
}
\usage{
hache(X, Y, weights, tol = 1E-6, echo = FALSE)
}  
\arguments{
  \item{X}{a matrix of ratios (contracts in lines, years in columns);}
  \item{Y}{a chosen design matrix;}
  \item{weights}{a matrix of weights corresponding to ratios;}
  \item{tol}{maximum relative error in the iterative procedure;}
  \item{echo}{logical; whether to echo iterative procedure or not.}
}
\details{
  The function is used to estimate the credibility premiums. More soon...
}
\value{
    For \code{hache}, an object of \code{\link[base]{class}}
  \code{"hache"}.

  An object of class \code{"hache"} is a list with the following
  components:
  \item{beta}{individual estimators of the intercept and the slope;}
  \item{betaAdj}{credibility adjusted individual estimators;}
  \item{betaTotal}{collective estimators of the intercept and the slope;}
  \item{cred}{an array containing  the credibility factors for each contract;}
  \item{s2}{estimator of the within contract heterogeneity parameter;}
  \item{a}{estimator of the between contract heterogeneity parameter;}
  \item{p}{a matrix containing the credibility premiums.}
}
\seealso{
  \code{\link{bstraub}}, \code{\link{cm}}
}
\author{
  Tommy Ouellet and Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(hachemeister)

X <- hachemeister[, 2:13]
W <- hachemeister[, 14:25]
Y <- matrix(c(rep(1,12),12:1), 12)

hache(X, Y, weights = W, tol = 0.000001)
}
\keyword{models}