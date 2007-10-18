\name{hist.grouped.data}
\alias{hist.grouped.data}
\title{Histogram for grouped data}
\description{
  This method for the generic function \code{\link{hist}} is mainly
  useful to plot the histogram of grouped data. If \code{plot = FALSE},
  the resulting object of class \code{"histogram"} is returned for
  compatibility with \code{\link{hist.default}}, but does not contain
  much information not already in \code{x}.
}
\usage{
\method{hist}{grouped.data}(x, freq = NULL, probability = !freq,
     density = NULL, angle = 45, col = NULL, border = NULL,
     main = paste("Histogram of" , xname),
     xlim = range(cj), ylim = NULL, xlab = xname, ylab,
     axes = TRUE, plot = TRUE, labels = FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"}; only the first
    column of frequencies is used.}
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a
    representation of frequencies, the \code{counts} component of
    the result; if \code{FALSE}, probability densities, component
    \code{density}, are plotted (so that the histogram has a total area
    of one). Defaults to \code{TRUE} \emph{iff} group boundaries are
    equidistant (and \code{probability} is not specified).}
  \item{probability}{an \emph{alias} for \code{!freq}, for S compatibility.}
  \item{density}{the density of shading lines, in lines per inch.
    The default value of \code{NULL} means that no shading lines
    are drawn. Non-positive values of \code{density} also inhibit the
    drawing of shading lines.}
  \item{angle}{the slope of shading lines, given as an angle in
    degrees (counter-clockwise).}
  \item{col}{a colour to be used to fill the bars.
    The default of \code{NULL} yields unfilled bars.}
  \item{border}{the color of the border around the bars. The default
    is to use the standard foreground color.}
  \item{main, xlab, ylab}{these arguments to \code{title} have useful
    defaults here.}
  \item{xlim, ylim}{the range of x and y values with sensible defaults.
    Note that \code{xlim} is \emph{not} used to define the histogram
    (breaks), but only for plotting (when \code{plot = TRUE}).}
  \item{axes}{logical. If \code{TRUE} (default), axes are draw if the
    plot is drawn.}
  \item{plot}{logical. If \code{TRUE} (default), a histogram is
    plotted, otherwise a list of breaks and counts is returned.}
  \item{labels}{logical or character. Additionally draw labels on top
    of bars, if not \code{FALSE}; see \code{\link{plot.histogram}}.}
  \item{\dots}{further graphical parameters passed to
    \code{\link{plot.histogram}} and their to \code{\link{title}} and
    \code{\link{axis}} (if \code{plot=TRUE}).}
}
\value{
  An object of class \code{"histogram"} which is a list with components:
  \item{breaks}{the \eqn{r + 1} group boundaries.}
  \item{counts}{\eqn{r} integers; the frequency within each group.}
  \item{density}{the relative frequencies within each group
    \eqn{n_j/n}{n[j]/n}, where \eqn{n_j}{n[j]} = \code{counts[j]}.}
  \item{intensities}{same as \code{density}. Deprecated, but retained
    for compatibility.}
  \item{mids}{the \eqn{r} group midpoints.}
  \item{xname}{a character string with the actual \code{x} argument name.}
  \item{equidist}{logical, indicating if the distances between
    \code{breaks} are all the same.}
}
\note{
  The resulting value does \emph{not} depend on the values of
  the arguments \code{freq} (or \code{probability})
  or \code{plot}. This is intentionally different from S.
}
\seealso{
  \code{\link{hist}} and \code{\link{hist.default}} for histograms of
  individual data and fancy examples.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\examples{
data(gdental)
hist(gdental)
}
\keyword{dplot}
\keyword{hplot}
\keyword{distribution}
