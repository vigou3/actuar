\name{NEWS}
\title{actuar News}
\encoding{UTF-8}

\section{CHANGES IN actuar VERSION 2.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{Support functions \code{[dpqr]logarithmic} and
	\code{[dpqr]zmlogarithmic} for the logarithmic (or log-series)
	and zero-modified logarithmic distributions.}
      \item{Support functions \code{[dpqr]ztpois} and
	\code{[dpqr]zmpois} for the zero-truncated and zero-modified
	Poisson distributions.}
      \item{Support functions \code{[dpqr]ztnbinom} and
	\code{[dpqr]zmnbinom} for the zero-truncated and zero-modified
	negative binomial distributions.}
      \item{Support functions \code{[dpqr]ztgeom} and
	\code{[dpqr]zmgeom} for the zero-truncated and zero-modified
	geometric distributions.}
      \item{Support functions \code{[dpqr]ztbinom} and
	\code{[dpqr]zmbinom} for the zero-truncated and zero-modified
	binomial distributions.}
      \item{\code{aggregateDist} now accepts
        \code{"zero-truncated binomial"},
	\code{"zero-truncated geometric"},
	\code{"zero-truncated negative binomial"},
	\code{"zero-truncated poisson"},
	\code{"zero-modified binomial"},
	\code{"zero-modified geometric"},
	\code{"zero-modified negative binomial"},
	\code{"zero-modified poisson"} and
	\code{"zero-modified logarithmic"}
	for argument \code{model.freq} with the \code{"recursive"}
	method.}
      \item{Extended range of admissible values for many limited
	expected value functions thanks to new C-level function
	\code{expint}, \code{betaint} and \code{gammaint}. These provide
	special integrals presented in the introduction of Appendix A of
	Klugman et al. (2012).

	Concerned functions are:
	\code{levtrbeta},
	\code{levgenpareto},
	\code{levburr},
	\code{levinvburr},
	\code{levpareto},
	\code{levinvpareto},
	\code{levllogis},
	\code{levparalogis},
	\code{levinvparalogis} in the Transformed Beta family, and
	\code{levinvtrgamma},
	\code{levinvgamma},
	\code{levinvweibull} in the Transformed Gamma family.}
      \item{Functions \code{expint}, \code{betaint} and \code{gammaint}
	to compute to special integrals mentioned above. These are
	merely convenience R interfaces to the C level functions. They
	are \emph{not} exported by the package.}
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item{In \code{aggregateDist} with the \code{"recursive"} method,
	if the length of \code{p0} is greater than one, only the first
	element is used, with a warning.}
      \item{\code{aggregateDist} with the \code{"recursive"} method and
	\code{model.freq = "logarithmic"} now uses the new
	\code{dlogarithmic} family of functions. Therefore,
	parametrization has changed from the one of Klugman et al.
	(2012) to the standard parametrization for the logarithmic
	distribution. Basically, any value of \code{prob} for the
	logarithmic parameter in previous versions of actuar should now
	be \code{1 - prob}.}
    }
  }
  \subsection{BUG FIX}{
    \itemize{
    }
  }
}

\note{
  News for actuar 1.2-2 and earlier can be found in file
  \file{NEWS.1.Rd}.
}
