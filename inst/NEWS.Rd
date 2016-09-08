\name{NEWS}
\title{actuar News}
\encoding{UTF-8}

\section{CHANGES IN actuar VERSION 2.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{Support functions \code{<d,p,q,r>logarithmic} for the
	logarithmic (or log-series) distribution.}
      \item{Support functions \code{<d,p,q,r>ztpois} for the zero
	truncated Poisson distribution.}
      \item{Support functions \code{<d,p,q,r>ztnbinom} for the zero
	truncated negative binomial distribution.}
      \item{Support functions \code{<d,p,q,r>ztgeom} for the zero
	truncated geometric distribution.}
      \item{Support functions \code{<d,p,q,r>ztbinom} for the zero
	truncated binomial distribution.}
      \item{\code{aggregateDist} now accepts
        \code{"zero truncated binomial"},
	\code{"zero truncated geometric"},
	\code{"zero truncated negative binomial"} and
	\code{"zero truncated poisson"}
	for argument \code{model.freq} with the \code{"recursive"}
	method.}
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

\section{OLDER NEWS}{
  News for actuar 1.2-2 and earlier can be found in file
  \file{NEWS.1.Rd}.
}
