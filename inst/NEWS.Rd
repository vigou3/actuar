\name{NEWS}
\title{actuar News}
\encoding{UTF-8}

\section{CHANGES IN actuar VERSION 2.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{New support functions \code{[dpqrm,lev,mgf]invgaussian} for
	the inverse gaussian distribution. The first three functions are
	C (read: faster) implementations of functions of the same name
	in package \pkg{statmod}.}
      \item{New support functions \code{[dpqrm]gumbel} for the Gumbel
	extreme value distribution.}
      \item{Extended range of admissible values for many limited
	expected value functions thanks to new C-level functions
	\code{expint}, \code{betaint} and \code{gammaint}. These provide
	special integrals presented in the introduction of Appendix A of
	Klugman et al. (2012); see also	\code{vignette("distributions")}.

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
      \item{New functions \code{expint}, \code{betaint} and
	\code{gammaint} to compute to the special integrals mentioned
	above. These are merely convenience R interfaces to the C level
	functions. They are \emph{not} exported by the package.}
      \item{New support functions \code{[dpqr]logarithmic} and
	\code{[dpqr]zmlogarithmic} for the logarithmic (or log-series)
	and zero-modified logarithmic distributions.}
      \item{New support functions \code{[dpqr]ztpois} and
	\code{[dpqr]zmpois} for the zero-truncated and zero-modified
	Poisson distributions.}
      \item{New support functions \code{[dpqr]ztnbinom} and
	\code{[dpqr]zmnbinom} for the zero-truncated and zero-modified
	negative binomial distributions.}
      \item{New support functions \code{[dpqr]ztgeom} and
	\code{[dpqr]zmgeom} for the zero-truncated and zero-modified
	geometric distributions.}
      \item{New support functions \code{[dpqr]ztbinom} and
	\code{[dpqr]zmbinom} for the zero-truncated and zero-modified
	binomial distributions.}
      \item{New vignette \code{"distributions"} reviews in great detail
	the continuous and discrete distributions provided in the
	package.}
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
      \item{New function \code{rmixture} to generate random variates
	from discrete mixtures, that is from random variables with
	densities of the form \eqn{f(x) = p_1 f_1(x) + ... + p_n
	f_n(x)}.}
      \item{New function \code{rcompound} to generate random variates
	from (non hierarchical) compound models of the form \eqn{S = X_1
	+ \dots + X_N}. Designed to be quite faster than \code{simul}
	for non hierarchical models.}
    \item{New function code{rcomppois} that is a simplified version of
	\code{rcompound} for the very common compound Poisson case.}
      \item{\code{simul} now accepts a basic (named or not) vector for
	argument \code{nodes} when simulating from a non hierarchical
	compound model. But really, one should use \code{rcompound} for
	such cases.}
      \item{New alias \code{rcomphierarc} for \code{simul} that better
	fits within the usual naming scheme of random generation
	functions.}
      \item{Functions \code{grouped.data} and \code{ogive} now accept
	individual data in argument. The former will group the data
	using \code{hist} (therefore, all the algorithms to
	compute the number of breakpoints available  in \code{hist} are
	also available in \code{grouped.data}). \code{ogive} will first
	create a grouped data object and then compute the ogive.

	While the two functions are not guaranteed to be backward
	compatible (the number and position of the arguments have
	changed), standard calls should not be affected.}
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item{The  material on probability laws of vignette \code{"lossdist"}
	has been moved to the new vignette \code{"distributions"} (see
	the previous section).}
      \item{The first argument of the \code{mgffoo} functions has
	changed from \code{x} to \code{t}. This is a more common
	notation for moment generating functions.}
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
      \item{The aim of vignette \code{"simulation"} is changed from
	\dQuote{simulation of compound hierarchical models} to
	\dQuote{simulation of insurance data with \pkg{actuar}} as it
	also covers the new functions \code{rmixture} and
	\code{rcompound}.}
    }
  }
  \subsection{DEPRECATED}{
    \itemize{
      \item{Functions \code{[m,lev,mgf]invGauss} that complemented
	functions \code{[dpqr]invGauss} of package \pkg{SuppDists} are
	deprecated in favor of the new complete set of functions
	\code{[dpqrm,lev,mgf]invgauss}.}
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
