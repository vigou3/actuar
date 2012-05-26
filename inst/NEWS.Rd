\name{NEWS}
\title{actuar News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.1-4}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The function resulting from elev() for individual data is
      now faster for a large number of limits. (Thanks to Frank Zhan
      \email{FrankZhan@donegalgroup.com} for the catch and report.)
    }
  }
}

\section{CHANGES IN VERSION 1.1-3}{
  \subsection{BUG FIX}{
    \itemize{
      \item Resolved symbol clash at C level tickled by package
      GeneralizedHyperbolic on Solaris.
      \item Wrong result given by levinvGauss() because the upper tail
      of the normal distribution was used in the calculation instead of
      the lower tail. Thanks to Dan Murphy <chiefmurphy@gmail.com> for
      the heads up.
    }
  }
}

\section{CHANGES IN VERSION 1.1-2}{
  \subsection{BUG FIX}{
    \itemize{
      \item \code{discretize()} would return wrong results when argument
      \code{step} was omitted in favor of \code{by} \emph{and} the
      discretization method \code{unbiased} was used. (Thanks to
      Marie-Pier Côté \email{marie-pier.cote.11@ulaval.ca} for the
      catch.)
    }
  }
}

\section{CHANGES IN VERSION 1.1-1}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item CITATION file updated.
    }
  }

  \subsection{BUG FIX}{
    \itemize{
      \item \code{summary.cm()} could skip records in the output
      thinking they were duplicates.
    }
  }
}

\section{CHANGES IN VERSION 1.1-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New argument \code{convolve} in \code{aggregateDist()} to
      convolve the distribution obtained with the recursive method a
      number of times with itself. This is used for large portfolios
      where the expected number of claims is so large that recursions
      cannot start. Dividing the frequency parameter by \eqn{2^n} and
      convolving \eqn{n} times can solve the problem.

      \item New method of \code{diff()} for \code{"aggregateDist"}
      objects to return the probability mass function at the knots of
      the aggregate distribution. Valid (and defined) for
      \code{"recursive"}, \code{"exact"} and \code{"simulation"} methods
      only.

      \item Since the terminology Tail Value-at-Risk is often used
      instead of Conditional Tail Expectation, \code{TVaR()} is now an
      alias for \code{CTE()}.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Quantiles (and thus VaRs and CTEs) for
      \code{"aggregateDist"} objects where off by one knot of the
      distribution.

      \item \code{cm()} returned the internal classification codes
      instead of the original ones for hierarchical models. (Thanks to
      Zachary Martin for the heads up.)
    }
  }
}

\section{CHANGES IN VERSION 1.0-2}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Functions \code{m<foo>()} and \code{lev<foo>()} now return
      \code{Inf} instead of \code{NaN} for infinite moments. (Thanks to
      David Humke for the idea.)
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Non-ascii characters in one R source file prevented
      compilation of the package in a C locale (at least on OS X).

      \item For probability laws that have a strictly positive mode or a
      mode at zero depending on the value of one or more shape
      parameters, \code{d<foo>(0, ...)} did not handle correctly the
      case exactly at the boundary condition. (Thanks to Stephen L
      \email{bulls22eye@gmail.com} for the catch.)
    }
  }
}

\section{CHANGES IN VERSION 1.0-1}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{levinvpareto()} works for \code{order > -shape} and
      defaults to \code{order = 1}, like all other \code{lev<foo>()}
      functions.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Functions \code{d<foo>()} handle the case \code{x == 0}
      correctly.

      \item Functions \code{q<foo>()} return \code{NaN} instead of an
      error when argument \code{p} is outside \eqn{[0, 1]} (as in R).

      \item Functions \code{r<foo>()} for three parameter distributions
      (e.g. Burr) no longer wrongly display the \code{"NaNs produced"}
      warning message.

      \item The warning message \code{"NaNs produced"} was not (and
      could not be) translated.

      \item Function \code{levinvpareto()} computes limited moments for
      \code{order > -shape} using numerical integration.
    }
  }
}

\section{CHANGES IN VERSION 1.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Improved support for regression credibility models. There is
      now an option to make the computations with the intercept at the
      barycenter of time. This assures that the credibility adjusted
      regression line (or plane, or ...) lies between the individual and
      collective ones. In addition, contracts without data are now
      supported like in other credibility models.

      \item Argument \code{right} for \code{grouped.data()} to allow
      intervals closed on the right (default) or on the left.

      \item Method of \code{quantile()} for grouped data objects to
      compute the inverse of the ogive.
    }
  }

  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{cm()} no longer returns the values of the unbiased
      estimators when \code{method = "iterative"}.

      \item Specification of regression models in \code{cm()} has
      changed: one should now provide the regression model as a formula
      and the regressors in a separate matrix or data frame.

      \item Due to above change, \code{predict.cm()} now expects
      \code{newdata} to be a data frame as for
      \code{stats:::predict.lm()}.
    }
  }

  \subsection{DEFUNCT}{
    \itemize{
      \item Function \code{bstraub()} is no longer exported. Users are
      expected to use \code{cm()} as interface instead.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Functions \code{r<foo>()} are now more consistent in warning
      when \code{NA}s (specifically \code{NaN}s) are generated (as per
      the change in R 2.7.0).

      \item \code{frequency.portfolio()} was wrongly counting \code{NA}s.

      \item Domain of pdfs returned by \code{aggregateDist()} now
      restricted to \eqn{[0, 1]}.

      \item Quantiles are now computed correctly (and more efficiently)
      in 0 and 1 by \code{quantile.aggregateDist()}.

      \item \code{coverage()} no longer requires a cdf when it is not
      needed, namely when there is no deductible and no limit.
    }
  }
}

\section{CHANGES IN VERSION 0.9-7}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{plot} method for function objects returned by
      \code{ruin()}.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Calculation of the Bühlmann-Gisler and Ohlsson estimators
      was incorrect for hierarchical models with more than one level.

      \item Better display of first column for grouped data objects.

      \item Miscellaneous corrections to the vignettes.
    }
  }
}

\section{CHANGES IN VERSION 0.9-6}{
  \itemize{
    \item  Accented letters in comments removed to avoid compilation problems
    under MacOS X on CRAN (see thread starting at
    \url{https://stat.ethz.ch/pipermail/r-devel/2008-February/048391.html}).
  }
}

\section{CHANGES IN VERSION 0.9-5}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{simulation} vignette on usage of function
      \code{simul()}. Most of the material was previously in the
      \code{credibility} vignette.

      \item Examples of \code{ruin()} and \code{adjCoef()} added to the
      \code{risk} demo.
    }
  }

  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Following some negative comments on a function name VG had
      been using for years, function \code{simpf()} is renamed to
      \code{simul()} and the class of the output from \code{simpf} to
      \code{portfolio}.

      \item The components of the list returned by
      \code{severity.portfolio()} are renamed from \code{"first"} and
      \code{"last"} to \code{"main"} and \code{"split"}, respectively.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item \code{levinvgauss()} returned wrong results.

      \item Restructuring of the weights matrix in \code{simpf()} may
      fail with an incorrect number of columns.

      \item Fixed index entry of the credibility theory vignette.

      \item \code{adjCoef()} would only accept as argument \code{h} a
      function named \code{h}.

      \item \code{ruin()} built incorrect probability vector and
      intensity matrix for mixture of Erlangs.

      \item \code{CTE.aggregateDist()} sometimes gave values smaller
      than the VaR for recursive and simulation methods.
    }
  }
}

\section{CHANGES IN VERSION 0.9-4}{
  \itemize{
    \item Maintenance and new features release.
  }

  \subsection{NEW FEATURES -- LOSS DISTRIBUTIONS}{
    \itemize{
      \item Functions \code{mgffoo()} to compute the moment (or cumulant
      if \code{log = TRUE}) generating function of the following
      distributions: chi-square, exponential, gamma, inverse gaussian
      (from package \pkg{SuppDists}), inverse gamma, normal, uniform and
      phase-type (see below).

      \item Functions \code{mfoo()} to compute the raw moments of all
      the probability distributions supported in the package and the
      following of base R: chi-square, exponential, gamma, inverse
      gaussian (from package \pkg{SuppDists}), inverse gamma, normal,
      uniform.

      \item Functions \code{<d,p,mgf,m,r>phtype()} to compute the
      probability density function, cumulative distribution function,
      moment generating function, raw moments of, and to generate
      variates from, phase-type distributions.
    }
  }

  \subsection{NEW FEATURES -- RISK THEORY}{
    \itemize{
      \item Function \code{VaR()} with a method for objects of class
      \code{"aggregateDist"} to compute the Value at Risk of a
      distribution.

      \item Function \code{CTE()} with a method for objects of class
      \code{"aggregateDist"} to compute the Conditional Tail Expectation
      of a distribution.

      \item Function \code{adjCoef()} to compute the adjustment
      coefficient in ruin theory. If proportional or excess-of-loss
      reinsurance is included in the model, \code{adjCoef()} returns a
      function to compute the adjustment coefficient for given limits. A
      plot method is also included.

      \item Function \code{ruin()} returns a function to compute the
      infinite time probability of ruin for given initial surpluses in
      the Cramér-Lundberg and Sparre Andersen models. Most calculations
      are done using the cdf of phase-type distributions as per Asmussen
      and Rolski (1991).

      \item Calculations of the aggregate claim distribution using the
      recursive method much faster now that recursions are done in C.
    }
  }

  \subsection{NEW FEATURES -- CREDIBILITY THEORY}{
    \itemize{
      \item Modular rewrite of \code{cm()}: the function now calls
      internal functions to carry calculations for each supported
      credibility model. This is more efficient.

      \item Basic support for the regression model of Hachemeister in
      function \code{cm()}.

      \item For the hierarchical credibility model: support for the
      variance components estimators of Bühlmann and Gisler (2005) and
      Ohlsson (2005). Support remains for iterative pseudo-estimators.

      \item Calculations of iterative pseudo-estimators in hierarchical
      credibility are much faster now that they are done in C.
    }
  }

  \subsection{OTHER NEW FEATURES}{
    \itemize{
      \item Four new vignettes: introduction to the package and presentation
      of the features in loss distributions, risk theory and credibility
      theory.

      \item Portfolio simulation material of the \code{credibility} demo
      moved to demo \code{simulation}.
    }
  }

  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Argument \code{approx.lin} of
      \code{quantile.aggregateDist()} renamed \code{smooth}.

      \item Function \code{aggregateDist()} gains a \code{maxit}
      argument for the maximum number of recursions when using Panjer's
      algorithm. This is to avoid infinite recursion when the cumulative
      distribution function does not converge to 1.

      \item Function \code{cm()} gains a \code{maxit} argument for the
      maximum number of iterations in pseudo-estimators calculations.

      \item Methods of \code{aggregate()}, \code{frequency()},
      \code{severity()} and \code{weights()} for objects of class
      \code{"simpf"} gain two new arguments:

      \enumerate{
	\item \code{classification}; when \code{TRUE}, the columns
	giving the classification structure of the portfolio are
	excluded from the result. This eases calculation of loss ratios
	(aggregate claim amounts divided by the weights);

	\item \code{prefix}; specifies a prefix to use in column names,
        with sensible defaults to avoid name clashes for data and weight
        columns.
      }
    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item The way weights had to be specified for the
      \code{"chi-square"} method of \code{mde()} to give expected
      results was very unintuitive. The fix has no effect when using the
      default weights.

      \item The empirical step function returned by the
      \code{"recursive"} and \code{"convolution"} methods of
      \code{aggregateDist()} now correctly returns 1 when evaluated past
      its largest knot.
    }
  }

  \subsection{DEPRECATED}{
    \itemize{
      \item Direct usage of \code{bstraub()} is now deprecated in favor
      of \code{cm()}. The function will remain in the package since it
      is used internally by \code{cm()}, but it will not be exported in
      future releases of the package. The current format of the results
      is also deprecated.
    }
  }
}

\section{CHANGES IN VERSION 0.9-3}{
  \subsection{DEPRECATED, DEFUNCT OR NO BACKWARD COMPATIBILITY}{
    \itemize{
      \item The user interface of \code{coverage()} has changed. Instead
      of taking in argument the name of a probability law (say
      \code{foo}) and require that functions \code{dfoo()} and
      \code{pfoo()} exist, \code{coverage()} now requires a function
      name or function object to compute the cdf of the unmodified
      random variable and a function name or function object to compute
      the pdf. If both functions are provided, \code{coverage()} returns
      a function to compute the pdf of the modified random variable; if
      only the cdf is provided, \code{coverage()} returns the cdf of the
      modified random variable. Hence, argument \code{cdf} is no longer
      a boolean. The new interface is more in line with other functions
      of the package.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Methods of \code{summary()} and \code{print.summary()} for
      objects of class \code{"cm"} were not declared in the NAMESPACE
      file.

      \item Various fixes to the demo files.
    }
  }
}

\section{CHANGES IN VERSION 0.9-2}{
  \itemize{
    Major official update. This version is not backward compatible with
    the 0.1-x series. Features of the package can be split in the following
    categories: loss distributions modeling, risk theory, credibility
    theory.
  }

  \subsection{NEW FEATURES -- LOSS DISTRIBUTIONS}{
    \itemize{
      \item Functions \code{<d,p,q,r>foo()} to compute the density
      function, cumulative distribution function, quantile function of,
      and to generate variates from, all probability distributions of
      Appendix A of Klugman et al. (2004), \emph{Loss Models, Second
      Edition} (except the inverse gaussian and log-t) not already in R.
      Namely, this adds the following distributions (the root is what
      follows the \code{d}, \code{p}, \code{q} or \code{r} in function
      names):

      \tabular{ll}{
       DISTRIBUTION NAME         \tab ROOT \cr
       Burr                      \tab \code{burr} \cr
       Generalized beta          \tab \code{genbeta} \cr
       Generalized Pareto        \tab \code{genpareto} \cr
       Inverse Burr              \tab \code{invburr} \cr
       Inverse exponential       \tab \code{invexp} \cr
       Inverse gamma	         \tab \code{invgamma} \cr
       Inverse Pareto            \tab \code{invpareto} \cr
       Inverse paralogistic      \tab \code{invparalogis} \cr
       Inverse transformed gamma \tab \code{invtrgamma} \cr
       Inverse Weibull           \tab \code{invweibull} \cr
       Loggamma                  \tab \code{loggamma} \cr
       Loglogistic               \tab \code{llogis} \cr
       Paralogistic              \tab \code{paralogis} \cr
       Pareto                    \tab \code{pareto} \cr
       Single parameter Pareto   \tab \code{pareto1} \cr
       Transformed beta          \tab \code{trbeta} \cr
       Transformed gamma         \tab \code{trgamma}
     }

     All functions are coded in C for efficiency purposes and should
     behave exactly like the functions in base R. For all distributions
     that have a scale parameter, the corresponding functions have
     \code{rate = 1} and \code{scale = 1/rate} arguments.

     \item Functions \code{<m,lev>foo()} to compute the \eqn{k}-th raw
     (non-central) moment and \eqn{k}-th limited moment for all the
     probability distributions mentioned above, plus the following ones
     of base R: beta, exponential, gamma, lognormal and Weibull.

     \item Facilities to store and manipulate grouped data (stored in an
     interval-frequency fashion). Function \code{grouped.data()} creates
     a grouped data object similar to a data frame. Methods of
     \code{"["}, \code{"[<-"}, \code{mean()} and \code{hist()} created
     for objects of class \code{"grouped.data"}.

     \item Function \code{ogive()} --- with appropriate methods of
     \code{knots()}, \code{plot()}, \code{print()} and \code{summary()}
     --- to compute the ogive of grouped data. Usage is in every respect
     similar to \code{stats:::ecdf()}.

     \item Function \code{elev()} to compute the empirical limited
     expected value of a sample of individual or grouped data.

     \item Function emm() to compute the k-th empirical raw (non-central)
     moment of a sample of individual or grouped data.

     \item Function \code{mde()} to compute minimum distance estimators
     from a sample of individual or grouped data using one of three
     distance measures: Cramer-von Mises (CvM), chi-square, layer
     average severity (LAS). Usage is similar to \code{fitdistr()} of
     package \pkg{MASS}.

     \item Function \code{coverage()} to obtain the pdf or cdf of the
     payment per payment or payment per loss random variable under any
     combination of the following coverage modifications: ordinary of
     franchise deductible, policy limit, coinsurance, inflation. The
     result is a function that can be used in fitting models to data
     subject to such coverage modifications.

     \item Individual dental claims data set \code{dental} and grouped dental
     claims data set \code{gdental} of Klugman et al. (2004), \emph{Loss Models,
     Second Edition}.
    }
  }

  \subsection{NEW FEATURES -- RISK THEORY}{
    \itemize{
      \item Function \code{aggregateDist()} returns a function to
      compute the cumulative distribution function of the total amount
      of claims random variable for an insurance portfolio using any of
      the following five methods:

      \enumerate{
	\item exact calculation by convolutions (using function convolve() of
	package \pkg{stats};
	\item recursive calculation using Panjer's algorithm;
	\item normal approximation;
	\item normal power approximation;
	\item simulation.
      }

      The modular conception of \code{aggregateDist()} allows for easy
      inclusion of additional methods. There are special methods of
      \code{print()}, \code{summary()}, \code{quantile()} and
      \code{mean()} for objects of class \code{"aggregateDist"}. The objects
      otherwise inherit from classes \code{"ecdf"} (for methods 1, 2 and 3) and
      \code{"function"}.

      See also the DEPRECATED, DEFUNCT OR NO BACKWARD COMPATIBILITY
      section below.

      \item Function \code{discretize()} to discretize a continuous
      distribution using any of the following four methods:

      \enumerate{
	\item upper discretization, where the discretized cdf is always above
	the true cdf;
	\item lower discretization, where the discretized cdf is always under
	the true cdf;
	\item rounding, where the true cdf passes through the midpoints of the
	intervals of the discretized cdf;
	\item first moment matching of the discretized and true
	distributions.
      }

      Usage is similar to \code{curve()} of package \pkg{graphics}.
      Again, the modular conception allows for easy inclusion of
      additional discretization methods.
    }
  }

  \subsection{NEW FEATURES -- CREDIBILITY THEORY}{
    \itemize{

      \item Function \code{simpf()} can now simulate data for
      hierarchical portfolios of any number of levels. Model
      specification changed completely; see the DEPRECATED, DEFUNCT OR
      NO BACKWARD COMPATIBILITY below. The function is also
      significantly (\eqn{\sim 10\times}{~10x}) faster than the previous
      version.

      \item Generic function \code{severity()} defined mostly to provide
      a method for objects of class \code{"simpf"}; see below.

      \item Methods of \code{aggregate()}, \code{frequency()},
      \code{severity()} and \code{weights()} to extract information from
      objects of class \code{"simpf"}:

      \enumerate{
	\item \code{aggregate()} returns the matrix of aggregate claim
	amounts per node;
	\item \code{frequency()} returns the matrix of the number of
	claims per node;
	\item \code{severity()} returns the matrix of individual claim
	amounts per node;
	\item \code{weights()} returns the matrix of weights
	corresponding to the data.
      }

      Summaries can be done in various ways; see \code{?simpf.summaries}

      \item Function \code{cm()} (for \emph{c}redibility \emph{m}odel)
      to compute structure parameters estimators for hierarchical
      credibility models, including the Bühlmann and Bühlmann-Straub
      models. Usage is similar to \code{lm()} of package \pkg{stats} in
      that the hierarchical structure is specified by means of a formula
      object and data is extracted from a matrix or data frame. There
      are special methods of \code{print()}, \code{summary()} for
      objects of class \code{"cm"}. Credibility premiums are computed
      using a method of \code{predict()}; see below.

      For simple Bühlmann and Bühlmann-Straub models, \code{bstraub()}
      remains simpler to use and faster.

      \item Function \code{bstraub()} now returns an object of class
      \code{"bstraub"} for which there exist print and summary
      methods. The function no longer computes the credibility
      premiums; see the DEPRECATED, DEFUNCT OR NO BACKWARD
      COMPATIBILITY section below.

      \item Methods of \code{predict()} for objects of class \code{"cm"}
      and \code{"bstraub"} created to actually compute the credibility
      premiums of credibility models. Function \code{predict.cm()} can
      return the premiums for specific levels of a hierarchical
      portfolio only.
    }
  }

  \subsection{OTHER NEW FEATURES}{
    \itemize{
      \item Function \code{unroll()} to unlist a list with a
      \code{"dim"} attribute of length 0, 1 or 2 (that is, a vector or
      matrix of vectors) according to a specific dimension. Currently
      identical to \code{severity.default()} by lack of a better usage
      of the default method of \code{severity()}.

      \item Three new demos corresponding to the three main fields of
      actuarial science covered by the package.

      \item French translations of the error and warning messages.

      \item The package now has a name space.
    }
  }

  \subsection{DEPRECATED, DEFUNCT OR NO BACKWARD COMPATIBILITY}{
    \itemize{
      \item Function \code{panjer()}, although still present in the
      package, should no longer be used directly. Recursive calculation
      of the aggregate claim amount should be done with
      \code{aggregateDist()}. Further, the function is not backward
      compatible: model specification has changed, discretization of the
      claim amount distribution should now be done with
      \code{discretize()}, and the function now returns a function to
      compute the cdf instead of a simple vector of probabilities.

      \item Model specification for \code{simpf()} changed completely
      and is not backward compatible with previous versions of the
      package. The new scheme allows for much more general models.

      \item Function \code{rearrangepf()} is defunct and has been
      replaced by methods of \code{aggregate()}, \code{frequency()} and
      \code{severity()}.

      \item Function \code{bstraub()} no longer computes the credibility
      premiums. One should now instead use \code{predict()} for this.

      \item The data set \code{hachemeister} is no longer a list but rather a
      matrix with a state specification.
    }
  }
}


\section{CHANGES IN VERSION 0.1-3}{
  \itemize{
    \item Fixed the dependency on R >= 2.1.0 since the package uses
    function \code{isTRUE()}.
  }
}

\section{CHANGES IN VERSION 0.1-2}{
  \itemize{
    \item First public release.

    \item Fixed an important bug in \code{bstraub()}: when calculating
    the range of the weights matrix, \code{NA}s were not excluded.

    \item Miscellaneous documentation corrections.
  }
}


\section{CHANGES IN VERSION 0.1-1}{
  \itemize{
    \item Initial release.

    \item Contains functions \code{bstraub()}, \code{simpf()},
    \code{rearrangepf()} and \code{panjer()}, and the dataset
    \code{hachemeister}.
  }
}
