\name{NEWS}
\title{actuar News}
\encoding{UTF-8}

\section{LATER NEWS}{
  This file covers NEWS for the 1.x series. News for actuar 2.0-0 and
  later can be found in file \file{NEWS.Rd}.
}

\section{CHANGES IN VERSION 1.2-2}{
  \subsection{BUG FIX}{
    \itemize{
      \item dpareto() did not handle the case x == 0 correctly.
    }
  }
}

\section{CHANGES IN VERSION 1.2-1}{
  \subsection{(MORE OR LESS) USER-VISIBLE CHANGES}{
    \itemize{
      \item The package now depends on R >= 3.3.0 since it uses
      chkDots() in a few methods that do not use the content of their
      '...' argument.
      \item ogive() lost its argument '...' as it was unused anyway.
      \item severity.portfolio() calls unroll() directly instead of
      relying on the default method to be identical to unroll().
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Deleted an unwanted debugging message ("local") printed by
      CTE() at every execution.
      \item predict.cm() and summary.cm() now treat the '...' argument
      as advertised in the help file.
      \item Fixed bad examples in a few probability law help files that
      returned unintended results such as Inf or NaN.
    }
  }
  \subsection{MAINTENANCE}{
    \itemize{
      \item C-level function log1pexp(...) used in a few places in lieu
      of log1p(exp(...)).
      \item Names of the internal utility macros defined in dpq.h
      changed from "R_<...>" to "ACT_<...>" to make it clearer that they
      are defined by the package (although they were copied from R
      sources).
    }
  }
}

\section{CHANGES IN VERSION 1.2-0}{
  \subsection{NEW FEATURE}{
    \itemize{
      \item In the computation of the CTE in the Normal Power case,
      numerical integration has been replaced by the explicit formula given
      in Castañer, A.; Claramunt, M.M.; Mármol, M. (2013). Tail value at
      risk. An analysis with the Normal-Power approximation. In Statistical
      and Soft Computing Approaches in Insurance Problems. Nova Science
      Publishers. ISBN 978-1-62618-506-7.
    }
  }
}

\section{CHANGES IN VERSION 1.1-10}{
  \subsection{BUG FIX}{
    \itemize{
      \item Results of 'cm' for hierarchical models would get
      incorrectly sorted when there were 10 nodes or more at a given
      level. Thanks to Dylan Wienke \email{dwienke2@gmail.com} for the
      catch.
    }
  }
}

\section{CHANGES IN VERSION 1.1-9}{
  \subsection{MAINTENANCE}{
    \itemize{
      \item Functions 'head' and 'tail' explicitly imported from package
      utils in NAMESPACE as per a new requirement of R 3.3.x.
    }
  }
}

\section{CHANGES IN VERSION 1.1-8}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Memory allocation problem at the C level in hierarc().
      Thanks to Prof. Ripley for identification of the problem and help
      solving it.
      \item Abusive use of abs() at the C level in a few places.
    }
  }
}

\section{CHANGES IN VERSION 1.1-7}{
  \subsection{BUG FIX}{
    \itemize{
      \item panjer() result was wrong for the "logarithmic" type of
      frequency distribution. Thanks to \email{mmclaramunt@ub.edu} for the
      catch.
    }
  }
}

\section{CHANGES IN VERSION 1.1-6}{
  \subsection{BUG FIX}{
    \itemize{
      \item Fixed a deprecated use of real().
    }
  }
}

\section{CHANGES IN VERSION 1.1-5}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Complete rewrite of coverage(); the function it creates no
      longer relies on ifelse() and, consequently, is much faster. The
      rewrite was motivated by a change in the way [dp]gamma() handle their
      arguments in R 2.15.1.
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item summary.ogive() no longer relies on length 'n' to be in the
      environment of a function created by approxfun(). Fix required by
      R >= 2.16.0.
    }
  }
}

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
      the lower tail. Thanks to Dan Murphy \email{chiefmurphy@gmail.com} for
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

\section{OLDER NEWS}{
  News for actuar 0.9.7 and earlier can be found in file \file{NEWS.0.Rd}.
}
