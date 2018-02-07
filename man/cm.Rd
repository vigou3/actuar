\name{cm}
\alias{cm}
\alias{print.cm}
\alias{predict.cm}
\alias{summary.cm}
\alias{print.summary.cm}
\title{Credibility Models}
\description{
  Fit the following credibility models: \enc{Bühlmann}{Buhlmann},
  \enc{Bühlmann}{Buhlmann}-Straub, hierarchical, regression
  (Hachemeister) or linear Bayes.
}
\usage{
cm(formula, data, ratios, weights, subset,
   regformula = NULL, regdata, adj.intercept = FALSE,
   method = c("Buhlmann-Gisler", "Ohlsson", "iterative"),
   likelihood, ...,
   tol = sqrt(.Machine$double.eps), maxit = 100, echo = FALSE)

\method{print}{cm}(x, \dots)

\method{predict}{cm}(object, levels = NULL, newdata, \dots)

\method{summary}{cm}(object, levels = NULL, newdata, \dots)

\method{print}{summary.cm}(x, \dots)
}
\arguments{
  \item{formula}{character string \code{"bayes"} or an object of
    class \code{"\link[stats]{formula}"}: a symbolic description of the
    model to be fit. The details of model specification are given
    below.}
  \item{data}{a matrix or a data frame containing the portfolio
    structure, the ratios or claim amounts and their associated weights,
    if any.}
  \item{ratios}{expression indicating the columns of \code{data}
    containing the ratios or claim amounts.}
  \item{weights}{expression indicating the columns of \code{data}
    containing the weights associated with \code{ratios}.}
  \item{subset}{an optional logical expression indicating a subset of
    observations to be used in the modeling process. All observations
    are included by default.}
  \item{regformula}{an object of class \code{"\link[stats]{formula}"}:
    symbolic description of the regression component (see
    \code{\link[stats]{lm}} for details). No left hand side is needed
    in the formula; if present it is ignored. If \code{NULL}, no
    regression is done on the data.}
  \item{regdata}{an optional data frame, list or environment (or object
    coercible by \code{\link[base]{as.data.frame}} to a data frame)
    containing the variables in the regression model.}
  \item{adj.intercept}{if \code{TRUE}, the intercept of the regression
    model is located at the barycenter of the regressor instead of the
    origin.}
  \item{method}{estimation method for the variance components of the
    model; see details below.}
  \item{likelihood}{a character string giving the name of the likelihood
    function in one of the supported linear Bayes cases; see details
    below.}
  \item{tol}{tolerance level for the stopping criteria for iterative
    estimation method.}
  \item{maxit}{maximum number of iterations in iterative estimation
    method.}
  \item{echo}{logical; whether to echo the iterative procedure or not}
  \item{x, object}{an object of class \code{"cm"}}
  \item{levels}{character vector indicating the levels to predict or to
    include in the summary; if \code{NULL} all levels are included.}
  \item{newdata}{data frame containing the variables used to predict
    credibility regression models.}
  \item{\dots}{parameters of the prior distribution for \code{cm};
    additional attributes to attach to the result for the
    \code{predict} and \code{summary} methods; further arguments to
    \code{\link[base]{format}} for the \code{print.summary} method;
    unused for the \code{print} method.}
}
\details{
  \code{cm} is the unified front end for credibility models fitting. The
  function supports hierarchical models with any number of levels (with
  \enc{Bühlmann}{Buhlmann} and \enc{Bühlmann}{Buhlmann}-Straub models as
  special cases) and the regression model of Hachemeister. Usage of
  \code{cm} is similar to \code{\link[stats]{lm}} for these cases.
  \code{cm} can also fit linear Bayes models, in which case usage is
  much simplified; see the section on linear Bayes below.

  When not \code{"bayes"}, the \code{formula} argument symbolically
  describes the structure of the portfolio in the form \eqn{~ terms}.
  Each term is an interaction between risk factors contributing to the
  total variance of the portfolio data. Terms are separated by \code{+}
  operators and interactions within each term by \code{:}. For a
  portfolio divided first into sectors, then units and finally
  contracts, \code{formula} would be \code{~ sector + sector:unit +
  sector:unit:contract}, where \code{sector}, \code{unit} and
  \code{contract} are column names in \code{data}. In general, the
  formula should be of the form \code{~ a + a:b + a:b:c + a:b:c:d +
  ...}.

  If argument \code{regformula} is not \code{NULL}, the regression model
  of Hachemeister is fit to the data. The response is usually time. By
  default, the intercept of the model is located at time origin. If
  argument \code{adj.intercept} is \code{TRUE}, the intercept is moved
  to the (collective) barycenter of time, by orthogonalization of the
  design matrix. Note that the regression coefficients may be difficult
  to interpret in this case.

  Arguments \code{ratios}, \code{weights} and \code{subset} are used
  like arguments \code{select}, \code{select} and \code{subset},
  respectively, of function \code{\link[base]{subset}}.

  Data does not have to be sorted by level. Nodes with no data (complete
  lines of \code{NA} except for the portfolio structure) are allowed,
  with the restriction mentioned above.
}
\section{Hierarchical models}{
  The credibility premium at one level is a convex combination between
  the linearly sufficient statistic of a node and the credibility
  premium of the level above. (For the first level, the complement of
  credibility is given to the collective premium.) The linearly
  sufficient statistic of a node is the credibility weighted average of
  the data of the node, except at the last level, where natural weights
  are used. The credibility factor of node \eqn{i} is equal to
  \deqn{\frac{w_i}{w_i + a/b},}{w[i]/(w[i] + a/b),}
  where \eqn{w_i}{w[i]} is the weight of the node used in the linearly
  sufficient statistic, \eqn{a} is the average within node variance and
  \eqn{b} is the average between node variance.
}
\section{Regression models}{
  The credibility premium of node \eqn{i} is equal to
  \deqn{y^\prime b_i^a,}{y' ba[i],}
  where \eqn{y} is a matrix created from \code{newdata} and
  \eqn{b_i^a}{ba[i]} is the vector of credibility adjusted regression
  coefficients of node \eqn{i}. The latter is given by
  \deqn{b_i^a = Z_i b_i + (I - Z_I) m,}{
    ba[i] = Z[i] b[i] + (I - Z[i]) m,}
  where \eqn{b_i}{b[i]} is the vector of regression coefficients based
  on data of node \eqn{i} only, \eqn{m} is the vector of collective
  regression coefficients, \eqn{Z_i}{Z[i]} is the credibility matrix and
  \eqn{I} is the identity matrix. The credibility matrix of node \eqn{i}
  is equal to
  \deqn{A^{-1} (A + s^2 S_i),}{A^(-1) (A + s2 S[i]),}
  where \eqn{S_i}{S[i]} is the unscaled regression covariance matrix of
  the node, \eqn{s^2}{s2} is the average within node variance and
  \eqn{A} is the within node covariance matrix.

  If the intercept is positioned at the barycenter of time, matrices
  \eqn{S_i}{S[i]} and \eqn{A} (and hence \eqn{Z_i}{Z[i]}) are diagonal.
  This amounts to use \enc{Bühlmann}{Buhlmann}-Straub models for each
  regression coefficient.

  Argument \code{newdata} provides the \dQuote{future} value of the
  regressors for prediction purposes. It should be given as specified in
  \code{\link[stats]{predict.lm}}.
}
\section{Variance components estimation}{
  For hierarchical models, two sets of estimators of the variance
  components (other than the within node variance) are available:
  unbiased estimators and iterative estimators.

  Unbiased estimators are based on sums of squares of the form
  \deqn{B_i = \sum_j w_{ij} (X_{ij} - \bar{X}_i)^2 - (J - 1) a}{%
    B[i] =  sum(j; w[ij] (X[ij] - Xb[i])^2 - (J - 1) a)}%
  and constants of the form
  \deqn{c_i = w_i - \sum_j \frac{w_{ij}^2}{w_i},}{%
    c[i] = w[i] - sum(j; w[ij]^2)/w[i],}%
  where \eqn{X_{ij}}{X[ij]} is the linearly sufficient statistic of
  level \eqn{(ij)}; \eqn{\bar{X_{i}}}{Xb[i]} is the weighted average of
  the latter using weights \eqn{w_{ij}}{w[ij]}; \eqn{w_i = \sum_j
  w_{ij}}{w[i] = sum(j; w[ij])}; \eqn{J} is the effective number of
  nodes at level \eqn{(ij)}; \eqn{a} is the within variance of this
  level. Weights \eqn{w_{ij}}{w[ij]} are the natural weights at the
  lowest level, the sum of the natural weights the next level and the
  sum of the credibility factors for all upper levels.

  The \enc{Bühlmann}{Buhlmann}-Gisler estimators (\code{method =
    "Buhlmann-Gisler"}) are given by%
  \deqn{b = \frac{1}{I} \sum_i \max \left( \frac{B_i}{c_i}, 0
    \right),}{%
    b = mean(max(B[i]/c[i], 0)),}%
  that is the average of the per node variance estimators truncated at
  0.

  The Ohlsson estimators (\code{method = "Ohlsson"}) are given by
  \deqn{b = \frac{\sum_i B_i}{\sum_i c_i},}{%
    b = sum(i; B[i]) / sum(i; c[i]),}%
  that is the weighted average of the per node variance estimators
  without any truncation. Note that negative estimates will be truncated
  to zero for credibility factor calculations.

  In the \enc{Bühlmann}{Buhlmann}-Straub model, these estimators are
  equivalent.

  Iterative estimators \code{method = "iterative"} are pseudo-estimators
  of the form
  \deqn{b = \frac{1}{d} \sum_i w_i (X_i - \bar{X})^2,}{%
    b = sum(i; w[i] * (X[i] - Xb)^2)/d,}
  where \eqn{X_i}{X[i]} is the linearly sufficient statistic of one
  level, \eqn{\bar{X}}{Xb} is the linearly sufficient statistic of
  the level above and \eqn{d} is the effective number of nodes at one
  level minus the effective number of nodes of the level above. The
  Ohlsson estimators are used as starting values.

  For regression models, with the intercept at time origin, only
  iterative estimators are available. If \code{method} is different from
  \code{"iterative"}, a warning is issued. With the intercept at the
  barycenter of time, the choice of estimators is the same as in the
  \enc{Bühlmann}{Buhlmann}-Straub model.
}
\section{Linear Bayes}{
  When \code{formula} is \code{"bayes"}, the function computes pure
  Bayesian premiums for the following combinations of distributions
  where they are linear credibility premiums:
  \itemize{
    \item{\eqn{X|\Theta = \theta \sim \text{Poisson}(\theta)}{X|\Theta ~ Poisson(\Theta)} and
      \eqn{\Theta \sim \text{Gamma}(\alpha, \lambda)}{\Theta ~ Gamma(\alpha, \lambda)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Exponential}(\theta)}{X|\Theta ~ Exponential(\Theta)} and
      \eqn{\Theta \sim  \text{Gamma}(\alpha, \lambda)}{\Theta ~ Gamma(\alpha, \lambda)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Gamma}(\tau, \theta)}{X|\Theta ~ Gamma(\tau, \Theta)} and
      \eqn{\Theta \sim  \text{Gamma}(\alpha, \lambda)}{\Theta ~ Gamma(\alpha, \lambda)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Normal}(\theta, \sigma_2^2)}{X|\Theta ~ Normal(\Theta, \sigma_2^2)} and
      \eqn{\Theta \sim \text{Normal}(\mu, \sigma_1^2)}{\Theta ~ Normal(\mu, \sigma_1^2)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Bernoulli}(\theta)}{X|\Theta ~ Bernoulli(\Theta)} and
      \eqn{\Theta \sim \text{Beta}(a, b)}{\Theta ~ Beta(a, b)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Binomial}(\nu, \theta)}{X|\Theta ~ Binomial(\nu, \Theta)} and
      \eqn{\Theta \sim \text{Beta}(a, b)}{\Theta ~ Beta(a, b)};}
    \item{\eqn{X|\Theta = \theta \sim \text{Geometric}(\theta)}{X|\Theta = \theta ~ Geometric(\theta)} and
      \eqn{\Theta \sim \text{Beta}(a, b)}{\Theta ~ Beta(a, b)}.}
    \item{\eqn{X|\Theta = \theta \sim \text{Negative Binomial}(r, \theta)}{X|\Theta ~ Negative Binomial(r, \Theta)} and
      \eqn{\Theta \sim \text{Beta}(a, b)}{\Theta ~ Beta(a, b)}.}}

  The following combination is also supported:
  \eqn{X|\Theta = \theta \sim \text{Single Parameter Pareto}(\theta)}{X|\Theta ~ Single Parameter Pareto(\Theta)}
  and \eqn{\Theta \sim \text{Gamma}(\alpha, \lambda)}{\Theta ~
  Gamma(\alpha, \lambda)}. In this case, the Bayesian estimator not of
  the risk premium, but rather of parameter \eqn{\theta} is linear with
  a \dQuote{credibility} factor that is not restricted to \eqn{(0, 1)}.  
      
  Argument \code{likelihood} identifies the distribution of \eqn{X|\Theta
    = \theta} as one of
  \code{"poisson"},
  \code{"exponential"},
  \code{"gamma"},
  \code{"normal"},
  \code{"bernoulli"},
  \code{"binomial"},
  \code{"geometric"},
  \code{"negative binomial"} or
  \code{"pareto"}.

  The parameters of the distributions of \eqn{X|\Theta = \theta} (when
  needed) and \eqn{\Theta} are set in \code{\dots} using the argument
  names (and default values) of \code{\link[stats]{dgamma}},
  \code{\link[stats]{dnorm}}, \code{\link[stats]{dbeta}},
  \code{\link[stats]{dbinom}}, \code{\link[stats]{dnbinom}} or
  \code{dpareto1}, as appropriate. For the Gamma/Gamma case, use
  \code{shape.lik} for the shape parameter \eqn{\tau} of the Gamma
  likelihood. For the Normal/Normal case, use \code{sd.lik} for the
  standard error \eqn{\sigma_2} of the Normal likelihood.

  Data for the linear Bayes case may be a matrix or data frame as usual;
  an atomic vector to fit the model to a single contract; missing or
  \code{NULL} to fit the prior model. Arguments \code{ratios},
  \code{weights} and \code{subset} are ignored.  
}
\value{
  Function \code{cm} computes the structure parameters estimators of the
  model specified in \code{formula}. The value returned is an object of
  class \code{cm}.

  An object of class \code{"cm"} is a list with at least the following
  components:
  \item{means}{a list containing, for each level, the vector of linearly
    sufficient statistics.}
  \item{weights}{a list containing, for each level, the vector of total
    weights.}
  \item{unbiased}{a vector containing the unbiased variance components
    estimators, or \code{NULL}.}
  \item{iterative}{a vector containing the iterative variance components
    estimators, or \code{NULL}.}
  \item{cred}{for multi-level hierarchical models: a list containing,
    the vector of credibility factors for each level. For one-level
    models: an array or vector of credibility factors.}
  \item{nodes}{a list containing, for each level, the vector of the
    number of nodes in the level.}
  \item{classification}{the columns of \code{data} containing the
    portfolio classification structure.}
  \item{ordering}{a list containing, for each level, the affiliation of
    a node to the node of the level above.}

  Regression fits have in addition the following components:
  \item{adj.models}{a list containing, for each node, the credibility
    adjusted regression model as obtained with
    \code{\link[stats]{lm.fit}} or \code{\link[stats]{lm.wfit}}.}
  \item{transition}{if \code{adj.intercept} is \code{TRUE}, a transition
    matrix from the basis of the orthogonal design matrix to the basis
    of the original design matrix.}
  \item{terms}{the \code{\link[stats]{terms}} object used.}

  The method of \code{predict} for objects of class \code{"cm"} computes
  the credibility premiums for the nodes of every level included in
  argument \code{levels} (all by default). Result is a list the same
  length as \code{levels} or the number of levels in \code{formula}, or
  an atomic vector for one-level models.
}
\references{
  \enc{Bühlmann}{Buhlmann}, H. and Gisler, A. (2005), \emph{A Course in
    Credibility Theory and its Applications}, Springer.

  Belhadj, H., Goulet, V. and Ouellet, T. (2009), On parameter
  estimation in hierarchical credibility, \emph{Astin Bulletin} \bold{39}.

  Goulet, V. (1998), Principles and application of credibility
  theory, \emph{Journal of Actuarial Practice} \bold{6}, ISSN
  1064-6647.

  Goovaerts, M. J. and Hoogstad, W. J. (1987), \emph{Credibility Theory},
  Surveys of Actuarial Studies, No. 4, Nationale-Nederlanden N.V.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Xavier Milhaud, Tommy Ouellet, Louis-Philippe Pouliot
}
\seealso{
  \code{\link[base]{subset}}, \code{\link[stats]{formula}},
  \code{\link[stats]{lm}}, \code{\link[stats]{predict.lm}}.
}
\examples{
data(hachemeister)

## Buhlmann-Straub model
fit <- cm(~state, hachemeister,
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit				# print method
predict(fit)			# credibility premiums
summary(fit)			# more details

## Two-level hierarchical model. Notice that data does not have
## to be sorted by level
X <- data.frame(unit = c("A", "B", "A", "B", "B"), hachemeister)
fit <- cm(~unit + unit:state, X, ratio.1:ratio.12, weight.1:weight.12)
predict(fit)
predict(fit, levels = "unit")	# unit credibility premiums only
summary(fit)
summary(fit, levels = "unit")	# unit summaries only

## Regression model with intercept at time origin
fit <- cm(~state, hachemeister,
          regformula = ~time, regdata = data.frame(time = 12:1),
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit
predict(fit, newdata = data.frame(time = 0))
summary(fit, newdata = data.frame(time = 0))

## Same regression model, with intercept at barycenter of time
fit <- cm(~state, hachemeister, adj.intercept = TRUE,
          regformula = ~time, regdata = data.frame(time = 12:1),
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
fit
predict(fit, newdata = data.frame(time = 0))
summary(fit, newdata = data.frame(time = 0))

## Poisson/Gamma pure Bayesian model
fit <- cm("bayes", data = c(5, 3, 0, 1, 1),
          likelihood = "poisson", shape = 3, rate = 3)
fit
predict(fit)
summary(fit)

## Normal/Normal pure Bayesian model
cm("bayes", data = c(5, 3, 0, 1, 1),
   likelihood = "normal", sd.lik = 2,
   mean = 2, sd = 1)
}
\keyword{models}
