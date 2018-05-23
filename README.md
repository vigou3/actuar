# actuar
[![Travis-CI Build Status](https://travis-ci.org/vigou3/actuar.svg?branch=master)](https://travis-ci.org/vigou3/actuar) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/actuar)](https://cran.r-project.org/package=actuar) ![downloads](http://cranlogs.r-pkg.org/badges/grand-total/actuar)

**actuar** is a package providing additional actuarial science
functionality to the [R](https://r-project.org) statistical system.
The project was officially launched in 2005 and is under active
development.

The current feature set of the package can be split into five main
categories: 

1. Additional probability distributions to model insurance loss
   amounts and loss frequency (19 continuous heavy tailed
   distributions, see the list [below](#distributions); the
   Poisson-inverse Gaussian discrete distribution; zero-truncated and
   zero-modified extensions of the standard discrete distributions;
   phase-type distributions);
2. Loss distributions modeling (extensive support for grouped data;
   empirical raw and limited moments; minimum distance estimation);
3. Risk and ruin theory (discretization of the claim amount
   distribution and computation of the aggregate claim amount
   distribution; computation of the adjustment coefficient and ruin
   probabilities); 
4. Simulation of discrete mixtures, compound models and compound
   hierarchical models;
5. Credibility theory (Bühlmann, Bühlmann-Straub, hierarchical,
   regression and linear Bayes models).

The package includes extensive documentation in the form of package
*vignettes*. Each vignette focuses on a feature set of the package. To
get the list of available vignettes, enter at the R command prompt:

```R
vignette(package = "actuar")
```

## Installation

You can install the stable version of the package from the 
[Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/package=actuar):

```R
install.packages("actuar")
```

## Citation

To cite package **actuar** in publications see the output of

```R
citation(package = "actuar")
```

## License

**actuar** is free software licensed under the [GNU General Public
License (GPL)](https://www.gnu.org/copyleft/gpl.html).

## Philosophy

As much as possible, the developers have tried to keep the user
interface of the various functions of the package consistent.
Moreover, the package follows the general R philosophy of working with
model objects. This means that instead of merely returning, say, a
vector of probabilities, many functions will return an object
containing, among other things, the said probabilities. The object can
then be manipulated at one's will using various extraction, summary or
plotting functions.

## <a name="distributions"></a> Additional continuous distributions

**actuar** provides support functions for all the probability
distributions found in Appendix&nbsp;A of 
[*Loss Models: From Data to Decisions*, 4th Edition](https://www.wiley.com/en-us/Loss+Models%3A+From+Data+to+Decisions%2C+4th+Edition-p-9781118411650)
and not already present in base R, excluding the log-t, but
including the loggamma distribution. These distributions mostly fall
under the umbrella of extreme value or heavy tailed distributions.

The list of distributions supported by **actuar** is as follows, using
the nomenclature of *Loss Models*.

###  Transformed beta family

- Transformed beta
- Burr
- Loglogistic
- Paralogistic
- Generalized Pareto
- Pareto
- Inverse Burr
- Inverse Pareto
- Inverse paralogistic

### Transformed gamma family

- Transformed gamma
- Inverse transformed gamma
- Inverse gamma
- Inverse Weibull
- Inverse exponential

###  Other

- Loggamma
- Gumbel
- Inverse Gaussian
- Single parameter Pareto
- Generalized beta
