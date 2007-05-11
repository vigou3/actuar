### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the credibility theory facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)


###
### PORTFOLIO SIMULATION
###

### Function 'simpf' can be used to simulate portfolios of data for
### one or multi level hierarchical models. It is possible to simulate
### claim frequencies and claim amounts separately. Weights can be
### used in the models to allow for varying exposure. Examples below
### range from the simple to the more elaborate.
###
### The result of 'simpf' is a two dimension list. Four functions are
### available to extract useful information from this type of object:
### "simpf" class methods of functions 'aggregate', 'frequency',
### 'severity' and 'weights'.

## A simple Compound Poisson model: S_t = C_1 + ... + C_{N_t}, with
##
##   N_t ~ Poisson(10)
##     C ~ Lognormal(log(1500) - 1, 1)
##
## for t = 1, ..., 10. The name of the components below serves no
## purpose. 'severity' can extract individual claim amounts separately
## for columns specified in argument 'splitcol'.
pf <- simpf(list(y = 10), model.freq = expression(y = rpois(10)),
            model.sev = expression(y = rlnorm(log(1500) - 1, 1)))
pf                                      # print method
aggregate(pf)                           # aggregate claim amounts
frequency(pf)                           # frequencies
severity(pf)                            # individual claim amounts
severity(pf, splitcol = 10)             # last period separate


## Simple (continuous) mixture of models:
##
##   S_t|Theta ~ Poisson(Theta)
##       Theta ~ Gamma(2, 1).
##
## Any names can be used in the model.
pf <- simpf(list(Theta = 1, S = 10),
            model.freq = expression(Theta = rgamma(2, 1), S = rpois(Theta)))
aggregate(pf)                           # actual data
frequency(pf)                           # same, here

## Model with with mixtures for both frequency and severity: S_{it} =
## C_{i1} + ... + C_{it N_{it}}, with
##
##   N_{it}|Lambda_i ~ Poisson(Lambda_i)
##          Lambda_i ~ Gamma(2, 1)
##   C_{itu}|Theta_i ~ Lognormal(Theta_i, 1)
##           Theta_i ~ N(5, 1)
##
## This is considered a two-level hierarchical model.
pf <- simpf(list(entity = 10, year = 5),
            model.freq = expression(entity = rgamma(2, 1),
                year = rpois(entity)),
            model.sev = expression(entity = rnorm(5, 1),
                year = rlnorm(entity, 1)))
pf
aggregate(pf)
frequency(pf)

## Same model as above, but with weights incorporated into the model:
##
##   N_{it}|Lambda_i ~ Poisson(w_{it} * Lambda_i)
##          Lambda_i ~ Gamma(2, 1)
##   C_{itu}|Theta_i ~ Lognormal(Theta_i, 1)
##           Theta_i ~ N(5, 1)
##
## The string "weights" should appear in the model specification
## wherever weights are to be used. Argument 'weights' of 'simpf'
## should be a vector listing the weights in lexicographic order, that
## is all weights of entity 1, then all weights of entity 2, and so
## on.
wit <- runif(10, 2, 10)
( wit <- runif(50, rep(0.5 * wit, each = 5), rep(1.5 * wit, each = 5)) )
( pf <- simpf(list(entity = 10, year = 5),
              model.freq = expression(entity = rgamma(2, 1),
                  year = rpois(weights * entity)),
              model.sev = expression(entity = rnorm(5, 1),
                  year = rlnorm(entity, 1)),
              weights = wit) )
weights(pf)                             # extraction of weights

## Four-level hierarchical model (sector, unit, contract, years
## [data]). Claim severity varies only by sector and unit:
##
##  N_{ijkt}|Lambda_k ~ Poisson(w_{ijkt} * Lambda_k)
##     Lambda_k|Phi_j ~ Gamma(Phi_j, 1)
##        Phi_j|Psi_i ~ Gamma(Psi_i, 0.1)
##              Psi_i ~ Exponential(2)
##  C_{ijktu}|Theta_j ~ Lognormal(Theta_j, 1)
##       Theta_j|Nu_i ~ Normal(Nu_i, 1)
##               Nu_i ~ Normal(2, 0.1)
##
## The number of "nodes" at each level is different:
##
##     Sectors: 2
##       Units: 3 in sector 1
##              4 in sector 2
##   Contracts: 10 in unit 1 of sector 1
##               5 in unit 2 of sector 1
##               8 in unit 3 of sector 1
##               5 in unit 1 of sector 2
##               7 in unit 2 of sector 2
##              11 in unit 3 of sector 2
##               4 in unit 4 of sector 2
##       Years: 6 everywhere
##
## Hence, a total of (10 + 5 + 8 + 5 + 7 + 11 + 4) * 6 = 300 nodes at
## the last level.
wijkt <- runif(50, 2, 10)
wijkt <- runif(300, rep(0.5 * wijkt, each = 6), rep(1.5 * wijkt, each = 6))
nodes <- list(sector = 2, unit = c(3, 4),
              contract = c(10, 5, 8, 5, 7, 11, 4), year = 6)
mf <- expression(sector = rexp(2),
                 unit = rgamma(sector, 0.1),
                 contract = rgamma(unit, 1),
                 year = rpois(weights * contract))
ms <- expression(sector = rnorm(2, sqrt(0.1)),
                 unit = rnorm(sector, 1),
                 contract = NULL,
                 year = rlnorm(unit, 1))
pf <- simpf(nodes, model.freq = mf, model.sev = ms, weights = wijkt)
frequency(pf)
weights(pf)


###
### CREDIBILITY CALCULATIONS
###

### Currently, two functions are available to carry credibility
### calculations: 'bstraub' and 'cm'. The first one is limited to
### calculations for the Bühlmann and Bühlmann-Straub models, but its
### usage is simpler and it is much faster. The second one can in
### addition do calculations for hierarchical models.

## The package provides the famous data set of Hachemeister (1975) as
## a matrix of 5 lines (one for each state) and 25 columns (the state
## number, 12 periods of ratios, 12 periods of corresponding weights).
data(hachemeister)
hachemeister

## Calculations for a Bühlmann model using 'bstraub'. By default, the
## iterative Bishel-Straub estimator is used for the between variance
## structure parameter. 'bstraub' returns the estimators of the
## structure parameters only. To obtain the credibility premiums, use
## 'predict'.
fit <- bstraub(hachemeister[, 2:13])
fit                          # print method
summary(fit)                 # more information
fit$individual               # individual (weighted) averages
fit$weights                  # total weights
fit$collective               # estimator of the collective mean
fit$s2                       # estimator of the within variance
fit$unbiased                 # unbiased estimator of between variance
fit$iterative                # iterative estimator of between variance
predict(fit)                 # credibility premiums

## Calculations for a Bühlmann-Straub model require a matrix of
## weights.
fit <- bstraub(hachemeister[, 2:13], hachemeister[, 14:25])
summary(fit)
predict(fit)

## For simple models, function 'cm' is somewhat unduly complicated to
## use. It is designed to be used similarly to 'lm' of package
## stats. One must specify the (hierarchical) structure of the data
## using a formula object, a matrix or data frame, and which columns
## contain the ratios and the weights.
##
## Same calculations as above. (Note that the unbiased between
## variance estimator is not available.)
fit <- cm(~state, hachemeister,
          ratios = ratio.1:ratio.12, weights = weight.1:weight.12)
summary(fit)

## Calculations for the four level hierarchical portfolio simulated
## above.
DF <- data.frame(aggregate(pf), weight = weights(pf)[, -(1:3)])

fit <- cm(~sector + sector:unit + sector:unit:contract,
          data = DF, ratios = year.1:year.6,
          weights = weight.year.1:weight.year.6)
fit
predict(fit)                          # credibility premiums
predict(fit, levels = "unit")         # unit credibility premiums only
summary(fit)                          # portfolio summary
summary(fit, levels = "unit")         # unit portfolio summary only
