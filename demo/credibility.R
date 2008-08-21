### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the credibility theory facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)


## The package provides the famous data set of Hachemeister (1975) as
## a matrix of 5 lines (one for each state) and 25 columns (the state
## number, 12 periods of ratios, 12 periods of corresponding weights).
data(hachemeister)
hachemeister

## Fitting of a Buhlmann model to the Hachemeister data set using
## function 'cm'. The interface of the function is similar to 'lm'.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12)
fit                          # print method
summary(fit)                 # more information
fit$means                    # (weighted) averages
fit$weights                  # total weights
fit$unbiased                 # unbiased variance estimators
predict(fit)                 # credibility premiums

## Fitting of a Buhlmann-Straub model require weights. Here, iterative
## estimators of the variance components are used.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12, method = "iterative")
summary(fit)
predict(fit)
## plot the contract number 2, prediction and heterogeneity of the portfolio
plot.cm(fit, levels = NULL, node = 2, which = 1:2, from = NULL, to = NULL)

## Fitting of a Hachemeister regression model. This requires to
## specify a vector or matrix of regressors with argument 'xreg'.
## 'xreg' must be a matrix and represents the regressor = time here.
## The boolean adjust is set at TRUE so as to adjust intercept
## and slope in the regression model. Iterative estimators.
fit <- cm(~state, hachemeister, ratios = ratio.1:ratio.12,
          weights = weight.1:weight.12, regformula = ~time,
          regdata = data.frame(time=12:1), adj.intercept = TRUE, method = "iterative")
summary(fit, newdata = data.frame(time = 0))     # 'newdata' is the future value of regressor
predict(fit, newdata = data.frame(time = 0))
## graph of the Hachemeister regression lines, contract No 4
plot.cm(fit, levels = NULL, node = 4, which = 1:2, from = 1, to = 12)

## Simulation of a three level hierarchical portfolio.
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
wijkt <- runif(50, 2, 10)
wijkt <- runif(300, rep(0.5 * wijkt, each = 6), rep(1.5 * wijkt, each = 6))
pf <- simul(nodes, model.freq = mf, model.sev = ms, weights = wijkt)

## Fitting of a hierarchical model to the portfolio simulated above.
DB <- cbind(weights(pf, prefix = "weight."),
            aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE))

fit <- cm(~sector + sector:unit + sector:unit:contract,
          data = DB, ratios = year.1:year.6,
          weights = weight.year.1:weight.year.6)
fit
predict(fit)                          # credibility premiums
predict(fit, levels = "unit")         # unit credibility premiums only
summary(fit)                          # portfolio summary
summary(fit, levels = "unit")         # unit portfolio summary only
## plot the premium lines and prediction, then the heterogeneity
## at levels 'sector', for the sector 2.
plot.cm(fit, levels = "sector", node = 2, which = 1:2, from = NULL, to = NULL)
## the same at levels 'unit', for the unit 3 in the sector 2.
plot.cm(fit, levels = "unit", node = c(2,3), which = 1:2, from = NULL, to = NULL)
## the same at levels 'contract', for the contract 1 in the
## unit 4 of the sector 2.
plot.cm(fit, levels = "contract", node = c(2,4,1), which = 1:2, from = NULL, to = NULL)


## Simulation of a 4 levels hierarchical portfolio.
nodes <- list(sector = 2, unit = c(2, 2), interm = c(3, 2, 4, 1),
              contract = c(4,3,2,5,7,2,5,1,3,6), year = 6)
mf <- expression(sector = rexp(2),
                 unit = rgamma(sector, 0.1),
                 interm = rgamma(unit, 0.1),
                 contract = rgamma(interm, 1),
                 year = rpois(weights * contract))
ms <- expression(sector = rnorm(2, sqrt(0.1)),
                 unit = rnorm(sector, 1),
                 interm = rnorm(unit, 1),
                 contract = NULL,
                 year = rlnorm(unit, 1))
wijkt <- runif(38, 2, 10)
wijkt <- runif(228, rep(0.5 * wijkt, each = 6), rep(1.5 * wijkt, each = 6))
pf <- simul(nodes, model.freq = mf, model.sev = ms, weights = wijkt)

## Fitting of a hierarchical model to the portfolio simulated above.
DB <- cbind(weights(pf, prefix = "weight."),
            aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE))

fit <- cm(~sector + sector:unit + sector:unit:interm + sector:unit:interm:contract,
          data = DB, ratios = year.1:year.6,
          weights = weight.year.1:weight.year.6)
predict(fit)                          # credibility premiums
predict(fit, levels = "interm")         # unit credibility premiums only
summary(fit)                          # portfolio summary
summary(fit, levels = "interm")         # unit portfolio summary only
## plot the premium lines and prediction, then the heterogeneity
## at levels 'sector', for the sector 2.
plot.cm(fit, levels = "sector", node = 2, which = 1:2, from = NULL, to = NULL)
## the same at levels 'unit', for the unit 1 in the sector 2.
plot.cm(fit, levels = "unit", node = c(2,1), which = 1:2, from = NULL, to = NULL)
## the same at levels 'interm', for the interm 3 in the
## unit 1 of the sector 2.
plot.cm(fit, levels = "interm", node = c(2,1,4), which = 1:2, from = NULL, to = NULL)
## the same at levels 'contract', for the contract 2 in the
## interm 3 of the unit 1 of the sector 2.
plot.cm(fit, levels = "contract", node = c(2,1,4,2), which = 1:2, from = NULL, to = NULL)
