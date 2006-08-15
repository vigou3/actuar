
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

### Method of knots() for objects of class 'ogive'. Identical to
### stats::knots.stepfun().
knots.ogive <- stats:::knots.stepfun

