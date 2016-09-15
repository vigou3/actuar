### ===== actuar: An R Package for Actuarial Science =====
###
### The exponential integral and two integrals related to the
### incomplete beta function and the incomplete gamma function. Used
### only at the C level in the package. The functions therein provide
### R interfaces just in case these could be useful.
###
### The functions are *not* exported.
###
### See Appendix A of Klugman, Panjer and Willmot (2012), Loss Models,
### Fourth Edition, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## see src/betaint.c
betaint <- function(x, a, b)
    .External("actuar_do_dpq", "betaint", x, a, b, FALSE)

## see src/expint.c
expint <- function(x)
    .External("actuar_do_dpq", "expint", x, FALSE)

## see src/gammaint.c
gammaint <- function(x, a)
    .External("actuar_do_dpq", "gammaint", x, a, FALSE)

