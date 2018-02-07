### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### The "beta integral"
###
###    B(a, b; x) = Gamma(a + b) int_0^x t^(a-1) (1 - t)^(b-1) dt
###
### a > 0, b != -1, -2, ..., 0 < x < 1. This mathematicacal is only
### used at the C level in the package. The function therein provides
### an R interface just in case it could be useful.
###
### The function is *not* exported.
###
### See Appendix A of Klugman, Panjer and Willmot (2012), Loss Models,
### Fourth Edition, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## see src/betaint.c
betaint <- function(x, a, b)
    .External(C_actuar_do_dpq, "betaint", x, a, b, FALSE)

