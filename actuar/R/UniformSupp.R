### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}exp functions to compute raw and limited
### moments, and the moment generating function
### for the Uniform distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter ?? of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>


mgfunif <- function(x, min = 0, max = 1, log = FALSE)
    .External("do_dpq", "mgfunif", x, min, max, log)
