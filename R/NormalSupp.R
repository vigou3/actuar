### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,mgf}norm functions to compute raw and limited
### moments, and the moment generating function
### for the Normal distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter 13 of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mnorm <- function(order, mean = 0, sd = 1)
{
    if( as.integer(order) == order )
        .External("do_dpq", "mnorm", order, mean, sd, FALSE)
    else
        stop("non integer order is not supported")
}
    
mgfnorm <- function(x, mean = 0, sd = 1, log = FALSE)
    .External("do_dpq", "mgfnorm", x, mean, sd, log)
