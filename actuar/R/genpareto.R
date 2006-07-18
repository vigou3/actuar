### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}generalized pareto functions. 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dgenpareto <- function(x, shape1, scale, shape2, log = FALSE)
     .External("do_dpq", "dgenpareto", x, shape1, scale, shape2, log)

pgenpareto <- function(q, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "pgenpareto", q, shape1, scale, shape2, lower.tail, log.p) 

qgenpareto <- function(p, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qgenpareto", p, shape1, scale, shape2, lower.tail, log.p) 

rgenpareto <- function(n, shape1, scale, shape2)
     .External("do_random", "rgenpareto", n, shape1, scale, shape2)

mgenpareto <- function(k, shape1, scale, shape2, log = FALSE)
     .External("do_dpq", "mgenpareto", k, shape1, scale, shape2, log)

levgenpareto <- function(d, shape1, scale, shape2, order = 1, log = FALSE)
     .External("do_dpq", "levgenpareto", d, shape1, scale, shape2, order, log)
