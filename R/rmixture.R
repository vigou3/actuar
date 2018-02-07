### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Simulation of discrete mixtures
###
###    f(x) = p_1 f_1(x) + ... + p_n f_n(x).
###
### Uses the syntax of simul() for model specfification.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

rmixture <- function(n, probs, models)
{
    ## Number of models in the mixture.
    m <- max(length(probs), length(models))

    ## Number of variates from each model. Note that 'rmultinom' will
    ## normalize probabilities to sum 1.
    x <- rmultinom(1, n, prob = rep_len(probs, m))

    ## Auxiliary function to generate 'n' variates from the model
    ## given in 'expr'.
    f <- function(n, expr)
    {
        expr$n <- n
        eval(expr)
    }

    ## Simulate from each model the appropriate number of times and
    ## return result as an atomic vector.
    unlist(mapply(f, n = x, expr = rep_len(models, m),
                  SIMPLIFY = FALSE))
}


