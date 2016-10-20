### ===== actuar: An R Package for Actuarial Science =====
###
### Simulation of standard, non hierarchical, compound models. Uses a
### simplified version of the syntax of simul() for model
### specfification.
###
### Where simul() was developed for flexibility, this one is aimed at
### execution speed. Various algorithms where tested. No argument
### validity checks.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

rcompound <- function(n, model.freq, model.sev)
{
    ## Convert model expressions into language objects.
    cl.freq <- substitute(model.freq)
    cl.sev <- substitute(model.sev)

    ## If model expressions are wrapped into 'expression' as in
    ## 'simul', get rid of the call.
    if (cl.freq[[1]] == "expression")
        cl.freq <- cl.freq[[-1]]
    if (cl.sev[[1]] == "expression")
        cl.sev <- cl.sev[[-1]]

    ## Initialize the output vector. We will use the fact that 'res'
    ## is filled with zeros later.
    res <- numeric(n)

    ## Add the number of variates to the 'model.freq' call.
    cl.freq$n <- n

    ## Generate frequencies.
    N <- eval(cl.freq)

    ## Add the number of variates to the 'model.sev' call.
    cl.sev$n <- sum(N)

    ## Generate all severities.
    x <- eval(cl.sev)

    ## Create a vector that will be used as a factor to regroup
    ## severities for the computation of aggregate values. Idea:
    ## assign one integer to each frequency and repeat that integer a
    ## number of times equal to the frequency. For example, if the
    ## frequencies are (2, 0, 1, 3), then the vector will be (1, 1, 3,
    ## 4, 4, 4).
    f <- rep.int(seq_len(n), N)

    ## Compute aggregate values and put them in the appropriate
    ## positions in the output vector. The positions corresponding to
    ## zero frequencies are already initialized with zeros.
    res[which(N != 0)] <- tapply(x, f, sum)

    res
}
