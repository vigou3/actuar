### ===== actuar: an R package for Actuarial Science =====
###
### Simulation of a aggregate claim amounts
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

simS <- function(n, model.freq, model.sev)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()

    ## Prepare the call to simpf() by building up 'nodes'.
    level.names <- names(if (is.null(model.freq)) model.sev else model.freq)
    nlevels <- length(level.names)
    nodes <- as.list(c(rep(1, nlevels - 1), n))
    names(nodes) <- level.names

    ## Get sample
    x <- aggregate(simpf(nodes = nodes,
                         model.freq = model.freq,
                         model.sev = model.sev))[-1]

    FUN <- ecdf(x)
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("Call", Call, env = environment(FUN))
    comment(FUN) <- "Approximation by simulation"
    FUN
}
