### ===== actuar: an R package for Actuarial Science =====
###
### Simulation of a aggregate claim amounts
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

simS <- function(n, model.freq, model.sev)
{
    call <- match.call()

    x <- drop(aggregate(simpf(contracts = 1, years = n,
                              model.freq = model.freq,
                              model.sev = model.sev)))

    FUN <- ecdf(x)
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("call", call, env = environment(FUN))
    comment(FUN) <- "Approximation by simulation"
    FUN
}
