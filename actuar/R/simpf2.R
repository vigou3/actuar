simpf2 <- function(nodes, model.freq = NULL, model.sev = NULL, weights = NULL)
{
    ## Sanity checks: level names (and consequently the number of
    ## levels) should be the same everywhere.
    if (!(identical(names(nodes), names(model.freq)) &
          identical(names(nodes), names(model.sev))))
        stop("level names are different in 'nodes', 'model.freq' and 'model.sev'")

    ## Level names
    level.names <- names(nodes)

    ## Number of levels
    nlevels <- length(nodes)

    ## Number of nodes at the current level. At first, there is only
    ## the portfolio.
    n.current <- 1

    ## Simulate risk (or mixing) parameters for each level
    ## (e.g. class, contract). The last level contains the actual
    ## data.
    for (i in seq_len(nlevels))
    {
        ## Going down one level in the model. Remember the *total*
        ## number of nodes at this level.
        n.above <- sum(n.current)

        ## Number of nodes at the new level, recycled if necessary.
        n.current <- rep(nodes[[i]], length = n.above)

        ## Extract simulation model for the level.
        Call <- model.freq[[i]]

        ## Correctly repeat the mixing parameters of all levels above
        ## the current one. Normally, only the immediately above
        ## mixing parameter will be used in the model for a level, but
        ## the code here allows for all preceeding parameters to be
        ## used.
        for (param in intersect(all.vars(Call), level.names))
        {
            eval(substitute(x <- rep(x, n.current),
                            list(x = as.name(param))))
            cat(param, "\n", eval(as.name(param)), fill = TRUE)
        }

        ## Add the number of variates to the call.
        Call$n <- sum(n.current)
        print(Call)

        ## Simulation of the mixing parameters or the data. In the
        ## latter case, store the results in a fixed variable name.
        assign(ifelse(i < nlevels, level.names[[i]], "frequencies"),
               eval(Call))
        cat(eval(as.name(ifelse(i < nlevels, level.names[[i]], "frequencies"))), "\n\n", fill = TRUE)
    }

    ## Repeat the same procedure for the severity model, with two
    ## differences:
    ##
    ## 1. the nodes at the last level is not the years of
    ##    observations, but rather the frequencies simulated above;
    ## 2. since there can be many claim amounts per contract, the
    ##    simulation of data (the last level) is postponed to after
    ##    the loop, when the simulation will be done with 'lapply'.
    n.current <- 1                      # get back to portfolio level

    for (i in seq_len(nlevels))
    {
        n.above <- sum(n.current)
        n.current <- rep(nodes[[i]], length = n.above)
        Call <- model.sev[[i]]
        for (param in intersect(all.vars(Call), level.names))
        {
            eval(substitute(x <- rep(x, n.current),
                            list(x = as.name(param))))
            cat(param, "\n", eval(as.name(param)), fill = TRUE)
        }

        ## The rest of the procedure differs depending if we are still
        ## simulating mixing parameters or claim amounts.
        if (i < nlevels)
        {
            ## Simulation of mixing parameters is identical to the
            ## simulation of frequencies.
            Call$n <- sum(n.current)
            print(Call)
            assign(level.names[[i]], eval(Call))
        }
        else
        {
            ## For the simulation of claim amounts, the number of
            ## variates is rather given by the 'frequencies'
            ## object. Furthermore, the mixing parameters must be
            ## replicated once more to match the vector of
            ## frequencies.
            for (param in intersect(all.vars(Call), level.names))
            {
                eval(substitute(x <- rep(x, frequencies),
                                list(x = as.name(param))))
                cat(param, "\n", eval(as.name(param)), fill = TRUE)
            }
            Call$n <- sum(frequencies)
            print(Call)
            severities <-eval(Call)
        }
        cat(eval(as.name(ifelse(i < nlevels, level.names[[i]], "severities"))), "\n\n", fill = TRUE)
    }
    severities
}

nodes <- list(group = 3,
              contract = c(3, 4, 2),
              year = c(5, 4, 3, 4, 3, 5, 4, 4, 5))

model.freq <- expression(group = rexp(0.5),
    contract = rgamma(group, 1),
    year = rpois(contract))

model.sev <- expression(group = rnorm(6, 0.1),
    contract = rnorm(group, 1),
    year = rlnorm(contract, 1))

simpf2(nodes, model.freq, model.sev)
