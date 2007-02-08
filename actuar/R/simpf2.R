simpf2 <- function(nodes, model, weights)
{
    ## Level names. Used to store variables later.
    nnodes <- names(nodes)

    ## Sanity check: level names should be the same everywhere.
    if (!identical(names(nodes), names(model)))
        stop("level names are different in 'nodes' and 'model'")

    ## Number of nodes one level above the current level. At first,
    ## there is only the portfolio.
    n.above <- 1

    for (i in seq_along(nodes))
    {
        ## Number of nodes at the current level, recycled if
        ## necessary
        n.current <- rep(nodes[[i]], length = n.above)

        ## Extract simulation model for the level
        Call <- model[[i]]

        ## Correctly repeat the mixing parameters of all levels above
        ## the current one. Normally, only the immediately above
        ## mixing parameter will be used in the model for a level, but
        ## the code here allows for other preceeding parameters to be
        ## used.
        for (param in intersect(all.vars(Call), nnodes))
        {
            eval(substitute(x <- rep(x, n.current),
                            list(x = as.name(param))))
            cat(param, "\n", eval(as.name(param)), fill = TRUE)
        }

        ## Add the number of variates to the model.
        Call$n <- sum(n.current)
        print(Call)

        ## Simulation
        assign(nnodes[[i]], eval(Call))
        cat(eval(as.name(nnodes[[i]])), "\n\n", fill = TRUE)

        ## Remember the total number of nodes at the current level for
        ## the next iteration.
        n.above <- sum(n.current)
    }

    list(group, contract, year)

}

nodes <- list(group = 3,
              contract = c(3, 4, 2),
              year = c(5, 4, 3, 4, 3, 5, 4, 4, 5))

model <- expression(group = rgamma(2, 3),
    contract = rgamma(group, 4),
    year = rpois(contract))

simpf2(nodes = nodes, model = model)

model <- expression(group = rgamma(2, 3),
    contract = rgamma(group, 4),
    year = rpois(contract * group))

simpf2(nodes = nodes, model = model)
