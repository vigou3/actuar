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

    #for (i in seq_along(nodes))
    for (i in 1)
    {
        n.current <- rep(nodes[[i]], length = n.above)
        Call <- substitute(model[[i]])
        print(Call)
        Call$n <- n.current
        assign(nnodes[[i]], eval(Call))
        n.above <- n.current
    }

}

nodes <- list(class = 3,
              contract = c(3, 4, 2),
              year = c(5, 4, 3, 4, 3, 5, 4, 4, 5))
model <- list(class = expression(rgamma(2, 3)),
              contract = expression(rgamma(class, 4)),
              year = expression(rpois(contract)))

model <-

simpf2(nodes = nodes, model = list(class = rgamma(2, 3),
              contract = rgamma(class, 4),
              year = rpois(contract)))
