simpf2 <- function(nodes, model.freq = NULL, model.sev = NULL, weights = NULL)
{
    ## Sanity checks: at least either of 'model.freq' or 'model.sev'
    ## should be specified; level names (and consequently the number
    ## of levels) should be the same everywhere.
    hasfreq <- !is.null(model.freq)     # frequency model present?
    hassev  <- !is.null(model.sev)      # severity model present?
    if (!hasfreq & !hassev)
        stop("one of 'model.freq' or 'model.sev' must be non-NULL")
    if ((hasfreq & !identical(names(nodes), names(model.freq))) |
        (hassev  & !identical(names(nodes), names(model.sev))))
        stop("level names are different in 'nodes', 'model.freq' and 'model.sev'")

    ## Level names
    level.names <- names(nodes)

    ## Number of levels
    nlevels <- length(nodes)

    ## Recycling of the number of nodes (if needed) must be done
    ## "manually". We do it here once and for all since in any case
    ## below we will need to know the total number of nodes in the
    ## portfolio. Furthermore, the recycled 'nodes' will be returned
    ## by the function.
    ##
    ## Number of nodes at the current level. At first, there is only
    ## the portfolio.
    n.current <- 1

    for (i in seq_len(nlevels))
    {
        ## Going down one level in the model. Remember the *total*
        ## number of nodes at this level.
        n.above <- sum(n.current)

        ## Number of nodes at the new level. This is where recycling
        ## is done.
        nodes[[i]] <- n.current <- rep(nodes[[i]], length = n.above)
    }

    ## Simulation of the risk (or mixing) parameters for each level
    ## (e.g. class, contract) and, at the last level, the actual
    ## frequencies. If 'model.freq' is NULL, this is equivalent to
    ## having one claim per node.
    if (hasfreq)
    {
        for (i in seq_len(nlevels))
        {
            ## Number of nodes at the current level
            n.current <- nodes[[i]]

            ## Extract simulation model for the level.
            Call <- model.freq[[i]]

            ## Correctly repeat the mixing parameters of all levels
            ## above the current one. Normally, only the immediately
            ## above mixing parameter will be used in the model for a
            ## level, but the code here allows for all preceeding
            ## parameters to be used.
            for (param in intersect(all.vars(Call), level.names))
                eval(substitute(x <- rep.int(x, n.current),
                                list(x = as.name(param))))

            ## Add the number of variates to the call.
            Call$n <- sum(n.current)

            ## Simulation of the mixing parameters or the data. In the
            ## latter case, store the results in a fixed variable
            ## name.
            assign(ifelse(i < nlevels, level.names[[i]], "frequencies"),
                   eval(Call))
        }
    }
    else
        frequencies <- rep(1, sum(nodes[[nlevels]]))

    ## Simulation of the claim amounts. If 'model.sev' is NULL, this
    ## is equivalent to simulating frequencies only.
    if (hassev)
    {
        ## Repeat the same procedure as for the frequency model, with
        ## one difference: when reaching the last level (claim
        ## amounts), the number of variates to simulate is not given
        ## by the number of nodes but rather by the number of claims,
        ## as simulated above.
        for (i in seq_len(nlevels))
        {
            n.current <- nodes[[i]]         # no need to recycle this time
            Call <- model.sev[[i]]
            for (param in intersect(all.vars(Call), level.names))
                eval(substitute(x <- rep.int(x, n.current),
                                list(x = as.name(param))))

            ## The rest of the procedure differs depending if we are
            ## still simulating mixing parameters or claim amounts.
            if (i < nlevels)
            {
                ## Simulation of mixing parameters is identical to the
                ## simulation of frequencies.
                Call$n <- sum(n.current)
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
                    eval(substitute(x <- rep.int(x, frequencies),
                                    list(x = as.name(param))))
                Call$n <- sum(frequencies)
                severities <-eval(Call)
            }
        }
    }
    else
        severities <- rep(1, sum(frequencies))

    ## Vector 'severities' contains the claim amounts. We must now
    ## distribute these claim amounts to the appropriate nodes. This
    ## is complicated by the possibility to have different number of
    ## nodes (years of experience) for each "contract". The result
    ## must be a matrix with the number of columns equal to the
    ## maximum number of last level nodes.
    ##
    ## The number of nodes (years of observation) per "contract" is
    ## given by 'n.current' since we reached the last level in (either
    ## one of) the above loops.
    ##
    ## Assign a unique id to each node, leaving gaps for nodes without
    ## observations.
    ind <- unlist(mapply(seq,
                         from = seq(by = max(n.current), along = n.current),
                         length = n.current))

    ## Assign a node to each claim amount, then split claim amounts by
    ## node. The result is a list with one element per node where
    ## there are claims.
    f <- rep.int(ind, frequencies)
    severities <- split(severities, f)

    ## Nodes with frequency equal to 0, which is different from having
    ## no observation (NA).
    freq0 <- ind[which(frequencies == 0)]

    ## Finally, rearrange the list of claim amounts in a matrix;
    ##
    ##      number of rows: number of nodes at the penultimate level
    ##                      (number of contracts)
    ##   number of columns: maximum number of nodes at the last level
    ##                      (number of years of observation).
    ##
    ## Moreover, assign a value of 'numeric(0)' to nodes with a
    ## frequency of 0.
    nrow <- length(n.current)
    ncol <- max(n.current)
    data <- as.list(rep.int(NA, nrow * ncol))
    data[unique(f)] <- severities
    data[freq0] <- lapply(rep.int(0, length(freq0)), numeric)
    data <- matrix(data, nrow, ncol, byrow = TRUE,
                   dimnames = list(NULL,
                   paste(level.names[nlevels], seq_len(ncol), sep = ".")))

    ## Finally, create a matrix where each row contains the series of
    ## identifiers for a "contract" in the portfolio, e.g. if the data
    ## is denoted X_{ijkt}, one line of the matrix will contain
    ## subscripts i, j and k. As we move from right to left in the
    ## columns of 'm', the subcripts are increasingly repeated.
    ncol <- nlevels - 1
    m <- matrix(NA, nrow, ncol,
                dimnames = list(NULL, head(level.names, ncol)))
    for (i in seq_len(ncol - 1))        # all but the last column
    {
        ## Vector 'x' will originally contain all subscripts for one
        ## level. These subscripts are then repeated as needed to give
        ## the desired result. To avoid another explicit loop, I use a
        ## 'lapply' with a direct assignment in the current
        ## frame. Somewhat unusual, but this is the simplest procedure
        ## I came up with.
        x <- unlist(lapply(nodes[[i]], seq))
        lapply(nodes[(i + 1):(nlevels - 1)],
               function(v) assign("x", rep(x, v), envir = parent.frame(2)))
        m[, i] <- x
    }
    m[, ncol] <- unlist(lapply(nodes[[ncol]], seq)) # last column

    structure(list(freq = frequencies,
                   data = data,
                   classification = m,
                   nodes = nodes,
                   model.freq = model.freq,
                   model.sev = model.sev),
              class = "simpf")
}


print.simpf <- function(x, ...)
{
    cat("\nPortfolio of claim amounts \n\n")
    nn <- names(x$nodes)
    nc <- max(nchar(nn))
    if (!is.null(x$model.freq))
    {
        cat("  Frequency model\n")
        cat(paste("    ",
                  format(nn, width = nc),
                  " ~ ",
                  x$model.freq,
                  "\n", sep = ""), sep = "")
    }
    if (!is.null(x$model.sev))
    {
        cat("  Severity model\n")
        cat(paste("    ",
                  format(nn, width = nc),
                  " ~ ",
                  x$model.freq,
                  "\n", sep = ""), sep = "")
    }
    cat("\n  Number of claims per node: \n\n")
    compfreq <- function(x) if (identical(x, NA)) NA else length(x)
    print(cbind(x$classification,
                array(sapply(x$data, compfreq), dim(x$data),
                      dimnames = dimnames(x$data))))
    invisible(x)
}

### TODO

## * Vérifier si l'un des modèles est NULL
## * Vérifier la prise en compte des poids

### TESTS
nodes <- list(group = 3,
              contract = c(3, 4, 2),
              year = c(5, 4, 3, 4, 3, 5, 4, 4, 5))

model.freq <- expression(group = rexp(0.5),
    contract = rgamma(group, 1),
    year = rpois(contract))

model.sev <- expression(group = rnorm(6, 0.1),
    contract = rnorm(group, 1),
    year = rlnorm(contract, 1))

pf <- simpf2(nodes, model.freq, model.sev)
pf <- simpf2(nodes, model.freq, NULL)
pf

###

nodes <- list(sector = 2,
              unit = c(3, 4),
              contract = c(3, 4, 3, 4, 2, 3, 4),
              year = rep.int(5, 23))

model.freq <- expression(sector = rexp(1),
    unit = rexp(sector),
    contract = rgamma(unit, 1),
    year = rpois(contract))

model.sev <- expression(sector = rnorm(6, 0.1),
    unit = rnorm(sector, 1),
    contract = rnorm(unit, 1),
    year = rlnorm(contract, 1))

pf <- simpf2(nodes, model.freq, model.sev)


modelfreq <- list(dist1 = "pois",
                  par1 = list(lambda = quote(Lambda * weights)),
                  dist2 = "gamma",
                  par2 = c(shape = 2, rate = 1))
modelsev<-list(dist1 = "lnorm",
               par1 = list(meanlog = quote(Theta), sdlog = 1),
               dist2 = "norm",
               par2 = c(mean = 5, sd = 1))
system.time(simpf(100, 12, modelfreq, modelsev))

model.freq <- expression(contrat = rgamma(2, 1),
    year = rpois(contrat))
model.sev <- expression(contrat = rnorm(5, 1),
    year = rlnorm(contrat, 1))
simpf2(nodes = list(contrat = 100, year = 5), model.freq, model.sev)
system.time(simpf2(nodes = list(contrat = 100, year = 5), model.freq, model.sev))
