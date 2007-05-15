### ===== actuar: an R package for Actuarial Science =====
###
### Simulation of hierarchical portfolios of data. Claim number and
### claim amounts in any given node are simulated independently. Both
### frequency and severity models can be mixtures of distributions.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair, Louis-Philippe Pouliot and Tommy Ouellet

simpf <- function(nodes, model.freq = NULL, model.sev = NULL, weights = NULL, model = "hierarchical")
{
    ## Model validation : hierarchical or crossed model.
    if ((model != "hierarchical") && (model != "crossed"))
        stop("model must be either a hierarchical model or a crossed model")
    
    ## Sanity checks: 'nodes' must be a named list; at least either of
    ## 'model.freq' or 'model.sev' should be non-NULL; level names
    ## (and consequently the number of levels) should be the same
    ## everywhere.
    if (missing(nodes) || !is.list(nodes))
        stop("'nodes' must be a named list")

    has.freq <- !(is.null(model.freq) || # frequency model present?
                 all(sapply(model.freq, is.null)))
    has.sev  <- !(is.null(model.sev) ||  # severity model present?
                 all(sapply(model.sev, is.null)))
    if (!has.freq && !has.sev)
        stop("one of 'model.freq' or 'model.sev' must be non-NULL")
    if ((has.freq && !identical(names(nodes), names(model.freq))) ||
        (has.sev  && !identical(names(nodes), names(model.sev))))
        stop("level names different in 'nodes', 'model.freq' and 'model.sev'")

    ## The function is written for models with at least two levels
    ## (entity and year). If there is only one, add a dummy level to
    ## avoid scattering the code with conditions.
    if (length(nodes) < 2)
    {
        nodes <- c(node = 1, nodes)
        model.freq <-
            if (has.freq) c(expression(node = NULL), model.freq) else NULL
        model.sev <-
            if (has.sev) c(expression(node = NULL), model.sev) else NULL
    }

    ## Frequently used quantities
    level.names <- names(nodes)         # level names
    nlevels <- length(nodes)            # number of levels

    ## If it is a crossed model, we must ensure that every level has
    ## only one single parameter. Also, the penultimate level must equal
    ## the multiplication of all preceding levels.
    if (model == "crossed")
    {
        for (i in 1:(nlevels - 1))
            if (length(nodes[[i]]) != 1)
                stop("invalid parameters for a crossed model; all node levels must have a length of one")
        prec <- 1
        for (i in 1:(nlevels-2)) prec <- prec * nodes[[i]]
        if (nodes[nlevels-1] != prec)
            stop("the penultimate level must equal the multiplication of all preceding levels")
    }
    
    ## Recycling of the number of nodes (if needed) must be done
    ## "manually". We do it here once and for all since in any case
    ## below we will need to know the total number of nodes in the
    ## portfolio. Furthermore, the recycled list 'nodes' will be
    ## returned by the function.
    if (model == "hierarchical")
        for (i in 2:nlevels)       # first node doesn't need recycling
            nodes[[i]] <- rep(nodes[[i]], length = sum(nodes[[i - 1]]))

    ## Simulation of the frequency risk (or mixing) parameters for
    ## each level (e.g. class, contract) and, at the last level, the
    ## actual frequencies. If 'model.freq' is NULL, this is equivalent
    ## to having one claim per node.
    if (has.freq)
    {
        ## Normally, only the immediately above mixing parameter will
        ## be used in the model for a level, but the code here allows
        ## for more general schemes. For this to work, all mixing
        ## parameters have to be correctly recycled at each level
        ## where they *could* be used. Since the model at any level
        ## could be NULL, 'params' will keep track of the mixing
        ## parameters that were simulated in previous iteration of the
        ## forthcoming loop.
        params <- character(0)

        for (i in seq_len(nlevels))
        {
            ## Number of nodes at the current level
            n.current <- nodes[[i]]

            ## Extract simulation model for the level.
            Call <- model.freq[[i]]

            if (model == "hierarchical")
            ## Repeat the mixing parameters of all levels above the
            ## current one that were simulated in the past.
                for (j in seq_along(params))
                    eval(substitute(x <- rep.int(x, n.current),
                                    list(x = as.name(params[[j]]))))

            if (model == "crossed")
            {
                if (i == nlevels - 1)
                {
                    for (j in seq_len(i - 1))
                    {
                        temp <- 1
                        if (j == 1) before <- 1
                        else before <- for (k in 1:(j - 1)) temp <- temp * nodes[[k]]
                        temp <- 1
                        if (j == i - 1) after <- 1
                        else after <- for(k in (j + 1):(i - 1)) temp <- temp * nodes[[k]]
                        eval(substitute(x <- rep(x, before, each=after),
                                        list(x = as.name(params[[j]])))) 
                    }       
                }
                
                else if (i == nlevels)
                    eval(substitute(x <- rep(x, each = nodes[[i]]),
                                    list(x = as.name(params[[i - 1]])))) 
            }
            
            ## Simulate data only if there is a model at the current
            ## level.
            if (!is.null(Call))
            {
                ## Add the number of variates to the call.
                if (model == "hierarchical") Call$n <- sum(n.current)
                if (model == "crossed")
                {
                    if (i == nlevels) Call$n <- nodes[[i-1]] * nodes[[i]]
                    else Call$n <- nodes[[i]]
                }

                ## Simulation of the mixing parameters or the data. In
                ## the latter case, store the results in a fixed
                ## variable name.
                if (i < nlevels)
                {
                    assign(level.names[[i]], eval(Call))
                    params[i] <- level.names[[i]] # remember the parameter
                }
                else
                    frequencies <- eval(Call)
               
            }
        }
    }
    else
        frequencies <- rep.int(1, sum(nodes[[nlevels]]))

    ## Simulation of the claim amounts. If 'model.sev' is NULL, this
    ## is equivalent to simulating frequencies only.
    if (has.sev)
    {
        ## Repeat the same procedure as for the frequency model, with
        ## one difference: when reaching the last level (claim
        ## amounts), the number of variates to simulate is not given
        ## by the number of nodes but rather by the number of claims
        ## as found in 'frequencies'.
        params <- character(0)

        for (i in seq_len(nlevels))
        {
            n.current <- nodes[[i]]
            Call <- model.sev[[i]]

            if (model == "hierarchical")
                for (j in seq_along(params))
                    eval(substitute(x <- rep.int(x, n.current),
                                    list(x = as.name(params[[j]]))))

            if (model == "crossed")
            {
                if (i == nlevels - 1)
                {
                    for (j in seq_len(i - 1))
                    {
                        temp <- 1
                        if (j == 1) before <- 1
                        else before <- for (k in 1:(j - 1)) temp <- temp * nodes[[k]]
                        temp <- 1
                        if (j == i - 1) after <- 1
                        else after <- for(k in (j + 1):(i - 1)) temp <- temp * nodes[[k]]
                        eval(substitute(x <- rep(x, before, each=after),
                                        list(x = as.name(params[[j]])))) 
                    }       
                }
                
                else if (i == nlevels)
                    eval(substitute(x <- rep(x, each = nodes[[i]]),
                                    list(x = as.name(params[[i - 1]])))) 
            }

            if (!is.null(Call))
            {
                ## The rest of the procedure differs depending if we
                ## are still simulating mixing parameters or claim
                ## amounts.
                if (i < nlevels)
                {
                    ## Simulation of mixing parameters is identical to the
                    ## simulation of frequencies.
                    if (model == "hierarchical") Call$n <- sum(n.current)
                    if (model == "crossed") Call$n <- nodes[[i]]
                    assign(level.names[[i]], eval(Call))
                    params[i] <- level.names[[i]]
                }
                else
                {
                    if (model == "hierarchical")
                    {
                        ## For the simulation of claim amounts, the number
                        ## of variates is rather given by the
                        ## 'frequencies' object. Furthermore, the mixing
                        ## parameters must be recycled once more to match
                        ## the vector of frequencies.
                        for (p in intersect(all.vars(Call), params))
                            eval(substitute(x <- rep.int(x, frequencies),
                                            list(x = as.name(p))))
                        Call$n <- sum(frequencies)
                        severities <- eval(Call)
                    }
                    if (model == "crossed")
                    {
                        eval(substitute(x <- rep.int(x, frequencies),
                                        list(x = as.name(params[[i - 1]]))))
                        Call$n <- sum(frequencies)
                        severities <- eval(Call)
                    }
                }
            }
        }
    }
    else
        severities <- rep.int(1, sum(frequencies))

    ## We must now distribute the claim amounts in vector 'severities'
    ## to the appropriate nodes. This is complicated by the
    ## possibility to have different number of nodes (years of
    ## observation) for each "entity". The result must be a matrix
    ## with the number of columns equal to the maximum number of last
    ## level nodes.
    ##
    ## The number of nodes (years of observation) per "entity" is
    ## given by 'n.current' since we reached the last level in (either
    ## one of) the above loops.
    ##
    ## Assign a unique ID to each node, leaving gaps for nodes without
    ## observations.
    if (model == "hierarchical")
        ind <- unlist(mapply(seq,
                             from = seq(by = max(n.current), along = n.current),
                             length = n.current))
    if (model == "crossed")
        ind <- unlist(mapply(seq,
                             from = seq(by = max(n.current),
                             along = seq_len(nodes[[nlevels - 1]])),
                             length = n.current))

    ## Repeating the vector of IDs according to the frequencies
    ## effectively assigns a node ID to each claim amount. The vector
    ## of claim amounts is then split by node, yielding a list where
    ## each element corresponds to a node with claims.
    f <- rep.int(ind, frequencies)
    severities <- split(severities, f)
    
    ## Identify nodes with frequency equal to 0, which is different
    ## from having no observation (NA).
    freq0 <- ind[which(frequencies == 0)]

    ## Rearrange the list of claim amounts in a matrix;
    ##
    ##      number of rows: number of nodes at the penultimate level
    ##                      (number of entities)
    ##   number of columns: maximum number of nodes at the last level
    ##                      (number of years of observation).
    ##
    ## Moreover, assign a value of 'numeric(0)' to nodes with a
    ## frequency of 0.
    if (model == "hierarchical") nrow <- length(n.current)     # number of entities
    if (model == "crossed") nrow <- nodes[[nlevels - 1]]       # number of entities
    ncol <- max(n.current)                                     # number of years
    res <- as.list(rep.int(NA, nrow * ncol))
    res[unique(f)] <- severities
    res[freq0] <- lapply(rep.int(0, length(freq0)), numeric)
    res <- matrix(res, nrow, ncol, byrow = TRUE,
                  dimnames = list(NULL,
                  paste(level.names[nlevels], seq_len(ncol), sep = ".")))

    ## Finally, create a matrix where each row contains the series of
    ## identifiers for an "entity" in the portfolio, e.g. if the data
    ## is denoted X_{ijkt}, one line of the matrix will contain
    ## subscripts i, j and k. As we move from right to left in the
    ## columns of 'm', the subcripts are increasingly repeated.
    if (model == "hierarchical") ncol <- nlevels - 1
    if (model == "crossed") ncol <- nlevels - 2
    
    m <- matrix(1, nrow, ncol,
                dimnames = list(NULL, head(level.names, ncol)))

    
    if (model == "hierarchical")
    {
        for (i in seq_len(ncol - 1))    # all but the last column
        {
            ## Vector 'x' will originally contain all subscripts for one
            ## level. These subscripts are then repeated as needed to give
            ## the desired result. To avoid another explicit loop, I use a
            ## 'lapply' with a direct assignment in the current
            ## frame. Somewhat unusual, but this is the simplest procedure
            ## I managed to come up with.
            x <- unlist(lapply(nodes[[i]], seq))
            lapply(nodes[(i + 1):(nlevels - 1)],
                   function(v) assign("x", rep.int(x, v), envir = parent.frame(2)))
            m[, i] <- x
        }
        m[, ncol] <- unlist(lapply(nodes[[ncol]], seq)) # last column
    }

    if (model == "crossed")
    {
        for (i in seq_len(ncol))
        {
            temp <- 1
            if (i == 1) before <- 1
            else before <- for (k in 1:(i - 1)) temp <- temp * nodes[[k]]
            temp <- 1
            if (i == ncol) after <- 1
            else after <- for(k in (i + 1):ncol) temp <- temp * nodes[[k]]
            x <- rep(seq_len(nodes[[i]]), before, each=after)
            m[, i] <- x
        }
    }
    
    ## Reshape weights as a matrix, if necessary.
    weights <- if (is.null(weights))
        NULL
    else
    {
        ## Integrate NAs into the weights matrix as appropriate.
        w <- rep.int(NA, nrow * ncol)
        w[ind] <- weights
        matrix(w, nrow = nrow, byrow = TRUE, dimnames = dimnames(res))
    }
    
    ## Return object of class 'simpf'
    structure(list(data = res,
                   weights = weights,
                   classification = m,
                   nodes = nodes,
                   model.freq = model.freq,
                   model.sev = model.sev,
                   model = model),
              class = "simpf")
}

### 'print' method for 'simpf' objects
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
                  x$model.sev,
                  "\n", sep = ""), sep = "")
    }
    cat("\n  Number of claims per node: \n\n")
    print(frequency(x))
    invisible(x)
}
