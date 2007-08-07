### ===== actuar: an R package for Actuarial Science =====
###
### Credibility Models
###
### Fit a credibility model in the formulation of variance components
### as described in Dannenburg, Kaas and Goovaerts (1996). Models
### supported are part of a generalized hierarchical credibility
### theory as introduced in Dannenburg (1995).
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,

cm <- function(formula, data, ratios, weights, subset,
               TOL = 1E-6, echo = FALSE)
{
    Call <- match.call()

    ## === MODEL ANALYSIS ===
    ##
    ## Decompose the formula giving the portfolio structure. Attribute
    ## "order" gives the interaction level of each term in the
    ## formula. In hierarchical structures, each term should represent
    ## a different level, hence there should not be any duplicates in
    ## this attribute. The column names in 'data' containing the
    ## portfolio structure can be obtained from the rownames of the
    ## matrix in attribute "factors".
    ##
    ## Note that the very last level, the data, is not taken into
    ## account here.
    tf <- terms(formula)
    level.numbers <- attr(tf, "order")           # level IDs
    level.names <- rownames(attr(tf, "factors")) # level names
    nlevels <- length(level.names)               # number of levels
    nlevels1p <- nlevels + 1                     # frequently used

    ## Sanity checks
    ##
    ## 1. only hierarchical interactions are allowed in 'formula';
    ## 2. if 'ratios' is missing, all columns of 'data' are taken to
    ##    be ratios, so 'weights' should also be missing.
    if (any(duplicated(level.numbers)))
        stop("unsupported interactions in 'formula'")
    if (missing(ratios) & !missing(weights))
        stop("ratios have to be supplied if weights are")

    ## === DATA EXTRACTION ===
    ##
    ## 'data' is split into three matrices: one for the portfolio
    ## structure, one for the ratios and one for the weights. They are
    ## obtained via calls to subset() built from this function's
    ## call. That way, arguments 'ratios', 'weights' and 'subset' are
    ## not evaluated before being passed to subset(). Argument
    ## matching is as follows:
    ##
    ##   Argument of cm()     Argument of subset()
    ##   ================     ====================
    ##     data                 x
    ##     ratios               select
    ##     weights              select
    ##     subset               subset
    ##
    ## Positions of the arguments that will be needed.
    m <- match(c("data", "ratios", "weights", "subset"), names(Call), 0)

    ## Extraction of the portfolio structure. Arguments 'data' and
    ## 'subset' are passed to subset().
    cl <- Call[c(1, m[c(1, 4)])]        # use data and subset only
    cl[[1]] <- as.name("subset")        # change function name
    names(cl)[2] <- "x"                 # argument matching
    cl$select <- level.names            # add argument 'select'
    levs <- eval(cl, parent.frame())    # extraction

    ## Object 'levs' is a data frame or matrix with as many colums as
    ## there are levels in the model (still notwithstanding the data
    ## level). Rows contain nodes identifiers which can be
    ## anything. For calculations, these identifiers are converted
    ## into simple subscripts (i, j, k, ...) as used in mathematical
    ## notation.
    ##
    ## Note that 'apply' will coerce to a matrix.
    ilevs <- apply(levs, 2, function(x) as.integer(factor(x)))

    ## To symmetrize further calculations, bind a column of ones
    ## representing the affiliation to the global portfolio.
    ilevs <- cbind(pf = 1, ilevs)

    ## Extraction of the ratios. If argument 'ratios' is missing, then
    ## use all columns of 'data' except those of the portfolio
    ## structure.
    cl$select <-
        if (missing(ratios))
            setdiff(colnames(data), level.names)
        else
            Call[[m[2]]]
    ratios <- as.matrix(eval(cl, parent.frame())) # ratios as matrix

    ## Creation of a weight matrix. Extract from data if argument
    ## 'weights' is specified, otherwise create a matrix of ones. For
    ## extraction, the only change from ratio extraction is the
    ## content of element "select" of the call.
    weights <-
        if (missing(weights))
        {
            if (any(is.na(ratios)))
                stop("missing ratios not allowed when weights are not supplied")
            array(1, dim(ratios))       # matrix of ones
        }
        else
        {
            cl$select <- Call[[m[3]]]
            as.matrix(eval(cl, parent.frame())) # weights as matrix
        }

    ## Sanity check if weights and ratios correspond.
    if(!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in weights and in ratios")

    ## === NUMBER OF NODES AND SPLITTING FACTORS ===
    ##
    ## Future computation of per level summaries will require a set of
    ## factors based on the number of nodes in each level. An example
    ## will best explain what is achieved here: suppose there are two
    ## sectors; sector 1 has 3 units; sector 2 has 2 units. To make
    ## per sector summaries, the following factors can be used to
    ## split the unit data: 1 1 1 2 2.
    ##
    ## Generating such factors first requires to know the number of
    ## nodes at each level in a format identical to the 'nodes'
    ## argument of simpf(). [In the previous example, the number of
    ## nodes would be 'list(2, c(3, 2))'.] Then, the factors are
    ## obtained by repeating a sequence the same length as the number
    ## of nodes at one level [2] according to the number of nodes at
    ## the level below [c(3, 2)].
    ##
    ## 1. Calculation of the number of nodes: the main idea is to
    ## create a unique ID for each node by pasting together the
    ## elements in the rows of 'ilevs'. This is not required for the
    ## lowest level (the entities), though, since they are known to
    ## all be different.
    fx <- vector("list", nlevels1p)
    fx[[nlevels1p]] <- factor(ilevs[, nlevels1p]) # entity level

    fnodes <- nnodes <- vector("list", nlevels)

    for (i in nlevels:1)
    {
        fx[[i]] <- factor(apply(ilevs[, seq.int(i), drop = FALSE], 1,
                                paste, collapse = ""))
        ## 'as.vector' below is used to get rid of names
        nnodes[[i]] <- as.vector(sapply(split(fx[[i + 1]], fx[[i]]),
                                        function(x) length(unique(x))))
    }

    ## 2. Generation of the factors. Following the rule described
    ## above, this could simply be
    ##
    ##   fnodes <- lapply(nnodes, function(x) rep(seq_along(x), x))
    ##
    ## However, this will not work if rows of data are not sorted per
    ## level. (In the example above, if the data of unit 3 of sector 1
    ## is at the bottom of the matrix of data, then the factors need
    ## to be 1 1 2 2 1.)
    ##
    ## The solution is actually simple: converting the entity level
    ## factors ('fx[[nlevels]]') to integers will assure that any
    ## summary made using these factors will be sorted. This done, it
    ## is possible to use the command above for the upper levels.
    fnodes[[nlevels]] <- as.integer(fx[[nlevels]])
    fnodes[-nlevels] <-
        lapply(nnodes[-nlevels], function(x) rep(seq_along(x), x))

    ## === PER ENTITY SUMMARIES ===
    ##
    ## Individual weighted averages. It could happen that an entity
    ## has no observations, for example when applying the model on
    ## claim amounts. In such a situation, put the total weight of the
    ## entity and the weighted average both equal to zero. That way,
    ## the premium will be equal to the credibility weighted average,
    ## as it should, but the entity will otherwise have no
    ## contribution in the calculations.
    weights.s <- rowSums(weights, na.rm = TRUE)
    ratios.w <- ifelse(weights.s > 0, rowSums(weights * ratios, na.rm = TRUE) / weights.s, 0)

    ## === EFFECTIVE NUMBER OF NODES ===
    ##
    ## Given the possibility to have whole levels with no data, as
    ## explained above, it is necessary to count the *effective*
    ## number of nodes in each level, that is the number of nodes with
    ## data. This comes this late since it relies on 'weights.s'.
    ##
    ## Where 'nnodes' is a list where each element is a vector of the
    ## number of nodes for each classification of the level above,
    ## object 'eff.nnodes' will be a simple vector of the total number
    ## of non "empty" nodes at each level. That is, if it weren't for
    ## entities with no data, the effective number of nodes would
    ## simply be 'sapply(nnodes, sum)'.
    eff.nnodes <- numeric(nlevels)
    w <- weights.s
    for (i in nlevels:1)
    {
        eff.nnodes[i] <- sum(w > 0)      # effective number of nodes
        w <- tapply(w, fnodes[[i]], sum) # running totals
    }

    ## === DENOMINATORS OF VARIANCE ESTIMATORS ===
    ##
    ## The denominators for all the variance estimators never
    ## change. The denominator at one level is equal to the total
    ## number of nodes at that level minus the total number of nodes
    ## at the level above. At the lowest level (the denominator of
    ## s^2), this is
    ##
    ##   number of (non NA) ratios - (effective) number of entities.
    ##
    ## The number of (non missing) ratios is not included in
    ## 'eff.nnodes'.  For the portfolio level, the denominator is
    ##
    ##   (effective) number of "sectors" - 1
    ##
    ## The 1 neither is included in 'eff.nnodes'.
    denoms <- diff(c(1, eff.nnodes, sum(!is.na(ratios))))

    ## Final sanity checks
    if (any(!denoms))
        stop("there must be at least two nodes at every level")
    if (ncol(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")

    ## === ESTIMATION OF s^2 ===
    s2 <-  sum(weights * (ratios - ratios.w)^2, na.rm = TRUE) /
        denoms[nlevels1p]

    ## === ESTIMATION OF THE OTHER STRUCTURE PARAMETERS ===
    ##
    ## Create vectors to hold values to be computed at each level
    ## (from portfolio to entity), namely: the total node weights, the
    ## node weighted averages, the between variance and the node
    ## credibility factors.
    ##
    ## Only credibility factors are not computed for the portfolio
    ## level, hence this list is one shorter than the others.
    ##
    ## At the entity level: total weight is the sum of the natural
    ## weights, weighted average uses the natural weights and between
    ## variance is s^2.
    ##
    ## All upper levels: total weight is the sum of the credibility
    ## factors of the level below, weighted average uses credibility
    ## factors, between variance estimated recursively and credibility
    ## factor use total weight of the level, between variance of the
    ## level below (hence the within variance) and between variance of
    ## the current level.
    ##
    ## The vector of variance is initialized with s^2 to provide
    ## starting values for the iterative estimation procedure.
    tweights <- vector("list", nlevels1p)       # total level weights
    wmeans <- vector("list", nlevels1p)         # weighted averages
    b <- rep(s2, nlevels1p)                     # variance estimators
    cred <- vector("list", nlevels)             # credibility factors

    ## Values already computed at the entity level.
    tweights[[nlevels1p]] <- as.vector(weights.s); # rm(weights.s)
    wmeans[[nlevels1p]] <- as.vector(ratios.w);    # rm(ratios.w)

    ## Iterative estimation of the structure parameters
    .External("cm", cred, tweights, wmeans, fnodes, denoms, b, TOL, echo)

    ## Final credibility factors and weighted averages (computed with
    ## the latest structure parameters). If a variance estimator is equal
    ## to (or tends toward) zero then we estimate the means with the
    ## total level weights instead of the credibility factors.
    for (i in nlevels:1)
    {
        cred[[i]] <- 1/(1 + ifelse(is.na(b[i + 1]), 0, b[i + 1])/(b[i] * tweights[[i + 1]]))
        tweights[[i]] <- as.vector(tapply(ifelse(cred[[i]] > 0, cred[[i]], tweights[[i + 1]]),
                                          fnodes[[i]],
                                          sum))
        wmeans[[i]] <- ifelse(tweights[[i]] > 0,
                              as.vector(tapply(ifelse(cred[[i]] > 0, cred[[i]], tweights[[i + 1]]) * wmeans[[i + 1]],
                                               fnodes[[i]],
                                               sum) / tweights[[i]]),
                              0)
        if (b[i] < TOL^2) b[i] <- NA
    }
    
    ## Transfer level names to lists
    names(tweights) <- names(wmeans) <- names(b) <-
        c("portfolio", level.names)
    names(cred) <- names(nnodes) <- names(fnodes) <-
        level.names

    ## Results
    structure(list(means = wmeans,
                   weights = tweights,
                   variances = b,
                   cred = cred,
                   levels = levs,
                   nodes = nnodes,
                   ordering = fnodes,
                   call = Call),
              class = "cm")
}

print.cm <- function(x, ...)
{
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:", x$means[[1]], "\n")
    for (i in seq.int(nlevels))
    {
        if (i == 1)
            cat("  Between", level.names[i], "variance:",
                x$variances[i], "\n")
        else
            cat("  Within ", level.names[i - 1],
                "/Between ", level.names[i], " variance: ",
                x$variances[i], "\n", sep = "")
    }
    cat("  Within", level.names[nlevels], "variance:",
        x$variances[nlevels + 1],"\n", fill = TRUE)
}

predict.cm <- function(object, levels = NULL, ...)
{
    ## The credibility premium of a node at one level is equal to
    ##
    ##   p + z * (m - p)
    ##
    ## where 'p' is the credibility premium of the level above (or the
    ## collective premium for the portfolio), 'z' is the credibility
    ## factor of the node, and 'm' is the weighted average of the
    ## node.
    fnodes <- object$ordering
    cred <- object$cred
    means <- object$means
    level.names <- names(object$nodes)

    plevs <-
        if (is.null(levels))
            seq_along(level.names)
        else
            pmatch(levels, level.names)
    if (any(is.na(plevs)))
        stop("invalid level name")
    n <- max(plevs)

    res <- vector("list", n)

    ## First level credibility premiums
    res[[1]] <- means[[1]] + cred[[1]] * (means[[2]] - means[[1]])

    for (i in seq(2, length.out = n - 1))
    {
        p <- res[[i - 1]][fnodes[[i]]]
        res[[i]] <- p + cred[[i]] * (means[[i + 1]] - p)
    }

    structure(res[plevs], names = level.names[plevs], ...)
}

summary.cm <- function(object, levels = NULL, ...)
{
    level.names <- names(object$nodes)

    plevs <-
        if (is.null(levels))
            seq_along(level.names)
        else
            pmatch(levels, level.names)
    if (any(is.na(plevs)))
        stop("invalid level name")

    structure(list(means = object$means[c(1, plevs + 1)],
                   weights = object$weights[c(1, plevs + 1)],
                   variances = object$variances[sort(unique(c(plevs, plevs + 1)))],
                   cred = object$cred[plevs],
                   levels = object$levels[, seq.int(max(plevs)), drop = FALSE],
                   nodes = object$nodes[plevs],
                   call = object$call,
                   premiums = predict(object, levels)),
              class = c("summary.cm", class(object)), ...)
}

print.summary.cm <- function(x, ...)
{
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:", x$means[[1]], "\n")
    for (i in seq.int(nlevels))
    {
        if (i == 1)
            cat("  Between", level.names[i], "variance:",
                x$variances[i], "\n")
        else
            cat("  Within ", level.names[i - 1],
                "/Between ", level.names[i], " variance: ",
                x$variances[i], "\n", sep = "")
    }
    cat("  Within", level.names[nlevels], "variance:",
        x$variances[nlevels + 1],"\n", fill = TRUE)
    cat("Detailed premiums\n\n")
    for (i in seq.int(nlevels))
    {
        cat("  Level:", level.names[i], "\n")
        level.id <- match(level.names[i], colnames(x$levels))
        levs <- x$levels[, seq.int(level.id), drop = FALSE]
        m <- duplicated(apply(levs, 1, paste, collapse = ""))
        y <- cbind(as.matrix(levs[!m, , drop = FALSE]),
                   format(x$means[[i + 1]], ...),
                   format(x$weights[[i + 1]], ...),
                   format(x$cred[[i]], ...),
                   format(x$premiums[[i]], ...))
        colnames(y) <- c(colnames(levs),
                         "Ind. premium", "Weight",
                         "Cred. factor", "Cred. premium")
        rownames(y) <- rep("   ", nrow(y))
        print(y, quote = FALSE, right = TRUE, ...)
        cat("\n")
    }
    invisible(x)
}
