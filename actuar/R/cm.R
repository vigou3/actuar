### ===== actuar: an R package for Actuarial Science =====
###
### Main interface to credibility model fitting functions.
###
### AUTHORS: Louis-Philippe Pouliot, Tommy Ouellet,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>.

cm <- function(formula, data, ratios, weights, subset, design, ...)
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

    ## Sanity checks
    ##
    ## 1. only hierarchical interactions are allowed in 'formula';
    ## 2. hierarchical regression models are not supported.
    ## 3. if 'ratios' is missing, all columns of 'data' are taken to
    ##    be ratios, so 'weights' should also be missing;
    ##
    if (any(duplicated(level.numbers)))
        stop("unsupported interactions in 'formula'")
    if (nlevels > 1 && !missing(design))
        stop("hierarchical regression models are not supported")
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

    ## Dispatch to appropriate calculation function
    if (nlevels < 2)
    {
        if (missing(design))
        {
            ## *** The 'old.format = FALSE' argument is necessary in
            ## *** the deprecation phase of the output format of
            ## *** bstraub(). Delete once phase is over.
            res <- bstraub(ratios, weights, ..., old.format = FALSE)
            res$classification <- levs
            res$ordering <- list(seq_along(levs))
        }
        else
            res <- hache(ratios, weights, design, ...)
    }
    else
        res <- hierarc(ratios, weights, ilevs, ...)

    ## Transfer level names to lists
    names(res$weights) <- names(res$means) <- names(res$unbiased) <-
        c("portfolio", level.names)
    names(res$iterative) <- if (!is.null(res$iterative)) names(res$unbiased)
    names(res$nodes) <- names(res$ordering) <- level.names
    if (is.list(res$cred))
        names(res$cred) <- level.names

    ## Results
    class(res) <- c("cm", class(res))
    attr(res, "call") <- Call
    res
}

print.cm <- function(x, ...)
{
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    b <- if (is.null(x$iterative)) x$unbiased else x$iterative

    cat("Call:\n")
    print(attr(x, "call"))
    cat("\n")

    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:", x$means[[1]], "\n\n")
    for (i in seq.int(nlevels))
    {
        if (i == 1)
            cat("  Between", level.names[i], "variance:",
                b[i], "\n")
        else
            cat("  Within ", level.names[i - 1],
                "/Between ", level.names[i], " variance: ",
                b[i], "\n", sep = "")
    }
    cat("  Within", level.names[nlevels], "variance:",
        b[nlevels + 1],"\n", fill = TRUE)
}

predict.cm <- function(object, levels = NULL, ...)
{
    ## Convert the character 'levels' argument into numeric and defer
    ## to next method.
    level.names <- names(object$nodes)

    levels <-
        if (is.null(levels))
            seq_along(level.names)
        else
            pmatch(levels, level.names)
    if (any(is.na(levels)))
        stop("invalid level name")
    NextMethod()
}

summary.cm <- function(object, levels = NULL, ...)
{
    level.names <- names(object$nodes)

    if (length(level.names) == 1)
    {
        ## Single level cases (Bühlmann-Straub and Hachemeister):
        ## return the object with the following modifications: put
        ## credibility factors into a list and add a list of the
        ## credibility premiums.
        object$premiums <- list(predict(object))
        object$cred <- list(object$cred)
        class(object) = c("summary.cm", class(object))
    }
    else
    {
        ## Multi-level case (hierarchical): select result of the
        ## appropriate level(s).
        plevs <-
            if (is.null(levels))
                seq_along(level.names)
            else
                pmatch(levels, level.names)
        if (any(is.na(plevs)))
            stop("invalid level name")

        object$premiums <- predict(object, levels) # new element
        object$means <- object$means[c(1, plevs + 1)]
        object$weights <- object$weights[c(1, plevs + 1)]
        object$unbiased <- object$unbiased[sort(unique(c(plevs, plevs + 1)))]
        object$iterative <- object$iterative[sort(unique(c(plevs, plevs + 1)))]
        object$cred <- object$cred[plevs]
        object$classification <- object$classification[, seq.int(max(plevs)), drop = FALSE]
        object$nodes <- object$nodes[plevs]
        class(object) <- c("summary.cm", class(object))
    }
    object
}

print.summary.cm <- function(x, ...)
{
    nlevels <- length(x$nodes)
    level.names <- names(x$nodes)
    NextMethod()                        # print.cm()
    cat("Detailed premiums\n\n")
    for (i in seq.int(nlevels))
    {
        cat("  Level:", level.names[i], "\n")
        level.id <- match(level.names[i], colnames(x$classification))
        levs <- x$classification[, seq.int(level.id), drop = FALSE]
        m <- duplicated(apply(levs, 1, paste, collapse = ""))
        y <- cbind(as.matrix(levs[!m, , drop = FALSE]),
                   format(x$means[[i + 1]], ...),
                   format(x$weights[[i + 1]], ...),
                   format(x$cred[[i]], ...),
                   format(x$premiums[[i]], ...))
        colnames(y) <- c(colnames(levs),
                         "Indiv. mean", "Weight",
                         "Cred. factor", "Cred. premium")
        rownames(y) <- rep("   ", nrow(y))
        print(y, quote = FALSE, right = TRUE, ...)
        cat("\n")
    }
    invisible(x)
}
